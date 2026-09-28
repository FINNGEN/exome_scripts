version 1.0

workflow gnomad_wes_finns_chrom {
  input {
    File vcf_list
    String root_name                # output filename root, e.g. "gnomAD_v4_Finns_subset"
    String genotype_filter
    String variant_filter
    Int cpu_count
    Int? test_sample_count
    File norm_fasta
    File denials                  # sample IDs to remove, one per line (see ParallelPreFilterByRegion)
    File aliases                  # tab-delimited alias groups (one group per line, same members
                                   # share a group) — used to expand `denials` to every alias of
                                   # each denied ID before exclusion, since the ID actually present
                                   # in a given VCF's header may be an alias rather than the ID
                                   # recorded in the denials list (see ExpandDenials)
   }

  Array[String] vcf_files = read_lines(vcf_list)

  # Per-chromosome positions to exclude from PreFilter, keyed by chrom name.
  # chr19:5787204 is a repeat-expansion locus — 2,026 of 2,357 records in a
  # 40bp window are stacked at that single position. Reading through it
  # (whether via one unbroken stream or any chunk whose region spans across
  # it) deterministically triggers a BGZF read failure. An independent
  # byte-level integrity check of the raw compressed bytes at that exact
  # location (walking 1,406 consecutive BGZF blocks, each with its CRC32 and
  # size recomputed via plain Python zlib/struct, entirely outside
  # bcftools/htslib) found no structural or checksum errors — the data
  # itself is not corrupted; see CHR19_PREFILTER_INVESTIGATION.md for the
  # full record. Add more entries here if other loci are identified.
  Map[String, Array[String]] prefilter_exclude_positions = {"chr19": ["5787204"]}

  # Expand denials to cover alias variants before any per-chromosome PreFilter call uses it.
  # Dataset-independent, so this runs once for the whole workflow rather than once per chromosome.
  call ExpandDenials {
    input:
      denials = denials,
      aliases = aliases
  }

  scatter (vcf in vcf_files) {
    if (defined(test_sample_count)) {
      call SubsetSamples {
        input:
          input_vcf = vcf,
          sample_count = select_first([test_sample_count])
      }
    }

    String vcf_to_filter = if defined(SubsetSamples.subset_vcf)
                            then select_first([SubsetSamples.subset_vcf]) + ""
                            else vcf

    # Chrom name pulled from the filename (e.g. ..._chr19.vcf.gz -> "chr19"), used
    # only to look up prefilter_exclude_positions above — not a File/task input, so
    # this stays a cheap workflow-level string expression. Only indexes the map
    # once the key is confirmed present, since WDL map indexing throws on a
    # missing key.
    String vcf_chrom = sub(basename(vcf), ".*(chr[0-9XY]+).*", "$1")
    Array[String] chrom_exclude_positions = if vcf_chrom == "chr19" then prefilter_exclude_positions["chr19"] else []

    call ComputeStats as OriginalStats {
      input:
        input_vcf = vcf_to_filter,
        cpu_count = cpu_count
    }

    # size() here is a workflow-level expression: Cromwell resolves it with a metadata
    # lookup against the GCS object, not a call input, so it never triggers localization
    # the way a File-typed task input/declaration would.
    # disk_gb here is 2x the source size + a flat buffer, not 1x + buffer: unlike
    # the old single-stream PreFilter (one input, one output, nothing else on
    # disk), this task's concat step has every chunk_*.vcf.gz coexisting with
    # the growing concatenated output right up until cleanup runs at the very
    # end — peak local disk usage is close to 2x the final output size, not 1x.
    # A run across many chromosomes hit "Write failed, wrote -1 instead of N
    # bytes" partway through concat — a classic ENOSPC signature — because the
    # old formula only budgeted 1x + 20GB.
    call ParallelPreFilterByRegion {
      input:
        input_vcf = vcf_to_filter,
        positions = OriginalStats.positions,
        denials = ExpandDenials.expanded_denials,
        exclude_positions = write_lines(chrom_exclude_positions),
        cpu_count = cpu_count,
        disk_gb = ceil(size(vcf_to_filter, "GB")) * 2 + 20
    }

    call ParallelFilterByRegion {
      input:
        input_vcf = ParallelPreFilterByRegion.prefiltered_vcf + "",
        positions = OriginalStats.positions,
        genotype_filter = genotype_filter,
        variant_filter = variant_filter,
        cpu_count = cpu_count,
        norm_fasta = norm_fasta,
        disk_gb = ceil(size(ParallelPreFilterByRegion.prefiltered_vcf, "GB")) * 2 + 20
    }

    call ComputeStats as FilteredStats {
      input:
        input_vcf = ParallelFilterByRegion.filtered_vcf + "",
        cpu_count = cpu_count
    }
    
    call ValidateFiltering {
      input:
        original_sample_vcf = OriginalStats.sample_vcf,
        original_sample_vcf_tbi = OriginalStats.sample_vcf_tbi,
        original_stats = OriginalStats.stats,
        filtered_sample_vcf = FilteredStats.sample_vcf,
        filtered_sample_vcf_tbi = FilteredStats.sample_vcf_tbi,
        filtered_stats = FilteredStats.stats,
        genotype_filter = genotype_filter,
        variant_filter = variant_filter
    }
  }


  call SummaryStats {
    input:
      vcf_file_names = vcf_files,
      original_stats = OriginalStats.stats,
      filtered_stats = FilteredStats.stats
  }

  call ConcatVcfs {
    input:
      input_vcfs = ParallelFilterByRegion.filtered_vcf,   # Array[File] coerced to Array[String] — no localisation
      summary_report = SummaryStats.report,
      root_name = root_name,   # explicit output root — not derived from the raw input VCF's name
      disk_gb = ceil(size(ParallelFilterByRegion.filtered_vcf, "GB")) + 20
  }

  output {
    Array[File] filtered_vcfs = ParallelFilterByRegion.filtered_vcf
    Array[File] filtered_vcf_tbis = ParallelFilterByRegion.filtered_vcf_tbi
    Array[File] original_stats = OriginalStats.stats
    Array[File] filtered_stats = FilteredStats.stats
    Array[File] validation_reports = ValidateFiltering.report
    File report = ConcatVcfs.report
    File concatenated_vcf = ConcatVcfs.concatenated_vcf
    File concatenated_vcf_tbi = ConcatVcfs.concatenated_vcf_tbi
  }
}

task ComputeStats {
  input {
    String input_vcf
    Int cpu_count = 8
    Int disk_gb = 20
  }

  command <<<
  set -euo

  echo "=== Computing statistics and creating sample ==="

  CHUNKS=~{cpu_count}
  fuse_vcf=$(echo "~{input_vcf}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')

  # Get chromosome from first variant
  chrom=$(bcftools query -f '%CHROM\n' "$fuse_vcf" | head -n 1)
  echo "Chromosome: $chrom"

  # Get chromosome and contig length from index
  read chrom_idx contig_len < <(bcftools index -s "$fuse_vcf" | awk '{print $1, $2}')
  echo "Contig length: $contig_len"

  # Get first position
  first_pos=$(bcftools view -H "$fuse_vcf" | head -n 1 | cut -f2)
  echo "First variant position: $first_pos"

  # Binary search for last variant position
  echo "Starting binary search for last variant..."
  low=$first_pos
  high=$contig_len
  while (( low <= high )); do
    mid=$(( (low + high) / 2 ))
    if bcftools view -H -r "$chrom:$mid-$high" "$fuse_vcf" 2>/dev/null | head -n 1 | grep -q .; then
      low=$(( mid + 1 ))
    else
      high=$(( mid - 1 ))
    fi
  done

  # Get exact last position from narrow window
  search_start=$(( high > 10000 ? high - 10000 : first_pos ))
  last_pos=$(bcftools view -H -r "$chrom:$search_start-$high" "$fuse_vcf" | tail -n 1 | cut -f2)
  echo "Last variant position: $last_pos"
  echo ""

  # Split interval equally using Python linspace
  echo "Splitting into $CHUNKS equal regions for parallel processing..."
  python3 -c "import numpy as np; [open(f'region_chunk_{i:02d}','w').write(f'$chrom\t{int(s)}\t{int(e)}\n') for i,(s,e) in enumerate(zip(np.linspace($first_pos,$last_pos,$CHUNKS+1)[:-1], np.linspace($first_pos,$last_pos,$CHUNKS+1)[1:]))]"

  # Create processing script
  cat > extract_chunk.sh << 'SCRIPT_EOF'
  #!/bin/bash
  input_file="$1"
  region_file="$2"

  bcftools query -f '%POS\n' "$input_file" -R "$region_file" > "${region_file}.positions"
  SCRIPT_EOF
  chmod +x extract_chunk.sh

  echo "Extracting positions in parallel..."
  ls region_chunk_* | sort -V | parallel -j $CHUNKS "./extract_chunk.sh '$fuse_vcf' {}"

  echo "Concatenating and sorting position files..."
  cat region_chunk_*.positions | sort -n -u > positions.txt

  # Count total variants — via bcftools index -s (true record count), NOT via
  # positions.txt: that file is deduplicated (sort -n -u) for region-chunk
  # boundary purposes only, so wc -l on it silently undercounts any position
  # with more than one record (e.g. multiallelic sites pre-normalisation).
  variant_count=$(bcftools index -s "$fuse_vcf" | awk '{sum+=$3} END {print sum}')
  echo "Total variants: $variant_count"

  # Count samples
  n_samples=$(bcftools query -l "$fuse_vcf" | wc -l)
  echo "Total samples: $n_samples"

  # Cleanup intermediate files
  rm -f region_chunk_*.positions extract_chunk.sh region_chunk_*
  
  # Create a small sample VCF (100 variants) for validation
  echo "Creating sample VCF for validation (100 variants)..."
  bcftools view -h "$fuse_vcf" | bgzip -c > sample.vcf.gz
  bcftools view -H "$fuse_vcf" | head -n 100 | bgzip -c >> sample.vcf.gz
  tabix -p vcf sample.vcf.gz
  echo "Sample VCF created: $(bcftools view -H sample.vcf.gz | wc -l) variants"

  # Create stats file
  printf "variant_count\t%s\nn_samples\t%s\n" "$variant_count" "$n_samples" > stats.txt

  echo "=== Complete ==="
  >>>

  output {
    File positions = "positions.txt"
    File sample_vcf = "sample.vcf.gz"
    File sample_vcf_tbi = "sample.vcf.gz.tbi"
    File stats = "stats.txt"
  }

  runtime {
    memory: "4G"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu_count
    preemptible: 1
  }
}

# ---------------------------------------------------------------------------
# Expands a raw denials list (one ID per line) to include every alias of each
# denied ID, per the alias groups in `aliases` (same tab-delimited group-per-
# line format as resolve_mapping.py's DEFAULT_ALIASES / load_aliases). Needed
# because the ID recorded in the denials list and the ID actually present in a
# given exome VCF's header can be different aliases of the same participant —
# bcftools view -S ^denials only matches literal strings, so PreFilter would
# otherwise fail to exclude (or hard-error on) a denied sample whose header ID
# is an alias rather than the one on the denials list.
# ---------------------------------------------------------------------------
task ExpandDenials {
  input {
    File denials    # plain list, one ID per line
    File aliases    # tab-delimited alias groups, one group per line
  }

  command <<<
  set -euo pipefail
  python3 << 'PY'
denied = {line.strip() for line in open("~{denials}") if line.strip()}
expanded = set(denied)
with open("~{aliases}") as fh:
    for line in fh:
        ids = [x.strip() for x in line.strip().split("\t") if x.strip()]
        if len(ids) < 2:
            continue
        if denied & set(ids):
            expanded.update(ids)
with open("expanded_denials.txt", "w") as out:
    for id_ in sorted(expanded):
        out.write(id_ + "\n")
print(f"Expanded {len(denied)} denied IDs to {len(expanded)} IDs (incl. aliases)")
PY
  >>>

  output {
    File expanded_denials = "expanded_denials.txt"
  }

  runtime {
    memory: "2G"
    disks: "local-disk 10 HDD"
    cpu: 1
    preemptible: 1
  }
}

# ---------------------------------------------------------------------------
# ParallelPreFilterByRegion — same commands.txt / `parallel` / concat pattern
# as ParallelFilterByRegion below, reusing OriginalStats' positions for
# chunking, but running PreFilter's own logic (remove denied samples, AC==0,
# annotate IDs) per-chunk instead of as one unbroken full-chromosome stream.
#
# Replaces the old monolithic PreFilter task: that task's single
# uninterrupted `bcftools view` over the full remote/gcsfuse-mounted VCF was
# empirically confirmed (across 3 independent Cromwell runs, at 8GB and 32GB
# memory) to hit a deterministic BGZF read failure on chr19 — memory bumps
# and retries never fixed it, since the failure recurred at the identical
# byte offset regardless of memory allocation. A local prototype
# (scripts/tests/chunked_prefilter_prototype.sh) validated that bounding the
# same read into small regions — including the exact regions where the
# monolithic read always failed — lets every chunk complete cleanly.
# ---------------------------------------------------------------------------
task ParallelPreFilterByRegion {
  input {
    String input_vcf
    File positions
    File denials    # sample IDs to remove, one per line
    File exclude_positions    # positions to skip entirely, one per line (may be empty)
    Int cpu_count
    Int disk_gb
  }

  String base_name = basename(basename(basename(input_vcf, ".vcf.gz"), ".vcf.bgz"), ".bcf")
  Int memory_gb = 64

  command <<<
  set -euo pipefail

  input_file=$(echo "~{input_vcf}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')
  CHUNKS=~{cpu_count}

  echo "=== Parallel Pre-Filter by Region: remove denied samples, AC==0, and annotate IDs ==="
  echo "Input: $input_file"
  echo "Removing denied samples listed in ~{denials}"
  echo "(--force-samples: denials list may include IDs not present in this VCF's header)"
  echo "CPU cores: $CHUNKS"
  echo ""

  # Get chromosome from first variant
  # process substitution, not a direct pipe: head's early exit would otherwise
  # SIGPIPE bcftools and (under pipefail) kill the script with 141
  chrom=$(head -n 1 <(bcftools query -f '%CHROM\n' "$input_file"))
  echo "Chromosome: $chrom"
  echo ""

  # Drop any excluded positions from the position list first (for -T target
  # correctness — a literal excluded position should never appear in any
  # chunk's target list).
  n_excluded=$(wc -l < ~{exclude_positions})
  echo "Excluded positions: $n_excluded"
  if [[ "$n_excluded" -gt 0 ]]; then
    grep -vxFf ~{exclude_positions} ~{positions} > positions_filtered.txt
  else
    cp ~{positions} positions_filtered.txt
  fi

  # Normal N-way split, unchanged from before exclusion support existed.
  echo "Splitting positions into $CHUNKS chunks..."
  split -n l/$CHUNKS -d positions_filtered.txt pos_chunk_tmp_

  # For each resulting chunk, check whether any excluded position's value
  # falls strictly inside that chunk's [first_pos, last_pos] span — even
  # though the excluded position itself was already dropped above, its
  # NEIGHBORS on either side can easily still land in the same chunk, so the
  # chunk's numeric span can still cross over it. If so, sub-split just that
  # one chunk's position list into a before/after pair around the excluded
  # position, so no chunk's eventual -r region ever spans across it — -r does
  # index-based seeking followed by a SEQUENTIAL decompress through every
  # BGZF block in that span, independent of what's in -T, so filtering the
  # position out of -T alone would not be enough. Total chunk count grows by
  # one for each chunk that needed splitting; that's fine, per-chunk runtime
  # stays roughly constant either way.
  chunk_num=0
  for tmp_file in pos_chunk_tmp_*; do
    first_pos=$(head -n 1 "$tmp_file")
    last_pos=$(tail -n 1 "$tmp_file")

    split_files=("$tmp_file")
    if [[ "$n_excluded" -gt 0 ]]; then
      for excl_pos in $(sort -n ~{exclude_positions}); do
        if (( excl_pos > first_pos && excl_pos < last_pos )); then
          echo "  chunk spanning $first_pos-$last_pos crosses excluded position $excl_pos — splitting into 2"
          awk -v hi="$excl_pos" '$1<hi' "$tmp_file" > "${tmp_file}.before"
          awk -v lo="$excl_pos" '$1>lo' "$tmp_file" > "${tmp_file}.after"
          split_files=("${tmp_file}.before" "${tmp_file}.after")
          break
        fi
      done
    fi

    for f in "${split_files[@]}"; do
      [[ -s "$f" ]] || continue   # a before/after piece can be empty if the excluded position was near an edge
      variant_count=$(wc -l < "$f")
      f_first=$(head -n 1 "$f")
      f_last=$(tail -n 1 "$f")

      chunk_id=$(printf '%02d' $chunk_num)
      pos_name="pos_chunk_${chunk_id}"
      region_name="region_chunk_${chunk_id}"

      # Create position file with CHROM POS format (required by -T)
      awk -v chr="$chrom" '{print chr "\t" $1}' "$f" > "$pos_name"

      # Create region file for efficient seeking
      printf "%s:%d-%d\n" "$chrom" "$f_first" "$f_last" > "$region_name"

      echo "  chunk_${chunk_id}: $chrom:$f_first-$f_last ($variant_count variants)"
      chunk_num=$((chunk_num + 1))
    done
  done
  rm -f pos_chunk_tmp_*
  echo "Total chunks: $chunk_num"
  echo ""

  # Create commands file with each command fully self-contained
  echo "Creating commands for parallel processing..."
  rm -f commands.txt
  chunk_num=0
  for pos_file in pos_chunk_*; do
    chunk_id=$(printf '%02d' $chunk_num)
    region=$(cat "region_chunk_${chunk_id}")
    cat >> commands.txt << EOF
set -o pipefail && echo "Processing: chunk_${chunk_id} (region: $region)" && bcftools view "$input_file" -r "$region" -T "pos_chunk_${chunk_id}" -S ^~{denials} --force-samples -Ou | bcftools +fill-tags -Ou -- -t AC | bcftools view -i 'AC>0' -Ou | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o "chunk_${chunk_id}.vcf.gz" && echo "  ✓ Done: chunk_${chunk_id}"
EOF
    chunk_num=$((chunk_num + 1))
  done

  echo "Generated commands (first 3):"
  head -3 commands.txt
  echo ""

  echo "Processing chunks in parallel..."
  cat commands.txt | parallel -j $CHUNKS

  echo ""
  echo "Concatenating chunks..."
  ls chunk_*.vcf.gz | sort -V > chunk_list.txt
  bcftools concat -n -f chunk_list.txt -Oz -o ~{base_name}.prefiltered.vcf.gz

  echo "Indexing output..."
  tabix --threads $(nproc) -p vcf ~{base_name}.prefiltered.vcf.gz

  # Cleanup
  echo "Cleaning up temporary files..."
  rm -f chunk_*.vcf.gz chunk_list.txt pos_chunk_* region_chunk_* commands.txt positions_filtered.txt

  echo "=== Complete! ==="
  >>>

  output {
    File prefiltered_vcf = "~{base_name}.prefiltered.vcf.gz"
    File prefiltered_vcf_tbi = "~{base_name}.prefiltered.vcf.gz.tbi"
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu_count
    preemptible: 1
  }
}

task ParallelFilterByRegion {
  input {
    String input_vcf
    File positions
    String genotype_filter
    String variant_filter
    Int cpu_count
    File norm_fasta
    Int disk_gb
  }

  File norm_fasta_fai = norm_fasta + ".fai"
  String base_name = basename(basename(basename(input_vcf, ".vcf.gz"), ".vcf.bgz"), ".bcf")
  Int memory_gb = 64

  command <<<
  set -euo pipefail

  input_file=$(echo "~{input_vcf}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')
  CHUNKS=~{cpu_count}
  touch ~{norm_fasta_fai}

  echo "=== Parallel Filter by Region ==="
  echo "Input: $input_file"
  echo "Genotype filter: ~{genotype_filter}"
  echo "Variant filter: ~{variant_filter}"
  echo "CPU cores: $CHUNKS"
  echo ""
  
  # Get chromosome from first variant
  # process substitution, not a direct pipe: head's early exit would otherwise
  # SIGPIPE bcftools and (under pipefail) kill the script with 141
  chrom=$(head -n 1 <(bcftools query -f '%CHROM\n' "$input_file"))
  echo "Chromosome: $chrom"
  echo ""
  
  # Split positions into equal chunks locally
  echo "Splitting positions into $CHUNKS equal chunks..."
  split -n l/$CHUNKS -d ~{positions} pos_chunk_tmp_
  
  # Create properly formatted position and region files
  chunk_num=0
  for tmp_file in pos_chunk_tmp_*; do
    variant_count=$(wc -l < "$tmp_file")
    first_pos=$(head -n 1 "$tmp_file")
    last_pos=$(tail -n 1 "$tmp_file")
    
    chunk_id=$(printf '%02d' $chunk_num)
    pos_name="pos_chunk_${chunk_id}"
    region_name="region_chunk_${chunk_id}"
    
    # Create position file with CHROM POS format (required by -T)
    awk -v chr="$chrom" '{print chr "\t" $1}' "$tmp_file" > "$pos_name"
    
    # Create region file for efficient seeking
    printf "%s:%d-%d\n" "$chrom" "$first_pos" "$last_pos" > "$region_name"
    
    echo "  chunk_${chunk_id}: $chrom:$first_pos-$last_pos ($variant_count variants)"
    chunk_num=$((chunk_num + 1))
  done
  rm -f pos_chunk_tmp_*
  echo ""

  # Create commands file with each command fully self-contained
  echo "Creating commands for parallel processing..."
  rm -f commands.txt
  chunk_num=0
  for pos_file in pos_chunk_*; do
    chunk_id=$(printf '%02d' $chunk_num)
    region=$(cat "region_chunk_${chunk_id}")
    cat >> commands.txt << EOF
set -o pipefail && echo "Processing: chunk_${chunk_id} (region: $region)" && bcftools view "$input_file" -r "$region" -T "pos_chunk_${chunk_id}" -Ou | bcftools norm -f '~{norm_fasta}' -m -any -c x -Ou | bcftools +setGT -Ou -- -t q -n . -i '~{genotype_filter}' | bcftools +fill-tags -Ou -- -t AC | bcftools view -i '~{variant_filter}' -Oz -o "chunk_${chunk_id}.vcf.gz" && echo "  ✓ Done: chunk_${chunk_id}"
EOF
    chunk_num=$((chunk_num + 1))
  done
  
  echo "Generated commands (first 3):"
  head -3 commands.txt
  echo ""

  echo "Processing chunks in parallel..."
  cat commands.txt | parallel -j $CHUNKS

  echo ""
  echo "Concatenating chunks..."
  ls chunk_*.vcf.gz | sort -V > chunk_list.txt
  bcftools concat -n -f chunk_list.txt -Oz -o ~{base_name}.filtered.vcf.gz
  
  echo "Indexing output..."
  tabix --threads $(nproc) -p vcf ~{base_name}.filtered.vcf.gz
  
  # Cleanup
  echo "Cleaning up temporary files..."
  rm -f chunk_*.vcf.gz chunk_list.txt pos_chunk_* region_chunk_* commands.txt
  
  echo "=== Complete! ==="
  >>>

  output {
    File filtered_vcf = "~{base_name}.filtered.vcf.gz"
    File filtered_vcf_tbi = "~{base_name}.filtered.vcf.gz.tbi"
  }

  runtime {
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu_count
    preemptible: 1
  }
}


task ValidateFiltering {
  input {
    File original_sample_vcf
    File original_sample_vcf_tbi
    File original_stats
    File filtered_sample_vcf
    File filtered_sample_vcf_tbi
    File filtered_stats
    String genotype_filter
    String variant_filter
  }

  Int disk_size = 10

  command <<<
  set -euo

  echo "=== Validating VCF Filtering ===" > report.txt
  echo "Sample-based validation (using up to 10k variants)" >> report.txt
  echo "" >> report.txt

  # Touch indices
  touch ~{original_sample_vcf_tbi}
  touch ~{filtered_sample_vcf_tbi}

  # Read variant counts
  orig_count=$(grep "variant_count" ~{original_stats} | cut -f2)
  filt_count=$(grep "variant_count" ~{filtered_stats} | cut -f2)
  
  echo "Original variant count: $orig_count" >> report.txt
  echo "Filtered variant count: $filt_count" >> report.txt
  
  if [[ "$orig_count" -gt 0 ]]; then
    pct_kept=$(awk "BEGIN {printf \"%.1f\", ($filt_count / $orig_count) * 100}")
    pct_dropped=$(awk "BEGIN {printf \"%.1f\", (($orig_count - $filt_count) / $orig_count) * 100}")
  else
    pct_kept="0.0"
    pct_dropped="0.0"
  fi
  
  echo "Kept: ${pct_kept}%" >> report.txt
  echo "Dropped: ${pct_dropped}%" >> report.txt
  echo "" >> report.txt

  # Test 1: Check genotypes were set to missing
  echo "Test 1: Checking genotypes were set to missing (~{genotype_filter})" >> report.txt
  echo "--------------------------------------------------------------" >> report.txt

  found=0
  bcftools query -i 'GT!="mis"' -f '%CHROM\t%POS[\t%SAMPLE\t%GT\t%DP\t%GQ]\n' "~{original_sample_vcf}" 2>/dev/null | \
    awk '($4 < 10 || $5 < 20)' | head -100 | \
    while read chrom pos sample gt dp gq; do
      if bcftools view -H -r "$chrom:$pos" "~{filtered_sample_vcf}" 2>/dev/null | grep -q .; then
        after_gt=$(bcftools query -s "$sample" -r "$chrom:$pos" -f '[\t%GT]\n' "~{filtered_sample_vcf}" 2>/dev/null | tr -d '\t' | tr -d ' ')
        
        if [[ "$after_gt" == "./." ]] || [[ "$after_gt" == "." ]]; then
          echo "✓ PASS: $chrom:$pos sample=$sample GT: $gt → $after_gt" >> report.txt
        else
          echo "✗ FAIL: $chrom:$pos sample=$sample GT: $gt → $after_gt (not missing)" >> report.txt
        fi
        found=1
        break
      fi
    done

  if [[ $found -eq 0 ]]; then
    echo "⚠ NOTE: All low-quality genotypes resulted in variants being removed (AC=0)" >> report.txt
  fi

  echo "" >> report.txt

  # Test 2: Check AC field exists
  echo "Test 2: Checking AC field is present in filtered VCF" >> report.txt
  echo "-----------------------------------------------------" >> report.txt

  first_var=$(bcftools view -H "~{filtered_sample_vcf}" | head -1 | awk '{print $1":"$2}')

  if [[ -n "$first_var" ]]; then
    ac_value=$(bcftools query -r "$first_var" -f '%AC\n' "~{filtered_sample_vcf}" 2>/dev/null)
    
    if [[ -n "$ac_value" ]]; then
      echo "✓ PASS: AC field present (example: $first_var AC=$ac_value)" >> report.txt
    else
      echo "✗ FAIL: AC field missing in filtered VCF" >> report.txt
    fi
  else
    echo "⚠ NOTE: No variants in filtered file" >> report.txt
  fi

  echo "" >> report.txt

  # Test 3: Check variant IDs
  echo "Test 3: Checking variant IDs (CHROM_POS_REF_ALT format)" >> report.txt
  echo "--------------------------------------------------------" >> report.txt

  sample_id=$(bcftools view -H "~{filtered_sample_vcf}" 2>/dev/null | head -1 | awk '{print $3}')
  if [[ "$sample_id" =~ ^[^_]+_[0-9]+_.+_.+$ ]]; then
    echo "✓ PASS: IDs formatted as CHROM_POS_REF_ALT (example: $sample_id)" >> report.txt
  else
    echo "✗ FAIL: ID not in expected format (got: $sample_id)" >> report.txt
  fi

  echo "" >> report.txt
  echo "=== Validation Complete ===" >> report.txt
  >>>

  output {
    File report = "report.txt"
  }

  runtime {
    disks: "local-disk ~{disk_size} HDD"
  }
}

task SubsetSamples {
  input {
    String input_vcf
    Int sample_count
    Int disk_gb = 50
  }

  String base_name = basename(basename(basename(input_vcf, ".vcf.gz"), ".vcf.bgz"), ".bcf")
  String output_vcf = base_name + ".subset_" + sample_count + "samples.vcf.gz"
  String output_tbi = base_name + ".subset_" + sample_count + "samples.vcf.gz.tbi"

  command <<<
  set -euo

  THREADS=$(nproc)
  NCOLS=$((9 + ~{sample_count}))
  fuse_vcf=$(echo "~{input_vcf}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')

  echo "=== Subsetting VCF to first ~{sample_count} samples (columns 1-$NCOLS) ==="

  # Extract metadata lines (##) - keep intact
  echo "Processing metadata..."
  bcftools view -h "$fuse_vcf" | grep "^##" | bgzip -@ $THREADS -c > ~{output_vcf}

  # Extract and cut column header line (#CHROM)
  echo "Processing column header..."
  bcftools view -h "$fuse_vcf" | grep "^#CHROM" | cut -f 1-$NCOLS | bgzip -@ $THREADS -c >> ~{output_vcf}

  # Extract and cut body (skip all header lines starting with #)
  echo "Processing body..."
  zcat "$fuse_vcf" | grep -v "^#" | cut -f 1-$NCOLS | bgzip -@ $THREADS -c >> ~{output_vcf}

  echo "Indexing..."
  tabix --threads $THREADS -p vcf ~{output_vcf}
  echo "=== Complete ==="
  >>>

  output {
    File subset_vcf = output_vcf
    File subset_vcf_tbi = output_tbi
  }

  runtime {
    memory: "8G"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: 4
    preemptible: 1
  }
}

task ConcatVcfs {
  input {
    Array[String] input_vcfs   # Array[File] coerced to Array[String] at call site — no localisation
    File summary_report
    String root_name
    Int disk_gb
  }

  command <<<
  set -euo
  THREADS=$(nproc)

  echo "=== Concatenating VCF shards ==="
  sed 's|gs://[^/]*/|/mnt/disks/gcs/|' ~{write_lines(input_vcfs)} | \
    awk -F'/' '{print $NF"\t"$0}' | sort -V | cut -f2- > sorted_vcf_list.txt
  bcftools concat --threads $THREADS -f sorted_vcf_list.txt -Oz -o ~{root_name}.QC_ANNOTATED.vcf.gz

  echo "Indexing..."
  tabix --threads $THREADS -p vcf ~{root_name}.QC_ANNOTATED.vcf.gz

  cp ~{summary_report} ~{root_name}.QC_ANNOTATED.report.txt
  echo "=== Complete ==="
  >>>

  output {
    File concatenated_vcf = "~{root_name}.QC_ANNOTATED.vcf.gz"
    File concatenated_vcf_tbi = "~{root_name}.QC_ANNOTATED.vcf.gz.tbi"
    File report = "~{root_name}.QC_ANNOTATED.report.txt"
  }

  runtime {
    memory: "8 GB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: 16
  }
}

task SummaryStats {
  # Builds a per-chromosome stats table. Also guesses a root name from the input
  # VCF filenames (root_name.txt / root_name output below), but that guess is
  # informational only — the actual output filename root comes from the
  # workflow-level `root_name` input, not from this task.
  input {
    Array[String] vcf_file_names
    Array[File] original_stats
    Array[File] filtered_stats
  }

  command <<<
  set -euo

  echo "=== Creating Summary Statistics ==="

  # Derive root name from first VCF: strip _chrN.vcf.gz (or .vcf.bgz/.bcf)
  first_vcf=$(head -1 ~{write_lines(vcf_file_names)})
  filename=$(basename "$first_vcf")
  root=$(echo "$filename" | sed -E 's/_chr[0-9XYxy]+\.(vcf\.gz|vcf\.bgz|bcf)$//')
  echo "$root" > root_name.txt
  echo "Root name: $root"

  # Create header
  echo -e "chromosome\toriginal_variants\tfiltered_variants\tpercent_dropped\toriginal_samples\tfiltered_samples\tsamples_removed" > summary.report.txt
  
  # Process each VCF file
  idx=0
  while IFS= read -r vcf_file; do
    # Extract chromosome from filename (e.g., chr1.vcf.gz -> chr1)
    filename=$(basename "$vcf_file")
    chrom=$(echo "$filename" | sed -E 's/.*_(chr[0-9XY]+).*\.vcf\.(gz|bgz)$/\1/' | sed -E 's/^.*chr([0-9XY]+).*$/chr\1/' | sed -E 's/.*[^a-zA-Z](chr[0-9XY]+)[^a-zA-Z].*$/\1/' | head -1)
    
    # Fallback: try to extract just the chromosome number/letter
    if [[ -z "$chrom" ]] || [[ "$chrom" == "$filename" ]]; then
      chrom=$(echo "$filename" | grep -oE '(chr)?[0-9XY]+' | head -1)
    fi
    
    # If still empty, use filename
    if [[ -z "$chrom" ]]; then
      chrom=$(echo "$filename" | sed 's/\..*//')
    fi
    
    # Read original and filtered counts from stats files (not VCFs)
    orig_stats_file=$(echo "~{sep=' ' original_stats}" | cut -d' ' -f$((idx+1)))
    filt_stats_file=$(echo "~{sep=' ' filtered_stats}" | cut -d' ' -f$((idx+1)))
    
    orig_count=$(grep "variant_count" "$orig_stats_file" | cut -f2)
    filt_count=$(grep "variant_count" "$filt_stats_file" | cut -f2)
    
    # Calculate percent dropped
    if [[ "$orig_count" -gt 0 ]]; then
      pct_dropped=$(awk "BEGIN {printf \"%.2f\", (($orig_count - $filt_count) / $orig_count) * 100}")
    else
      pct_dropped="0.00"
    fi

    # Sample count is the same cohort across every chromosome of this dataset —
    # tracked per-chromosome as a sanity check, not summed into the TOTAL row below.
    orig_samples=$(grep "n_samples" "$orig_stats_file" | cut -f2)
    filt_samples=$(grep "n_samples" "$filt_stats_file" | cut -f2)
    samples_removed=$((orig_samples - filt_samples))

    echo -e "${chrom}\t${orig_count}\t${filt_count}\t${pct_dropped}\t${orig_samples}\t${filt_samples}\t${samples_removed}" >> summary.report.txt
    idx=$((idx + 1))
  done < <(cat << 'EOF'
~{sep='\n' vcf_file_names}
EOF
)

  # Add totals row
  echo "" >> summary.report.txt
  total_orig=$(grep "variant_count" ~{sep=' ' original_stats} | cut -f2 | awk '{sum+=$1} END {print sum}')
  total_filt=$(grep "variant_count" ~{sep=' ' filtered_stats} | cut -f2 | awk '{sum+=$1} END {print sum}')
  total_pct=$(awk "BEGIN {printf \"%.2f\", (($total_orig - $total_filt) / $total_orig) * 100}")

  echo -e "TOTAL\t${total_orig}\t${total_filt}\t${total_pct}\t${orig_samples}\t${filt_samples}\t${samples_removed}" >> summary.report.txt
  
  echo ""
  echo "Summary Table:"
  column -t summary.report.txt
  
  echo "=== Complete ==="
  >>>

  output {
    String root_name = read_string("root_name.txt")
    File report = "summary.report.txt"
  }

  runtime {
    memory: "2G"
    disks: "local-disk 10 HDD"
    cpu: 1
    preemptible: 1
  }
}

