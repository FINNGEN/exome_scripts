version 1.0

# exome_qc.wdl — unified QC pipeline reconciling wes_chrom.wdl and single_file_qc.wdl.
#
# Each dataset maps a name to a manifest File listing its VCF path(s), one per
# line (same manifest-file convention as wes_chrom.json/daly_qc.json's
# vcf_list input):
#   - manifest has exactly one line  -> treated as a "block" file (all
#                                       chromosomes in one VCF, single_file_qc's
#                                       use case) -> auto-split by SplitByChromosome
#   - manifest has more than one line -> already split per chromosome
#                                       (wes_chrom's use case) -> used as-is,
#                                       no split step run
# Both cases converge on the same per-chromosome QC chain (ComputeStats /
# PreFilter / ParallelFilterByRegion / ValidateFiltering — copied verbatim
# from wes_chrom.wdl on this fuse-test branch).
#
# Nested scatters aren't used (they don't work reliably on this Cromwell
# setup) — instead there are three flat scatters plus a gather task:
#   1. scatter over datasets: optionally split block files into per-chromosome
#      pieces; flatten() combines every dataset's chromosome files into one
#      flat list. Each dataset also exposes dataset_name + chrom_count (a
#      free expression, no task call) so GatherByDataset can later work out
#      which flat result belongs to which dataset.
#   2. scatter over that flat list: the per-chromosome QC chain, dataset-
#      agnostic — identical to wes_chrom.wdl.
#   3. GatherByDataset expands dataset_name/chrom_count into a per-chromosome
#      tag column and regroups the flat per-chromosome results back by
#      dataset (plain paste/awk grouping, writing one manifest per dataset)
#      — this replaces the per-dataset grouping a nested scatter would
#      otherwise have given for free.
#   4. scatter over the distinct datasets from step 3: SummaryStats +
#      ConcatVcfs once per dataset, reading that dataset's own manifests.
#
# Large VCF reads use the GCS FUSE mount (/mnt/disks/gcs/...) instead of
# Cromwell File localisation, matching wes_chrom.wdl/daly_qc.wdl on this same
# branch. NOTE: this assumes the execution VM has gcsfuse-mounted whichever
# bucket(s) these paths live in — unverified for production paths as of this
# branch; see conversation notes before running for real.

workflow exome_qc {
  input {
    Array[Array[String]] datasets   # [dataset_name, manifest_file] rows — manifest lists that dataset's VCF path(s), one per line
    String genotype_filter
    String variant_filter
    Int    cpu_count
    File   norm_fasta
    File?  denials                # optional: sample IDs to remove, one per line (see PreFilter)
    Int?   test_sample_count
  }

  # ── Phase 1: optional per-dataset split, flat scatter (no nesting) ────────
  # Map isn't usable here (this Cromwell setup rejects both scattering over a
  # Map and the keys() function) — datasets is Array[Array[String]] instead
  # ([name, manifest_file] rows), indexed via range(length(...)), the same
  # pattern already proven working in exome_duplicates.wdl's vcf_pairs input.
  scatter (i in range(length(datasets))) {
    String      dataset_name  = datasets[i][0]
    Array[File] dataset_files = read_lines(datasets[i][1])
    Boolean     needs_split   = length(dataset_files) == 1

    if (needs_split) {
      call SplitByChromosome {
        input:
          input_vcf = dataset_files[0],
          cpu_count = cpu_count
      }
    }

    # Unifies the "split a block file" and "already pre-split" cases into one
    # Array[File] — select_first unwraps SplitByChromosome's optional output
    # (it only ran when needs_split), so both branches stay File-typed with
    # no File/String coercion needed at this level.
    Array[File] chrom_vcfs = if needs_split
                              then select_first([SplitByChromosome.chrom_vcfs])
                              else dataset_files

    # How many chromosome files this dataset contributes to all_chrom_vcfs
    # below — a free expression (no task call). Combined with dataset_name,
    # this is all GatherByDataset needs to reconstruct which flat result
    # belongs to which dataset, without a dedicated per-dataset tagging task.
    Int chrom_count = length(chrom_vcfs)
  }

  Array[File] all_chrom_vcfs = flatten(chrom_vcfs)

  # ── Phase 2: flat per-chromosome QC chain, dataset-agnostic ───────────────
  scatter (vcf in all_chrom_vcfs) {
    if (defined(test_sample_count)) {
      call SubsetSamples {
        input:
          input_vcf = vcf,
          sample_count = select_first([test_sample_count])
      }
    }

    File vcf_to_filter = if defined(SubsetSamples.subset_vcf)
                          then select_first([SubsetSamples.subset_vcf])
                          else vcf

    call ComputeStats as OriginalStats {
      input:
        input_vcf = vcf_to_filter,
        cpu_count = cpu_count
    }

    call PreFilter {
      input:
        input_vcf = vcf_to_filter,
        cpu_count = cpu_count,
        denials = denials
    }

    call ParallelFilterByRegion {
      input:
        input_vcf = PreFilter.prefiltered_vcf + "",
        positions = OriginalStats.positions,
        genotype_filter = genotype_filter,
        variant_filter = variant_filter,
        cpu_count = cpu_count,
        norm_fasta = norm_fasta
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

  # ── Phase 3: regroup flat per-chromosome results back by dataset ─────────
  call GatherByDataset {
    input:
      dataset_names  = dataset_name,
      chrom_counts   = chrom_count,
      vcfs           = ParallelFilterByRegion.filtered_vcf,
      original_stats = OriginalStats.stats,
      filtered_stats = FilteredStats.stats
  }

  # ── Phase 4: one SummaryStats + ConcatVcfs per dataset ───────────────────
  scatter (i in range(length(GatherByDataset.distinct_datasets))) {
    String      ds_name      = GatherByDataset.distinct_datasets[i]
    Array[File] ds_vcfs      = read_lines(GatherByDataset.vcf_manifests[i])
    Array[File] ds_origstats = read_lines(GatherByDataset.origstats_manifests[i])
    Array[File] ds_filtstats = read_lines(GatherByDataset.filtstats_manifests[i])

    call SummaryStats {
      input:
        vcf_file_names = ds_vcfs,
        original_stats = ds_origstats,
        filtered_stats = ds_filtstats
    }

    call ConcatVcfs {
      input:
        input_vcfs = ds_vcfs,   # Array[File] coerced to Array[String] — no localisation
        summary_report = SummaryStats.report,
        root_name = ds_name
    }
  }

  output {
    Array[String] dataset_names         = GatherByDataset.distinct_datasets
    Array[File]   concatenated_vcfs     = ConcatVcfs.concatenated_vcf
    Array[File]   concatenated_vcf_tbis = ConcatVcfs.concatenated_vcf_tbi
    Array[File]   reports               = ConcatVcfs.report
    Array[String] chrom_dataset_tags    = GatherByDataset.chrom_dataset_tags
    Array[File]   validation_reports    = ValidateFiltering.report
  }
}


# ---------------------------------------------------------------------------
# Split a single "block" VCF (all chromosomes in one file) into one VCF per
# chromosome, so it can join the same per-chromosome QC chain as cohorts that
# already arrive pre-split. Chromosome discovery + chr-prefix renaming ported
# from single_file_qc.wdl's FilterByChromosome task — the QC filtering itself
# (norm/genotype-mask/AC recompute/variant filter) happens later, per
# chromosome, via PreFilter + ParallelFilterByRegion below, not here.
# ---------------------------------------------------------------------------
task SplitByChromosome {
  input {
    String input_vcf
    Int    cpu_count = 8
    Int    disk_gb   = 100
  }

  command <<<
  set -euo pipefail
  fuse_vcf=$(echo "~{input_vcf}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')
  CHUNKS=~{cpu_count}

  echo "=== Splitting block VCF by chromosome ==="

  mapfile -t chromosomes < <(bcftools index -s "$fuse_vcf" | awk '{print $1}')
  num_chroms=${#chromosomes[@]}
  echo "Found $num_chroms chromosomes: ${chromosomes[*]}"
  if [[ $num_chroms -eq 0 ]]; then
    echo "Error: No chromosomes found in VCF" >&2
    exit 1
  fi

  # Detect chromosome naming: if VCF uses non-chr names, rename to chr prefix in output
  FIRST_CHROM="${chromosomes[0]}"
  RENAME_TO_CHR=""
  if [[ "$FIRST_CHROM" != chr* ]]; then
    echo "VCF uses non-chr chromosome names — renaming to chr prefix in output"
    RENAME_TO_CHR=$(mktemp)
    for i in $(seq 1 22) X Y MT M; do
      echo "$i chr$i" >> "$RENAME_TO_CHR"
    done
  fi
  RENAME_STEP="cat"
  [[ -n "$RENAME_TO_CHR" ]] && RENAME_STEP="bcftools annotate --rename-chrs ${RENAME_TO_CHR} -Ou"

  # Generate per-chromosome scripts and run in parallel.
  # Variables expand at generation time — no quoting or function-export issues.
  SCRIPT_DIR=$(mktemp -d)
  for chrom in "${chromosomes[@]}"; do
    safe_chrom=$(echo "$chrom" | sed 's/[*:\/]/_/g')
    output="split_${safe_chrom}.vcf.gz"
    cat > "${SCRIPT_DIR}/run_${safe_chrom}.sh" << SCRIPT
#!/bin/bash
set -euo pipefail
echo "Splitting chromosome: ${chrom}"
printf '${chrom}\t0\t9999999999\n' > "${output}.region.bed"
bcftools view -R "${output}.region.bed" "${fuse_vcf}" | \\
    tr -d '\0' | \\
    ${RENAME_STEP} | \\
    bcftools view -Oz -o "${output}"
tabix -p vcf "${output}"
rm -f "${output}.region.bed"
echo "Done: ${output}"
SCRIPT
  done

  echo "Splitting $num_chroms chromosomes across $CHUNKS parallel jobs..."
  ls "${SCRIPT_DIR}"/run_*.sh | parallel --bar -j "$CHUNKS" 'bash {}'
  rm -rf "$SCRIPT_DIR"

  echo "Split into $(ls split_*.vcf.gz | wc -l) / $num_chroms per-chromosome files"
  >>>

  output {
    Array[File] chrom_vcfs     = glob("split_*.vcf.gz")
    Array[File] chrom_vcf_tbis = glob("split_*.vcf.gz.tbi")
  }

  runtime {
    memory: "8G"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu_count
    preemptible: 1
  }
}


# ---------------------------------------------------------------------------
# Regroups the flat, dataset-agnostic per-chromosome QC results back by
# dataset, writing one manifest file per distinct dataset for each of
# vcfs/original_stats/filtered_stats. Plain paste + awk grouping — stands in
# for what a nested scatter would otherwise give for free.
#
# dataset_names/chrom_counts (one entry per dataset, from Phase 1 — no extra
# task needed to produce them) are first expanded into a per-chromosome tag
# column matching the order flatten() produced all_chrom_vcfs/the Phase 2
# results in, then that tag column is used to group vcfs/original_stats/
# filtered_stats by dataset.
# ---------------------------------------------------------------------------
task GatherByDataset {
  input {
    Array[String] dataset_names    # one per dataset
    Array[Int]    chrom_counts     # length(chrom_vcfs) per dataset, same order as dataset_names
    Array[String] vcfs             # Array[File] coerced to Array[String] — paths only, not read here
    Array[String] original_stats
    Array[String] filtered_stats
  }

  command <<<
  set -euo pipefail

  # Expand dataset_names by chrom_counts into a per-chromosome tag column —
  # e.g. names=[WES,BOTNIA] counts=[24,1] -> WES x24 then BOTNIA x1, matching
  # the order flatten() used to build all_chrom_vcfs in the workflow.
  paste ~{write_lines(dataset_names)} ~{write_lines(chrom_counts)} | \
    while IFS=$'\t' read -r name count; do
      for i in $(seq 1 "$count"); do echo "$name"; done
    done > tags.txt

  paste tags.txt ~{write_lines(vcfs)} ~{write_lines(original_stats)} ~{write_lines(filtered_stats)} > combined.tsv

  # Distinct dataset names, in first-seen order
  awk -F'\t' '!seen[$1]++ {print $1}' combined.tsv > distinct_datasets.txt

  n=0
  while IFS= read -r ds; do
    idx=$(printf '%03d' $n)
    awk -F'\t' -v ds="$ds" '$1==ds {print $2}' combined.tsv > "vcf_manifest_${idx}.txt"
    awk -F'\t' -v ds="$ds" '$1==ds {print $3}' combined.tsv > "origstats_manifest_${idx}.txt"
    awk -F'\t' -v ds="$ds" '$1==ds {print $4}' combined.tsv > "filtstats_manifest_${idx}.txt"
    n=$((n + 1))
  done < distinct_datasets.txt

  echo "Grouped $(wc -l < combined.tsv) chromosome results into $n dataset(s)"
  >>>

  output {
    Array[String] distinct_datasets   = read_lines("distinct_datasets.txt")
    Array[String] chrom_dataset_tags  = read_lines("tags.txt")
    Array[File]   vcf_manifests       = glob("vcf_manifest_*.txt")
    Array[File]   origstats_manifests = glob("origstats_manifest_*.txt")
    Array[File]   filtstats_manifests = glob("filtstats_manifest_*.txt")
  }

  runtime {
    memory: "2G"
    disks: "local-disk 10 HDD"
    cpu: 1
    preemptible: 1
  }
}


# ---------------------------------------------------------------------------
# Below: identical to wes_chrom.wdl's per-chromosome QC chain on this same
# fuse-test branch (String inputs, GCS-FUSE reads instead of File
# localisation). Copied in rather than imported so exome_qc.wdl stays a
# single self-contained WDL, matching the rest of this repo's convention —
# no cross-file imports are used anywhere else in this codebase.
# ---------------------------------------------------------------------------

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
  ls region_chunk_* | sort -V | parallel -j $CHUNKS './extract_chunk.sh "$fuse_vcf" {}'

  echo "Concatenating and sorting position files..."
  cat region_chunk_*.positions | sort -n -u > positions.txt

  # Count total variants
  variant_count=$(wc -l < positions.txt)
  echo "Total variants: $variant_count"

  # Cleanup intermediate files
  rm -f region_chunk_*.positions extract_chunk.sh region_chunk_*

  # Create a small sample VCF (100 variants) for validation
  echo "Creating sample VCF for validation (100 variants)..."
  bcftools view -h "$fuse_vcf" | bgzip -c > sample.vcf.gz
  bcftools view -H "$fuse_vcf" | head -n 100 | bgzip -c >> sample.vcf.gz
  tabix -p vcf sample.vcf.gz
  echo "Sample VCF created: $(bcftools view -H sample.vcf.gz | wc -l) variants"

  # Create stats file
  echo -e "variant_count\t$variant_count" > stats.txt

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

task PreFilter {
  input {
    String input_vcf
    Int    cpu_count
    File?  denials    # optional: sample IDs to remove, one per line
    Int    disk_gb = 50
  }

  String base_name = basename(basename(basename(input_vcf, ".vcf.gz"), ".vcf.bgz"), ".bcf")

  command <<<
  set -euo
  THREADS=$(nproc)
  fuse_vcf=$(echo "~{input_vcf}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')

  echo "=== Pre-Filter: remove denied samples, AC==0, and annotate IDs ==="
  echo "Input: ~{input_vcf}"

  DENIALS_FLAG=""
  ~{if defined(denials) then "DENIALS_FLAG='-S ^" + denials + "'" else ""}
  if [[ -n "$DENIALS_FLAG" ]]; then
    echo "Removing denied samples listed in ~{if defined(denials) then denials else ""}"
  fi

  TARGET_SIZE=$(stat -c%s "$fuse_vcf")
  bcftools view --threads $THREADS $DENIALS_FLAG "$fuse_vcf" -Ou | \
    bcftools +fill-tags --threads $THREADS -Ou -- -t AC | \
    bcftools view --threads $THREADS -i 'AC>0' -Ou | \
    bcftools annotate --threads $THREADS --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Oz | \
    pv -s $TARGET_SIZE -N "prefilter" -i 60 > ~{base_name}.prefiltered.vcf.gz

  tabix -p vcf ~{base_name}.prefiltered.vcf.gz
  echo "=== Complete! ==="
  >>>

  output {
    File prefiltered_vcf = "~{base_name}.prefiltered.vcf.gz"
    File prefiltered_vcf_tbi = "~{base_name}.prefiltered.vcf.gz.tbi"
  }

  runtime {
    memory: "8 GB"
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
    Int disk_gb = 100
  }

  File norm_fasta_fai = norm_fasta + ".fai"
  String base_name = basename(basename(basename(input_vcf, ".vcf.gz"), ".vcf.bgz"), ".bcf")
  Int memory_gb = 64

  command <<<
  set -euo

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
  chrom=$(bcftools query -f '%CHROM\n' "$input_file" | head -n 1)
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
echo "Processing: chunk_${chunk_id} (region: $region)" && bcftools view "$input_file" -r "$region" -T "pos_chunk_${chunk_id}" -Ou | bcftools norm -f '~{norm_fasta}' -m -any -c x -Ou | bcftools +setGT -Ou -- -t q -n . -i '~{genotype_filter}' | bcftools +fill-tags -Ou -- -t AC | bcftools view -i '~{variant_filter}' -Oz -o "chunk_${chunk_id}.vcf.gz" && echo "  ✓ Done: chunk_${chunk_id}"
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
  tabix -p vcf ~{base_name}.filtered.vcf.gz

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
  tabix -p vcf ~{output_vcf}
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
    Int disk_gb = 100
  }

  command <<<
  set -euo
  THREADS=$(nproc)

  echo "=== Concatenating VCF shards ==="
  sed 's|gs://[^/]*/|/mnt/disks/gcs/|' ~{write_lines(input_vcfs)} | \
    awk -F'/' '{print $NF"\t"$0}' | sort -V | cut -f2- > sorted_vcf_list.txt
  bcftools concat --threads $THREADS -f sorted_vcf_list.txt -Oz -o ~{root_name}.QC_ANNOTATED.vcf.gz

  echo "Indexing..."
  tabix -p vcf ~{root_name}.QC_ANNOTATED.vcf.gz

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
  echo -e "chromosome\toriginal_variants\tfiltered_variants\tpercent_dropped" > summary.report.txt

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

    echo -e "${chrom}\t${orig_count}\t${filt_count}\t${pct_dropped}" >> summary.report.txt
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

  echo -e "TOTAL\t${total_orig}\t${total_filt}\t${total_pct}" >> summary.report.txt

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
