version 1.0

workflow bge_qc {
  input {
    File          vcf_list                                    # one GCS path per line
    File          rename_file                                 # TSV: FINNGENID_finngen(1) FINNGENID_biobank(2) SAMPLE_ID(3) ...
    String        root_name                                   # output filename root, e.g. "Blended_Genome_Exome_scizophrenia_bipolar_controls"
    String        filter_expression = "FILTER~'NO_HQ_GENOTYPES'"
    Int           cpu_count        = 8
    Int           vcf_max_gb       = 25  # size of the largest VCF; drives disk allocation
    Int           chunk_multiplier = 3   # chunks = cpu_count × multiplier; tune to stay within auth token window
    File          denials                # sample IDs to remove, one per line (see ParallelFilter)
    File          aliases                # tab-delimited alias groups (one group per line, same members
                                          # share a group) — used to expand `denials` to every alias of
                                          # each denied ID before exclusion, since the ID actually present
                                          # in a given VCF's header (post-rename) may be an alias rather
                                          # than the ID recorded in the denials list (see ExpandDenials)
  }

  Array[String] vcfs = read_lines(vcf_list)

  # Expand denials to cover alias variants before any per-chromosome ParallelFilter call uses it.
  # Dataset-independent, so this runs once for the whole workflow rather than once per chromosome.
  call ExpandDenials {
    input:
      denials = denials,
      aliases = aliases
  }

  scatter (vcf in vcfs) {

    call ComputeStats {
      input: vcf = vcf
    }

    call AnnotateAndRename {
      input:
        vcf         = vcf,
        rename_file = rename_file,
        cpu_count   = cpu_count,
        vcf_max_gb  = vcf_max_gb
    }

    call ParallelFilter {
      input:
        input_vcf          = AnnotateAndRename.output_vcf + "",
        filter_expression  = filter_expression,
        cpu_count          = cpu_count,
        vcf_max_gb         = vcf_max_gb,
        chunk_multiplier   = chunk_multiplier,
        denials            = ExpandDenials.expanded_denials
    }

    call ComputeStats as FilteredStats {
      input: vcf = ParallelFilter.filtered_vcf + ""
    }

    call ValidateFiltering {
      input:
        original_vcf      = vcf,
        filtered_vcf      = ParallelFilter.filtered_vcf + "",
        original_stats    = ComputeStats.stats,
        filtered_stats    = FilteredStats.stats,
        rename_report     = AnnotateAndRename.rename_report,
        filter_expression = filter_expression
    }
  }

  call SummaryStats {
    input:
      vcf_file_names = vcfs,
      original_stats = ComputeStats.stats,
      filtered_stats = FilteredStats.stats
  }

  call SortAndMerge {
    input:
      vcf_files      = ParallelFilter.filtered_vcf,   # Array[File] coerced to Array[String] — no localisation
      summary_report = SummaryStats.report,
      root_name      = root_name,                     # explicit output root — not derived from the raw input VCF's name
      cpu_count      = cpu_count
  }

  output {
    Array[File] original_stats      = ComputeStats.stats
    Array[File] filtered_stats      = FilteredStats.stats
    Array[File] validation_reports  = ValidateFiltering.report
    File        merged_vcf          = SortAndMerge.merged_vcf
    File        merged_vcf_tbi      = SortAndMerge.merged_vcf_tbi
    File        report              = SortAndMerge.report
  }
}


# ── ComputeStats ─────────────────────────────────────────────────────────────
# Reads only the remote index — no VCF download.

task ComputeStats {
  input {
    String vcf
  }

  command <<<
  set -euo pipefail
  fuse_vcf=$(echo "~{vcf}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')

  chrom=$(bcftools index -s "$fuse_vcf" | head -n 1 | cut -f1)
  variant_count=$(bcftools index -s "$fuse_vcf" | awk '{sum+=$3} END {print sum}')
  n_samples=$(bcftools query -l "$fuse_vcf" | wc -l)

  echo "Chromosome:    $chrom"
  echo "Variants:      $variant_count"
  echo "Samples:       $n_samples"

  printf "chrom\t%s\nvariant_count\t%s\nn_samples\t%s\n" \
    "$chrom" "$variant_count" "$n_samples" > stats.txt
  >>>

  output {
    File stats = "stats.txt"
  }

  runtime {
    memory:      "2G"
    disks:       "local-disk 10 HDD"
    cpu:         1
    preemptible: 1
  }
}


# ── AnnotateAndRename ─────────────────────────────────────────────────────────
# 1. Fetches VCF header from remote GCS.
# 2. Adds any missing FILTER definitions.
# 3. Builds a SAMPLE_ID → FINNGENID_finngen rename map for samples in this VCF.
# 4. Applies header annotation + sample rename in a single bcftools reheader call.
#    bcftools reheader rewrites only the header blocks and copies BGZF body
#    verbatim — fast regardless of VCF size.
# 5. Reports how many samples were renamed.

task AnnotateAndRename {
  input {
    String vcf
    File   rename_file
    Int    cpu_count  = 8
    Int    vcf_max_gb = 25
  }

  Int    disk_gb  = vcf_max_gb + 10   # output ≈ input (reheader copies body verbatim)
  String out_name = basename(vcf)

  command <<<
  set -euo pipefail
  VCF=$(echo "~{vcf}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')
  OUT="~{out_name}"

  # ── 1. Fetch remote header ────────────────────────────────────────────
  bcftools view -h "$VCF" > old_header.txt

  # ── 2. Identify and append missing FILTER definitions ────────────────
  declare -A FILTER_DESC
  FILTER_DESC["NO_HQ_GENOTYPES"]="Site has no high quality variant genotypes (GQ>=20, DP>=10, AB>=0.2 for het). If only one genotype exists the filter is not applied."
  FILTER_DESC["ExcessHet"]="Excess heterozygosity (z-score < -4.5, phred 54.69). Suggests mapping errors or contamination."
  FILTER_DESC["LowQual"]="QUALapprox too low (< 60 for SNPs, < 69 for Indels)."
  FILTER_DESC["EXCESS_ALLELES"]="Excess alternate alleles above threshold (e.g. >100)."
  FILTER_DESC["OUTSIDE_OF_TARGETS"]="Exome only. Site is outside the target intervals of the assay."

  mapfile -t header_filters < <(awk -F'[=,]' '/^##FILTER=/{print $3}' old_header.txt)

  > missing_filters.txt
  for filter in "${!FILTER_DESC[@]}"; do
      found=0
      for h in "${header_filters[@]}"; do [[ "$h" == "$filter" ]] && found=1 && break; done
      [[ $found -eq 0 ]] && \
          echo "##FILTER=<ID=$filter,Description=\"${FILTER_DESC[$filter]}\">" \
          >> missing_filters.txt
  done

  if [[ -s missing_filters.txt ]]; then
      echo "Adding $(wc -l < missing_filters.txt) missing FILTER header(s)"
      (head -n -1 old_header.txt; cat missing_filters.txt; tail -n 1 old_header.txt) \
          > new_header.txt
  else
      echo "No missing FILTER headers"
      cp old_header.txt new_header.txt
  fi

  # ── 3. Build rename map: SAMPLE_ID (col 3) → FINNGENID_finngen (col 1) ─
  bcftools query -l "$VCF" > current_samples.txt
  n_samples=$(wc -l < current_samples.txt)

  # Skip header row; extract col3=SAMPLE_ID, col1=FINNGENID_finngen
  awk 'NR>1{print $3"\t"$1}' ~{rename_file} > all_mappings.txt

  > rename.txt
  while IFS= read -r sample; do
      new=$(awk -F'\t' -v s="$sample" '$1==s{print $2; exit}' all_mappings.txt)
      [[ -n "$new" ]] && echo -e "$sample\t$new" >> rename.txt
  done < current_samples.txt

  # ── Resolve duplicate new IDs ────────────────────────────────────────
  # Two old IDs mapping to the same new ID → add _dup / _dup2 etc. suffix.
  awk 'BEGIN{OFS="\t"} {
      count[$2]++
      if      (count[$2] == 1) print $1, $2
      else if (count[$2] == 2) print $1, $2 "_dup"
      else                     print $1, $2 "_dup" (count[$2]-1)
  }' rename.txt > rename_dedup.txt
  mv rename_dedup.txt rename.txt

  n_dups=$(awk '{print $2}' rename.txt | sort | uniq -d | wc -l)
  if [[ $n_dups -gt 0 ]]; then
      echo "WARNING: $n_dups duplicate new ID(s) resolved with _dup suffix:"
      awk '{print $2}' rename.txt | sort | uniq -d | while read -r dup_id; do
          grep -w "$dup_id\|${dup_id}_dup" rename.txt
      done
  fi

  n_renamed=$(wc -l < rename.txt)
  echo "Samples in VCF:     $n_samples"
  echo "Samples renamed:    $n_renamed / $n_samples"
  echo "Samples kept as-is: $(( n_samples - n_renamed ))"
  echo "Duplicate new IDs:  $n_dups"
  printf "n_samples\t%s\nn_renamed\t%s\nn_dups\t%s\n" \
      "$n_samples" "$n_renamed" "$n_dups" > rename_report.txt

  # ── 4. Apply header + rename in one bcftools reheader pass ────────────
  if [[ -s rename.txt ]]; then
      bcftools reheader -h new_header.txt -s rename.txt -o "$OUT" "$VCF"
  else
      echo "No samples to rename — applying header annotation only"
      bcftools reheader -h new_header.txt -o "$OUT" "$VCF"
  fi

  tabix --threads $(nproc) -p vcf "$OUT"
  echo "Done: $OUT"
  >>>

  output {
    File output_vcf     = out_name
    File output_vcf_tbi = out_name + ".tbi"
    File rename_report  = "rename_report.txt"
  }

  runtime {
    memory:      "4G"
    disks:       "local-disk ~{disk_gb} HDD"
    cpu:         cpu_count
    preemptible: 1
  }
}


# ── ParallelFilter ────────────────────────────────────────────────────────────
# Splits the chromosome into cpu_count equal position windows, filters each
# chunk in parallel with GNU parallel, then concatenates.
# Input VCF is local (File) since it was produced by AnnotateAndRename.

# ── ExpandDenials ─────────────────────────────────────────────────────────────
# Expands a raw denials list (one ID per line) to include every alias of each
# denied ID, per the alias groups in `aliases` (same tab-delimited group-per-
# line format as resolve_mapping.py's DEFAULT_ALIASES / load_aliases). Needed
# because the ID recorded in the denials list and the ID actually present in a
# given VCF's header (post-AnnotateAndRename) can be different aliases of the
# same participant — bcftools view -S ^denials only matches literal strings,
# so ParallelFilter would otherwise fail to exclude (or hard-error on) a
# denied sample whose header ID is an alias rather than the one on the list.

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
    memory:      "2G"
    disks:       "local-disk 10 HDD"
    cpu:         1
    preemptible: 1
  }
}


task ParallelFilter {
  input {
    String input_vcf
    String filter_expression
    File   denials    # sample IDs to remove, one per line
    Int    cpu_count        = 8
    Int    vcf_max_gb       = 25
    Int    chunk_multiplier = 3
  }

  Int disk_gb = vcf_max_gb * 2 + 20  # N chunk outputs + filtered VCF; input read via GCS FUSE

  command <<<
  set -euo pipefail
  fuse_vcf=$(echo "~{input_vcf}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')
  CHUNKS=$(( ~{cpu_count} * ~{chunk_multiplier} ))   # chunk count now purely for parallelism (FUSE has no auth-token expiry window)

  # ── 1. Get chromosome and contig length from index ───────────────────
  chrom=$(bcftools index -s "$fuse_vcf" | awk '{print $1; exit}')
  contig_len=$(bcftools index -s "$fuse_vcf" | awk '{print $2; exit}')
  echo "Contig: $chrom  length: $contig_len  chunks: $CHUNKS"

  # ── 2. Build chunk BED files spanning the full contig ────────────────
  # Splitting on contig_len works for both WGS and exome VCFs — chunks with
  # no variants are handled gracefully (tabix seeks return immediately, empty
  # chunks are excluded before concat).
  chunk_size=$(( (contig_len + CHUNKS - 1) / CHUNKS ))   # ceiling division
  for i in $(seq 0 $(( CHUNKS - 1 ))); do
      start=$(( i * chunk_size + 1 ))
      # Last chunk end is 10× the contig length — bcftools clamps to the last variant.
      end=$(( i == CHUNKS - 1 ? contig_len * 10 : start + chunk_size - 1 ))
      printf "%s\t%d\t%d\n" "$chrom" "$start" "$end" > "chunk_$(printf '%02d' $i).bed"
  done

  # ── 3. Build one command per line, run all chunks in parallel ────────
  echo "Removing denied samples listed in ~{denials}"
  echo "(--force-samples: denials list may include IDs not present in this VCF's header)"
  for i in $(seq 0 $(( CHUNKS - 1 ))); do
      bed="chunk_$(printf '%02d' $i).bed"
      out="chunk_$(printf '%02d' $i).vcf.gz"
      echo "bcftools view -S ^~{denials} --force-samples '$fuse_vcf' -R $bed -Ou | bcftools view -e \"~{filter_expression}\" -Ou | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o $out && echo \"done chunk $i / $((CHUNKS-1))\""
  done > chunks.sh

  parallel -j $(nproc) < chunks.sh

  # ── 4. Concatenate all chunks and index ──────────────────────────────
  ls chunk_*.vcf.gz | sort -V > chunk_list.txt
  echo "Concatenating $(wc -l < chunk_list.txt) chunks..."
  bcftools concat -n -f chunk_list.txt -Oz -o filtered.vcf.gz
  rm -f chunk_*.vcf.gz chunk_*.vcf.gz.tbi chunk_*.bed

  tabix --threads $(nproc) -p vcf filtered.vcf.gz

  variant_count=$(bcftools index -s filtered.vcf.gz | awk '{sum+=$3} END {print sum}')
  echo "Variants after filtering: $variant_count"
  echo "$variant_count" > variant_count.txt
  >>>

  output {
    File filtered_vcf     = "filtered.vcf.gz"
    File filtered_vcf_tbi = "filtered.vcf.gz.tbi"
    Int  variant_count    = read_int("variant_count.txt")
  }

  runtime {
    memory:      "8G"
    disks:       "local-disk ~{disk_gb} HDD"
    cpu:         cpu_count
    preemptible: 1
  }
}


# ── ValidateFiltering ─────────────────────────────────────────────────────────
# Queries original and filtered VCFs remotely (no localisation).
# Includes rename summary from AnnotateAndRename.

task ValidateFiltering {
  input {
    String original_vcf
    String filtered_vcf
    File   original_stats
    File   filtered_stats
    File   rename_report
    String filter_expression
  }

  command <<<
  set -euo pipefail
  fuse_filtered_vcf=$(echo "~{filtered_vcf}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')

  echo "=== Validating VCF Filtering ===" > report.txt
  echo "" >> report.txt

  # ── Variant counts ────────────────────────────────────────────────────
  orig_count=$(grep "variant_count" ~{original_stats} | cut -f2)
  filt_count=$(grep "variant_count" ~{filtered_stats} | cut -f2)
  pct_kept=$(awk "BEGIN {printf \"%.1f\", ($filt_count / $orig_count) * 100}")
  pct_dropped=$(awk "BEGIN {printf \"%.1f\", (($orig_count - $filt_count) / $orig_count) * 100}")

  echo "Original variants: $orig_count" >> report.txt
  echo "Filtered variants: $filt_count  (kept ${pct_kept}%  dropped ${pct_dropped}%)" >> report.txt
  echo "" >> report.txt

  # ── Test 1: no variants matching filter expression remain ─────────────
  echo "Test 1: Filter expression check (~{filter_expression})" >> report.txt
  bcftools view -H -i '~{filter_expression}' "$fuse_filtered_vcf" 2>/dev/null \
      | head -n 1 > check_filter.txt 2>/dev/null || true
  if [[ -s check_filter.txt ]]; then
      echo "✗ FAIL: Found variants matching filter expression (should be 0)" >> report.txt
  else
      echo "✓ PASS: No variants matching filter expression in filtered VCF" >> report.txt
  fi
  echo "" >> report.txt

  # ── Test 2: variant ID format (CHROM_POS_REF_ALT) ────────────────────
  echo "Test 2: Variant ID format check" >> report.txt
  bcftools view -H "$fuse_filtered_vcf" 2>/dev/null \
      | head -n 1 > check_id.txt 2>/dev/null || true
  sample_id=$(awk '{print $3}' check_id.txt)
  if [[ "$sample_id" =~ ^[^_]+_[0-9]+_.+_.+$ ]]; then
      echo "✓ PASS: IDs in CHROM_POS_REF_ALT format (e.g. $sample_id)" >> report.txt
  else
      echo "✗ FAIL: Unexpected ID format (got: $sample_id)" >> report.txt
  fi
  echo "" >> report.txt

  # ── Sample rename summary ─────────────────────────────────────────────
  echo "=== Sample Rename Summary ===" >> report.txt
  n_samples=$(grep "^n_samples" ~{rename_report} | cut -f2)
  n_renamed=$(grep "^n_renamed" ~{rename_report} | cut -f2)
  n_dups=$(grep "^n_dups"    ~{rename_report} | cut -f2)
  echo "Total samples:              $n_samples" >> report.txt
  echo "Samples renamed:            $n_renamed / $n_samples" >> report.txt
  echo "Duplicate new IDs resolved: $n_dups"    >> report.txt
  echo "" >> report.txt
  echo "=== Validation Complete ===" >> report.txt
  >>>

  output {
    File report = "report.txt"
  }

  runtime {
    memory:      "2G"
    disks:       "local-disk 10 HDD"
    cpu:         1
    preemptible: 1
  }
}


# ── SummaryStats ──────────────────────────────────────────────────────────────
# Builds a per-chromosome stats table. Also guesses a root name from the input
# VCF filenames (root_name.txt / root_name output below), but that guess is
# informational only — the actual output filename root comes from the
# workflow-level `root_name` input, not from this task.

task SummaryStats {
  input {
    Array[String] vcf_file_names
    Array[File]   original_stats
    Array[File]   filtered_stats
  }

  command <<<
  set -euo

  first_vcf=$(head -1 ~{write_lines(vcf_file_names)})
  root=$(basename "$first_vcf" | sed -E 's/_chr[0-9XYxy]+\.(vcf\.gz|vcf\.bgz|bcf)$//')
  echo "$root" > root_name.txt
  echo "Root name: $root"

  echo -e "chromosome\toriginal_variants\tfiltered_variants\tpercent_dropped\toriginal_samples\tfiltered_samples\tsamples_removed" > summary.report.txt

  idx=0
  while IFS= read -r vcf_file; do
      filename=$(basename "$vcf_file")
      chrom=$(echo "$filename" | grep -oE 'chr[0-9XY]+' | head -1)
      [[ -z "$chrom" ]] && chrom=$(echo "$filename" | sed 's/\..*//')

      orig_stats_file=$(echo "~{sep=' ' original_stats}" | cut -d' ' -f$((idx+1)))
      filt_stats_file=$(echo "~{sep=' ' filtered_stats}" | cut -d' ' -f$((idx+1)))

      orig_count=$(grep "variant_count" "$orig_stats_file" | cut -f2)
      filt_count=$(grep "variant_count" "$filt_stats_file" | cut -f2)
      pct_dropped=$(awk "BEGIN {printf \"%.2f\", (($orig_count - $filt_count) / $orig_count) * 100}")

      # Sample count is the same cohort across every chromosome of this dataset —
      # tracked per-chromosome as a sanity check, not summed into the TOTAL row below.
      orig_samples=$(grep "n_samples" "$orig_stats_file" | cut -f2)
      filt_samples=$(grep "n_samples" "$filt_stats_file" | cut -f2)
      samples_removed=$((orig_samples - filt_samples))

      echo -e "${chrom}\t${orig_count}\t${filt_count}\t${pct_dropped}\t${orig_samples}\t${filt_samples}\t${samples_removed}" >> summary.report.txt
      idx=$((idx + 1))
  done < ~{write_lines(vcf_file_names)}

  total_orig=$(grep "variant_count" ~{sep=' ' original_stats} | cut -f2 | awk '{sum+=$1} END {print sum}')
  total_filt=$(grep "variant_count" ~{sep=' ' filtered_stats} | cut -f2 | awk '{sum+=$1} END {print sum}')
  total_pct=$(awk "BEGIN {printf \"%.2f\", (($total_orig - $total_filt) / $total_orig) * 100}")
  echo -e "TOTAL\t${total_orig}\t${total_filt}\t${total_pct}\t${orig_samples}\t${filt_samples}\t${samples_removed}" >> summary.report.txt

  echo ""; echo "Summary:"; column -t summary.report.txt
  >>>

  output {
    String root_name = read_string("root_name.txt")
    File   report    = "summary.report.txt"
  }

  runtime {
    memory:      "2G"
    disks:       "local-disk 10 HDD"
    cpu:         1
    preemptible: 1
  }
}


# ── SortAndMerge ──────────────────────────────────────────────────────────────
# Concatenates per-chromosome filtered VCFs in WDL scatter order (chromosome
# order is guaranteed by the input vcf_list order — no filename sort needed).

task SortAndMerge {
  input {
    Array[String] vcf_files   # Array[File] coerced to Array[String] at call site — no localisation
    File          summary_report
    String        root_name
    Int           cpu_count = 8
    Int           disk_gb   = 100
  }

  command <<<
  set -euo pipefail

  # Use WDL scatter order — chromosome order is guaranteed by vcf_list input order.
  sed 's|gs://[^/]*/|/mnt/disks/gcs/|' ~{write_lines(vcf_files)} > vcf_list.txt
  echo "Merging ~{length(vcf_files)} VCFs in scatter order"

  bcftools concat -n -f vcf_list.txt -Oz -o ~{root_name}.QC_ANNOTATED.vcf.gz
  tabix --threads $(nproc) -p vcf ~{root_name}.QC_ANNOTATED.vcf.gz
  cp ~{summary_report} ~{root_name}.QC_ANNOTATED.report.txt

  echo "Done: ~{root_name}.QC_ANNOTATED.vcf.gz"
  >>>

  output {
    File merged_vcf     = "~{root_name}.QC_ANNOTATED.vcf.gz"
    File merged_vcf_tbi = "~{root_name}.QC_ANNOTATED.vcf.gz.tbi"
    File report         = "~{root_name}.QC_ANNOTATED.report.txt"
  }

  runtime {
    memory:      "8G"
    disks:       "local-disk ~{disk_gb} HDD"
    cpu:         cpu_count
    preemptible: 1
  }
}
