#!/usr/bin/env python3
"""sync_data_release.py — mirror the curated subset of data/ into release/.

The two trees are not a straight mirror: release/ only carries a fixed set of
deliverables, some renamed, some split by file type. Mapping (see data/README.txt):

  data/finngen_R14_exome_id_mapping.tsv  -> release/data/finngen_R14_exome_id_mapping.tsv
  data/ld/finngen_R14_exome.ld.tsv.gz            -> release/data/finngen_R14_exome.ld.tsv.gz
  data/ld/finngen_R14_exome.ld_annotated.tsv.gz  -> release/data/finngen_R14_exome.ld_annotated.tsv.gz
  data/QC/*.vcf.gz(.tbi)                 -> release/data/qc_vcf_full/
  data/QC/*.report.txt                   -> release/documentation/
  data/renamed_vcf_chr/*                 -> release/data/renamed_vcf_chr/
  data/plink_fg_merged_chr/*             -> release/data/plink_fg_merged_chr/
  data/finngen_R14_exome_id_mapping_flowchart.png -> release/documentation/FG_EXOME_resolved_flowchart.png

The public readme is not read from data/ at all — it's this repo's own
data/finngen_R14_exome_readme.md, pushed to release/finngen_R14_exome_readme
(no extension, matching the existing object name).

Everything else under data/ (annotate/, ld/'s other intermediate files, README.txt)
is internal-only and is never copied. Each mapped destination is kept in sync with
`gsutil rsync -d`, so files removed from the data/ side (e.g. denied/duplicate
samples dropped during a rerun) are also removed from release/ — nothing stale
is left behind. Single files whose source doesn't exist yet (e.g. the annotated
LD output, pending a rerun) are skipped with a warning rather than failing the
whole sync.

Usage:
  python sync_data_release.py gs://fg-3/exome_v2/data gs://fg-3/exome_v2/release [--dry-run] [--skip-validate]
"""
import argparse, pathlib, subprocess, sys

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent

# (label, local path relative to repo root, dst_rel under release/)
LOCAL_FILES = [
    ("readme", "data/finngen_R14_exome_readme.md", "finngen_R14_exome_readme"),
    ("id_mapping_flowchart", "data/finngen_R14_exome_id_mapping_flowchart.png",
     "documentation/FG_EXOME_resolved_flowchart.png"),
]

CATEGORIES = [
    # (label, src_rel, dst_rel, rsync_exclude_regex_or_None)
    ("QC vcf/tbi",          "QC",               "data/qc_vcf_full",       r".*\.report\.txt$"),
    ("QC reports",          "QC",               "documentation",          r".*\.vcf\.gz(\.tbi)?$"),
    ("renamed_vcf_chr",     "renamed_vcf_chr",   "data/renamed_vcf_chr",   None),
    ("plink_fg_merged_chr", "plink_fg_merged_chr","data/plink_fg_merged_chr", None),
]

SINGLE_FILES = [
    ("id_mapping.tsv",        "finngen_R14_exome_id_mapping.tsv",       "data/finngen_R14_exome_id_mapping.tsv"),
    ("ld raw",                "ld/finngen_R14_exome.ld.tsv.gz",         "data/finngen_R14_exome.ld.tsv.gz"),
    ("ld annotated",          "ld/finngen_R14_exome.ld_annotated.tsv.gz", "data/finngen_R14_exome.ld_annotated.tsv.gz"),
]


def exists(gs_path):
    return subprocess.run(["gsutil", "-q", "stat", gs_path]).returncode == 0


def run(cmd, dry_run):
    print(f"$ {' '.join(cmd)}", flush=True)
    if dry_run:
        return
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(f"ERROR: command failed ({result.returncode}): {' '.join(cmd)}")


def rsync(src, dst, exclude, dry_run):
    cmd = ["gsutil", "-m", "rsync", "-r", "-d"]
    if exclude:
        cmd += ["-x", exclude]
    if dry_run:
        cmd += ["-n"]
    cmd += [src, dst]
    run(cmd, dry_run=False)  # -n already makes gsutil itself a no-op preview


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("data_path", help="gs://.../data root")
    p.add_argument("release_path", help="gs://.../release root")
    p.add_argument("--dry-run", action="store_true", help="Preview with gsutil rsync -n; no changes made")
    p.add_argument("--skip-validate", action="store_true",
                   help="Skip running validate_samples_denials.py against release_path afterward")
    args = p.parse_args()

    data = args.data_path.rstrip("/")
    release = args.release_path.rstrip("/")

    for label, src_rel, dst_rel, exclude in CATEGORIES:
        print(f"\n== {label} ==", flush=True)
        rsync(f"{data}/{src_rel}", f"{release}/{dst_rel}", exclude, args.dry_run)

    for label, src_rel, dst_rel in SINGLE_FILES:
        print(f"\n== {label} ==", flush=True)
        src = f"{data}/{src_rel}"
        if not exists(src):
            print(f"WARNING: {src} not found — skipping (not yet generated?)", flush=True)
            continue
        run(["gsutil", "cp", src, f"{release}/{dst_rel}"], args.dry_run)

    for label, local_rel, dst_rel in LOCAL_FILES:
        print(f"\n== {label} ==", flush=True)
        local = REPO_ROOT / local_rel
        if not local.exists():
            sys.exit(f"ERROR: {local} not found")
        run(["gsutil", "cp", str(local), f"{release}/{dst_rel}"], args.dry_run)

    if args.skip_validate:
        return
    if args.dry_run:
        print("\n(dry run — skipping validate_samples_denials.py; nothing was actually copied)")
        return

    print(f"\n== validate_samples_denials.py {release} ==")
    validate_script = pathlib.Path(__file__).resolve().parent / "validate_samples_denials.py"
    result = subprocess.run([sys.executable, str(validate_script), release])
    if result.returncode != 0:
        sys.exit(f"ERROR: validate_samples_denials.py failed against {release}")


if __name__ == "__main__":
    main()
