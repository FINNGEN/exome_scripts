#!/usr/bin/env python3
"""
Build a leads table of coding-unexplained credible sets from FinnGen group reports.

Produces two output files from a given PREFIX:
  PREFIX.leads.tsv.gz   — one row per credible set, all coding-unexplained CSs,
                          columns: phenotype, locus_id, lead_mlogp, lead_beta, good_cs
                          (cached: if file exists it is read directly on next run)
  PREFIX.hits.tsv.gz    — one row per unique locus_id, with a hits JSON column
                          {phenotype: [mlogp, beta], ...} after applying --mlogp
                          and --good-cs filters
  PREFIX.annotated.tsv.gz — LD file filtered to non-coding FG / coding exome pairs,
                          joined to the hits table on FG_SNP == locus_id;
                          columns: EXOME_SNP, FG_SNP, R2, hits

Usage:
    build_leads_table.py --credible_groups FILE_LIST --prefix PREFIX [OPTIONS]

    build_leads_table.py --credible_groups groups.txt --prefix r14

    # With filters for the hits table and LD annotation:
    build_leads_table.py --credible_groups groups.txt --prefix r14 \\
        --ld_file finngen_R14_exome_r0.6_ld.tsv.gz --mlogp 6.0 --good-cs

    # Test with first N files (default N=5):
    build_leads_table.py --credible_groups groups.txt --prefix r14 \\
        --ld_file finngen_R14_exome_r0.6_ld.tsv.gz --test
"""

import argparse
import json
import os
import sys
import tempfile

import pandas as pd
from tqdm import tqdm

TEST_N = 5

KEEP_COLS = ["phenotype_abbreviation", "locus_id", "lead_mlogp", "lead_beta", "good_cs"]


def read_file_list(path: str) -> list[str]:
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def concat_to_tmp(paths: list[str]) -> str:
    """Concatenate all group_report files into one temp TSV, keeping a single header."""
    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", dir="/tmp", delete=False, prefix="group_reports_"
    )
    header_written = False
    for path in tqdm(paths, desc="Concatenating", unit="file", file=sys.stderr):
        try:
            with open(path) as f:
                header = f.readline()
                if not header_written:
                    tmp.write(header)
                    header_written = True
                tmp.write(f.read())
        except Exception as e:
            print(f"  WARNING: could not read {path}: {e}", file=sys.stderr)
    tmp.close()
    return tmp.name


def build_leads(tmp_path: str) -> pd.DataFrame:
    print(f"Reading concatenated file {tmp_path}", file=sys.stderr)
    df = pd.read_csv(tmp_path, sep="\t", low_memory=False)
    mask = df["best_coding_var"].isna() | (df["best_coding_var"].astype(str).str.strip() == "NA")
    return df[mask][KEEP_COLS].copy()


def load_ld(ld_file: str, ld_cache_path: str) -> pd.DataFrame:
    if os.path.exists(ld_cache_path):
        print(f"Loading cached LD table from {ld_cache_path}", file=sys.stderr)
        return pd.read_csv(ld_cache_path, sep="\t")
    print(f"Reading LD file {ld_file}", file=sys.stderr)
    ld = pd.read_csv(ld_file, sep="\t")
    ld = ld[(ld["is_fg_coding"] == False) & (ld["is_ex_coding"] == True)]
    print(f"  {len(ld)} non-coding FG / coding exome LD pairs", file=sys.stderr)
    ld.to_csv(ld_cache_path, sep="\t", index=False, compression="gzip")
    print(f"LD cache written to {ld_cache_path}", file=sys.stderr)
    return ld


def build_annotated(ld: pd.DataFrame, hits: pd.DataFrame) -> pd.DataFrame:
    merged = ld.merge(hits, left_on="FG_SNP", right_on="locus_id", how="inner")
    return merged[["EXOME_SNP", "FG_SNP", "R2", "hits"]]


def build_hits(leads: pd.DataFrame) -> pd.DataFrame:
    hits: dict[str, dict] = {}
    for _, row in leads.iterrows():
        locus = row["locus_id"]
        if locus not in hits:
            hits[locus] = {}
        hits[locus][row["phenotype_abbreviation"]] = [round(row["lead_mlogp"], 4), round(row["lead_beta"], 4)]

    return pd.DataFrame([
        {"locus_id": locus, "hits": json.dumps(d)}
        for locus, d in hits.items()
    ])


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--credible_groups", required=True, metavar="FILE",
                        help="Text file with one group_report path per line")
    parser.add_argument("--ld_file", required=True, metavar="FILE",
                        help="Exome LD TSV (gzipped); columns: FG_SNP, EXOME_SNP, R2, is_fg_coding, is_ex_coding")
    parser.add_argument("--prefix", required=True, metavar="PREFIX",
                        help="Output prefix; generates PREFIX.leads.tsv.gz, PREFIX.hits.tsv.gz, PREFIX.annotated.tsv.gz")
    parser.add_argument(
        "--mlogp", type=float, default=None, metavar="THRESH",
        help="Minimum lead_mlogp for the hits table (default: no filter)",
    )
    parser.add_argument(
        "--good-cs", action="store_true",
        help="Restrict hits table to good_cs == True",
    )
    parser.add_argument(
        "--test", type=int, default=None, metavar="N",
        help=f"Process only the first N files (default when flag omitted: {TEST_N})",
        nargs="?", const=TEST_N,
    )
    args = parser.parse_args()

    leads_path = f"{args.prefix}.leads.tsv.gz"
    hits_path = f"{args.prefix}.hits.tsv.gz"
    ld_cache_path = f"{args.prefix}.ld_noncoding.tsv.gz"
    annotated_path = f"{args.prefix}.annotated.tsv.gz"

    # --- Step 1: leads table (cached) ---
    if os.path.exists(leads_path):
        print(f"Loading cached leads table from {leads_path}", file=sys.stderr)
        leads = pd.read_csv(leads_path, sep="\t")
    else:
        paths = read_file_list(args.credible_groups)
        if args.test is not None:
            paths = paths[:args.test]
            print(f"TEST MODE: {len(paths)} files", file=sys.stderr)
        else:
            print(f"{len(paths)} group_report files to process", file=sys.stderr)

        tmp_path = concat_to_tmp(paths)
        try:
            leads = build_leads(tmp_path)
        finally:
            os.unlink(tmp_path)

        if leads.empty:
            print("No rows passed filters — exiting.", file=sys.stderr)
            return

        leads.to_csv(leads_path, sep="\t", index=False, compression="gzip")
        print(f"Leads table written to {leads_path}", file=sys.stderr)

    print(
        f"  {len(leads)} coding-unexplained CSs across {leads['phenotype_abbreviation'].nunique()} phenotypes",
        file=sys.stderr,
    )

    # --- Step 2: hits table (cached, filtered) ---
    if os.path.exists(hits_path):
        print(f"Loading cached hits table from {hits_path}", file=sys.stderr)
        hits = pd.read_csv(hits_path, sep="\t")
    else:
        filtered = leads.copy()
        if args.good_cs:
            filtered = filtered[filtered["good_cs"] == True]
            print(f"  After good_cs filter: {len(filtered)} rows", file=sys.stderr)
        if args.mlogp is not None:
            filtered = filtered[filtered["lead_mlogp"] >= args.mlogp]
            print(f"  After mlogp >= {args.mlogp} filter: {len(filtered)} rows", file=sys.stderr)

        hits = build_hits(filtered)
        hits.to_csv(hits_path, sep="\t", index=False, compression="gzip")
        print(f"Hits table written to {hits_path}  ({len(hits)} unique loci)", file=sys.stderr)

    # --- Step 3: annotated LD table ---
    ld = load_ld(args.ld_file, ld_cache_path)
    annotated = build_annotated(ld, hits)
    annotated.to_csv(annotated_path, sep="\t", index=False, compression="gzip")
    print(
        f"Annotated table written to {annotated_path}  "
        f"({len(annotated)} rows, {annotated['EXOME_SNP'].nunique()} unique exome variants)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
