#!/usr/bin/env python3
"""
Build a leads table from FinnGen group reports and annotate with exome LD.

Produces:
  PREFIX.leads.tsv.gz     — one row per credible set with cs_type annotation
  PREFIX.leads.pkl        — same, as pickle (cache)
  PREFIX.annotated.tsv.gz — LD file joined to leads on FG_SNP == locus_id

Usage:
    build_leads_table.py --credible_groups groups.txt --prefix r14 \\
        --ld_file exome_finngen_ld.tsv.gz

    # Test with first N files (default N=5):
    build_leads_table.py --credible_groups groups.txt --prefix r14 \\
        --ld_file exome_finngen_ld.tsv.gz --test
"""

import argparse
import gzip
import sys

import pandas as pd
from tqdm import tqdm

TEST_N = 5

READ_COLS = ["phenotype_abbreviation", "locus_id", "lead_mlogp", "lead_beta",
             "lead_af_alt", "good_cs", "best_coding_var", "functional_variants_relaxed"]
KEEP_COLS = ["PHENO", "locus_id", "lead_mlogp", "lead_beta", "lead_af_alt", "good_cs"]


def read_file_list(path: str) -> list[str]:
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def concat_to_df(paths: list[str]) -> pd.DataFrame:
    dfs = []
    for path in tqdm(paths, desc="Reading", unit="file", file=sys.stderr):
        try:
            dfs.append(pd.read_csv(path, sep="\t", usecols=READ_COLS, low_memory=False))
        except Exception as e:
            print(f"  WARNING: could not read {path}: {e}", file=sys.stderr)
    return pd.concat(dfs, ignore_index=True)


def _parse_functional_first(series: pd.Series) -> pd.DataFrame:
    """Extract variant name and r2 from the first entry of functional_variants_relaxed.

    Format: variant|consequence|gene|r2;...
    """
    def parse(val):
        if pd.isna(val):
            return None, None
        first = val.split(";")[0]
        parts = first.split("|")
        return parts[0], float(parts[3])

    parsed = series.map(parse)
    return pd.DataFrame(parsed.tolist(), index=series.index, columns=["functional_var", "functional_var_r2"])


def build_leads(df: pd.DataFrame) -> pd.DataFrame:
    df = df.rename(columns={"phenotype_abbreviation": "PHENO"})

    has_coding     = df["best_coding_var"].notna()
    has_functional = df["functional_variants_relaxed"].notna()

    df["cs_type"] = "NA"
    df.loc[has_functional & ~has_coding, "cs_type"] = "functional_relaxed"
    df.loc[has_coding,                   "cs_type"] = "coding"

    parsed = _parse_functional_first(df["functional_variants_relaxed"])
    df["functional_var"]    = parsed["functional_var"]
    df["functional_var_r2"] = parsed["functional_var_r2"]

    return df[KEEP_COLS + ["cs_type", "functional_var", "functional_var_r2"]].copy()


def build_annotated(ld_file: str, leads: pd.DataFrame, out_path: str, chunksize: int = 500_000) -> tuple[int, int]:
    """Read LD file in chunks, join each to leads, write directly to gzipped output."""
    n_rows = 0
    unique_exome: set = set()
    header_written = False

    print(f"Reading LD file {ld_file}", file=sys.stderr)
    with gzip.open(out_path, "wt") as out:
        for chunk in tqdm(pd.read_csv(ld_file, sep="\t", chunksize=chunksize),
                          desc="Joining LD", unit="chunk", file=sys.stderr):
            merged = chunk.merge(leads, left_on="FG_SNP", right_on="locus_id", how="inner").drop(columns="locus_id")
            if not merged.empty:
                merged.to_csv(out, sep="\t", index=False, header=not header_written, na_rep="NA")
                header_written = True
                n_rows += len(merged)
                unique_exome.update(merged["EXOME_SNP"])

    return n_rows, len(unique_exome)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--credible_groups", required=True, metavar="FILE",
                        help="Text file with one group_report path per line")
    parser.add_argument("--ld_file", required=True, metavar="FILE",
                        help="Exome LD TSV (gzipped); columns: FG_SNP, EXOME_SNP, R2, exome_consequence, EXOME_AF")
    parser.add_argument("--prefix", required=True, metavar="PREFIX",
                        help="Output prefix; generates PREFIX.leads.tsv.gz and PREFIX.annotated.tsv.gz")
    parser.add_argument(
        "--test", type=int, default=None, metavar="N",
        help=f"Process only the first N files (default when flag omitted: {TEST_N})",
        nargs="?", const=TEST_N,
    )
    args = parser.parse_args()

    leads_pkl      = f"{args.prefix}.leads.pkl"
    leads_path     = f"{args.prefix}.leads.tsv.gz"
    annotated_path = f"{args.prefix}.annotated.tsv.gz"

    # --- Step 1: leads table (cached) ---
    if os.path.exists(leads_pkl):
        print(f"Loading cached leads from {leads_pkl}", file=sys.stderr)
        leads = pd.read_pickle(leads_pkl)
    else:
        paths = read_file_list(args.credible_groups)
        if args.test is not None:
            paths = paths[:args.test]
            print(f"TEST MODE: {len(paths)} files", file=sys.stderr)
        else:
            print(f"{len(paths)} group_report files to process", file=sys.stderr)

        raw = concat_to_df(paths)
        leads = build_leads(raw)

        if leads.empty:
            print("No rows — exiting.", file=sys.stderr)
            return

        leads.to_pickle(leads_pkl)
        print(f"Leads pickle written to {leads_pkl}", file=sys.stderr)
        leads.to_csv(leads_path, sep="\t", index=False, compression="gzip", na_rep="NA")
        print(f"Leads TSV written to {leads_path}", file=sys.stderr)

    print(f"  {len(leads)} CSs across {leads['PHENO'].nunique()} phenotypes", file=sys.stderr)

    # --- Step 2: annotated LD table ---
    n_rows, n_exome = build_annotated(args.ld_file, leads, annotated_path)
    print(f"Annotated table written to {annotated_path}  ({n_rows} rows, {n_exome} unique exome variants)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
