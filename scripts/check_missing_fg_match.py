#!/usr/bin/env python3
"""check_missing_fg_match.py — for exome samples with no KING genetic match
(STATUS MISSING in resolve_mapping.py's id_mapping.tsv), check whether the
exome ID nonetheless corresponds to a real FinnGen participant — directly or
via a known alias — and report the matched FG ID and cohort. IDs absent from
the COV file are also checked against a denials list (samples removed from
analysis upstream); a denials match reports cohort DENIAL.

Output columns: ID  FG_ID  cohort  DATASET  (FG_ID/cohort are NA if no match is found)
Output contains real sample IDs — never commit it to the repo.

Usage:
  python check_missing_fg_match.py id_mapping.tsv --cov R14_COV_V0.FID.txt.gz \\
      [--aliases FILE] [--denials FILE] [--out FILE]
  python check_missing_fg_match.py --test
"""
import argparse, os, pathlib, sys
import pandas as pd

from resolve_mapping import load_aliases, DEFAULT_ALIASES

FG_ID_LEN = 10  # canonical FinnGen ID = "FG" + 8 characters

_SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
TEST_DIR     = _SCRIPT_DIR.parent / "test"
TEST_MAPPING = TEST_DIR / "missing_id_mapping.tsv"
TEST_COV     = TEST_DIR / "missing_cov.txt.gz"
TEST_ALIASES = TEST_DIR / "missing_aliases.tsv"
TEST_DENIALS = TEST_DIR / "missing_denials.txt"


def load_cov(path):
    df = pd.read_csv(path, sep="\t", compression="infer", usecols=["FID", "cohort"], dtype=str)
    return dict(zip(df["FID"], df["cohort"]))


def load_denials(path):
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def resolve_one(exome_id, fid_cohort, ag, denials):
    """Return (fg_id, cohort, kind) with kind in {DIRECT, ALIAS, DENIAL, NO}."""
    direct = exome_id[:FG_ID_LEN]
    if direct in fid_cohort:
        return direct, fid_cohort[direct], "DIRECT"

    group = ag.get(exome_id) or ag.get(direct) or frozenset()
    hits = {m[:FG_ID_LEN] for m in group} & fid_cohort.keys()
    if len(hits) == 1:
        fg_id = next(iter(hits))
        return fg_id, fid_cohort[fg_id], "ALIAS"
    if len(hits) > 1:
        print(f"WARNING: {exome_id} has {len(hits)} alias candidates in COV file "
              f"(expected at most 1) — treating as unresolved", file=sys.stderr)

    if direct in denials:
        return direct, "DENIAL", "DENIAL"
    den_hits = {m[:FG_ID_LEN] for m in group} & denials
    if len(den_hits) == 1:
        return next(iter(den_hits)), "DENIAL", "DENIAL"
    if len(den_hits) > 1:
        print(f"WARNING: {exome_id} has {len(den_hits)} alias candidates in denials list "
              f"(expected at most 1) — treating as unresolved", file=sys.stderr)

    return "NA", "NA", "NO"


def main():
    p = argparse.ArgumentParser(
        description="Check whether MISSING exome samples plausibly match a real FinnGen ID.")
    p.add_argument("mapping", nargs="?",
                   help="id_mapping.tsv from resolve_mapping.py (required unless --test)")
    p.add_argument("--cov", help="R14_COV_*.FID.txt.gz file with FID + cohort columns "
                                 "(required unless --test)")
    p.add_argument("--aliases", default=DEFAULT_ALIASES,
                   help="Tab-delimited alias file (one group per line)")
    p.add_argument("--denials", default=None,
                   help="Plain list (one ID per line) of samples removed from analysis upstream")
    p.add_argument("--out", default=None, help="Output TSV path (default: <mapping-stem>_fg_match.tsv)")
    p.add_argument("--test", action="store_true", help="Run on built-in synthetic test dataset")
    args = p.parse_args()

    if args.test:
        for f in (TEST_MAPPING, TEST_COV, TEST_ALIASES, TEST_DENIALS):
            if not f.exists():
                sys.exit(f"ERROR: test file not found: {f}")
        mapping_path, cov_path, aliases_path, denials_path = TEST_MAPPING, TEST_COV, TEST_ALIASES, TEST_DENIALS
        out_path = TEST_DIR / "missing_fg_match_output.tsv"
    else:
        if not args.mapping or not args.cov:
            p.error("mapping and --cov are required unless --test is given")
        mapping_path, cov_path, aliases_path, denials_path = args.mapping, args.cov, args.aliases, args.denials
        stem = os.path.splitext(args.mapping)[0]
        out_path = args.out or f"{stem}_fg_match.tsv"

    df = pd.read_csv(mapping_path, sep="\t", dtype=str)
    if miss := {"QUERY", "STATUS", "DATASET"} - set(df.columns):
        sys.exit(f"ERROR: missing columns in mapping file: {miss}")
    missing = df[df["STATUS"].str.startswith("MISSING")]

    fid_cohort = load_cov(cov_path)
    ag = load_aliases(str(aliases_path)) if aliases_path else {}
    denials = load_denials(denials_path) if denials_path else set()

    rows = []
    kinds = {"DIRECT": 0, "ALIAS": 0, "DENIAL": 0, "NO": 0}
    kind_by_dataset = {}
    for _, r in missing.iterrows():
        exome_id, dataset = r["QUERY"], r["DATASET"]
        fg_id, cohort, kind = resolve_one(exome_id, fid_cohort, ag, denials)
        rows.append({"ID": exome_id, "FG_ID": fg_id, "cohort": cohort, "DATASET": dataset})
        kinds[kind] += 1
        kind_by_dataset.setdefault(dataset, {"DIRECT": 0, "ALIAS": 0, "DENIAL": 0, "NO": 0})[kind] += 1

    out = pd.DataFrame(rows, columns=["ID", "FG_ID", "cohort", "DATASET"])
    out.to_csv(out_path, sep="\t", index=False)

    n = len(out)
    print(f"\n{mapping_path} -> {len(df)} rows, {n} MISSING\n")
    print(f"  {'DIRECT':<8} {kinds['DIRECT']:>6}")
    print(f"  {'ALIAS':<8} {kinds['ALIAS']:>6}")
    print(f"  {'DENIAL':<8} {kinds['DENIAL']:>6}")
    print(f"  {'NO':<8} {kinds['NO']:>6}")
    print("\nBy origin dataset:")
    for ds, c in sorted(kind_by_dataset.items()):
        print(f"  {ds:<10} DIRECT={c['DIRECT']:>5}  ALIAS={c['ALIAS']:>5}  "
              f"DENIAL={c['DENIAL']:>5}  NO={c['NO']:>5}")
    print(f"\nWritten -> {out_path}  (contains real IDs — do not commit)")


if __name__ == "__main__":
    main()
