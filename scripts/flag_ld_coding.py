#!/usr/bin/env python3
"""
Filter a raw plink2 .vcor file to FG→exome pairs and optionally flag coding status.

Inputs:
  --vcor        raw plink2 .vcor[.gz]             (#CHROM_A POS_A ID_A CHROM_B POS_B ID_B UNPHASED_R2)
  --fg_bim      FG plink BIM file                 (used to identify which IDs are FG variants)
  --annot       VEP annotation TSV[.gz] or .pkl   (optional)

Without --annot: outputs FG_SNP, EXOME_SNP, R2 for FG→exome pairs only.
With --annot:    also adds is_fg_coding and is_ex_coding columns.

--annot accepts:
  - a raw TSV (plain or .gz): read, cached as .pkl in the working directory
  - a .pkl file: loaded directly, skipping the TSV read

Variant IDs are normalized (strip 'chr', underscores → colons) to match the
annotation 'variant' column format (e.g. 21:5121014:C:A).
"""

import argparse
import csv
import gzip
import io
import os
import pickle
import sys
from contextlib import contextmanager
from pathlib import Path
import pandas as pd
from tqdm import tqdm


def open_maybe_gz(path):
    with open(path, "rb") as f:
        magic = f.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


@contextmanager
def open_with_progress(path: str):
    """Open a plain or gzipped file with a byte-level tqdm progress bar."""
    with open(path, "rb") as f:
        is_gz = f.read(2) == b"\x1f\x8b"
    file_size = os.path.getsize(path)
    pbar = tqdm(total=file_size, unit="B", unit_scale=True, unit_divisor=1024,
                desc=Path(path).name, dynamic_ncols=True)
    raw = open(path, "rb")

    class _Tracked:
        def read(self, size=-1):
            chunk = raw.read(size)
            pbar.update(len(chunk))
            return chunk

    try:
        if is_gz:
            yield io.TextIOWrapper(gzip.GzipFile(fileobj=_Tracked()))
        else:
            yield io.TextIOWrapper(raw)
    finally:
        raw.close()
        pbar.close()



CODING_CONSEQUENCES = {
    "missense_variant",
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "stop_lost",
    "inframe_insertion",
    "inframe_deletion",
}


def strip_suffixes(name: str) -> str:
    """Strip .gz then .vcor/.ld/.tsv from a filename stem. foo.vcor.gz → foo"""
    p = Path(name)
    if p.suffix == ".gz":
        p = p.with_suffix("")
    if p.suffix in (".vcor", ".ld", ".tsv"):
        p = p.with_suffix("")
    return str(p)


def normalize(variant_id: str) -> str:
    """chr21_5121014_C_A  →  21:5121014:C:A"""
    return variant_id.removeprefix("chr").replace("_", ":")


def flip(variant_id: str) -> str:
    """21:5121014:C:A  →  21:5121014:A:C (swap REF and ALT)"""
    parts = variant_id.split(":")
    if len(parts) == 4:
        parts[2], parts[3] = parts[3], parts[2]
    return ":".join(parts)


def load_fg_ids(bim_path: str) -> set:
    """Read FG variant IDs from BIM column 2 into a set."""
    fg_ids = set()
    with open(bim_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                fg_ids.add(parts[1])
    print(f"Loaded {len(fg_ids)} FG variant IDs from {bim_path}", file=sys.stderr)
    return fg_ids


def load_consequence(annot_path: str, id_col: str) -> dict:
    if Path(annot_path).suffix == ".pkl":
        print(f"Loading annotation from pickle {annot_path}", file=sys.stderr)
        with open(annot_path, "rb") as f:
            return pickle.load(f)

    pickle_path = Path(Path(annot_path).name).with_suffix(".pkl")
    if pickle_path.exists():
        print(f"Loading annotation from cache {pickle_path}", file=sys.stderr)
        with open(pickle_path, "rb") as f:
            return pickle.load(f)

    print(f"Reading annotation from {annot_path}", file=sys.stderr)
    annot = pd.read_csv(open_maybe_gz(annot_path), sep="\t", usecols=[id_col, "most_severe"])
    consequence = dict(zip(annot[id_col], annot["most_severe"]))

    with open(pickle_path, "wb") as f:
        pickle.dump(consequence, f)
    print(f"Cached annotation to {pickle_path}", file=sys.stderr)

    return consequence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcor",        required=True,  help="plink2 .vcor file, plain or .gz")
    parser.add_argument("--fg_bim",      required=True,  help="FG plink BIM file")
    parser.add_argument("--annot",       default=None,   help="VEP annotation TSV (plain/.gz) or .pkl — optional")
    parser.add_argument("--annot_id_col",default="rsid", help="Variant ID column in annotation (default: rsid)")
    parser.add_argument("--out",         default=None,   help="Output file (default: <workdir>/<vcor_stem>_fg_exome_ld.tsv)")
    parser.add_argument("--test",        nargs="?", const=1000, default=None, type=int,
                        help="Only process first N data lines (default N=1000 if flag given without value)")
    args = parser.parse_args()

    fg_ids      = load_fg_ids(args.fg_bim)
    consequence = load_consequence(args.annot, args.annot_id_col) if args.annot else None

    stem = strip_suffixes(Path(args.vcor).name)
    out_path = args.out or str(Path.cwd() / f"{stem}_fg_exome_ld.tsv")

    out_fieldnames = ["FG_SNP", "EXOME_SNP", "R2"]
    if consequence:
        out_fieldnames += ["is_fg_coding", "is_ex_coding"]

    missing_fg, missing_ex = set(), set()
    n_written = 0

    with open_with_progress(args.vcor) as fin, open(out_path, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = reader.fieldnames
        if not fieldnames:
            raise ValueError(f"vcor file has no header or is empty: {args.vcor}")

        writer = csv.DictWriter(fout, fieldnames=out_fieldnames, delimiter="\t")
        writer.writeheader()

        for i, row in enumerate(reader):
            if args.test is not None and i >= args.test:
                break

            fg_id  = row["ID_A"]
            ex_id  = row["ID_B"]

            # keep only FG→exome pairs (ID_B not a FG variant)
            if ex_id in fg_ids:
                continue

            out_row = {"FG_SNP": fg_id, "EXOME_SNP": ex_id, "R2": row["UNPHASED_R2"]}

            if consequence:
                fg_cons = consequence.get(fg_id)
                ex_cons = consequence.get(ex_id)

                if fg_cons is None:
                    missing_fg.add(fg_id)
                if ex_cons is None:
                    missing_ex.add(ex_id)

                out_row["is_fg_coding"] = fg_cons in CODING_CONSEQUENCES
                out_row["is_ex_coding"] = ex_cons in CODING_CONSEQUENCES

            writer.writerow(out_row)
            n_written += 1

    print(f"Written {n_written} FG-exome pairs to {out_path}", file=sys.stderr)

    if missing_fg or missing_ex:
        err_path = out_path.replace(".tsv", "_missing.txt")
        with open(err_path, "w") as ferr:
            for v in sorted(missing_fg):
                ferr.write(f"FG\t{v}\n")
            for v in sorted(missing_ex):
                ferr.write(f"EXOME\t{v}\n")
        print(f"WARNING: {len(missing_fg)} FG and {len(missing_ex)} exome variants not found in annotation — see {err_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
