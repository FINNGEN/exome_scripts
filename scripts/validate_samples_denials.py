#!/usr/bin/env python3
"""validate_samples_denials.py — scan every VCF under a GCS path and confirm
none of their samples are on the denial list, directly or via a known alias.

Denial expansion mirrors the ExpandDenials task used by wdl/single_file_qc.wdl,
wdl/daly_qc.wdl and wdl/wes_chrom.wdl: the denials set is unioned with every
member of any alias group that intersects it, since the ID recorded in the
denial list and the ID actually present in a given VCF's header may be
different aliases of the same participant.

Sample lists are read remotely with `bcftools query -l` (no download) using
GCS_OAUTH_TOKEN for auth, the same approach the WDL tasks use.

Usage:
  python validate_samples_denials.py gs://bucket/path [--denials FILE] [--aliases FILE] [--workers N]
  python validate_samples_denials.py --test

Never commit this script's output — it lists real sample IDs.
"""
import argparse, os, pathlib, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed

from resolve_mapping import _parse_alias_lines

DEFAULT_DENIALS = "gs://thl-incoming-data/from_THL_registerteam/denial_lists/finngen_R14_finngenid_exclusion_list.txt"
DEFAULT_ALIASES = "gs://fg-3/exome_v2/inputs/finngen_R14_duplicate_list.txt"


def _localize(path, tmp_dir):
    """Copy a gs:// path to tmp_dir and return the local path; local paths pass through."""
    if not str(path).startswith("gs://"):
        return pathlib.Path(path)
    dest = tmp_dir / pathlib.Path(path).name
    result = subprocess.run(["gsutil", "cp", path, str(dest)], capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"ERROR: gsutil cp failed for {path}\n{result.stderr}")
    return dest


def load_expanded_denials(denials_path, aliases_path):
    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)
        local_denials = _localize(denials_path, tmp)
        local_aliases = _localize(aliases_path, tmp) if aliases_path else None

        denied = {line.strip() for line in local_denials.read_text().splitlines() if line.strip()}
        ag = _parse_alias_lines(local_aliases.read_text().splitlines(), label=str(aliases_path)) \
            if local_aliases else {}

    expanded = set(denied)
    for id_ in denied:
        group = ag.get(id_)
        if group:
            expanded.update(group)

    print(f"Loaded {len(denied)} denied IDs, expanded to {len(expanded)} incl. aliases", file=sys.stderr)
    return denied, expanded


def list_vcfs(gs_path):
    target = gs_path.rstrip("/") + "/**/*.vcf.gz"
    result = subprocess.run(["gsutil", "ls", target], capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"ERROR: gsutil ls failed for {target}\n{result.stderr}")
    return sorted(
        line.strip() for line in result.stdout.splitlines()
        if line.strip().startswith("gs://")
    )


def get_gcs_token():
    result = subprocess.run(
        ["gcloud", "auth", "application-default", "print-access-token"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        sys.exit(f"ERROR: could not obtain a GCS access token\n{result.stderr}")
    return result.stdout.strip()


def query_samples(vcf_path, env):
    result = subprocess.run(["bcftools", "query", "-l", vcf_path],
                             capture_output=True, text=True, env=env)
    if result.returncode != 0:
        return vcf_path, None, result.stderr.strip()
    return vcf_path, result.stdout.splitlines(), None


def classify(sample_id, denied, expanded):
    if sample_id in denied:
        return "DIRECT"
    if sample_id in expanded:
        return "ALIAS"
    return None


def main():
    p = argparse.ArgumentParser(
        description="Validate that no VCF under a GCS path contains a denied sample (directly or via alias).")
    p.add_argument("gs_path", nargs="?", help="GCS path to scan recursively for *.vcf.gz files")
    p.add_argument("--denials", default=DEFAULT_DENIALS,
                   help="Plain list (one ID per line), local or gs://; one sample ID per line")
    p.add_argument("--aliases", default=DEFAULT_ALIASES,
                   help="Tab-delimited alias groups (one group per line), local or gs://")
    p.add_argument("--workers", type=int, default=8, help="Parallel bcftools query -l calls")
    p.add_argument("--test", action="store_true", help="Run built-in synthetic self-test")
    args = p.parse_args()

    if args.test:
        sys.exit(run_test())

    if not args.gs_path:
        p.error("gs_path is required unless --test is given")

    denied, expanded = load_expanded_denials(args.denials, args.aliases)

    vcfs = list_vcfs(args.gs_path)
    if not vcfs:
        sys.exit(f"ERROR: no *.vcf.gz files found under {args.gs_path}")
    print(f"Found {len(vcfs)} VCF(s) under {args.gs_path}", file=sys.stderr)

    env = os.environ.copy()
    env["GCS_OAUTH_TOKEN"] = get_gcs_token()

    violations, errors = {}, {}
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futures = {ex.submit(query_samples, vcf, env): vcf for vcf in vcfs}
        for fut in as_completed(futures):
            vcf, samples, err = fut.result()
            if err:
                errors[vcf] = err
                continue
            hits = [(s, classify(s, denied, expanded)) for s in samples]
            hits = [(s, k) for s, k in hits if k]
            if hits:
                violations[vcf] = hits

    for vcf in vcfs:
        if vcf in errors:
            print(f"ERROR  {vcf}: {errors[vcf]}")
        elif vcf in violations:
            n_direct = sum(1 for _, k in violations[vcf] if k == "DIRECT")
            n_alias  = sum(1 for _, k in violations[vcf] if k == "ALIAS")
            print(f"FAIL   {vcf}: {len(violations[vcf])} denied sample(s) "
                  f"(direct={n_direct}, alias={n_alias})")
            for sample_id, kind in violations[vcf]:
                print(f"         {kind:<6} {sample_id}")
        else:
            print(f"OK     {vcf}")

    if errors:
        sys.exit(f"\n{len(errors)} VCF(s) failed to query — see ERROR lines above")
    if violations:
        n_hits = sum(len(v) for v in violations.values())
        sys.exit(f"\nFAILED: {n_hits} denied sample(s) found across {len(violations)} VCF(s)")

    print(f"\nPASSED: no denied samples found in {len(vcfs)} VCF(s)")
    return 0


def run_test():
    """Self-contained test: builds a synthetic local VCF + denials/aliases and
    checks that direct and alias-mediated denial hits are both caught, and
    clean samples are not flagged."""
    samples = ["SAMPLE_OK1", "SAMPLE_OK2", "SAMPLE_DENIED_DIRECT", "SAMPLE_ALIAS_OF_DENIED"]

    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)

        vcf = tmp / "test.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"
            "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t" + "\t".join(["0/1"] * len(samples)) + "\n"
        )
        vcf_gz = tmp / "test.vcf.gz"
        subprocess.run(["bcftools", "view", "-Oz", "-o", str(vcf_gz), str(vcf)], check=True)
        subprocess.run(["bcftools", "index", "-t", str(vcf_gz)], check=True)

        denials_path = tmp / "denials.txt"
        denials_path.write_text("SAMPLE_DENIED_DIRECT\nREF_DENIED_NOT_IN_VCF\n")

        aliases_path = tmp / "aliases.tsv"
        aliases_path.write_text("REF_DENIED_NOT_IN_VCF\tSAMPLE_ALIAS_OF_DENIED\n")

        denied, expanded = load_expanded_denials(str(denials_path), str(aliases_path))
        _, samples_out, err = query_samples(str(vcf_gz), os.environ.copy())
        if err:
            print(f"TEST FAILED: bcftools query -l errored: {err}")
            return 1

        hits = {s: classify(s, denied, expanded) for s in samples_out}
        expected = {
            "SAMPLE_OK1": None,
            "SAMPLE_OK2": None,
            "SAMPLE_DENIED_DIRECT": "DIRECT",
            "SAMPLE_ALIAS_OF_DENIED": "ALIAS",
        }

        ok = True
        for sample_id, exp_kind in expected.items():
            got = hits.get(sample_id)
            mark = "✓" if got == exp_kind else "✗"
            if got != exp_kind:
                ok = False
            print(f"  {mark}  {sample_id:<24} expected={exp_kind!r:<8} got={got!r}")

        print("\nTEST PASSED" if ok else "\nTEST FAILED")
        return 0 if ok else 1


if __name__ == "__main__":
    main()
