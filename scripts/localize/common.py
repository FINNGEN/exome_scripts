"""common.py — shared helpers for scripts/localize/*.py.

Each localize/<wdl>.py script maps that WDL's Cromwell outputs directly onto
the release/ file structure documented in release/finngen_R14_exome_readme
(the public release readme — "## File structure" section is the ground
truth for every path used here). Usage is the same for every script:

  python3 scripts/localize/<wdl>.py <workflow_id_or_metadata.json> <gs://bucket/.../release/>

The first argument can be a path to a metadata JSON file, or a bare Cromwell
workflow id — in which case it's looked up at
~/Dropbox/Projects/CromwellInteract/tmp/<id>.json (see resolve_metadata_path).

Prints a ready-to-review bash script of `gcloud storage cp` commands to
stdout; nothing is executed here.
"""
import json
import pathlib
import sys

METADATA_DIR = pathlib.Path("~/Dropbox/Projects/CromwellInteract/tmp").expanduser()


def resolve_metadata_path(arg):
    """Accept either a real path to a metadata JSON file, or a bare Cromwell
    workflow id, which is looked up as METADATA_DIR/<id>.json."""
    given = pathlib.Path(arg).expanduser()
    if given.exists():
        return given

    candidate = METADATA_DIR / f"{arg}.json"
    if candidate.exists():
        return candidate

    sys.exit(f"ERROR: no metadata file found for '{arg}' "
              f"(tried {given} and {candidate})")


def find_output(outputs, suffix):
    # keys are workflow-name-prefixed (e.g. "exome_reassign_ids.chrom_vcfs");
    # match by suffix so this doesn't break if the workflow is renamed/aliased.
    for key, value in outputs.items():
        if key.endswith("." + suffix):
            return value
    return None


def flatten(value):
    """Flatten one level of nesting (Array[Array[File]] -> Array[File])."""
    out = []
    for v in value:
        if isinstance(v, list):
            out.extend(x for x in v if x)
        elif v:
            out.append(v)
    return out


def load_metadata(metadata_path):
    with open(metadata_path) as f:
        return json.load(f)


def load_outputs(metadata_path):
    metadata = load_metadata(metadata_path)
    if metadata.get("status") != "Succeeded":
        print(f"WARNING: workflow status is '{metadata.get('status')}', not Succeeded — "
              f"outputs may be incomplete or missing", file=sys.stderr)
    outputs = metadata.get("outputs", {})
    if not outputs:
        sys.exit("ERROR: no 'outputs' key found in metadata JSON "
                  "(is this a completed top-level workflow metadata file?)")
    return outputs


def _sanity_check(dest_counts, dest_basenames):
    """Flags two classes of real mistake before any command is printed:
    1. different output arrays landing in the same destination dir with
       different lengths (a shard missing from one output — e.g. hwe_bims
       has 96 entries but hwe_fams only has 95);
    2. two different outputs (or two entries of the same output) resolving to
       the identical destination basename — gcloud storage cp would silently
       let one overwrite the other with no warning."""
    ok = True
    for dest_dir, counts in dest_counts.items():
        distinct = set(counts.values())
        if len(distinct) > 1:
            ok = False
            detail = ", ".join(f"{suffix}={n}" for suffix, n in counts.items())
            print(f"# SANITY WARNING: file counts differ in {dest_dir}: {detail} "
                  f"— a shard may be missing from one output", file=sys.stderr)

    for dest_dir, basenames in dest_basenames.items():
        for base, suffixes in basenames.items():
            if len(suffixes) > 1:
                ok = False
                print(f"# SANITY WARNING: {base} would be written to {dest_dir} by "
                      f"{len(suffixes)} source(s) ({', '.join(suffixes)}) — "
                      f"one will silently overwrite the other", file=sys.stderr)

    if ok:
        print("# sanity checks passed: no count mismatches or destination collisions", file=sys.stderr)


def build_commands(outputs, mapping, data_base):
    """mapping: list of (output_suffix, dest_subpath, fixed_filename_or_None, nested)
    fixed_filename_or_None: for a scalar output renamed to a fixed name on copy;
    None means an array output where each file's own basename is kept.
    nested: True if the output is Array[Array[File]] (e.g. exome_ld's merged_plink)."""
    if not data_base.endswith("/"):
        data_base += "/"

    commands = []
    dest_counts = {}     # dest_dir -> {suffix: file count}
    dest_basenames = {}  # dest_dir -> {basename: [suffixes writing it]}

    for suffix, subpath, fixed_name, nested in mapping:
        value = find_output(outputs, suffix)
        if value is None:
            print(f"# skip {suffix}: not present in outputs", file=sys.stderr)
            continue

        dest_dir = data_base + (subpath.rstrip("/") + "/" if subpath else "")

        if fixed_name is not None:
            if isinstance(value, list):
                print(f"# ERROR: {suffix} expected a scalar output but got an array — skipping", file=sys.stderr)
                continue
            commands.append(f"gcloud storage cp {value} {dest_dir}{fixed_name}")
            dest_basenames.setdefault(dest_dir, {}).setdefault(fixed_name, []).append(suffix)
        else:
            # count shards, not files: a nested output has N files per shard
            n_shards = len([s for s in value if s]) if isinstance(value, list) else 1
            if nested:
                value = flatten(value)
            elif not isinstance(value, list):
                value = [value]
            value = [v for v in value if v]  # drop nulls (e.g. optional array entries)
            if not value:
                print(f"# WARNING: {suffix} present but empty — nothing to copy", file=sys.stderr)
                continue

            dest_counts.setdefault(dest_dir, {})[suffix] = n_shards
            for v in value:
                base = v.rsplit("/", 1)[-1]
                dest_basenames.setdefault(dest_dir, {}).setdefault(base, []).append(suffix)

            srcs = " ".join(value)
            commands.append(f"gcloud storage cp {srcs} {dest_dir}")

    _sanity_check(dest_counts, dest_basenames)
    return commands


def main(mapping, doc):
    if len(sys.argv) != 3:
        print(doc)
        sys.exit(1)
    id_or_path, data_base = sys.argv[1], sys.argv[2]
    metadata_path = resolve_metadata_path(id_or_path)
    outputs = load_outputs(metadata_path)
    commands = build_commands(outputs, mapping, data_base)

    print("#!/usr/bin/env bash")
    print("set -euo pipefail")
    print()
    for cmd in commands:
        print(cmd)
        print()
