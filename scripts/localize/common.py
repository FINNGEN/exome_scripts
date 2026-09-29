"""common.py — shared helpers for scripts/localize/*.py.

Each localize/<wdl>.py script maps that WDL's Cromwell outputs directly onto
the release/ file structure documented in release/finngen_R14_exome_readme
(the public release readme — "## File structure" section is the ground
truth for every path used here). Usage is the same for every script:

  python3 scripts/localize/<wdl>.py <metadata.json> <gs://bucket/.../release/>

Prints a ready-to-review bash script of `gcloud storage cp` commands to
stdout; nothing is executed here.
"""
import json
import sys


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


def load_outputs(metadata_path):
    with open(metadata_path) as f:
        metadata = json.load(f)
    if metadata.get("status") != "Succeeded":
        print(f"WARNING: workflow status is '{metadata.get('status')}', not Succeeded — "
              f"outputs may be incomplete or missing", file=sys.stderr)
    outputs = metadata.get("outputs", {})
    if not outputs:
        sys.exit("ERROR: no 'outputs' key found in metadata JSON "
                  "(is this a completed top-level workflow metadata file?)")
    return outputs


def build_commands(outputs, mapping, data_base):
    """mapping: list of (output_suffix, dest_subpath, fixed_filename_or_None, nested)
    fixed_filename_or_None: for a scalar output renamed to a fixed name on copy;
    None means an array output where each file's own basename is kept.
    nested: True if the output is Array[Array[File]] (e.g. exome_ld's merged_plink)."""
    if not data_base.endswith("/"):
        data_base += "/"

    commands = []
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
        else:
            if nested:
                value = flatten(value)
            elif not isinstance(value, list):
                value = [value]
            value = [v for v in value if v]  # drop nulls (e.g. optional array entries)
            if not value:
                continue
            srcs = " ".join(value)
            commands.append(f"gcloud storage cp {srcs} {dest_dir}")

    return commands


def main(mapping, doc):
    if len(sys.argv) != 3:
        print(doc)
        sys.exit(1)
    metadata_path, data_base = sys.argv[1], sys.argv[2]
    outputs = load_outputs(metadata_path)
    commands = build_commands(outputs, mapping, data_base)

    print("#!/usr/bin/env bash")
    print("set -euo pipefail")
    print()
    for cmd in commands:
        print(cmd)
        print()
