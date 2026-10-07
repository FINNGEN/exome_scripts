#!/usr/bin/env python3
"""localize.py — build the `gcloud storage cp` commands that copy a completed
Cromwell run's outputs into release/, auto-detecting which WDL produced it
from the metadata JSON's workflowName (mapping tables live in
localize_config.py — see that file for what goes where and why).

Usage:
  python3 scripts/localize/localize.py <workflow_id_or_metadata.json> <gs://bucket/.../release/> [--workflow=NAME]

The first argument can be a path to a metadata JSON file, or a bare Cromwell
workflow id — in which case it's looked up at
~/Dropbox/Projects/CromwellInteract/tmp/<id>.json.

--workflow=NAME overrides auto-detection and forces a specific mapping from
localize_config.WORKFLOWS (useful if metadata lacks workflowName, or to
double-check a mapping choice).

Prints a ready-to-review bash script of `gcloud storage cp` commands to
stdout; nothing is executed here.
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import common
from localize_config import WORKFLOWS


def main():
    override = None
    args = []
    for a in sys.argv[1:]:
        if a.startswith("--workflow="):
            override = a.split("=", 1)[1]
        else:
            args.append(a)

    if len(args) != 2:
        print(__doc__)
        sys.exit(1)

    id_or_path, release_base = args
    metadata_path = common.resolve_metadata_path(id_or_path)
    metadata = common.load_metadata(metadata_path)
    workflow_name = override or metadata.get("workflowName")

    if workflow_name not in WORKFLOWS:
        sys.exit(f"ERROR: unrecognized workflowName '{workflow_name}' in {metadata_path} "
                  f"— known workflows: {', '.join(WORKFLOWS)}")

    print(f"# workflow: {workflow_name}{' (forced)' if override else ' (detected)'} — {metadata_path}",
          file=sys.stderr)
    outputs = common.load_outputs(metadata_path)
    commands = common.build_commands(outputs, WORKFLOWS[workflow_name], release_base)

    print("#!/usr/bin/env bash")
    print("set -euo pipefail")
    print()
    for cmd in commands:
        print(cmd)
        print()


if __name__ == "__main__":
    main()
