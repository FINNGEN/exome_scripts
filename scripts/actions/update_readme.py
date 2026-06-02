#!/usr/bin/env python3
"""update_readme.py — inject file contents into README.md between path markers.

Any block in README.md bounded by:
    <!-- BEGIN:relative/path/to/file -->
    ...
    <!-- END:relative/path/to/file -->
is replaced with the current contents of that file (path relative to repo root).

Usage:
    python scripts/update_readme.py
"""
import pathlib, re, sys

REPO   = pathlib.Path(__file__).resolve().parent.parent
README = REPO / "README.md"

MARKER = re.compile(
    r"<!-- BEGIN:(?P<path>[^>]+) -->\n(?:.*?\n)?<!-- END:(?P=path) -->",
    re.DOTALL,
)

def update(text):
    changed = 0
    def replace(m):
        nonlocal changed
        rel   = m.group("path").strip()
        src   = REPO / rel
        if not src.exists():
            print(f"  SKIP (not found): {rel}", file=sys.stderr)
            return m.group(0)
        content = src.read_text().strip()
        changed += 1
        print(f"  inject: {rel} ({len(content):,} chars)")
        return f"<!-- BEGIN:{rel} -->\n{content}\n<!-- END:{rel} -->"
    return MARKER.sub(replace, text), changed

text, n = update(README.read_text())
README.write_text(text)
print(f"Done — {n} block(s) updated in {README.relative_to(REPO)}")
