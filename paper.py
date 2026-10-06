#!/usr/bin/env python3
"""Single entry point for the FastaLake paper repository.

Commands
--------
``python paper.py list``
    List every unit (figure, table, text-statistics section, pipeline) with the state of
    its required files, the number of NUMBERS.tsv rows by status and OPEN_ITEMS entries.
``python paper.py build <unit>``
    Execute the unit's notebook(s) headless from ``source_data/`` (e.g. ``Fig2``,
    ``SFig04``, ``STab1``, ``methods``).
``python paper.py build-all``
    Execute every notebook.
``python paper.py check``
    Run the release scan, the number check and the self-tests of both.

All paths are relative to the repository root. No command reads outside the repository
except controlled inputs named by ``FASTALAKE_CONTROLLED_DIR`` and fonts named by
``FASTALAKE_FONT_DIR``.
"""

import argparse
import csv
import subprocess
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parent

#: Directories that hold units, and the files each kind of unit must contain.
UNIT_KINDS = {
    "figures": ["README.md", "{u}.ipynb", "panels.py", "NUMBERS.tsv", "STATS.tsv",
                "ASSEMBLY.md", "PROVENANCE.md", "OPEN_ITEMS.md", "UPSTREAM.json"],
    "tables": ["README.md", "{u}.ipynb", "NUMBERS.tsv", "PROVENANCE.md", "OPEN_ITEMS.md",
               "UPSTREAM.json"],
    "text_statistics": ["README.md", "{u}.ipynb", "stats.py", "NUMBERS.tsv", "STATS.tsv",
                        "OPEN_ITEMS.md", "UPSTREAM.json"],
    "pipelines": ["README.md", "RECEIPTS.md"],
}


def units():
    """Yield ``(kind, unit_name, path)`` for every unit directory."""
    for kind in UNIT_KINDS:
        base = REPO / kind
        for d in sorted(p for p in base.iterdir() if p.is_dir()) if base.is_dir() else []:
            yield kind, d.name, d


def find_unit(name):
    """Return the directory of unit ``name`` (case-insensitive), or exit with a message."""
    for _, u, p in units():
        if u.lower() == name.lower():
            return p
    sys.exit(f"unknown unit {name!r}; run 'python paper.py list'")


def number_status(path):
    """Count NUMBERS.tsv rows by status class (MATCH, v28f_change, OPEN, other)."""
    f = path / "NUMBERS.tsv"
    if not f.is_file():
        return Counter()
    with open(f, newline="", encoding="utf-8") as fh:
        return Counter((r.get("status") or "").split(":")[0] or "blank"
                       for r in csv.DictReader(fh, delimiter="\t"))


def open_items(path):
    """Count list entries (lines starting with '- ') in OPEN_ITEMS.md."""
    f = path / "OPEN_ITEMS.md"
    if not f.is_file():
        return 0
    return sum(1 for ln in f.read_text(encoding="utf-8").splitlines() if ln.lstrip().startswith("- "))


def cmd_list(_args):
    """Print one line per unit: kind, name, missing required files, number status, open items."""
    print(f"{'kind':16} {'unit':28} {'missing':>7}  numbers                         open")
    for kind, u, p in units():
        req = [r.format(u=u) for r in UNIT_KINDS[kind]]
        missing = [r for r in req if not (p / r).exists()]
        ns = number_status(p)
        nstr = " ".join(f"{k}={v}" for k, v in sorted(ns.items())) or "-"
        print(f"{kind:16} {u:28} {len(missing):>7}  {nstr:31} {open_items(p)}")
    return 0


def _py(*args):
    return subprocess.call([sys.executable, *map(str, args)], cwd=REPO)


def cmd_build(args):
    """Execute the notebook(s) of one unit."""
    p = find_unit(args.unit)
    return _py("tools/run_notebooks.py", p.relative_to(REPO))


def cmd_build_all(_args):
    """Execute every notebook in the repository."""
    return _py("tools/run_notebooks.py")


def cmd_check(_args):
    """Run release scan, number check and both self-tests; non-zero if any fails."""
    rc = 0
    for step in (["tools/scan_release.py", "--self-test"], ["tools/scan_release.py"],
                 ["tools/check_numbers.py", "--self-test"], ["tools/check_numbers.py"]):
        print(f"== {' '.join(step)}", flush=True)
        rc |= _py(*step)
    return rc


def main(argv=None):
    """Command-line entry point."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list", help="list units and their state").set_defaults(f=cmd_list)
    b = sub.add_parser("build", help="run one unit's notebook(s)")
    b.add_argument("unit")
    b.set_defaults(f=cmd_build)
    sub.add_parser("build-all", help="run every notebook").set_defaults(f=cmd_build_all)
    sub.add_parser("check", help="release scan and number check").set_defaults(f=cmd_check)
    a = ap.parse_args(argv)
    return a.f(a)


if __name__ == "__main__":
    sys.exit(main())
