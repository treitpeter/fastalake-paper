#!/usr/bin/env python3
"""Regenerate the reviewer index from the registered NUMBERS.tsv rows."""

import argparse
import collections
import csv
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
REF = "v42"


def escape(value):
    return str(value or "").replace("|", "\\|").replace("\n", " ").strip()


def paragraph_key(value):
    match = re.search(r"\d+", value or "")
    return int(match.group()) if match else 10**6


def render():
    records = []
    units = []
    for folder in ("figures", "tables", "text_statistics"):
        for path in sorted((REPO / folder).glob("*/NUMBERS.tsv")):
            with path.open(encoding="utf-8", newline="") as stream:
                reader = csv.DictReader(stream, delimiter="\t")
                required = {f"{REF}_{field}" for field in ("para", "quote", "printed_value", "status")}
                if not required.issubset(reader.fieldnames or []):
                    raise ValueError(f"{path.relative_to(REPO)}: missing {REF} reference columns")
                rows = [dict(row, unit=path.parent.name) for row in reader]
            units.append((path.parent.relative_to(REPO), rows))
            records.extend(rows)
    if not records:
        raise ValueError("No registered number rows found")
    counts = collections.Counter(row[f"{REF}_status"] for row in records)
    lines = [
        "# Reviewer guide", "",
        f"Find the code and source tables for {len(records):,} registered numerical claims below. "
        "The current manuscript reference is recorded in [MANUSCRIPT.json](MANUSCRIPT.json). "
        "See [change notes](../CHANGELOG.md) for the updates awaiting the next bioRxiv revision. "
        "Remaining gaps are listed in each unit's `OPEN_ITEMS.md`.", "",
        "Generated with `python tools/render_reviewer_guide.py`; edit the unit `NUMBERS.tsv` files, then regenerate. "
        "`python tools/render_reviewer_guide.py --check` detects stale documentation. No manuscript text is needed to build this index.", "",
        "Statuses: `MATCH` records agreement at printed precision; `OPEN` records an unresolved value or trace; "
        "`REMOVED_V30`, `REMOVED_V33`, `REMOVED_V36` and `REMOVED_V42` retain claims removed from those versions. "
        "Historical v28d statuses are shown separately, including `SUPERSEDED_V33`.", "",
        f"Totals: {len(records)} registered rows. {REF} status: "
        + ", ".join(f"{name} {count}" for name, count in sorted(counts.items())) + ".", "",
        "## Check a recorded claim", "",
        "Find the row, open its code and source table, then run `python paper.py build <unit>`. "
        f"To verify the printed value and its quote, run `python tools/check_numbers.py --unit <unit> --manuscript-text visible_{REF}.txt`. "
        "The manuscript text is supplied separately. Without it, only registered-value consistency is checked. "
        "The default check lists OPEN entries without failing; `--strict` also fails on them. "
        "Neither mode discovers manuscript numbers missing from the registered tables.", "",
        "## Per-unit index", "",
        f"| unit | folder | registered rows | {REF} MATCH | {REF} OPEN | removed |",
        "|---|---|---|---|---|---|",
    ]
    for path, rows in units:
        status = collections.Counter(row[f"{REF}_status"] for row in rows)
        removed = sum(count for name, count in status.items() if name.startswith("REMOVED_"))
        lines.append(f"| [{path.name}](../{path}/) | [{path}/](../{path}/) | {len(rows)} | {status['MATCH']} | {status['OPEN']} | {removed} |")
    lines += [
        "", "## Registered claims in reference manuscript order", "",
        f"| {REF} para | v28d line | claim quote | printed value | recomputed value | where computed | source_data file | {REF} status | v28d status | unit |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    records.sort(key=lambda row: (paragraph_key(row[f"{REF}_para"]), paragraph_key(row["v28d_line"]), row["unit"]))
    for row in records:
        values = [row[f"{REF}_para"], row["v28d_line"], row[f"{REF}_quote"], row[f"{REF}_printed_value"],
                  row[f"{REF}_recomputed_value"], row[f"{REF}_code"], row[f"{REF}_source_data"], row[f"{REF}_status"], row["status"], row["unit"]]
        lines.append("| " + " | ".join(escape(value) for value in values) + " |")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Fail if the committed guide is stale")
    args = parser.parse_args()
    path = REPO / "docs" / "REVIEWER_GUIDE.md"
    expected = render()
    if args.check:
        if not path.exists() or path.read_text(encoding="utf-8") != expected:
            print("STALE docs/REVIEWER_GUIDE.md; run python tools/render_reviewer_guide.py")
            return 1
        print("Reviewer guide matches all registered number rows")
    else:
        path.write_text(expected, encoding="utf-8")
        print("Wrote docs/REVIEWER_GUIDE.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
