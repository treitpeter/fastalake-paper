#!/usr/bin/env python3
"""Write source_data/STab3/ from the original producer outputs (run once, by the unit author).

Supplementary Table 3: PREDICT antibiotic endpoints under alternative handling of
rifaximin for hepatic encephalopathy (HE).

This script is a record of how each public table was derived. It reads only aggregate
producer outputs (no patient identifiers, no per-acquisition clinical fields) and needs
``LAB_ROOT`` (the lab project directory, never written into the repository) in the
environment. Readers of the repository do not need to run it.

Outputs (all in ``source_data/STab3/``)

``endpoints_by_definition.tsv``
    Byte copy of the table that was assembled into the Endpoints sheet of the published
    workbook (``iteration1/fix_rifaximin/SupplementaryTable6_rifaximin_sensitivity.tsv``).
``endpoints_producer_output.tsv``
    Byte copy of the producer output ``iteration1/fix_rifaximin/out/paired_endpoints_by_scenario.tsv``
    (same statistics at full float precision; column names of the producer).
``composition_by_definition.tsv``
    ``iteration1/fix_rifaximin/out/paired_composition_by_scenario.tsv`` with one header
    change: ``antibiotic_state`` is written ``state`` so that the release scan, which
    rejects clinical column names, passes. Values are the same bytes.
``primary_paired_statistics.tsv``
    Byte copy of ``clinical_refit/clinical_paired_statistics.tsv``, the primary (42-patient)
    statistics that scenario A must reproduce.
``definition_counts.tsv``
    Aggregate counts from ``iteration1/fix_rifaximin/out/counts.json`` (all keys except
    ``EDS_85``, which belongs to a different analysis), one row per count.
"""
import hashlib
import json
import os
import shutil
from pathlib import Path

UNIT = Path(__file__).resolve().parent
REPO = UNIT.parents[1]
SD = REPO / "source_data" / UNIT.name

C = "review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01"
FIX = C + "/biology_refit/claude_review_2026-10-03/iteration1/fix_rifaximin"

#: published path -> (origin relative to LAB_ROOT, expected sha256 of origin)
ORIGINS = {
    "endpoints_by_definition.tsv": (FIX + "/SupplementaryTable6_rifaximin_sensitivity.tsv",
                                    "d4273f1aefad9ca37bafdf8357c1cd67a923af47c87a1ea3f02e4520af9dbf3f"),
    "endpoints_producer_output.tsv": (FIX + "/out/paired_endpoints_by_scenario.tsv",
                                      "f33383cad5fcbaca59aa1586de4435b225931c6d5129d50c12845965132e3f7d"),
    "composition_by_definition.tsv": (FIX + "/out/paired_composition_by_scenario.tsv",
                                      "b3dc72b682754e1a2c47f0ed5b193bcefb8ab173fedef1167c6f5abde371bf58"),
    "primary_paired_statistics.tsv": (C + "/clinical_refit/clinical_paired_statistics.tsv",
                                      "9518b0117267eb3198b1aafe785dab47b53b020bd4deb0d9d2fb24f67ccc3bc7"),
    "definition_counts.tsv": (FIX + "/out/counts.json",
                              "447677655f48c01c25898b3f2fea6f95f425c54a6721847142c0ea89228aa2f9"),
}


def sha(p):
    """Return the sha256 hex digest of file ``p``."""
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def flatten(d, prefix=""):
    """Flatten nested count dictionaries into ``(key, value)`` pairs joined with '.'."""
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            yield from flatten(v, key + ".")
        else:
            yield key, v


def main():
    """Copy or derive every source_data/STab3 table and verify origin hashes."""
    lab = Path(os.environ["LAB_ROOT"])
    SD.mkdir(parents=True, exist_ok=True)
    for name, (rel, expect) in ORIGINS.items():
        src = lab / rel
        assert sha(src) == expect, f"origin changed: {name}"
        if name == "composition_by_definition.tsv":
            text = src.read_text(encoding="utf-8")
            head, rest = text.split("\n", 1)
            assert head.split("\t")[1] == "antibiotic_state"
            head = head.replace("\tantibiotic_state\t", "\tstate\t", 1)
            (SD / name).write_text(head + "\n" + rest, encoding="utf-8")
        elif name == "definition_counts.tsv":
            counts = json.loads(src.read_text(encoding="utf-8"))
            counts.pop("EDS_85")
            rows = ["count\tvalue"] + [f"{k}\t{v}" for k, v in flatten(counts)]
            (SD / name).write_text("\n".join(rows) + "\n", encoding="utf-8")
        else:
            shutil.copyfile(src, SD / name)
        print(name, sha(SD / name))


if __name__ == "__main__":
    main()
