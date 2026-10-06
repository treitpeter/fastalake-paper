#!/usr/bin/env python3
"""Copy the Figure 6 input tables from the lab project into ``source_data/Fig6/``.

This script documents how every file in ``source_data/Fig6/`` was made. It is not needed
to rebuild the figure; it needs the lab project directory, given as ``FASTALAKE_LAB_ROOT``
(the directory written ``LAB_ROOT`` in ``PROVENANCE.md``). ``FASTALAKE_PSEUDONYM_MAP``
names a file OUTSIDE the repository where the PREDICT pseudonym map is written; that map
is a controlled input and never enters the repository.

Transformations (all others are byte copies):

* ``capscan_terminal_participant_region.tsv`` and ``capscan_residue_participant_region.tsv``:
  rows of the A24 participant-region table with ``min_length == 9`` (the only length
  Figure 6 and its two test families use), decompressed and split by unit
  (``assigned_bare_sequence`` rows, used for residue composition, go to the second file).
  The full table (all lengths, 45 MB) exceeds the 20 MB file limit.
* ``capscan_histidine_primary_*.tsv``: written by ``panels.derive_histidine_inputs`` from
  the A28 CapScan family table, exactly as the original ``prep_histidine_inputs.py``.
  They are kept as files so that reviewers can open the plotted rows directly.
* ``predict_histidine_patient_enrichment.tsv``: the PREDICT patient table with the raw
  patient identifier replaced by a pseudonym ``P001``..``P093`` assigned in the order the
  original table lists patients (ascending raw identifier). No other column changes.

Usage
-----
    FASTALAKE_LAB_ROOT=... FASTALAKE_PSEUDONYM_MAP=.../predict_pseudonyms.tsv \
        python figures/Fig6/prepare_source_data.py
"""

import gzip
import hashlib
import os
import shutil
import sys
from pathlib import Path

import pandas as pd

UNIT = Path(__file__).resolve().parent
REPO = UNIT.parents[1]
SD = REPO / "source_data" / UNIT.name
sys.path.insert(0, str(UNIT))
import panels  # noqa: E402

#: Producer locations relative to LAB_ROOT.
CR = ("review/v15_integration/jess_integration_2026-09-29/analysis_followup/"
      "manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03")
C2 = CR + "/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6"
SOURCES = {
    "a24_participant_region": C2 + "/A24_fixed347/results/CAPSCAN/per_participant_region_endpoints.tsv.gz",
    "a24_families": C2 + "/A24_fixed347/results/CAPSCAN/complete_statistical_families.tsv",
    "a24_design": C2 + "/A24_fixed347/results/CAPSCAN/acquisition_design.tsv",
    "a28_families": C2 + "/A28_fixed347/results_revision4/full_CAPSCAN/complete_statistical_families.tsv",
    "a28_participant_region": C2 + "/A28_fixed347/results_revision4/full_CAPSCAN/participant_region_enrichment.tsv",
    "a28_participant_reciprocal": C2 + "/A28_fixed347/results_revision4/full_CAPSCAN/participant_reciprocal_controls.tsv",
    "predict_source": "MM_review_analysis/end_of_day_2026-09-15/histidine/PREDICT_primary_source_enrichment.tsv",
    "predict_reciprocal": "MM_review_analysis/end_of_day_2026-09-15/histidine/PREDICT_primary_reciprocal_tests.tsv",
    "predict_families": "MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/complete_statistical_families.tsv",
    "predict_patient": "MM_review_analysis/end_of_day_2026-09-15/histidine/PREDICT_patient_source_enrichment.tsv",
    "capscan_primary_source_n15": CR + "/final/fig6_n15/inputs_n15/histidine/CAPSCAN_primary_source_enrichment.tsv",
    "capscan_primary_reciprocal_n15": CR + "/final/fig6_n15/inputs_n15/histidine/CAPSCAN_primary_reciprocal_tests.tsv",
}


def sha256(path):
    """Return the hex sha256 of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    """Write every source_data/Fig6 table and print origin and output hashes."""
    lab = Path(os.environ["FASTALAKE_LAB_ROOT"])
    pmap = Path(os.environ["FASTALAKE_PSEUDONYM_MAP"])
    assert REPO not in pmap.resolve().parents, "pseudonym map must stay outside the repository"
    src = {k: lab / v for k, v in SOURCES.items()}
    for k, p in src.items():
        print("origin", k, sha256(p.resolve()))
    SD.mkdir(parents=True, exist_ok=True)

    with gzip.open(src["a24_participant_region"], "rt") as fh:
        pp = pd.read_csv(fh, sep="\t", dtype={"participant": str}, keep_default_na=False, na_values=[""])
    pp9 = pp[pp.min_length.eq(9)]
    bare = pp9.unit.eq("assigned_bare_sequence")
    pp9[~bare].to_csv(SD / panels.F_A24_PARTICIPANT, sep="\t", index=False)
    pp9[bare].to_csv(SD / panels.F_A24_RESIDUE, sep="\t", index=False)
    for key, name in [("a24_families", panels.F_A24_FAMILIES), ("a24_design", panels.F_A24_DESIGN),
                      ("a28_families", panels.F_A28_FAMILIES),
                      ("a28_participant_region", panels.F_A28_PARTICIPANT),
                      ("a28_participant_reciprocal", panels.F_A28_RECIPROCAL),
                      ("predict_source", panels.F_PREDICT_SOURCE),
                      ("predict_reciprocal", panels.F_PREDICT_RECIPROCAL),
                      ("predict_families", panels.F_PREDICT_FAMILIES)]:
        shutil.copyfile(src[key], SD / name)

    a28 = pd.read_csv(SD / panels.F_A28_FAMILIES, sep="\t")
    s, r = panels.derive_histidine_inputs(a28)
    s.to_csv(SD / panels.F_CAPSCAN_SOURCE, sep="\t", index=False)
    r.to_csv(SD / panels.F_CAPSCAN_RECIPROCAL, sep="\t", index=False)
    for mine, theirs in [(panels.F_CAPSCAN_SOURCE, "capscan_primary_source_n15"),
                         (panels.F_CAPSCAN_RECIPROCAL, "capscan_primary_reciprocal_n15")]:
        same = sha256(SD / mine) == sha256(src[theirs])
        print("derived", mine, "byte-identical to producer input:", same)
        assert same

    pt = pd.read_csv(src["predict_patient"], sep="\t", dtype={"participant": str})
    order = list(dict.fromkeys(pt.participant))
    assert order == sorted(order, key=int), "patients are not listed in ascending raw-id order"
    pseudo = {raw: f"P{i:03d}" for i, raw in enumerate(order, 1)}
    pmap.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"raw_patient_id": list(pseudo), "pseudonym": list(pseudo.values())}).to_csv(pmap, sep="\t", index=False)
    pt.assign(participant=pt.participant.map(pseudo)).rename(columns={"participant": "patient"}).to_csv(
        SD / panels.F_PREDICT_PATIENT, sep="\t", index=False)
    print("pseudonym map written outside the repository:", pmap.name, sha256(pmap))

    for f in sorted(SD.glob("*.tsv")):
        print("output", f.name, f.stat().st_size, sha256(f))


if __name__ == "__main__":
    main()
