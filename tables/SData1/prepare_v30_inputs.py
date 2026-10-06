#!/usr/bin/env python3
"""Copy the 347-roster inputs of the v30 Supplementary Data 1 into ``source_data/SData1/`` (run once, by the authors).

The v30 workbook re-indexes the five CapScan analyses that the v28f build left on the
331-acquisition roster (``panels.V28F_NOT_REINDEXED_IDS``). Their memberships come from
the participant-15 reruns, read from the lab project directory named by ``LAB_ROOT``:

* ``v30_capscan_crossomics_specimens.tsv``: specimens with a positive microbial-OG
  acquisition and a metabolite measurement on the 347 roster (``capscan_specimen_roster.tsv``
  of the 347 cross-omics run, upstream ``capscan_crossomics``; columns specimen,
  participant and region only, the device and covariate columns are dropped);
* ``v30_capscan_crossomics_participants.tsv``: per participant, whether protein and
  metabolite profiles cover all five sampling groups and whether DNA profiles do
  (``COMPLETE.json`` of the same run);
* ``v30_capscan_metabolite_specimens.tsv``: specimen identifiers with a row in the
  extended CapScan metabolite table used by that run (identifiers only, no values);
* ``v30_capscan_h_design.tsv``: acquisition, participant, region and specimen of the 347
  CapScan tasks of the fixed-347 A28 histidine run (its input manifest);
* ``v30_capscan_h_finite_n.tsv`` and ``v30_capscan_h_participants.tsv``: the histidine
  finite-n and finite-participant summaries of the fixed-347 A28 run (upstream
  ``terminal_histidine_inference``), computed with the code of
  ``pipelines/sample_rosters_workbook/build_histidine_scope.py`` (``summarise_histidine``);
* ``v30_inputs_receipt.tsv``: sha256 of every file read.

Usage
-----
    LAB_ROOT=<lab project dir> python tables/SData1/prepare_v30_inputs.py
"""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

UNIT = Path(__file__).resolve().parent
REPO = UNIT.parents[1]
SD = REPO / "source_data" / "SData1"
LAB_ROOT = Path(os.environ["LAB_ROOT"])
FINAL = ("review/v15_integration/jess_integration_2026-09-29/analysis_followup/"
         "manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03")
C2 = LAB_ROOT / FINAL / "coauthor_comment_check_v25_2026-10-05/compute/C2_participant15"
CROSS = C2 / "fig5g/crossomics_347"
H347 = C2 / "fig6/A28_fixed347/results_revision4/full_CAPSCAN"
METABOLITES = (LAB_ROOT / "review/v15_integration/jess_integration_2026-09-29/cohort_reconciliation/"
               "crossomics/metabolomics/CAPSCAN_metabolite_intensity.tsv")
H_MANIFEST = C2 / "fig6/A28_fixed347/inputs/manifest_revision4.json"
READ = []


def sha(p):
    """sha256 of a file, streamed."""
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def note(p):
    """Record a file read (path relative to LAB_ROOT and sha256) and return it."""
    READ.append(dict(file="LAB_ROOT/" + Path(p).relative_to(LAB_ROOT).as_posix(), sha256=sha(p)))
    return p


def summarise_histidine(folder):
    """Finite-n and finite-participant histidine summaries, as in build_histidine_scope.py.

    Parameters
    ----------
    folder : Path
        ``results_revision4/full_CAPSCAN`` of an A28 run.

    Returns
    -------
    tuple of DataFrame
        ``(finite_n, participants)``; unit: statistical test and endpoint cell.
    """
    t = pd.read_csv(note(folder / "complete_statistical_families.tsv"), sep="\t")
    s = t.groupby(["family", "endpoint"], dropna=False).agg(
        saved_tests=("n", "size"), minimum_saved_n=("n", "min"), maximum_saved_n=("n", "max"),
        saved_estimable_tests=("estimable", "sum")).reset_index()
    s["scope"] = "Saved test n, not independently replayed endpoint eligibility or statistical tests"
    s["source_file"] = "complete_statistical_families.tsv"
    rows = []
    for name, kind, keys, metric in [("participant_region_enrichment.tsv", "H_source_difference", ["route", "side", "origin", "region"], "difference"),
                                     ("participant_reciprocal_controls.tsv", "complementary_product_difference", ["route", "origin", "region"], "difference")]:
        d = pd.read_csv(note(folder / name), sep="\t", dtype={"participant": str})
        if "residue" in d:
            d = d[d.residue.eq("H")]
        for labels, g in d.groupby(keys, sort=False):
            finite = np.isfinite(pd.to_numeric(g[metric], errors="coerce"))
            ids = sorted(g.loc[finite, "participant"].unique(), key=int)
            rows.append(dict(endpoint_kind=kind, **dict(zip(keys, labels)), finite_participants=len(ids),
                             participant_ids=";".join(ids), source_rows=len(g), source_file=name))
    return s, pd.DataFrame(rows)


def main():
    """Write the v30 input tables into source_data/SData1."""
    spec = pd.read_csv(note(CROSS / "capscan_specimen_roster.tsv"), sep="\t", dtype=str)
    spec = spec.rename(columns={spec.columns[0]: "specimen_id_crossomics", "subject_id": "participant_id", "group": "region"})
    spec = spec[["specimen_id_crossomics", "participant_id", "region"]]
    spec.to_csv(SD / "v30_capscan_crossomics_specimens.tsv", sep="\t", index=False)
    done = json.loads(Path(note(CROSS / "COMPLETE.json")).read_text())
    assert done["status"] == "PASS"
    pm = {str(p) for p in done["protein_metabolite_participants"]}
    dna = {str(p) for p in done["dna_participants"]}
    pd.DataFrame([dict(participant_id=str(p), protein_metabolite_complete=str(p) in pm, dna_complete=str(p) in dna)
                  for p in range(1, 16)]).to_csv(SD / "v30_capscan_crossomics_participants.tsv", sep="\t", index=False)
    met = pd.read_csv(note(METABOLITES), sep="\t", index_col=0, usecols=[0])
    pd.DataFrame(dict(specimen_id_crossomics=sorted({str(x) for x in met.index}, key=lambda x: (len(x), x)))).to_csv(
        SD / "v30_capscan_metabolite_specimens.tsv", sep="\t", index=False)
    hdone = json.loads(Path(note(H347 / "COMPLETE.json")).read_text())
    assert hdone["status"] == "PASS" and int(hdone["acquisitions"]) == 347
    m = json.loads(Path(note(H_MANIFEST)).read_text())
    tasks = [t for t in m["tasks"] if t["cohort"] == "CAPSCAN"]
    pd.DataFrame([dict(acquisition_id=t["sample"], participant_id=str(t["design"]["participant"]),
                       region=t["design"]["region"], histidine_specimen_id=t["design"]["specimen_id"])
                  for t in tasks]).to_csv(SD / "v30_capscan_h_design.tsv", sep="\t", index=False)
    fin, par = summarise_histidine(H347)
    fin.to_csv(SD / "v30_capscan_h_finite_n.tsv", sep="\t", index=False)
    par.to_csv(SD / "v30_capscan_h_participants.tsv", sep="\t", index=False)
    pd.DataFrame(READ).to_csv(SD / "v30_inputs_receipt.tsv", sep="\t", index=False)
    print("written:", len(spec), "specimens;", len(tasks), "histidine acquisitions;", len(fin), "families;", len(par), "endpoint cells")


if __name__ == "__main__":
    main()
