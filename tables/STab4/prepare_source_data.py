"""Write source_data/STab4/ and expected/ from the published workbook (run once, by the authors).

This script documents how the Supplementary Table 4 inputs of this repository were
derived. It is not run by the notebook and needs the lab tree, named by the environment
variable ``FASTALAKE_LAB_ROOT``. Receipts go to ``FASTALAKE_WORK_DIR`` (outside the
repository).

Inputs (relative to LAB_ROOT, sha256 checked):

* ``PKG_ST4``: Supplementary_Table_4.xlsx of the Source Data package (the table as
  published with v28d).
* ``PROD_ST4``: the workbook written by the original producer (``tables_cr003.py``,
  variant verified50, which named it Supplementary_Table_3.xlsx before renumbering).
* ``FIVE``: the five-pair adjustment table of the Source Data package
  (Biology_refit/eds_refit/separate_adjustments/all_models.tsv); it carries the
  bootstrap and 84-patient values of v28d line 799 that are not in the workbook.

Steps:

1. Compare the published workbook with the producer workbook cell for cell. The
   package builder replaced the Read_me sheet, removed the producer's
   Current_EDS_update sheet and replaced participant codes by pseudonyms; every other
   cell must be equal. Writes ``package_vs_producer.json`` (counts only) to the work dir.
2. Write the per-observation and per-patient tables the workbook needs that are not
   already in ``source_data/SFig10/``. Clinical group, antibiotic state, acquisition
   and microbiome sample identifiers are withheld (CONVENTIONS section 8).
3. Write ``expected/Supplementary_Table_4_redacted.xlsx``: the published workbook with
   the withheld columns removed and sample identifiers replaced by observation index
   (259-observation sheets) or patient display position (85-patient sheets). Nothing
   else is changed.

Usage
-----
    FASTALAKE_LAB_ROOT=<lab root> FASTALAKE_WORK_DIR=<work dir> \\
        python tables/STab4/prepare_source_data.py
"""

import hashlib
import json
import os
import shutil
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

UNIT = Path(__file__).resolve().parent
REPO = UNIT.parents[1]
OUT = REPO / "source_data" / "STab4"
SF10 = REPO / "source_data" / "SFig10"
LAB = Path(os.environ["FASTALAKE_LAB_ROOT"])
WORK = Path(os.environ["FASTALAKE_WORK_DIR"])
BASE = ("review/v15_integration/jess_integration_2026-09-29/analysis_followup/"
        "manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/")
PKG = LAB / BASE / "final/figure_sd/source_data/package/FastaLake_Source_Data"
PKG_ST4 = PKG / "Supplementary_Tables/Supplementary_Table_4.xlsx"
PROD_ST4 = LAB / BASE / "iteration2/cr003_followup/tables/verified50/Supplementary_Table_3.xlsx"
FIVE = PKG / "Biology_refit/eds_refit/separate_adjustments/all_models.tsv"

INPUT_SHA = {
    PKG_ST4: "4fe89a55c7607b0e5932377cbd488c031bc1f0c7f99c255f6c9eb7d77747e403",
    PROD_ST4: "46ea54c41fb867a6810e672a391bffd0804aeb190f62ce4e05c48d7b554ccd42",
    FIVE: "1f6400e71348f3d87d5248f407f7b9f8f650f4334e9a5eddf21938d315ea1028",
}

#: Columns of the published workbook that are withheld from the repository.
WITHHELD = {"sample_id", "microbiome_id", "clinical_group", "antibiotic_state"}
#: Sheets whose columns (after feature_id) are the 259 observations.
OBS_MATRIX = {"All20_features_259_observations"}
#: Sheets whose columns (after feature_id) are the 85 patients.
PAT_MATRIX = {"85_patient_raw", "85_patient_z"}
#: Sheets with one row per observation (259) or per patient (85).
OBS_ROWS = {"Observation_metadata"}
PAT_ROWS = {"85_patient_order", "Focused_raw"}


def sha(p):
    """Return the hex sha256 of a file."""
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def cells(path):
    """Return {sheet: list of row tuples} with cached values."""
    b = openpyxl.load_workbook(path, data_only=True)
    return {s.title: [tuple(r) for r in s.iter_rows(values_only=True)] for s in b.worksheets}


def compare_package_producer(pkg, prod):
    """Count cell differences per sheet; only participant-code columns may differ."""
    report = {"package_sheets": list(pkg), "producer_sheets": list(prod), "sheets": {}}
    for name, rows in pkg.items():
        if name == "Read_me":
            report["sheets"][name] = "replaced by the Source Data package builder"
            continue
        other = prod[name]
        assert len(rows) == len(other), name
        header = rows[0]
        diff_cols = {}
        for r1, r2 in zip(rows[1:], other[1:]):
            for j, (x, y) in enumerate(zip(r1, r2)):
                if x != y and not (isinstance(x, float) and isinstance(y, float) and np.isnan(x) and np.isnan(y)):
                    diff_cols[header[j]] = diff_cols.get(header[j], 0) + 1
        assert header == other[0], name
        assert set(diff_cols) <= {"patient_id"}, (name, diff_cols)
        report["sheets"][name] = {"rows": len(rows), "cells_differing_by_column": diff_cols}
    report["producer_only_sheets"] = [s for s in prod if s not in pkg]
    return report


def redact(pkg, obs_ids, pat_pos):
    """Return {sheet: (rows, freeze)} of the published workbook with withheld columns removed."""
    out = {}
    for name, rows in pkg.items():
        header = list(rows[0])
        body = [list(r) for r in rows[1:]]
        if name in OBS_MATRIX or name in PAT_MATRIX:
            ids = obs_ids if name in OBS_MATRIX else pat_pos
            new_header = [header[0]] + [ids[s] for s in header[1:]]
            out[name] = [new_header] + body
            continue
        if name in OBS_ROWS or name in PAT_ROWS:
            key = "observation_index" if name in OBS_ROWS else "patient_position"
            ids = obs_ids if name in OBS_ROWS else pat_pos
            sj = header.index("sample_id")
            keep = [j for j, h in enumerate(header) if h not in WITHHELD]
            new_header = [key] + [("observation_index" if header[j] == "index" else header[j]) for j in keep]
            new_body = []
            for r in body:
                vals = [(r[j] + 1 if header[j] == "index" else r[j]) for j in keep]
                new_body.append([ids[r[sj]]] + vals)
            out[name] = [new_header] + new_body
            continue
        assert not (set(header) & WITHHELD), name
        out[name] = [header] + body
    return out


def main():
    """Derive source_data/STab4 and expected/ from the published workbook."""
    for p, h in INPUT_SHA.items():
        assert sha(p) == h, p
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    (UNIT / "expected").mkdir(exist_ok=True)

    pkg, prod = cells(PKG_ST4), cells(PROD_ST4)
    rep = compare_package_producer(pkg, prod)
    (WORK / "package_vs_producer.json").write_text(json.dumps(rep, indent=1) + "\n")

    # 259 observations, in the column order of the published 259-observation sheet
    om = pd.DataFrame(pkg["Observation_metadata"][1:], columns=pkg["Observation_metadata"][0])
    assert list(om.sample_id) == list(pkg["All20_features_259_observations"][0][1:])
    assert len(om) == 259
    obs_ids = {s: i + 1 for i, s in enumerate(om.sample_id)}
    obs = pd.DataFrame({"observation_index": np.arange(1, 260), "patient_id": om.patient_id.fillna(""),  # 3 observations have no patient id in the producer table either
                        "Human": om.Human.astype(float), "pg_parsimony": om.pg_parsimony.astype(int),
                        "dysbiosis_score": om.dysbiosis_score.astype(float)})

    po = pd.DataFrame(pkg["85_patient_order"][1:], columns=pkg["85_patient_order"][0])
    assert len(po) == 85 and list(po.sample_id) == list(pkg["85_patient_raw"][0][1:])
    pat_pos = {s: i + 1 for i, s in enumerate(po.sample_id)}
    assert all(om.sample_id.iloc[int(i)] == s for i, s in zip(po["index"], po.sample_id))
    link = pd.DataFrame({"patient_position": np.arange(1, 86),
                         "observation_index": po["index"].astype(int).to_numpy() + 1})
    # the 85-patient table of SFig10 must be the same patients in the same order
    sp = pd.read_csv(SF10 / "patients.tsv", sep="\t", float_precision="round_trip")
    assert list(sp.patient_pseudonym) == list(po.patient_id)
    assert np.array_equal(sp.dysbiosis_score.to_numpy(), po.dysbiosis_score.astype(float).to_numpy())

    readme = pd.DataFrame({"line": np.arange(1, len(pkg["Read_me"]) + 1),
                           "text": [r[0] if r[0] is not None else "" for r in pkg["Read_me"]]})

    for df, name in [(obs, "observations_259.tsv"), (link, "patients_85_observation_index.tsv"),
                     (readme, "read_me.tsv")]:
        df.to_csv(OUT / name, sep="\t", index=False)
        back = pd.read_csv(OUT / name, sep="\t", float_precision="round_trip", keep_default_na=False)
        for c in df.columns:
            if pd.api.types.is_float_dtype(df[c]):
                assert np.array_equal(df[c].to_numpy(), back[c].astype(float).to_numpy()), (name, c)
            else:
                assert df[c].astype(str).equals(back[c].astype(str)), (name, c)
    shutil.copy2(FIVE, OUT / "five_pair_all_models.tsv")

    red = redact(pkg, obs_ids, pat_pos)
    src = openpyxl.load_workbook(PKG_ST4)
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    for name, rows in red.items():
        sh = wb.create_sheet(name)
        for r in rows:
            sh.append(list(r))
        sh.freeze_panes = src[name].freeze_panes
    wb.save(UNIT / "expected" / "Supplementary_Table_4_redacted.xlsx")

    receipt = {p.name: sha(p) for p in sorted(OUT.glob("*.tsv"))}
    receipt["expected/Supplementary_Table_4_redacted.xlsx"] = sha(UNIT / "expected" / "Supplementary_Table_4_redacted.xlsx")
    (WORK / "prepare_receipt.json").write_text(json.dumps(dict(inputs={str(p.relative_to(LAB)): h for p, h in INPUT_SHA.items()},
                                                               outputs=receipt, package_vs_producer=rep), indent=1) + "\n")
    print(json.dumps(receipt, indent=1))
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
