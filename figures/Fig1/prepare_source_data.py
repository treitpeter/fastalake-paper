#!/usr/bin/env python3
"""Copy the Figure 1 source tables from the lab tree into ``source_data/Fig1/``.

This script documents how ``source_data/Fig1/`` was made. It is not needed to rebuild the
figure, and it only runs where the lab tree exists. ``LAB_ROOT`` (environment variable)
names the lab project directory; no absolute path is stored here.

Seven tables and one JSON are copied byte for byte (``COPIES``). Two JSON records are
reduced: absolute file paths and modification times are removed, every count and sha256
is kept unchanged (``REDUCED``). The script checks the sha256 of every origin file before
copying and refuses to run if any differs from the value recorded here.

Usage
-----
    LAB_ROOT=/path/to/lab/project python figures/Fig1/prepare_source_data.py
"""

import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "source_data" / "Fig1"

FINAL = ("review/v15_integration/jess_integration_2026-09-29/analysis_followup/"
         "manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final")
PKG = FINAL + "/figure_sd/source_data/package/FastaLake_Source_Data"

#: published name -> (origin relative to LAB_ROOT, sha256 of origin)
COPIES = {
    "acquisition_sizes_726.tsv": (FINAL + "/fig1bd_n15/acquisition_sizes_726.tsv",
        "fd9ade3bb4ebd2ee6c129f69819b9697cf20abd7c10c77a12644bef25396e386"),
    "acquisition_support_summary_726.tsv": (FINAL + "/fig1bd_n15/acquisition_support_summary_726.tsv",
        "6f81da0a03893bcdf2ed5406ddb0e52eb16791cfd46f55a3253b31c65ac7dec2"),
    "length_by_arm.tsv": ("MM_review_analysis/A02_4_2_lookup_specificity/results/length_by_arm.tsv",
        "63594d955d35bd564a32fe68317274537519376858d3bc0daee1b09cc3b1cc29"),
    "five_seed_summary.tsv": ("MM_review_analysis/A02_4_2_lookup_specificity/results/five_seed_summary.tsv",
        "ebb3de14cd36ab997a1e0e2245e882de30f822b3b06d78cbf42e8af79bf73cd6"),
    "donor_seed_means.tsv": ("MM_review_analysis/A11_4_2_control_reconciliation/results/donor_seed_means.tsv",
        "5c6b0af7c7224b7f1968b82a8a2c22ecae4d25679f691ef04c0957644c27ea66"),
    "donor_recount.tsv": ("MM_review_analysis/A11_4_2_control_reconciliation/results/donor_recount.tsv",
        "4e51dd75c37deed3f13c84ad2407032b01ae16d877efac9d1b60a174aafb5ea8"),
    "scramble_end_to_end.tsv": (PKG + "/Review_v16/Terminal_scramble/scramble_end_to_end.tsv",
        "844f92521718387c35fea524b5cc8502cac8847e0064dfea04c27637293ce865"),
    "selection_grid_acquisitions.tsv": (PKG + "/Supplementary_Figure_02/selection_grid_acquisitions.tsv",
        "8b472997023e29d4e24b5a9f7a06d9aceda6f45fc4c45391825c075d85aa3425"),
    "stage1_stats.json": (PKG + "/Supplementary_Figure_02/reference58/stage1_stats.json",
        "2be7eafbd02937b6df5ac05e1845c1a4d1d64203bfc62bd77ebb23739a21fcd8"),
}

#: published name -> (origin relative to LAB_ROOT, sha256 of origin); see reduce_record
REDUCED = {
    "human_reference.json": ("MM_review_analysis/A01_4_1_database_sizes/results/human_reference.json",
        "b705dfc4aba8718c83939a812e3de6614ce0a72eb9260c3e0c7b6ff632cf1079"),
    "worked_acquisition_prediction_support.json": (
        "MM_review_analysis/A09_4_1_prediction_support/results/HSTOOL/"
        "20181231_QX7_PhGe_LC95_SA_WorkflowReplicates_2/COMPLETE.json",
        "f4a34959914caa7dacbb3a3a9f85445d214cfbf6edf10f98134ba92e9c0fb066"),
}


def sha256(path):
    """Return the hex sha256 of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def reduce_record(obj):
    """Replace absolute ``path`` values by their basename and drop ``mtime_ns``.

    Parameters
    ----------
    obj : dict or list
        Parsed JSON record.

    Returns
    -------
    dict or list
        The same record without absolute paths or modification times.
    """
    if isinstance(obj, dict):
        return {k: (Path(v).name if k == "path" and isinstance(v, str) else reduce_record(v))
                for k, v in obj.items() if k != "mtime_ns"}
    if isinstance(obj, list):
        return [reduce_record(v) for v in obj]
    return obj


def main():
    """Copy and reduce the source tables; exit non-zero on any hash mismatch."""
    lab = os.environ.get("LAB_ROOT")
    if not lab:
        sys.exit("set LAB_ROOT to the lab project directory")
    lab = Path(lab)
    OUT.mkdir(parents=True, exist_ok=True)
    for name, (rel, want) in COPIES.items():
        src = lab / rel
        got = sha256(src)
        if got != want:
            sys.exit(f"{name}: origin sha256 {got} != recorded {want}")
        shutil.copyfile(src, OUT / name)
        print(f"copied {name} {got}")
    for name, (rel, want) in REDUCED.items():
        src = lab / rel
        if sha256(src) != want:
            sys.exit(f"{name}: origin sha256 differs from recorded {want}")
        rec = reduce_record(json.loads(src.read_text()))
        rec["_origin_sha256"] = sha256(src)
        (OUT / name).write_text(json.dumps(rec, indent=1) + "\n")
        print(f"reduced {name} from {rec['_origin_sha256']}")


if __name__ == "__main__":
    main()
