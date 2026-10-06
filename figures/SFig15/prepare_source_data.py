#!/usr/bin/env python3
"""Copy the Supplementary Figure 15 source tables from the lab tree into ``source_data/SFig15/``.

This script documents how ``source_data/SFig15/`` and ``figures/SFig15/expected/`` were
made. It is not needed to rebuild the figure, and it only runs where the lab tree exists.
``LAB_ROOT`` (environment variable) names the lab project directory; no absolute path is
stored here. Every origin file is checked against the sha256 recorded below before use.

Tables
------
compute_stage_summary.tsv
    Copied byte for byte from the 347-acquisition rerun of the producer's Computing block.
compute_task_records.tsv
    The task table that the stage summary is computed from (same rerun). One change: the
    ``sample`` column held raw acquisition file names; it is replaced by a pseudonymous
    ``acquisition`` id (``<cohort>_A<rank>``, rank of the file name within its cohort,
    sorted). Every other byte of every row is unchanged. The figure and all numbers use
    no file name.
campi_alphanovo_inference_sacct.tsv
    The Slurm accounting export for the CAMPI AlphaNovo inference arrays, written by the
    upstream compute-provenance step with ``sacct --parsable2 --noheader``. This script
    adds the header named by that step's ``--format`` and turns ``|`` into tabs; values
    are unchanged.

Expected renders (copied, not re-rendered)
------------------------------------------
panels_abcd.pdf   producer step 1 output (figure with in-image title)
SFig15_v28f.pdf   producer step 2 output (title redacted), the v28f vector figure
SFig15_v28f.png   producer step 2 raster (300 dpi), the v28f image

Usage
-----
    LAB_ROOT=/path/to/lab/project python figures/SFig15/prepare_source_data.py
"""

import hashlib
import os
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SD = REPO / "source_data" / "SFig15"
EXP = Path(__file__).resolve().parent / "expected"

FINAL = ("review/v15_integration/jess_integration_2026-09-29/analysis_followup/"
         "manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final")
RA = FINAL + "/reruns_347/A7_other/resource_accounting/fixed347/source_data"
REV = ("00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/"
       "nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/"
       "reappraisal_2026-09-24/supplements/revision_2026-09-24")
SF15 = FINAL + "/suppfigs_n15/SF15/new"

ORIGINS = {
    "summary": (RA + "/compute_stage_summary.tsv",
                "37c1d94d821d5ba9defb5ff5c8495905ce93209a136952dec8b600b74816498f"),
    "records": (RA + "/compute_task_records.tsv",
                "ea24e72fe85b7e4c0755a849d5822a5347d9d3685e73f765a641072c2d7551b8"),
    "sacct": (REV + "/source_data/compute/INFERENCE_SACCT.tsv", None),
    "panels_abcd.pdf": (SF15 + "/COMPUTE_RESOURCE_MAP.pdf", None),
    "SFig15_v28f.pdf": (SF15 + "/Supplementary_Figure_15.pdf",
                        "02d9c5ee79cf5e33a87912ec3e6b6612ffbfe01a65dbd67b523b5a91c6a11fdb"),
    "SFig15_v28f.png": (SF15 + "/Supplementary_Figure_15.png",
                        "b343e653aa38c78d90f6753f25225eec8ad4077f5ff2b2b6d85d31dfc286fc56"),
}

#: Column names of the sacct export, from the upstream ``--format`` argument.
SACCT_COLUMNS = ["JobID", "JobName", "State", "ElapsedRaw", "NCPUS", "TotalCPU", "MaxRSS",
                 "Start", "End", "NodeList"]


def sha256(p):
    """Return the hex sha256 of a file."""
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def origin(root, key):
    """Return the origin path for ``key`` after checking its recorded sha256."""
    rel, want = ORIGINS[key]
    p = root / rel
    got = sha256(p)
    if want is not None and got != want:
        sys.exit(f"{key}: sha256 {got} differs from recorded {want}")
    print(f"{key}\t{got}")
    return p


def pseudonymise_records(src, dst):
    """Write the task table with the ``sample`` column replaced by a pseudonymous id.

    Parameters
    ----------
    src, dst : pathlib.Path
        Origin and output TSV.

    Returns
    -------
    int
        Number of rows whose sample field was replaced.
    """
    lines = src.read_text().split("\n")
    head = lines[0].split("\t")
    i = head.index("sample")
    ci = head.index("cohort")
    rows = [ln.split("\t") for ln in lines[1:] if ln]
    names = {}
    for r in rows:
        if r[i]:
            names.setdefault(r[ci], set()).add(r[i])
    key = {}
    for co, s in names.items():
        for k, n in enumerate(sorted(s), start=1):
            key[(co, n)] = f"{co}_A{k:03d}"
    n = 0
    for r in rows:
        if r[i]:
            r[i] = key[(r[ci], r[i])]
            n += 1
    head[i] = "acquisition"
    out = ["\t".join(head)] + ["\t".join(r) for r in rows]
    dst.write_text("\n".join(out) + "\n")
    return n


def main():
    """Copy and transform every file; print sha256 of origins and outputs."""
    root = Path(os.environ["LAB_ROOT"])
    SD.mkdir(parents=True, exist_ok=True)
    EXP.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(origin(root, "summary"), SD / "compute_stage_summary.tsv")
    n = pseudonymise_records(origin(root, "records"), SD / "compute_task_records.tsv")
    print("pseudonymised rows", n)
    raw = origin(root, "sacct").read_text().splitlines()
    out = ["\t".join(SACCT_COLUMNS)] + [ln.replace("|", "\t") for ln in raw if ln]
    (SD / "campi_alphanovo_inference_sacct.tsv").write_text("\n".join(out) + "\n")
    for name in ("panels_abcd.pdf", "SFig15_v28f.pdf", "SFig15_v28f.png"):
        shutil.copyfile(origin(root, name), EXP / name)
    for f in sorted(SD.glob("*.tsv")) + sorted(EXP.glob("*")):
        print(f.relative_to(REPO), sha256(f))


if __name__ == "__main__":
    main()
