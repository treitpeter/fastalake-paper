#!/usr/bin/env python3
"""Make every file in ``source_data/SFig13/`` from the lab project directory.

This script documents how the Source Data of Supplementary Figure 13 were made. It is not
needed to rebuild the figure. It needs the lab project directory as ``FASTALAKE_LAB_ROOT``
(written ``LAB_ROOT`` in ``PROVENANCE.md``). Every origin sha256 is asserted.

Transformations
---------------
* ``capscan_participant_motif_differences.tsv.gz``: byte copy of the A28 (347 acquisitions,
  15 participants) participant-region motif table.
* ``predict_histidine_statistical_families.tsv``: byte copy of the A28 PREDICT families.
* ``capscan_histidine_motif_display.tsv``: the 840 displayed cells, written by
  ``panels.derive_display_table`` from ``source_data/Fig6/capscan_histidine_statistical_families.tsv``
  exactly as the producer wrote it; asserted byte-identical to the producer's Source Data.
* ``predict_histidine_acquisition_enrichment.tsv.gz``: the A28 PREDICT acquisition table with
  identifying columns replaced, as text, so that every numeric field is unchanged:
  ``participant`` (raw patient identifier) -> ``P001``..``P093`` in ascending raw-identifier
  order (the same assignment Figure 6 uses); ``region`` (PREDICT visit code) -> ``withheld``;
  ``specimen_id`` (contains patient and visit codes) -> ``P###_S##`` numbered in the sort
  order of (visit code, specimen id) within patient; ``sample`` (raw file name) ->
  ``P###_S##_A##`` numbered in raw file name order within specimen. Acquisitions without a
  patient (excluded by every analysis) keep an empty patient and visit and get ``U##``.
  Because the numbering preserves the sort order of the original keys, every group-by in the
  analysis code sees the groups in the same order. The mapping is not written anywhere.

Usage
-----
    FASTALAKE_LAB_ROOT=... python figures/SFig13/prepare_source_data.py
"""

import gzip
import io
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

CR = ("review/v15_integration/jess_integration_2026-09-29/analysis_followup/"
      "manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03")
C2 = CR + "/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6"
PRED = "MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT"
SOURCES = {
    "capscan_motif": (C2 + "/A28_fixed347/results_revision4/full_CAPSCAN/participant_motif_differences.tsv.gz",
                      "df84cc2720b0ed63448048ef0676ad33c44ef9fe46eb05414a1d98a3ecb5ae77"),
    "predict_acq": (PRED + "/acquisition_enrichment.tsv.gz",
                    "8e49c01bac8e2af1f8f2570c4718fc4bf7dea5d648af1ff62f6f5d88cc14c1ab"),
    "predict_families": (PRED + "/complete_statistical_families.tsv",
                         "85c774cb604a1d9e6c9526b04750b10d6601a07855fe1552357f66be016d28d7"),
    "producer_display": (CR + "/final/suppfigs_n15/SF13/source_data/CapScan_histidine_motif_display.tsv",
                         "269fa990cd52029a6dd4961e7bf21e347436c374c89bee4fa5c9631dd87364ea"),
}


def pseudonymise_predict(src, dst):
    """Rewrite the PREDICT acquisition table with pseudonymous keys (see module docstring)."""
    with gzip.open(src, "rt", newline="") as fh:
        lines = fh.read().split("\n")
    header = lines[0].split("\t")
    assert lines[-1] == ""
    rows = [ln.split("\t") for ln in lines[1:-1]]
    ip, ir, isp, isa = (header.index(c) for c in ["participant", "region", "specimen_id", "sample"])
    pts = sorted({r[ip] for r in rows if r[ip] != ""}, key=int)
    assert len(pts) == 93 and len({len(p) for p in pts}) == 1  # same order as text and as number
    pmap = {p: f"P{i:03d}" for i, p in enumerate(pts, 1)}
    spec = {}
    for p in pts:
        keys = sorted({(r[ir], r[isp]) for r in rows if r[ip] == p})
        for j, k in enumerate(keys, 1):
            spec[(p,) + k] = f"{pmap[p]}_S{j:02d}"
    assert len({k[2] for k in spec}) == len(spec)  # specimen ids unique across patients
    per_spec = {}
    for r in rows:
        if r[ip] != "":
            per_spec.setdefault((r[ip], r[ir], r[isp]), set()).add(r[isa])
    acq = {}
    for k, ps in spec.items():
        for j, a in enumerate(sorted(per_spec[k]), 1):
            assert a not in acq
            acq[a] = f"{ps}_A{j:02d}"
    un = sorted({r[isa] for r in rows if r[ip] == ""})
    for j, a in enumerate(un, 1):
        acq[a] = f"U{j:02d}"
    out = ["\t".join(header)]
    for r in rows:
        r = list(r)
        if r[ip] == "":
            assert r[ir] == ""
            r[isp] = acq[r[isa]]
        else:
            r[isp] = spec[(r[ip], r[ir], r[isp])]
            r[ip] = pmap[r[ip]]
            r[ir] = "withheld"
        r[isa] = acq[r[isa]]
        out.append("\t".join(r))
    data = ("\n".join(out) + "\n").encode()
    with open(dst, "wb") as fh, gzip.GzipFile(filename="", mode="wb", fileobj=fh, mtime=0, compresslevel=9) as gz:
        gz.write(data)
    return dict(patients=len(pts), specimens=len(spec), acquisitions=len(acq) - len(un), unassigned=len(un))


def main():
    """Write every source_data/SFig13 file and print origin and output hashes."""
    lab = Path(os.environ["FASTALAKE_LAB_ROOT"])
    src = {}
    for k, (rel, sha) in SOURCES.items():
        src[k] = lab / rel
        panels.check_sha256(src[k], sha)
        print("origin", k, sha)
    SD.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src["capscan_motif"], SD / panels.F_CAPSCAN_PARTICIPANT_MOTIF)
    shutil.copyfile(src["predict_families"], SD / panels.F_PREDICT_FAMILIES)

    fam = REPO / "source_data" / panels.F_CAPSCAN_FAMILIES[0] / panels.F_CAPSCAN_FAMILIES[1]
    panels.check_sha256(fam, panels.SHA256[panels.F_CAPSCAN_FAMILIES[1]])
    _, display, text = panels.derive_display_table(fam)
    buf = io.StringIO(); display.to_csv(buf, sep="\t", index=False)
    (SD / panels.F_CAPSCAN_DISPLAY).write_text(buf.getvalue())
    assert panels.sha256(SD / panels.F_CAPSCAN_DISPLAY) == panels.sha256(src["producer_display"])
    print("display table byte-identical to producer Source Data: True")

    print("pseudonymised:", pseudonymise_predict(src["predict_acq"], SD / panels.F_PREDICT_ACQUISITION))
    for f in sorted(SD.iterdir()):
        if f.is_file():
            print("output", f.name, f.stat().st_size, panels.sha256(f))


if __name__ == "__main__":
    main()
