#!/usr/bin/env python3
"""Write source_data/SFig04/ and expected/ from the original producer's inputs and outputs.

This script documents how every table in ``source_data/SFig04/`` was derived. It is run
once by the maintainers, on Slurm, with ``FASTALAKE_LAB_ROOT`` set to the lab project
directory (``LAB_ROOT`` in ``PROVENANCE.md``). Reviewers do not need to run it: the
notebook reads only ``source_data/SFig04/``.

Stages
------
``source``   copy or reshape the producer inputs into ``source_data/SFig04/*.tsv`` and
             write the reconciliation against the Source Data package workbook
             (``checks/source_data_vs_package.tsv``).
``expected`` copy the original producer's composite (PDF, 400 dpi PNG) and plotted values
             into ``expected/``, and crop one reference image per panel from the
             producer's PNG using the panel boxes computed by ``panels.panel_boxes``.

Usage
-----
    FASTALAKE_LAB_ROOT=... FASTALAKE_FONT_DIR=... python prepare_source_data.py source
    FASTALAKE_LAB_ROOT=... FASTALAKE_FONT_DIR=... python prepare_source_data.py expected
"""

import glob
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

UNIT = Path(__file__).resolve().parent
REPO = UNIT.parents[1]
SD = REPO / "source_data" / UNIT.name
LAB = Path(os.environ["FASTALAKE_LAB_ROOT"])

#: M09 minimum-length-7 sensitivity analysis (original 331-acquisition CapScan roster).
M09 = LAB / ("00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/"
             "nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/"
             "reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M09_short_peptides")
RUN = M09 / "run_2026-09-24_minlen7"
#: Final review tree holding the 347-acquisition extensions and the v28f producer.
FINAL = LAB / ("review/v15_integration/jess_integration_2026-09-29/analysis_followup/"
               "manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final")
A2 = FINAL / "reruns_347/A2_minlen7"
P347 = FINAL / "public347/minlen7_public347"
PRODUCER_OUT = FINAL / "suppfigs_n15/SF04/n347_panel_c"
PACKAGE_XLSX = FINAL / ("figure_sd/source_data/package/FastaLake_Source_Data/Supplementary_Figure_04/"
                        "LENGTH_SENSITIVITY_source_data.xlsx")

#: The q threshold of the fusion-control rows used in panel e (file POOLED_q01.tsv).
FUSION_Q = 0.01


def sha(p):
    """Return the hex sha256 of a file."""
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write(df, name):
    """Write one source-data table (tab-separated, no index) and return its path."""
    p = SD / name
    df.to_csv(p, sep="\t", index=False)
    return p


def stage_source():
    """Copy or reshape every table the unit needs and reconcile with the package workbook."""
    SD.mkdir(parents=True, exist_ok=True)
    origin = {}

    # a, b and the line-777 numbers: matched public FastaLake vs metagenome comparisons.
    src = RUN / "wave2/peptide_summary/matched_comparisons.tsv"
    # Filtered as text lines, so every value keeps its original decimal string (a pandas
    # write/read round trip can change the last digit of a float).
    lines = src.read_text().splitlines(keepends=True)
    hdr = lines[0].rstrip("\n").split("\t"); ia, ic = hdr.index("arm"), hdr.index("comparator")
    keep = [lines[0]] + [l for l in lines[1:] if l.split("\t")[ia] == "FL_public" and l.split("\t")[ic] == "MG"]
    (SD / "benchmark_matched_per_acquisition.tsv").write_text("".join(keep))
    bm = pd.read_csv(SD / "benchmark_matched_per_acquisition.tsv", sep="\t")
    origin["benchmark_matched_per_acquisition.tsv"] = (src, "text lines with arm == FL_public and comparator == MG; all five length populations; columns and values unchanged")

    # c: per-acquisition accepted-peptide counts by length population, 1,990 searches.
    src = RUN / "wave2/peptide_summary/per_acquisition_wide_claude.tsv"
    shutil.copyfile(src, SD / "per_acquisition_length_counts.tsv")
    origin["per_acquisition_length_counts.tsv"] = (src, "byte copy")
    # c: the 16 participant-15 CapScan FL_public searches added for the 347 roster.
    src = P347 / "per_acquisition_new16.tsv"
    shutil.copyfile(src, SD / "per_acquisition_length_counts_participant15.tsv")
    origin["per_acquisition_length_counts_participant15.tsv"] = (src, "byte copy")

    # d: CAMPI member-specific pooled peptide counts.
    src = RUN / "members/MEMBER_COUNTS_BY_LENGTH.tsv"
    shutil.copyfile(src, SD / "campi_member_counts.tsv")
    origin["campi_member_counts.tsv"] = (src, "byte copy")

    # e: pooled fusion-control diagnostic and the per-acquisition rows it pools.
    src = RUN / "fusion375/POOLED_q01.tsv"
    shutil.copyfile(src, SD / "campi_fusion_pooled_q01.tsv")
    origin["campi_fusion_pooled_q01.tsv"] = (src, "byte copy")
    files = sorted(glob.glob(str(RUN / "fusion375/receipts/summary_*.tsv")))
    assert len(files) == 30, len(files)
    out_lines = []
    for f in files:
        lines = Path(f).read_text().splitlines(keepends=True)
        hdr = lines[0].rstrip("\n").split("\t"); iq = hdr.index("q_threshold")
        out_lines = out_lines or [lines[0]]
        assert lines[0] == out_lines[0]
        out_lines += [l for l in lines[1:] if float(l.split("\t")[iq]) == FUSION_Q]
    assert len(out_lines) == 151, len(out_lines)
    (SD / "campi_fusion_per_acquisition_q01.tsv").write_text("".join(out_lines))
    origin["campi_fusion_per_acquisition_q01.tsv"] = (RUN / "fusion375/receipts/summary_00.tsv .. summary_29.tsv",
                                                     "30 files concatenated in file-name order as text lines; rows q_threshold == 0.01; values unchanged")

    # f: healthy-stool common-FASTA and joint searches.
    src = RUN / "wave2/healthy_strategies/per_acquisition.tsv"
    shutil.copyfile(src, SD / "healthy_strategies_per_acquisition.tsv")
    origin["healthy_strategies_per_acquisition.tsv"] = (src, "byte copy")

    # g: study-level protein groups; CapScan from the 347-acquisition rerun, others from M09.
    rows = []
    files = sorted(glob.glob(str(RUN / "grouping/COMPARE_*.json")))
    files = [str(A2 / "run/grouping/COMPARE_CAPSCAN_noclust.json") if f.endswith("COMPARE_CAPSCAN_noclust.json") else f
             for f in files]
    for f in files:
        d = json.load(open(f))
        c = d["comparison"]
        rows.append(dict(task=d["task"], roster="347" if "reruns_347" in f else "M09",
                         arms_used=";".join(f"{k}:{v}" for k, v in d["arms_used"].items()),
                         acquisitions_baseline=c["acquisitions"]["baseline"], acquisitions_fresh7=c["acquisitions"]["fresh7"],
                         baseline=c["inferred_study_groups"]["baseline"], fresh7=c["inferred_study_groups"]["fresh7"]))
        origin.setdefault("grouping_study_groups.tsv", ([], "one row per COMPARE_<task>.json: task, arms_used, comparison.acquisitions, comparison.inferred_study_groups"))[0].append(Path(f))
    write(pd.DataFrame(rows), "grouping_study_groups.tsv")

    # Search inventory behind the legend count (1,990 in v28d, 2,022 in v28f).
    inv = []
    for run, root, idx_ok in [("M09_331", RUN, None), ("A2_347", A2 / "run", set(int(x) for x in (A2 / "run/new_indices.txt").read_text().split())),
                              ("public347", P347, None)]:
        t = pd.read_csv(root / "tasks.tsv", sep="\t")
        if idx_ok is not None:
            t = t[t["index"].isin(idx_ok)]
        for r in t.itertuples():
            rc = json.load(open(root / "receipts" / f"receipt_{r.index:04d}.json"))
            assert (rc["cohort"], rc["arm"], rc["sample"]) == (r.cohort, r.arm, r.sample)
            inv.append(dict(run=run, task_index=r.index, cohort=r.cohort, arm=r.arm, sample=r.sample, lfq=bool(r.lfq),
                            sage_version=rc["sage_version"], sage_exit=rc["sage_exit"], receipt_pass=bool(rc["pass"])))
    inv = pd.DataFrame(inv)
    assert len(inv) == 2022 and inv.receipt_pass.all(), (len(inv), inv.receipt_pass.sum())
    write(inv, "search_inventory.tsv")
    origin["search_inventory.tsv"] = ([RUN / "tasks.tsv", A2 / "run/tasks.tsv", A2 / "run/new_indices.txt", P347 / "tasks.tsv"],
                                      "task rows (A2: new indices only) joined to receipts/receipt_<index>.json; config and path columns dropped")

    # Reconciliation with the Source Data package workbook (sheet by sheet).
    xl = pd.read_excel(PACKAGE_XLSX, sheet_name=None)
    rec = []
    ab = xl["a_b_benchmark_per_acquisition"]
    sel = bm[bm.population.isin(["baseline5", "fresh7", "nested9"])]
    rec.append(("benchmark_matched_per_acquisition.tsv", "a_b_benchmark_per_acquisition", len(ab), len(sel),
                bool(len(ab) == len(sel) and np.allclose(np.sort(ab.gain_pct), np.sort(sel.gain_pct))),
                "sorted gain_pct of populations baseline5, fresh7, nested9"))
    w = pd.read_csv(SD / "per_acquisition_length_counts.tsv", sep="\t")
    c_ = xl["c_per_arm_change"]
    arms = [("PREDICT", "FL_public", "PREDICT FL"), ("PREDICT", "MG", "PREDICT MG"), ("CAPSCAN", "FL_public", "CapScan FL"),
            ("CAPSCAN", "MG", "CapScan MG"), ("HSTOOL", "FL", "Healthy FL"), ("CAMPI", "FL", "CAMPI FL")]
    mine = [(g.nested7_vs_base5.median(), g.fresh7_vs_base5.median()) for co, arm, _ in arms for g in [w[(w.cohort == co) & (w.arm == arm)]]]
    rec.append(("per_acquisition_length_counts.tsv", "c_per_arm_change", len(c_), len(arms),
                bool(np.allclose(sorted(c_.filtered_to_7.tolist() + c_.new_search_7.tolist()), sorted(v for m in mine for v in m))),
                "medians on the 331 roster; the package sheet is the 331 version, v28f CapScan FL uses 347 (adds participant15 table)"))
    mem = pd.read_csv(SD / "campi_member_counts.tsv", sep="\t")
    d_ = xl["d_member_counts"]
    mm = mem[mem.arm.isin(["FL", "MG", "REF"]) & mem.member.isin(["Lactobacillus plantarum", "Clostridium butyricum", "Bifidobacterium longum"])
             & mem.population.isin(["baseline5", "fresh7", "nested9"])]
    rec.append(("campi_member_counts.tsv", "d_member_counts", len(d_), len(mm),
                bool(sorted(d_.pooled_specific.astype(int)) == sorted(mm.pooled_specific.astype(int))), "sorted pooled_specific"))
    e_ = xl["e_fusion"]
    fu = pd.read_csv(SD / "campi_fusion_pooled_q01.tsv", sep="\t")
    rec.append(("campi_fusion_pooled_q01.tsv", "e_fusion", len(e_), len(fu),
                bool(np.allclose(sorted(e_.pooled_diagnostic_pct), sorted(fu.pooled_diag_pct))), "sorted pooled diagnostic (%)"))
    f_ = xl["f_healthy"]
    hs = pd.read_csv(SD / "healthy_strategies_per_acquisition.tsv", sep="\t")
    med = [hs[(hs.strategy == s) & (hs.population == p)].count_change_pct.median() for s in ["common", "joint"] for p in ["baseline5", "fresh7"]]
    rec.append(("healthy_strategies_per_acquisition.tsv", "f_healthy", len(f_), 4,
                bool(np.allclose(sorted(f_.median_count_change_pct), sorted(med))), "sorted median count change (%)"))
    g_ = xl["g_grouping"]
    gr = pd.read_csv(SD / "grouping_study_groups.tsv", sep="\t")
    gr["change_pct"] = 100 * (gr.fresh7 / gr.baseline - 1)
    cap = gr.task.eq("CAPSCAN_noclust")
    same = np.allclose(sorted(g_.change_pct[~g_.astype(str).apply(lambda r: r.str.contains("CAPSCAN").any(), axis=1)]),
                       sorted(gr.change_pct[~cap]))
    rec.append(("grouping_study_groups.tsv", "g_grouping", len(g_), len(gr), bool(same),
                "non-CapScan rows equal; CapScan row differs by design (package 331: 174,673 -> 175,282; v28f 347: 177,717 -> 178,357)"))
    out = pd.DataFrame(rec, columns=["source_data_file", "package_sheet", "package_rows", "compared_rows", "values_agree", "comparison"])
    (UNIT / "checks").mkdir(exist_ok=True)
    out.to_csv(UNIT / "checks" / "source_data_vs_package.tsv", sep="\t", index=False)
    print(out.to_string())

    # origins for PROVENANCE.md
    lines = []
    for name, (srcs, how) in origin.items():
        srcs = srcs if isinstance(srcs, list) else [srcs]
        for s in srcs:
            s = Path(s)
            h = sha(s) if s.is_file() else "multiple files"
            lines.append(dict(file=name, origin="LAB_ROOT/" + str(s).replace(str(LAB) + "/", ""), sha256=h, transformation=how))
    lines.append(dict(file="(package workbook)", origin="LAB_ROOT/" + str(PACKAGE_XLSX).replace(str(LAB) + "/", ""), sha256=sha(PACKAGE_XLSX), transformation="compared, not copied"))
    pd.DataFrame(lines).to_csv(Path(os.environ["SFIG04_WORK"]) / "origins.tsv", sep="\t", index=False)
    for p in sorted(SD.glob("*.tsv")):
        print(p.name, sha(p), sum(1 for _ in open(p)) - 1)


def stage_expected():
    """Copy the producer's composite and plotted values, and crop per-panel reference images."""
    sys.path[:0] = [str(REPO), str(UNIT)]
    import panels
    from PIL import Image
    exp = UNIT / "expected"
    exp.mkdir(exist_ok=True)
    shutil.copyfile(PRODUCER_OUT / "out/SupplementaryFigure4_R5_400dpi.png", exp / "SFig04_v28f.png")
    shutil.copyfile(PRODUCER_OUT / "out/SupplementaryFigure4_R5.pdf", exp / "SFig04_v28f.pdf")
    shutil.copyfile(PRODUCER_OUT / "plotted_values/SupplementaryFigure4_plotted_values.tsv", exp / "SFig04_plotted_values.tsv")
    boxes = panels.panel_boxes(SD)
    boxes.to_csv(exp / "PANEL_BOXES.tsv", sep="\t", index=False)
    im = Image.open(exp / "SFig04_v28f.png")
    for b in boxes.itertuples():
        im.crop((b.x0, b.y0, b.x1, b.y1)).save(exp / f"panel_{b.panel}.png")
    for p in sorted(exp.iterdir()):
        print(p.name, sha(p))


if __name__ == "__main__":
    {"source": stage_source, "expected": stage_expected}[sys.argv[1]]()
