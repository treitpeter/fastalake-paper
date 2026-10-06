#!/usr/bin/env python3
"""Extract the Figure 4 source data from the lab analysis tree (run once, by the authors).

This script is not part of the notebook. It documents exactly how every table in
``source_data/Fig4/`` and every controlled file was cut from the original analysis
outputs. It needs read access to the lab project directory, given as the environment
variable ``FASTALAKE_LAB_ROOT`` (written ``LAB_ROOT`` in PROVENANCE.md). Paths below are
relative to it.

What it does
------------
* copies the upstream tables the notebook needs, keeping every value as stored;
* replaces PREDICT acquisition names by ``PREDICT_A###`` (sorted acquisition name) and
  PREDICT patient numbers by the pseudonyms of the published Source Data package
  (``PREDICT_P###``); CapScan participants keep their published numbers 1-15;
* drops clinical columns (antibiotic state, clinical group, raw patient numbers) from
  every public table; the two clinical inputs are written, unchanged, to the
  controlled directory given by ``--controlled-out``;
* optionally (``--verify-recovery``) re-derives the cross-fitted recovery from the
  accepted-peptide files of the cross-fit run and checks it against the stored table.

Usage (on Slurm)::

    FASTALAKE_LAB_ROOT=... python figures/Fig4/prepare_source_data.py \
        --out source_data/Fig4 --controlled-out <dir> [--verify-recovery]
"""

import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

#: Analysis roots, relative to LAB_ROOT.
MC = Path("review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01")
CR = MC / "biology_refit/claude_review_2026-10-03"
PKG = CR / "final/figure_sd/source_data/package/FastaLake_Source_Data"
V50 = CR / "iteration1/fix_eggnog_native/downstream/verified50/B"

#: Input files and the sha256 each must have (prefix is enough to fail on change).
INPUTS = {
    "crossfit": (CR / "iteration4/crossfit/run/PER_ACQUISITION.tsv", "f7616394a3ffd829"),
    "crossfit_summary": (CR / "iteration4/crossfit/run/SUMMARY.tsv", "abd33b20058b9cec"),
    "recovery": (CR / "final/review_0/recovery_per_acquisition.tsv", None),
    "bench0": (PKG / "Figure4/benchmark_0_per_acquisition.tsv", None),
    "bench1": (PKG / "Figure4/benchmark_1_per_acquisition.tsv", None),
    "audit": (PKG / "Figure4/scoring_counts_public_vs_specimen.tsv", "f005f7c000d0"),
    "clin_stats": (MC / "clinical_refit/clinical_paired_statistics.tsv", None),
    "clin_paired": (MC / "clinical_refit/paired_antibiotic_profiles.tsv", None),
    "comp_means": (CR / "iteration4/figures/shadow_d12/B/figures/clinical_composition_means.tsv", "d68bc857489e9954"),
    "hm_values": (V50 / "eds_refit/input_data/heatmap_85_values.npz", "11c171206d200b90"),
    "hm_features": (V50 / "eds_refit/input_data/heatmap_features.tsv", "b883e1659743f543"),
    "hm_order": (V50 / "eds_refit/input_data/heatmap_85_patient_order.tsv", "43dc72c4b31e4bcf"),
    "native_corr": (V50 / "cross_refit/native_correlations.tsv", None),
    "functional_corr": (V50 / "cross_refit/functional_correlations.tsv", None),
    "native_dist": (V50 / "cross_refit/native_distances.npz", None),
    "functional_dist": (V50 / "cross_refit/functional_distances.npz", None),
    "dna_yield": (PKG / "Figure4/capscan_assembly_dna_yield/capscan_assembly_dna_yield_gain_by_acquisition.tsv", None),
}

#: Native layers in the order the original analysis inserted them (sets the pair order).
NATIVE_LAYERS = ["DNA", "RNA", "Non-host protein", "Host protein", "Metabolites"]
#: Functional layers in the order the original analysis inserted them.
FUNCTIONAL_LAYERS = ["metaG", "metaT", "microbial_metaP", "host_metaP", "metabolites"]


def sha(path):
    """Return the hex sha256 of a file.

    Parameters
    ----------
    path : pathlib.Path
        File to hash.

    Returns
    -------
    str
        Hex digest.
    """
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write(df, path):
    """Write a tab-separated table without index, full float precision.

    Parameters
    ----------
    df : pandas.DataFrame
        Table to write.
    path : pathlib.Path
        Output file.

    Returns
    -------
    None
    """
    df.to_csv(path, sep="\t", index=False, float_format="%.17g")


def main():
    """Cut every Figure 4 table; see the module docstring.

    Parameters
    ----------
    None. Command-line arguments ``--out``, ``--controlled-out`` and
    ``--verify-recovery``; the lab root is read from FASTALAKE_LAB_ROOT.

    Returns
    -------
    None
    """
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--controlled-out", required=True)
    ap.add_argument("--verify-recovery", action="store_true")
    a = ap.parse_args()
    lab = Path(os.environ["FASTALAKE_LAB_ROOT"])
    out, ctl = Path(a.out), Path(a.controlled_out)
    out.mkdir(parents=True, exist_ok=True)
    ctl.mkdir(parents=True, exist_ok=True)
    src, origin = {}, {}
    for key, (rel, prefix) in INPUTS.items():
        p = lab / rel
        h = sha(p)
        assert prefix is None or h.startswith(prefix), (key, h)
        src[key], origin[key] = p, {"path": "LAB_ROOT/" + rel.as_posix(), "sha256": h}

    # ---------------- PREDICT pseudonyms ----------------
    bench = pd.concat([pd.read_csv(src["bench0"], sep="\t"), pd.read_csv(src["bench1"], sep="\t")], ignore_index=True)
    pred_all = sorted(bench.loc[bench.cohort.eq("PREDICT"), "sample"].unique())
    assert len(pred_all) == 306
    acq = {s: f"PREDICT_A{k + 1:03d}" for k, s in enumerate(pred_all)}
    subj = bench.drop_duplicates("sample").set_index("sample").subject_id

    def acq_id(cohort, sample):
        """Return the public acquisition id: the PREDICT_A### pseudonym for PREDICT, the CapScan name otherwise."""
        return acq[sample] if cohort == "PREDICT" else sample

    # ---------------- cross-fit per acquisition (Fig. 4b,c; text lines 123, 129, 133) ----------------
    x = pd.read_csv(src["crossfit"], sep="\t", dtype={"subject_id": str})
    rec = pd.read_csv(src["recovery"], sep="\t")
    assert len(x) == 578 and list(rec.i) == list(x.i)
    noid = sorted(s for s in x.subject_id.unique() if s.startswith("NOID_"))
    pid = []
    for r in x.itertuples():
        if r.cohort == "CAPSCAN":
            pid.append(str(r.subject_id))
        elif r.subject_id.startswith("NOID_"):
            pid.append(f"PREDICT_NOID_{noid.index(r.subject_id) + 1}")
        else:
            p = subj.loc[r.sample]
            assert isinstance(p, str) and p.startswith("PREDICT_P"), r.sample
            pid.append(p)
    # The bootstrap of the original analysis forms clusters in np.unique order of the raw
    # participant strings. participant_sort_rank keeps that order without the raw numbers.
    rank = {}
    for coh in ("PREDICT", "CAPSCAN"):
        u = np.unique(x.loc[x.cohort.eq(coh), "subject_id"].to_numpy())
        rank.update({(coh, s): k for k, s in enumerate(u)})
    cf = pd.DataFrame({
        "crossfit_index": x.i,
        "cohort": x.cohort,
        "acquisition_id": [acq_id(c, s) for c, s in zip(x.cohort, x["sample"])],
        "participant_id": pid,
        "participant_sort_rank": [rank[(c, s)] for c, s in zip(x.cohort, x.subject_id)],
        "sampling_site": np.where(x.cohort.eq("PREDICT"), "stool", x.group),
        "lda_fallback": x.lda_fallback.fillna(""),
    })
    for L in ("5", "9"):
        for col in ("FL_pooled", "MG_pooled", "frozen_fl", "frozen_mg"):
            cf[f"{col}_{L}"] = x[f"{col}_{L}"]
    cf["recovery_cf_5_pct"] = rec.rec_cf_5
    cf["recovery_cf_9_pct"] = rec.rec_cf_9
    write(cf, out / "crossfit_per_acquisition.tsv")

    # ---------------- same-spectrum benchmark rows (FL public vs MG) ----------------
    b = bench[bench.test_arm.eq("FL_public") & bench.reference_arm.eq("MG")].copy()
    assert b.groupby(["cohort", "min_length"]).size().to_dict() == {
        ("CAPSCAN", 5): 282, ("CAPSCAN", 9): 282, ("PREDICT", 5): 296, ("PREDICT", 9): 296}
    bb = pd.DataFrame({
        "cohort": b.cohort, "acquisition_id": [acq_id(c, s) for c, s in zip(b.cohort, b["sample"])],
        "participant_id": b.subject_id.fillna(""),
        "sampling_site": np.where(b.cohort.eq("PREDICT"), "stool", b.group),
        "min_length": b.min_length, "n_fl": b.n_test, "n_mg": b.n_reference, "overlap": b.overlap,
        "count_ratio": b.count_ratio, "mg_recovery_fraction": b.reference_recovery})
    write(bb, out / "same_spectrum_benchmark.tsv")

    # ---------------- PREDICT acquisition roster (306; 296 with a matched metagenome) ----------------
    matched = set(b.loc[b.cohort.eq("PREDICT"), "sample"])
    write(pd.DataFrame({"acquisition_id": [acq[s] for s in pred_all],
                        "participant_id": [subj.get(s, "") if isinstance(subj.get(s, ""), str) else "" for s in pred_all],
                        "matched_metagenome": [s in matched for s in pred_all]}),
          out / "predict_acquisition_roster.tsv")

    # ---------------- spectrum audit (Fig. 4d) ----------------
    s = pd.read_csv(src["audit"], sep="\t")
    assert len(s) == 296 and set(s["sample"]) == matched
    s.insert(0, "acquisition_id", [acq[v] for v in s["sample"]])
    write(s.drop(columns="sample").sort_values("acquisition_id"), out / "spectrum_audit_per_acquisition.tsv")

    # ---------------- clinical aggregates (public) ----------------
    shutil.copyfile(src["clin_stats"], out / "clinical_paired_statistics.tsv")
    m = pd.read_csv(src["comp_means"], sep="\t").rename(columns={"antibiotic_state": "condition"})
    write(m, out / "clinical_composition_means.tsv")

    # ---------------- heatmap (Fig. 4h): public values, controlled antibiotic track ----------------
    z = np.load(src["hm_values"], allow_pickle=True)
    feat = pd.read_csv(src["hm_features"], sep="\t")
    order = pd.read_csv(src["hm_order"], sep="\t")
    assert (order.sample_id.to_numpy() == z["samples"].astype(str)).all()
    assert (z["features"].astype(str) == feat.feature_id.to_numpy()).all()
    write(feat[["block", "feature_id", "label", "n_positive_85"]], out / "heatmap_features.tsv")
    zz = z["z"].astype(float)
    write(pd.DataFrame([dict(feature_id=f, display_position=j + 1, z=zz[i, j])
                        for i, f in enumerate(feat.feature_id) for j in range(zz.shape[1])]),
          out / "heatmap_z.tsv")
    write(pd.DataFrame({"display_position": np.arange(1, len(order) + 1), "eds_score": order.dysbiosis_score,
                        "human_fraction": order.Human}), out / "heatmap_tracks.tsv")
    shutil.copyfile(src["hm_order"], ctl / "heatmap_85_patient_order.tsv")
    shutil.copyfile(src["clin_paired"], ctl / "paired_antibiotic_profiles.tsv")

    # ---------------- cross-omics (Fig. 4i) ----------------
    shutil.copyfile(src["native_corr"], out / "crossomics_native_correlations.tsv")
    shutil.copyfile(src["functional_corr"], out / "crossomics_functional_correlations.tsv")
    nd, fd = np.load(src["native_dist"]), np.load(src["functional_dist"])
    cols = {f"native__{k}": nd[f"PREDICT__sensitivity_all_features__{k}"] for k in NATIVE_LAYERS}
    for lev in ("OG", "KO", "Pathway"):
        for k in FUNCTIONAL_LAYERS:
            cols[f"{lev}__{k}"] = fd[f"PREDICT__displayed77__{lev}__0.0__n77__{k}"]
    n = 77
    ia, ib = np.triu_indices(n, 1)
    d = pd.DataFrame({"patient_a": ia + 1, "patient_b": ib + 1})
    for k, v in cols.items():
        assert v.shape == (n, n) and np.allclose(v, v.T)
        d[k] = v[ia, ib]
    write(d, out / "crossomics_distances_predict77.tsv")

    # ---------------- CapScan DNA yield and database size (text line 129) ----------------
    shutil.copyfile(src["dna_yield"], out / "capscan_dna_yield_by_acquisition.tsv")

    if a.verify_recovery:
        base = lab / CR / "iteration4/crossfit/outputs"
        bad = 0
        for r, rr in zip(x.itertuples(), cf.itertuples()):
            p = base / r.cohort / r.sample
            S = {k: {t for t in (p / f"{k}_POOLED_ACCEPTED.txt").read_text().split() if 9 <= len(t) <= 50} for k in ("FL", "MG")}
            v = 100 * len(S["FL"] & S["MG"]) / len(S["MG"])
            bad += not np.isclose(v, rr.recovery_cf_9_pct, rtol=0, atol=1e-9)
        print(json.dumps({"recovery_rederived_from_accepted_files": 578, "mismatches": bad}))

    receipt = {"inputs": origin, "outputs": {p.name: sha(p) for p in sorted(out.glob("*.tsv"))},
               "controlled": {p.name: sha(p) for p in sorted(ctl.glob("*.tsv"))}}
    (out / "EXTRACTION_RECEIPT.json").write_text(json.dumps(receipt, indent=1) + "\n")
    print(json.dumps({k: len(v) for k, v in receipt.items()}))


if __name__ == "__main__":
    main()
