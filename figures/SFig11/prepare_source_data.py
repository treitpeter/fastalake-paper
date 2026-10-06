"""Copy the tables Supplementary Figure 11 needs from the lab archive into source_data/SFig11/.

This script documents how every file in ``source_data/SFig11/`` was made. It is not run by
the notebook and it cannot run outside the lab, because it reads the original producer
inputs and outputs. Set ``FASTALAKE_LAB_ROOT`` to the lab project directory (``LAB_ROOT`` in
PROVENANCE.md) to rerun it.

Files written, and the only transformations applied:

* ``functional_correlations_PREDICT.tsv``, ``native_correlations_PREDICT.tsv``: the PREDICT
  rows of the producer input tables ``SF11/inputs_n15/{functional,native}_correlations.tsv``,
  copied as text lines (header plus PREDICT lines, byte for byte). The CapScan rows of those
  producer inputs are byte-identical to ``source_data/Fig5/functional_correlations.tsv`` and
  ``source_data/Fig5/native_correlations.tsv`` (asserted below), so they are not duplicated.
* ``controlled_176_comparisons.tsv``: the producer input ``SF11/inputs_n15/controlled_176_comparisons.tsv``,
  copied byte for byte.
* ``distances_PREDICT_functional.tsv``, ``distances_PREDICT_native.tsv``,
  ``distances_PREDICT_controlled.tsv``, ``distances_CAPSCAN_native_prevalence10.tsv``,
  ``distances_CAPSCAN_controlled.tsv``: upper triangles (i < j, row-major order) of the
  Euclidean distance matrices saved by the upstream analyses (``*.npz``), in long format,
  floats written with ``repr`` (exact round trip, asserted). Patients are identified only by
  their 1-based position in the analysis roster (``profile_i``, ``profile_j``); CapScan
  profiles by participant number and sampling group, as in ``source_data/Fig5``. Only the
  matrices of the comparison families shown in the figure are exported; the
  ``wuyts_recipe`` sensitivity is not part of the 176 controlled estimates and is omitted.
* ``predict_native_retained_signal.tsv``: the PREDICT ``primary_prevalence10`` DNA and RNA
  rows of ``Crossomics/native/SOURCE_DATA_retained_signal_by_profile.tsv``; the patient
  identifier column ``profile`` is replaced by the 1-based row position within each layer.

Every input sha256 is asserted, so a changed archive file stops the script.
"""

import csv
import hashlib
import os
from pathlib import Path

import numpy as np
import pandas as pd

LAB = Path(os.environ["FASTALAKE_LAB_ROOT"])
REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "source_data" / "SFig11"

BR = "review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/"
CR = BR + "claude_review_2026-10-03/"
SF11 = CR + "final/suppfigs_n15/SF11/"
V50 = CR + "iteration1/fix_eggnog_native/downstream/verified50/B/cross_refit/"
C2 = CR + "coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig5g/crossomics_347/"
A4 = CR + "final/reruns_347/A4_controlled_distance/run347/"
XN = "review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/"

#: input (relative to LAB_ROOT) -> sha256
INPUTS = {
    SF11 + "inputs_n15/functional_correlations.tsv": "882c36b28646498e5a7d1516b2af6c7f80397e95f3881f0fd1bcef169d11b6b7",
    SF11 + "inputs_n15/native_correlations.tsv": "55d6a8c2aa0f79df9b7d7bbd0a0e845184a64ac2cf8ffc3e5f53992903abd5ea",
    SF11 + "inputs_n15/controlled_176_comparisons.tsv": "1705248887b442a01f3c52b36c584b3107db38d5f8e5c53067c32081d506ffca",
    V50 + "functional_distances.npz": "a8e0d5a9f7f2aa52f9dfbfad87a6ac720934a5458bb6f341bbabc49bd53648d0",
    V50 + "native_distances.npz": "8928b6d572273e1763cf3422fcd5fa28fd381fbbfce0e92ab5bd9d06efa65e51",
    C2 + "native_distances.npz": "d7b0bb2573b8f984b1e026ca5f678f74c597384789d556ecb2d5bd69e5b4c0e3",
    C2 + "functional_distances.npz": "7b9ba5caa5b78dd767ccb5c8087053d047379fbf5ea3760beabad5e694cb0508",
    BR + "distance_refit/PREDICT_distances.npz": "f596d528271810b8ad604b48c5a1f1b2f218870ba3b76d83d1c8661c2ed4c7ee",
    A4 + "CAPSCAN_distances.npz": "60c011d1de9e7851d2100ca1e7bf2001bd55091dba68f96af2ea749edf2322f6",
    XN + "SOURCE_DATA_retained_signal_by_profile.tsv": "761a6b1242d432ae9f5b6fb1cf30221b5cd9b45a8cc492996aeaac57b8125524",
}

#: Shared tables written by the Fig5 unit that hold the CapScan rows of the producer inputs.
FIG5 = {
    "functional_correlations.tsv": "4c3468fc11199d32c9dd99fbda823bb536f624d01102df87562b6e29906de3f8",
    "native_correlations.tsv": "0004f35d522422580e1d7cce3bbcbc2d2967a8fbc8af28668b73408d52305375",
}

REG = ["duodenum", "jejunum", "ileum", "colon", "stool"]
CAPSCAN_DONORS = [1, 2, 3, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15]  # crossomics_347 and run347 COMPLETE.json
PREDICT_LAYERS_FUNCTIONAL = ["metaG", "metaT", "microbial_metaP", "host_metaP", "metabolites"]
PREDICT_LAYERS_NATIVE = ["DNA", "RNA", "Non-host protein", "Host protein", "Metabolites"]
CAPSCAN_LAYERS_NATIVE = ["DNA", "Non-host protein", "Host protein", "Metabolites"]
CONTROLLED_SETS = [(f, d) for f in ["all_observed", "prevalence10"] for d in ["hellinger", "rank"]]


def sha256(path):
    """Return the hex sha256 of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_inputs():
    """Assert the sha256 of every archive input and of the shared Fig5 tables."""
    for rel, h in INPUTS.items():
        assert sha256(LAB / rel) == h, rel
    for name, h in FIG5.items():
        assert sha256(REPO / "source_data" / "Fig5" / name) == h, name


def split_rows(name):
    """Write the PREDICT lines of a producer input table; assert its CapScan lines equal Fig5's file."""
    lines = (LAB / SF11 / "inputs_n15" / name).read_bytes().decode("utf-8").splitlines(keepends=True)
    head, body = lines[0], lines[1:]
    pred = [ln for ln in body if ln.startswith("PREDICT\t")]
    cap = [ln for ln in body if ln.startswith("CAPSCAN\t")]
    assert len(pred) + len(cap) == len(body) and body == pred + cap
    assert (REPO / "source_data" / "Fig5" / name).read_bytes().decode("utf-8") == head + "".join(cap)
    out = OUT / name.replace(".tsv", "_PREDICT.tsv")
    out.write_bytes((head + "".join(pred)).encode("utf-8"))
    return out


def write_long(path, header, rows):
    """Write rows with ``repr`` floats; return the number of rows."""
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        n = 0
        for r in rows:
            w.writerow([repr(float(x)) if isinstance(x, float) else x for x in r])
            n += 1
    return n


def predict_rows(mats):
    """Yield (key..., profile_i, profile_j, distance) for a dict {key tuple: square matrix}."""
    for key, m in mats.items():
        n = len(m)
        iu = np.triu_indices(n, 1)
        for a, b, v in zip(iu[0], iu[1], m[iu]):
            yield (*key, int(a) + 1, int(b) + 1, float(v))


def capscan_rows(mats):
    """Yield (key..., participant_i, region_i, participant_j, region_j, distance)."""
    prof = [(d, r) for d in CAPSCAN_DONORS for r in REG]
    for key, m in mats.items():
        assert len(m) == len(prof)
        iu = np.triu_indices(len(m), 1)
        for a, b, v in zip(iu[0], iu[1], m[iu]):
            yield (*key, prof[a][0], prof[a][1], prof[b][0], prof[b][1], float(v))


def roundtrip(path, keycols, mats, capscan):
    """Assert that the written long table reproduces every exported matrix exactly."""
    t = pd.read_csv(path, sep="\t", float_precision="round_trip")
    for key, m in mats.items():
        sel = np.ones(len(t), bool)
        for c, v in zip(keycols, key):
            sel &= t[c].astype(str).to_numpy() == str(v)
        iu = np.triu_indices(len(m), 1)
        assert np.array_equal(t.loc[sel, "distance"].to_numpy(), m[iu]), (path.name, key)


def main():
    """Write every source_data/SFig11 file and print its sha256 and row count."""
    check_inputs()
    OUT.mkdir(parents=True, exist_ok=True)
    written = [split_rows("functional_correlations.tsv"), split_rows("native_correlations.tsv")]
    (OUT / "controlled_176_comparisons.tsv").write_bytes((LAB / SF11 / "inputs_n15/controlled_176_comparisons.tsv").read_bytes())
    written.append(OUT / "controlled_176_comparisons.tsv")

    fz = np.load(LAB / V50 / "functional_distances.npz")
    mats = {(lev, lay): fz[f"PREDICT__displayed77__{lev}__0.0__n77__{lay}"]
            for lev in ["OG", "KO", "Pathway"] for lay in PREDICT_LAYERS_FUNCTIONAL}
    p = OUT / "distances_PREDICT_functional.tsv"
    write_long(p, ["level", "layer", "profile_i", "profile_j", "distance"], predict_rows(mats))
    roundtrip(p, ["level", "layer"], mats, False)
    written.append(p)

    nz = np.load(LAB / V50 / "native_distances.npz")
    mats = {(an, lay): nz[f"PREDICT__{an}__{lay}"]
            for an in ["primary_prevalence10", "sensitivity_all_features"] for lay in PREDICT_LAYERS_NATIVE}
    p = OUT / "distances_PREDICT_native.tsv"
    write_long(p, ["analysis", "layer", "profile_i", "profile_j", "distance"], predict_rows(mats))
    roundtrip(p, ["analysis", "layer"], mats, False)
    written.append(p)

    dz = np.load(LAB / BR / "distance_refit/PREDICT_distances.npz")
    mats = {(f, d, "unadjusted", lay): dz[f"{f}__{d}__unadjusted__{lay}"]
            for f, d in CONTROLLED_SETS for lay in PREDICT_LAYERS_NATIVE}
    p = OUT / "distances_PREDICT_controlled.tsv"
    write_long(p, ["features", "distance_kind", "adjustment", "layer", "profile_i", "profile_j", "distance"], predict_rows(mats))
    roundtrip(p, ["features", "distance_kind", "adjustment", "layer"], mats, False)
    written.append(p)

    cz = np.load(LAB / C2 / "native_distances.npz")
    mats = {("primary_prevalence10", lay): cz[f"CAPSCAN__primary_prevalence10__{lay}"] for lay in CAPSCAN_LAYERS_NATIVE}
    p = OUT / "distances_CAPSCAN_native_prevalence10.tsv"
    write_long(p, ["analysis", "layer", "participant_i", "region_i", "participant_j", "region_j", "distance"], capscan_rows(mats))
    roundtrip(p, ["analysis", "layer"], mats, True)
    written.append(p)
    # The all-observed CapScan native matrices are in source_data/Fig5/native_distances_long.tsv.
    f5 = pd.read_csv(REPO / "source_data/Fig5/native_distances_long.tsv", sep="\t", float_precision="round_trip")
    for lay in CAPSCAN_LAYERS_NATIVE:
        m = cz[f"CAPSCAN__sensitivity_all_features__{lay}"]
        got = f5[(f5.analysis == "sensitivity_all_features") & (f5.layer == lay)].distance.to_numpy()
        assert np.array_equal(got, m[np.triu_indices(len(m), 1)]), ("Fig5 native", lay)
    f5 = pd.read_csv(REPO / "source_data/Fig5/functional_distances_long.tsv", sep="\t", float_precision="round_trip")
    fc = np.load(LAB / C2 / "functional_distances.npz")
    for lev in ["OG", "KO", "Pathway"]:
        for lay in ["metaG", "microbial_metaP", "host_metaP", "metabolites"]:
            m = fc[f"CAPSCAN__matched12__{lev}__0.0__n65__{lay}"]
            got = f5[(f5.scope == "matched12") & (f5.prevalence == 0) & (f5.level == lev) & (f5.layer == lay)].distance.to_numpy()
            assert np.array_equal(got, m[np.triu_indices(len(m), 1)]), ("Fig5 functional", lev, lay)

    az = np.load(LAB / A4 / "CAPSCAN_distances.npz")
    mats = {(f, d, adj, lay): az[f"{f}__{d}__{adj}__{lay}"]
            for f, d in CONTROLLED_SETS for adj in ["unadjusted", "region_adjusted"] for lay in CAPSCAN_LAYERS_NATIVE}
    p = OUT / "distances_CAPSCAN_controlled.tsv"
    write_long(p, ["features", "distance_kind", "adjustment", "layer", "participant_i", "region_i", "participant_j", "region_j", "distance"], capscan_rows(mats))
    roundtrip(p, ["features", "distance_kind", "adjustment", "layer"], mats, True)
    written.append(p)

    r = pd.read_csv(LAB / XN / "SOURCE_DATA_retained_signal_by_profile.tsv", sep="\t", float_precision="round_trip")
    r = r[(r.cohort == "PREDICT") & (r.analysis == "primary_prevalence10") & r.layer.isin(["DNA", "RNA"])].copy()
    r["profile_index"] = r.groupby("layer").cumcount() + 1
    p = OUT / "predict_native_retained_signal.tsv"
    write_long(p, ["cohort", "analysis", "layer", "profile_index", "retained_signal_fraction"],
               r[["cohort", "analysis", "layer", "profile_index", "retained_signal_fraction"]].itertuples(index=False, name=None))
    written.append(p)

    for p in written:
        n = sum(1 for _ in open(p, encoding="utf-8")) - 1
        print(p.name, sha256(p), n)


if __name__ == "__main__":
    main()
