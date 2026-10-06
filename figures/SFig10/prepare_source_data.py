"""Write source_data/SFig10/ from the original producer inputs (run once, by the authors).

This script documents how every table in ``source_data/SFig10/`` was derived. It is not
run by the notebook and needs the lab tree, named by the environment variable
``FASTALAKE_LAB_ROOT``.

Inputs (relative to LAB_ROOT, sha256 checked):

* ``EDS`` = the EDS_approved/source_data directory that the original Supplementary
  Fig. 10 builder read (see PROVENANCE.md).
* ``PKG`` = the Source Data package folder EDS_focused_20261003/source_data, used only
  for the pseudonymised patient identifiers (PREDICT_Pnnn).

Transformations:

* Per-patient table: display position 1-85 (row order of the producer's patient table,
  which is the EDS order), pseudonym, EDS, human-assigned LFQ share, protein-group
  count and the total KO signal per patient (sum over all KO columns of the KO profile
  the upstream analysis read, computed exactly as build_approved_eds_cr003.py does for
  the annotated-relative normalization). Clinical group, antibiotic state, raw patient identifiers and acquisition or
  microbiome sample identifiers are withheld (controlled input, CONVENTIONS section 8).
* NumPy archives are written as long tab-separated tables keyed by display position (85
  patients) or observation index (259 observations); sample identifiers are dropped.
  Floats are written with round-trip precision and re-read (pandas float_precision=
  "round_trip") to assert exact equality.
* Screen tables (augmented pair screen, EDS associations), feature tables and the
  producer's summary and adjustment tables are copied unchanged.
* ``producer_plotted_ab.tsv``: the header and the heatmap rows (panel ``a`` and ``b``) of
  the Source Data sheet Supplementary_Figure_10/Supplementary_Figure10_plotted_values.tsv,
  copied byte for byte. The ``a_track`` rows (EDS, human share and per-patient antibiotic
  state) are withheld.

Usage
-----
    FASTALAKE_LAB_ROOT=<lab root> python figures/SFig10/prepare_source_data.py
"""

import hashlib
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "source_data" / "SFig10"
LAB = Path(os.environ["FASTALAKE_LAB_ROOT"])
BASE = ("review/v15_integration/jess_integration_2026-09-29/analysis_followup/"
        "manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/")
EDS = LAB / BASE / "iteration2/cr003_followup/shadow/verified50/eds/EDS_approved/source_data"
KO = LAB / BASE / "iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/input_data/KO_profile.tsv.gz"
PKG = LAB / BASE / "final/figure_sd/source_data/package/FastaLake_Source_Data/EDS_focused_20261003/source_data"
PLOTTED = LAB / BASE / "final/figure_sd/source_data/package/FastaLake_Source_Data/Supplementary_Figure_10/Supplementary_Figure10_plotted_values.tsv"

#: sha256 of every input, as recorded in the producer's BUILD_VALIDATION.json.
INPUT_SHA = {
    EDS / "patient_order_and_covariates.tsv": "6d8b883529d90dbc682898c2025f01f62e781ba9a38c9ac072aac3a64ddbadc9",
    EDS / "focused_values.npz": "e5900ff12c7f7e540f31a0bb1796e0be1bb6341216c436e9039145dd752ddc31",
    EDS / "focused_features.tsv": "7945c40a934016cbfab1952c42b02cdf2854fb923b3815b166b308a81e352aae",
    EDS / "legacy_85_values.npz": "11c171206d200b9051c0f785041288bbc621be61695baae1bfa8eaffb88b7b37",
    EDS / "legacy_features.tsv": "b883e1659743f5437c3d1d3981abd07a36d1ae41ee5d87b6321596fe90968e51",
    EDS / "legacy_259_heatmap_values.npz": "b937bb50fc4c870173cd85638041805e8c5700ded180616be9cf91ef22c2f042",
    EDS / "augmented_pair_screen.tsv": "d7e59ab759d62aa5b2068ac5fbe0b774ece288ad3894bd633cfd219cac2a3967",
    EDS / "all_EDS_associations.tsv": "dcd15919dde69cce09d3fdd3e96c5517e35f1786b824bb5575ca7dc0131008dd",
    EDS / "selected_pair_summary.tsv": "42e89ff32dc0a4d3977feec112f3e1e827cb5940d82df4062d45df7f9cb8e8c1",
    EDS / "separate_adjustments_all_models.tsv": "fe72b674538a021d48f14fc66c50d8aaa9e454960ca636969276b99a2133917b",
    KO: "164fe6cde20cec1e6cc6803c7b2ad50a3f373f02e411e211769670b0f39cec02",
    PKG / "patient_order_and_covariates.tsv": "8a4844354b3dfcdc8087e19def097e25dcfebfc59fabf1f8bfbad918f06baa68",
    PLOTTED: "e08b1d79a049d1f0049c0d60ab2b5eea095723d0f518d14fb4232daac8c6755b",
}


def sha(p):
    """Return the hex sha256 of a file.

    Parameters
    ----------
    p : str or pathlib.Path
        File.

    Returns
    -------
    str
        Hex digest.
    """
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def wr(df, name):
    """Write ``df`` as TSV and assert that it re-reads to identical values.

    Parameters
    ----------
    df : pandas.DataFrame
        Table to write.
    name : str
        File name inside ``source_data/SFig10``.

    Returns
    -------
    None
    """
    df.to_csv(OUT / name, sep="\t", index=False)
    back = pd.read_csv(OUT / name, sep="\t", float_precision="round_trip")
    for c in df.columns:
        if pd.api.types.is_float_dtype(df[c]):
            assert np.array_equal(df[c].to_numpy(), back[c].to_numpy(), equal_nan=True), (name, c)
        else:
            assert df[c].astype(str).equals(back[c].astype(str)), (name, c)


def long(npz, key_name, keys, with_z=True):
    """Flatten a feature-by-column NumPy archive into a long table.

    Parameters
    ----------
    npz : numpy.lib.npyio.NpzFile
        Archive with ``raw`` (and ``z``) arrays of shape features x columns, and
        optionally ``features``.
    key_name : str
        Name of the column key, ``patient_position`` or ``observation_index``.
    keys : array_like
        Key value of each column.
    with_z : bool
        Also write the ``z`` array.

    Returns
    -------
    pandas.DataFrame
        One row per feature and column: feature, key, raw (and z).
    """
    rows = []
    for i, f in enumerate(npz["features"] if "features" in npz.files else range(npz["raw"].shape[0])):
        for j, k in enumerate(keys):
            r = {"feature": str(f), key_name: int(k), "raw": float(npz["raw"][i, j])}
            if with_z:
                r["z"] = float(npz["z"][i, j])
            rows.append(r)
    return pd.DataFrame(rows)


def heatmap_rows(src, dst):
    """Copy the header and the panel a and b rows of the plotted-values sheet, byte for byte.

    Parameters
    ----------
    src : pathlib.Path
        Source Data plotted-values sheet.
    dst : pathlib.Path
        Output file.

    Returns
    -------
    int
        Number of data rows written.
    """
    lines = src.read_bytes().split(b"\n")
    keep = [lines[0]] + [x for x in lines[1:] if x.split(b"\t", 1)[0] in (b"a", b"b")]
    dst.write_bytes(b"\n".join(keep) + b"\n")
    return len(keep) - 1


def main():
    """Derive every source_data/SFig10 table from the producer inputs.

    Returns
    -------
    None
    """
    for p, h in INPUT_SHA.items():
        assert sha(p) == h, p
    OUT.mkdir(parents=True, exist_ok=True)

    order = pd.read_csv(EDS / "patient_order_and_covariates.tsv", sep="\t")
    pseud = pd.read_csv(PKG / "patient_order_and_covariates.tsv", sep="\t")
    assert len(order) == 85 and (order.sample_id.to_numpy() == pseud.sample_id.to_numpy()).all()
    assert order.drop(columns="patient_id").equals(pseud.drop(columns="patient_id"))
    assert np.all(np.diff(order.dysbiosis_score.to_numpy()) >= 0)  # rows are in EDS order
    pos = np.arange(1, 86)
    ko = pd.read_csv(KO, sep="\t", index_col=0).drop(index="FunctionallyUnassigned", errors="ignore").loc[:, order.sample_id].T
    wr(pd.DataFrame({"patient_position": pos, "patient_pseudonym": pseud.patient_id,
                     "dysbiosis_score": order.dysbiosis_score, "Human": order.Human,
                     "pg_parsimony": order.pg_parsimony,
                     "ko_total_signal": ko.sum(axis=1).to_numpy(float)}), "patients.tsv")

    fv = np.load(EDS / "focused_values.npz", allow_pickle=False)
    assert list(fv["samples"]) == list(order.sample_id)
    wr(long(fv, "patient_position", pos), "focused_values.tsv")
    lv = np.load(EDS / "legacy_85_values.npz", allow_pickle=False)
    assert list(lv["samples"]) == list(order.sample_id)
    wr(long(lv, "patient_position", pos), "legacy_85_values.tsv")
    l259 = np.load(EDS / "legacy_259_heatmap_values.npz", allow_pickle=False)
    feats = pd.read_csv(EDS / "legacy_features.tsv", sep="\t").feature_id
    rows = [{"feature": f, "observation_index": j + 1, "raw": float(l259["raw"][i, j])}
            for i, f in enumerate(feats) for j in range(l259["raw"].shape[1])]
    wr(pd.DataFrame(rows), "legacy_259_values.tsv")

    for name in ["focused_features.tsv", "legacy_features.tsv", "augmented_pair_screen.tsv",
                 "all_EDS_associations.tsv", "selected_pair_summary.tsv",
                 "separate_adjustments_all_models.tsv"]:
        shutil.copy2(EDS / name, OUT / name)
    assert heatmap_rows(PLOTTED, OUT / "producer_plotted_ab.tsv") == 25 * 85
    for p in sorted(OUT.glob("*.tsv")):
        print(p.name, sha(p))


if __name__ == "__main__":
    main()
