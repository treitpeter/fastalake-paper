"""Supplementary Table 4: Enterotype Dysbiosis Score (EDS) displays in the PREDICT cohort.

Regenerates the eleven sheets of the workbook and recomputes the statistics it reports.

Sheets and the source_data files they are built from:

* Read_me: ``source_data/STab4/read_me.tsv``
* All20_features_259_observations: ``source_data/SFig10/legacy_259_values.tsv``,
  ``source_data/SFig10/legacy_features.tsv``
* Fixed_feature_key: ``source_data/SFig10/legacy_features.tsv``
* Observation_metadata: ``source_data/STab4/observations_259.tsv``
* 85_patient_order: ``source_data/SFig10/patients.tsv``,
  ``source_data/STab4/patients_85_observation_index.tsv``
* 85_patient_raw, 85_patient_z: ``source_data/SFig10/legacy_85_values.tsv``
* Focused_raw: ``source_data/SFig10/focused_values.tsv`` and the two patient tables
* Focused_features, Focused_pairs, Focused_adjustments:
  ``source_data/SFig10/focused_features.tsv``, ``selected_pair_summary.tsv``,
  ``separate_adjustments_all_models.tsv``

Statistics read in addition ``source_data/SFig10/augmented_pair_screen.tsv``,
``source_data/SFig10/all_EDS_associations.tsv`` and
``source_data/STab4/five_pair_all_models.tsv``.

The repository version withholds acquisition and microbiome sample identifiers, clinical
group and antibiotic state (CONVENTIONS section 8). Sample identifiers are replaced by
observation index (1-259) or patient display position (1-85, EDS order). With
``FASTALAKE_CONTROLLED_DIR`` set, the complete published layout is rebuilt and the
statistics that need antibiotic state and clinical group are recomputed.

Statistical logic is ported from the original producers: ``corr`` and ``bh`` from the
EDS screen (analyze.py), ``estimate`` and ``blocks`` from build_approved_eds_cr003.py.
"""

import hashlib
import itertools
import os
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
from scipy.stats import rankdata, spearmanr, t

N_FEATURES_ORIGINAL = 20        # v28d line 481: earlier 20-feature selection
N_OBSERVATIONS = 259            # v28d line 481: 259 matched observations incl. repeated visits
N_PATIENTS = 85                 # v28d lines 137, 481, 799: earliest matched observation per patient
MIN_POSITIVE_N = 43             # v28d line 799: at least 43 patients with positive values for a test
Q_CUTOFF = 0.05                 # v28d line 799: q<0.05
PAIRS = [                       # v28d line 799: the three post hoc pairs (metabolite, KO, labels, reaction)
    ("HILIC_FT23003", "K01812", "Glucuronic acid", "Glucuronate isomerase", "R01482"),
    ("HILIC_FT24070", "K01667", "Tryptophan", "Tryptophanase", "R00673"),
    ("C18_F2423", "K01442", "Glycochenodeoxycholic acid", "Bile salt hydrolase", "R03975"),
]
COVARIATES = ["host", "depth", "antibiotics", "clinical"]   # v28d line 799: 16 combinations
NORMALIZATIONS = ["native", "annotated_relative"]          # v28d line 799
SCOPES = ["pair_complete", "common_84"]                    # v28d line 799: pair-specific and common subsets
JOINT_MODEL = "host_depth_antibiotic_clinical"             # v28d lines 137, 799: joint adjustment
DEPTH_TRANSFORM = "log2"        # v28d line 799: log2 protein-group depth, ranked
FIVE_PAIR_BOOTSTRAP_DRAWS = 1000        # v28d line 799; producer run_analysis.py bootstrap_n
FIVE_PAIR_SEED_BASE = 2026093032        # producer run_analysis.py: default_rng(2026093032 + pair index); not stated in v28d
FIVE_PAIR_INTERVAL = (0.025, 0.975)     # v28d line 799: pointwise 95% percentile intervals

#: Published columns withheld from the repository (controlled input).
WITHHELD = ["sample_id", "microbiome_id", "clinical_group", "antibiotic_state"]
#: Controlled files (controlled/CONTROLLED_INPUTS.md): pseudonymised Source Data copies.
CONTROLLED = {
    "patient_order_and_covariates.tsv": "8a4844354b3dfcdc8087e19def097e25dcfebfc59fabf1f8bfbad918f06baa68",
    "legacy_259_heatmap_columns.tsv": "a5c5bf922943e00937d77bf3fbd8b5cf2973d68fdb1d846cc2ceeceecc134740",
    "Supplementary_Table_4.xlsx": "4fe89a55c7607b0e5932377cbd488c031bc1f0c7f99c255f6c9eb7d77747e403",
}
SHEETS = ["Read_me", "All20_features_259_observations", "Fixed_feature_key", "Observation_metadata",
          "85_patient_order", "85_patient_raw", "85_patient_z", "Focused_raw", "Focused_features",
          "Focused_pairs", "Focused_adjustments"]
FREEZE = {"Focused_raw": "A2", "Focused_features": "A2", "Focused_pairs": "A2", "Focused_adjustments": "A2"}
WORDS = {3: "three", 5: "five", 6: "six"}


def _read(path, producer_parse=False):
    """Read a source_data table.

    ``producer_parse=True`` uses pandas' default float parser, as the original producer
    did when it read its own TSV files into the workbook (tables_cr003.py ``rd``); that
    parser can differ from the exact value in the last binary digit, and the published
    cells carry its result. All statistics use the exact (round-trip) values.
    """
    if producer_parse:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path, sep="\t", float_precision="round_trip")


def _matrix(long, features, key, n, col="raw"):
    """Feature-by-position matrix from a long table (rows in ``features`` order)."""
    m = long.pivot(index="feature", columns=key, values=col).reindex(index=features, columns=range(1, n + 1))
    return m.to_numpy(float)


def load(repo):
    """Read every source_data table the unit needs.

    Parameters
    ----------
    repo : pathlib.Path
        Repository root.

    Returns
    -------
    dict
        Tables (pandas) and feature-by-patient matrices (numpy).
    """
    s10, s4 = repo / "source_data" / "SFig10", repo / "source_data" / "STab4"
    d = {k: _read(s10 / f"{k}.tsv") for k in ["patients", "focused_features", "legacy_features", "selected_pair_summary",
                                               "separate_adjustments_all_models", "augmented_pair_screen", "all_EDS_associations"]}
    # workbook copies of the tables the producer read from TSV (see _read)
    d["sheet"] = {k: _read(s10 / f"{k}.tsv", producer_parse=True)
                  for k in ["legacy_features", "focused_features", "selected_pair_summary", "separate_adjustments_all_models"]}
    d["observations"] = _read(s4 / "observations_259.tsv").fillna({"patient_id": ""})
    d["link"] = _read(s4 / "patients_85_observation_index.tsv")
    d["read_me"] = pd.read_csv(s4 / "read_me.tsv", sep="\t", keep_default_na=False)
    d["five"] = _read(s4 / "five_pair_all_models.tsv")
    fv, lv, l259 = (_read(s10 / f) for f in ["focused_values.tsv", "legacy_85_values.tsv", "legacy_259_values.tsv"])
    ff, lf = list(d["focused_features"].feature), list(d["legacy_features"].feature_id)
    d["focused_raw"] = _matrix(fv, ff, "patient_position", N_PATIENTS)
    d["focused_raw_sheet"] = _matrix(_read(s10 / "focused_values.tsv", producer_parse=True), ff, "patient_position", N_PATIENTS)
    d["focused_z"] = _matrix(fv, ff, "patient_position", N_PATIENTS, "z")
    d["legacy_raw"] = _matrix(lv, lf, "patient_position", N_PATIENTS)
    d["legacy_z"] = _matrix(lv, lf, "patient_position", N_PATIENTS, "z")
    d["legacy_259"] = _matrix(l259, lf, "observation_index", N_OBSERVATIONS)
    assert (d["link"].patient_position.to_numpy() == d["patients"].patient_position.to_numpy()).all()
    return d


def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_controlled(d):
    """Read the controlled covariate files and the published workbook, if available.

    Returns
    -------
    dict or None
        ``pat`` (85 rows, display order), ``obs`` (259 rows), ``published`` (path); None
        with a printed message when ``FASTALAKE_CONTROLLED_DIR`` or a file is missing.
    """
    root = os.environ.get("FASTALAKE_CONTROLLED_DIR")
    paths = {k: Path(root) / k for k in CONTROLLED} if root else {}
    missing = [k for k in CONTROLLED if not paths or not paths[k].is_file()]
    if missing:
        print("SKIPPED: needs controlled input " + ", ".join(missing) + " (see controlled/CONTROLLED_INPUTS.md)")
        return None
    for k, h in CONTROLLED.items():
        assert _sha(paths[k]) == h, f"controlled input {k} has an unexpected sha256"
    pat = _read(paths["patient_order_and_covariates.tsv"])
    obs = _read(paths["legacy_259_heatmap_columns.tsv"])
    p, o = d["patients"], d["observations"]
    assert len(pat) == N_PATIENTS and len(obs) == N_OBSERVATIONS
    for a, b in [(pat, p), (obs, o)]:   # rows must be the same patients / observations in the same order
        assert np.array_equal(a.Human.to_numpy(), b.Human.to_numpy())
        assert np.array_equal(a.pg_parsimony.to_numpy(), b.pg_parsimony.to_numpy())
        assert np.array_equal(a.dysbiosis_score.to_numpy(), b.dysbiosis_score.to_numpy())
    assert (pat["index"].to_numpy() + 1 == d["link"].observation_index.to_numpy()).all()
    return {"pat": pat, "obs": obs, "published": paths["Supplementary_Table_4.xlsx"]}


# ---- workbook --------------------------------------------------------------------------
def _cell(v):
    """Python value for an xlsx cell; NaN becomes an empty cell."""
    if v is None:
        return None
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return None if np.isnan(v) else float(v)
    if isinstance(v, (np.bool_,)):
        return bool(v)
    return v


def _frame_rows(df):
    return [list(df.columns)] + [[_cell(v) for v in r] for r in df.itertuples(index=False, name=None)]


def _matrix_rows(features, ids, values):
    return [["feature_id"] + list(ids)] + [[f] + [_cell(v) for v in row] for f, row in zip(features, values)]


def build_sheets(d, ctrl=None):
    """Return {sheet: list of rows} for Supplementary Table 4.

    Parameters
    ----------
    d : dict
        Output of ``load``.
    ctrl : dict or None
        Output of ``load_controlled``. None gives the repository layout (withheld columns
        removed, sample identifiers replaced by observation index or patient position);
        otherwise the complete published layout.

    Returns
    -------
    dict
        Sheet name to list of rows (first row is the header, except Read_me).
    """
    p, o, link = d["patients"], d["observations"], d["link"]
    lf = list(d["legacy_features"].feature_id)
    pid = p.patient_pseudonym.to_numpy()
    if ctrl is None:
        obs_ids, pat_ids = list(o.observation_index), list(p.patient_position)
        meta = o.assign(patient_id=o.patient_id.where(o.patient_id != "", None))[["observation_index", "patient_id", "Human", "pg_parsimony", "dysbiosis_score"]]
        order = pd.DataFrame({"patient_position": p.patient_position, "observation_index": link.observation_index,
                              "patient_id": pid, "Human": p.Human, "pg_parsimony": p.pg_parsimony,
                              "dysbiosis_score": p.dysbiosis_score})
    else:
        c, co = ctrl["pat"], ctrl["obs"]
        obs_ids, pat_ids = list(co.sample_id), list(c.sample_id)
        meta = co.copy()
        order = c[["microbiome_id", "index", "sample_id", "patient_id", "clinical_group", "antibiotic_state",
                   "Human", "pg_parsimony", "dysbiosis_score"]].copy()
    focused = order.copy()
    if ctrl is not None:   # published Focused_raw column order
        focused = c[["sample_id", "microbiome_id", "index", "patient_id", "clinical_group", "antibiotic_state",
                     "Human", "pg_parsimony", "dysbiosis_score"]].copy()
    for f, v in zip(d["focused_features"].feature, d["focused_raw_sheet"]):
        focused[f] = v
    sheets = {
        "Read_me": [[t if t != "" else None] for t in d["read_me"].text],
        "All20_features_259_observations": _matrix_rows(lf, obs_ids, d["legacy_259"]),
        "Fixed_feature_key": _frame_rows(d["sheet"]["legacy_features"]),
        "Observation_metadata": _frame_rows(meta),
        "85_patient_order": _frame_rows(order),
        "85_patient_raw": _matrix_rows(lf, pat_ids, d["legacy_raw"]),
        "85_patient_z": _matrix_rows(lf, pat_ids, d["legacy_z"]),
        "Focused_raw": _frame_rows(focused),
        "Focused_features": _frame_rows(d["sheet"]["focused_features"]),
        "Focused_pairs": _frame_rows(d["sheet"]["selected_pair_summary"]),
        "Focused_adjustments": _frame_rows(d["sheet"]["separate_adjustments_all_models"]),
    }
    return {k: sheets[k] for k in SHEETS}


def write_workbook(sheets, path):
    """Write ``sheets`` to an xlsx file (one sheet per entry, freeze panes as published)."""
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    for name, rows in sheets.items():
        sh = wb.create_sheet(name)
        for r in rows:
            sh.append(r)
        if name in FREEZE:
            sh.freeze_panes = FREEZE[name]
    wb.save(path)
    return path


def read_workbook(path):
    """Return {sheet: list of row tuples} of cached cell values."""
    b = openpyxl.load_workbook(path, data_only=True, read_only=True)
    out = {s.title: [tuple(r) for r in s.iter_rows(values_only=True)] for s in b.worksheets}
    b.close()
    return out


def _same(x, y):
    if x == y:
        return True
    if x is None or y is None:
        return False
    return isinstance(x, float) and isinstance(y, float) and np.isnan(x) and np.isnan(y)


def compare_workbooks(built, reference):
    """Compare two workbooks cell for cell (exact equality; empty equals empty).

    Returns
    -------
    pandas.DataFrame
        One row per sheet: rows, columns, cells compared, cells differing, first difference.
    """
    a, b = read_workbook(built), read_workbook(reference)
    assert list(a) == list(b), (list(a), list(b))
    rows = []
    for name in a:
        ra, rb = a[name], b[name]
        ncol = max(max(len(r) for r in ra), max(len(r) for r in rb))
        diff, first, cells = 0, "", 0
        for i in range(max(len(ra), len(rb))):
            x = list(ra[i]) if i < len(ra) else []
            y = list(rb[i]) if i < len(rb) else []
            x += [None] * (ncol - len(x)); y += [None] * (ncol - len(y))
            for j in range(ncol):
                cells += 1
                if not _same(x[j], y[j]):
                    diff += 1
                    first = first or f"row {i + 1} col {j + 1}: {x[j]!r} vs {y[j]!r}"
        rows.append(dict(sheet=name, rows=len(ra), reference_rows=len(rb), columns=ncol, cells=cells,
                         differing=diff, first_difference=first))
    return pd.DataFrame(rows)


# ---- statistics ------------------------------------------------------------------------
def bh(p):
    """Benjamini-Hochberg q values (``bh`` of the EDS screen, analyze.py); NaN stays NaN."""
    p = np.asarray(p, float)
    v = np.full(len(p), np.nan)
    idx = np.where(np.isfinite(p))[0]
    s = idx[np.argsort(p[idx])]
    v[s] = np.minimum.accumulate((p[s] * len(s) / np.arange(1, len(s) + 1))[::-1])[::-1].clip(0, 1)
    return v


def corr(x, y, c=None):
    """Spearman or partial rank correlation with a two-sided t test (``corr`` of analyze.py).

    Ranks on the complete subset, regresses both ranks on [1, covariates], correlates the
    residuals; df = n - 2 - k with k the covariate rank.

    Returns
    -------
    tuple
        (n, rho, P, df, t statistic); NaN when n < ``MIN_POSITIVE_N``.
    """
    ok = np.isfinite(x) & np.isfinite(y)
    if c is not None:
        ok &= np.isfinite(c).all(axis=1)
    n = int(ok.sum())
    if n < MIN_POSITIVE_N or np.unique(x[ok]).size < 3 or np.unique(y[ok]).size < 3:
        return n, np.nan, np.nan, np.nan, np.nan
    a, b, k = rankdata(x[ok]), rankdata(y[ok]), 0
    if c is not None:
        C = np.column_stack([np.ones(n), c[ok]])
        k = np.linalg.matrix_rank(C) - 1
        a = a - C @ np.linalg.lstsq(C, a, rcond=None)[0]
        b = b - C @ np.linalg.lstsq(C, b, rcond=None)[0]
    r = np.corrcoef(a, b)[0, 1]
    df = n - k - 2
    tt = r * np.sqrt(df / max(1e-16, 1 - r * r))
    return n, float(r), float(2 * t.sf(abs(tt), df)), int(df), float(tt)


def _positive(v):
    return np.where(np.isfinite(v) & (v > 0), v, np.nan)


def _pair(d, m, k):
    f = list(d["focused_features"].feature)
    return d["focused_raw"][f.index(k)].astype(float), d["focused_raw"][f.index(m)].astype(float)


def _blocks(d, ctrl, idx):
    """Covariate blocks ranked within the analysis subset (``blocks`` of build_approved_eds_cr003.py)."""
    p = d["patients"]
    b = {"host": rankdata(p.Human.to_numpy()[idx])[:, None],
         "depth": rankdata(np.log2(p.pg_parsimony.to_numpy()[idx]))[:, None]}
    if ctrl is not None:
        c = ctrl["pat"]
        b["antibiotics"] = pd.get_dummies(c["antibiotic_state"], drop_first=True).to_numpy(float)[idx]
        b["clinical"] = pd.get_dummies(c["clinical_group"], drop_first=True).to_numpy(float)[idx]
    return b


def _estimate(x, y, block, subset):
    """Partial rank correlation (``estimate`` of build_approved_eds_cr003.py)."""
    Z = np.column_stack([np.ones(len(x))] + [block[k] for k in subset])
    ranks = np.column_stack([rankdata(x), rankdata(y)])
    res = ranks - Z @ np.linalg.lstsq(Z, ranks, rcond=None)[0]
    return float(np.corrcoef(res.T)[0, 1])


def stat_counts(d, sheets):
    """Counts printed in the legend (v28d line 481) and Methods (line 799).

    Returns
    -------
    dict
        Original features and observations, patients, focused pairs and measurements,
        overview features with a positive value, adjustment estimates in the sheet.
    """
    return dict(original_features=len(sheets["Fixed_feature_key"]) - 1,
                original_observations=len(sheets["All20_features_259_observations"][0]) - 1,
                observation_rows=len(sheets["Observation_metadata"]) - 1,
                patients=len(sheets["85_patient_order"]) - 1,
                focused_pairs=len(sheets["Focused_pairs"]) - 1,
                focused_measurements=len(sheets["Focused_features"]) - 1,
                adjustment_estimates=len(sheets["Focused_adjustments"]) - 1,
                overview_features_positive=int((np.nan_to_num(d["legacy_raw"]) > 0).any(axis=1).sum()))


def stat_focused_unadjusted(d):
    """Spearman correlation between each displayed metabolite and its KO signal.

    Statistical unit: patient (earliest matched observation); n = patients with positive
    finite values for both features (85, 84, 85). Two-sided t test, df = n - 2.

    Returns
    -------
    pandas.DataFrame
        One row per pair: n, rho, t, df, P, and scipy spearmanr rho for comparison.
    """
    rows = []
    for m, k, ml, kl, _ in PAIRS:
        x, y = _pair(d, m, k)
        n, r, p, df, tt = corr(_positive(x), _positive(y))
        ok = np.isfinite(_positive(x)) & np.isfinite(_positive(y))
        rows.append(dict(metabolite=m, KO=k, label=f"{ml}-{kl}", n=n, rho=r, t=tt, df=df, P=p,
                         scipy_rho=float(spearmanr(x[ok], y[ok]).statistic)))
    return pd.DataFrame(rows)


def stat_focused_features(d):
    """Spearman correlation of each displayed feature with EDS over its positive values.

    Statistical unit: patient; n = patients with a positive value (85 for all six).
    q is recomputed by Benjamini-Hochberg over the complete native unadjusted EDS family
    of the feature's layer (saved P values of all_EDS_associations.tsv).

    Returns
    -------
    pandas.DataFrame
        feature, n, rho, t, df, P, recomputed family q, family size.
    """
    e = d["all_EDS_associations"].copy()
    e["q_re"] = e.groupby(["layer", "normalization", "model"]).p.transform(lambda v: bh(v))
    eds = d["patients"].dysbiosis_score.to_numpy(float)
    rows = []
    for f, layer, v in zip(d["focused_features"].feature, d["focused_features"].layer, d["focused_raw"]):
        n, r, p, df, tt = corr(_positive(v.astype(float)), eds)
        fam = e[(e.layer == layer) & (e.normalization == "native") & (e.model == "unadjusted")]
        row = fam[fam.feature == f].iloc[0]
        rows.append(dict(feature=f, layer=layer, n=n, rho=r, t=tt, df=df, P=p, q=float(row.q_re),
                         family_size=int(fam.p.notna().sum()), saved_q=float(row.q)))
    return pd.DataFrame(rows)


def stat_family_q(d):
    """Benjamini-Hochberg q of the three pairs, recomputed in their complete screen families.

    Family: all tested pairs (finite P) of one (layer, normalization, model) of the
    augmented linked-pair screen; here layer KO, native, unadjusted or joint-adjusted.

    Returns
    -------
    tuple
        (per-pair table: model, family size, n, rho, P, q; largest absolute difference
        between recomputed and saved q_augmented over the whole screen).
    """
    s = d["augmented_pair_screen"].copy()
    s["q_re"] = s.groupby(["layer", "normalization", "model"]).p.transform(lambda v: bh(v))
    diff = float(np.nanmax(np.abs(s.q_re - s.q_augmented)))
    rows = []
    for m, k, *_ in PAIRS:
        for model in ["unadjusted", JOINT_MODEL]:
            fam = s[(s.layer == "KO") & (s.normalization == "native") & (s.model == model)]
            r = fam[(fam.metabolite == m) & (fam.protein_feature == k)].iloc[0]
            rows.append(dict(metabolite=m, KO=k, model=model, family_size=int(fam.p.notna().sum()),
                             n=int(r.n), rho=float(r.rho), P=float(r.p), q=float(r.q_re)))
    return pd.DataFrame(rows), diff


def stat_adjustment_grid(d, ctrl=None):
    """All 16 covariate combinations x 2 normalizations x 2 subsets x 3 pairs (192 estimates).

    Statistical unit: patient; n = 85, 84, 85 (pair-complete subsets) and 84 (common to
    the three pairs). Continuous covariates are ranked within each subset. Point
    estimates only (v28d line 799). Estimates that need antibiotic state or clinical group
    are NaN without the controlled input.

    Returns
    -------
    tuple
        (table with the columns of Focused_adjustments, size of the common subset).
    """
    p = d["patients"]
    ko_total = p.ko_total_signal.to_numpy(float)
    subsets = [c for n in range(len(COVARIATES) + 1) for c in itertools.combinations(COVARIATES, n)]
    masks = {}
    for m, k, *_ in PAIRS:
        x, y = _pair(d, m, k)
        masks[m] = np.isfinite(_positive(x)) & np.isfinite(_positive(y))
    common = np.logical_and.reduce(list(masks.values()))
    rows = []
    for m, k, *_ in PAIRS:
        for scope, mask in [("pair_complete", masks[m]), ("common_84", common)]:
            idx = np.flatnonzero(mask)
            for norm in NORMALIZATIONS:
                x, y = _pair(d, m, k)
                if norm == "annotated_relative":
                    x = x / ko_total
                xx, yy = x[idx], y[idx]
                block = _blocks(d, ctrl, idx)
                est = [(_estimate(xx, yy, block, s) if all(c in block for c in s) else np.nan) for s in subsets]
                for s, r in zip(subsets, est):
                    rows.append(dict(metabolite=m, KO=k, scope=scope, normalization=norm,
                                     model="+".join(s) if s else "none", n=len(idx), rho=r,
                                     delta_from_unadjusted=r - est[0]))
    return pd.DataFrame(rows), int(common.sum())


def stat_joint_adjusted(d, ctrl):
    """Partial rank correlation of each pair adjusted jointly for host share, depth, antibiotics and group.

    Screen covariate matrix: ranks of human share and log2 depth over all 85 patients,
    dummy-coded antibiotic state and clinical group. Statistical unit: patient; n = 85,
    84, 85; two-sided t test with residual df. Needs the controlled input.

    Returns
    -------
    pandas.DataFrame or None
        n, rho, t, df, P per pair; None without the controlled input.
    """
    if ctrl is None:
        return None
    p, c = d["patients"], ctrl["pat"]
    C = np.column_stack([rankdata(p.Human), rankdata(np.log2(p.pg_parsimony)),
                         pd.get_dummies(c[["antibiotic_state", "clinical_group"]], drop_first=True).to_numpy(float)])
    rows = []
    for m, k, *_ in PAIRS:
        x, y = _pair(d, m, k)
        n, r, pv, df, tt = corr(_positive(x), _positive(y), C)
        rows.append(dict(metabolite=m, KO=k, n=n, rho=r, t=tt, df=df, P=pv))
    return pd.DataFrame(rows)


def stat_five_pair(d):
    """Counts of the five-pair adjustment decomposition (v28d line 799; Source Data only).

    Read from the saved model table; the bootstrap itself is not rerun here.

    Returns
    -------
    dict
        pairs, covariate combinations per pair/subset/normalization, size of the common
        subset, bootstrap draws recorded.
    """
    f = d["five"]
    per = f.groupby(["id", "scope", "normalization"]).model.nunique()
    common = f[f.scope != "pair_complete"]
    return dict(pairs=int(f.id.nunique()), combinations=int(per.unique()[0]) if per.nunique() == 1 else -1,
                common_n=int(common.n.unique()[0]) if common.n.nunique() == 1 else -1,
                draws=int(f.bootstrap_draws.dropna().unique()[0]) if f.bootstrap_draws.dropna().nunique() == 1 else -1,
                rows=len(f), scopes=sorted(f.scope.unique()))


def word(n):
    """Number word as printed in the manuscript (three, five, six)."""
    return WORDS[int(n)]


def stat_overview_features(d):
    """Recompute Fixed_feature_key: positive patients and Spearman rho with EDS (85 patients).

    Statistical unit: patient; rho over patients with a positive value (computed when at
    least 3), as in the producer check of analyze.py.

    Returns
    -------
    pandas.DataFrame
        feature_id, n_positive, rho and the largest absolute difference to the sheet.
    """
    eds = d["patients"].dysbiosis_score.to_numpy(float)
    lf = d["legacy_features"]
    rows = []
    for f, v, nref, rref in zip(lf.feature_id, d["legacy_raw"], lf.n_positive_85, lf.rho_EDS_85):
        pos = np.isfinite(v) & (v > 0)
        r = float(spearmanr(v[pos], eds[pos]).statistic) if pos.sum() >= 3 else np.nan
        rows.append(dict(feature_id=f, n_positive=int(pos.sum()), n_sheet=int(nref), rho=r, rho_sheet=rref,
                         abs_diff=abs(r - rref) if np.isfinite(r) else np.nan))
    return pd.DataFrame(rows)


def stat_z(d):
    """Largest absolute difference between stored and recomputed within-feature z scores.

    z = (value - mean) / population SD over the positive finite values of a feature across
    the 85 patients (v28d line 801); non-detections stay empty.

    Returns
    -------
    dict
        max_abs_difference over the overview and focused matrices, and whether the empty
        cells coincide.
    """
    def z(raw):
        pos = np.where(np.isfinite(raw) & (raw > 0), raw, np.nan)
        with np.errstate(invalid="ignore", divide="ignore"):
            mu = np.array([np.nanmean(r) if np.isfinite(r).any() else np.nan for r in pos])
            sd = np.array([np.nanstd(r) if np.isfinite(r).any() else np.nan for r in pos])
            return (pos - mu[:, None]) / sd[:, None]
    zl, zf = z(d["legacy_raw"]), z(d["focused_raw"])
    diff = max(float(np.nanmax(np.abs(zl - d["legacy_z"]))), float(np.nanmax(np.abs(zf - d["focused_z"]))))
    same = bool(np.array_equal(np.isnan(zl), np.isnan(d["legacy_z"])) and np.array_equal(np.isnan(zf), np.isnan(d["focused_z"])))
    return dict(max_abs_difference=diff, empty_cells_identical=same)
