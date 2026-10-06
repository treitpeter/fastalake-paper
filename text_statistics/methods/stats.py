"""Statistics printed only in the main Methods (v28d lines 197-307).

Which numbers
-------------
* line 215: reservoir sizes and the number of assembly FASTAs pooled per cohort reservoir;
* line 227: records and contaminants of the supplied CAMPI reference and MG FASTAs;
* line 233 and 303: the five PREDICT acquisitions searched with reconverted mzML files;
* line 267: participant bootstrap of the spectrum-split comparison, the 578 cross-fitted
  acquisitions, the three acquisitions with a singular pooled fit, and the interpreter;
* line 269: exact-sequence unions of the selected FASTAs;
* line 275: patient bootstrap and sign-flip test of the antibiotic comparison;
* line 299: compute times and hardware.

Source data
-----------
``source_data/methods/``: reservoir_sources.tsv, campi_comparator_records.tsv,
union_contributions.tsv, crossfit_lda_fallback_347.tsv,
split_spectra_mg_per_acquisition.tsv, split_spectra_mg_summary_original.tsv,
predict_reconversion.tsv, compute_nodes.tsv, crossfit_interpreter.tsv.
Shared tables of other units (named, not copied): ``source_data/Fig1/acquisition_sizes_726.tsv``,
``source_data/Fig4/crossfit_per_acquisition.tsv``,
``source_data/Fig4/clinical_paired_statistics.tsv``,
``source_data/SFig15/compute_task_records.tsv``, ``source_data/SFig15/compute_stage_summary.tsv``,
``source_data/SFig15/campi_alphanovo_inference_sacct.tsv``.
Controlled: ``paired_antibiotic_profiles.tsv`` from ``FASTALAKE_CONTROLLED_DIR``.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd

UNIT_DIR = Path(__file__).resolve().parent
REPO = UNIT_DIR.parents[1]
SD = REPO / "source_data"

# ---------------------------------------------------------------- choices (named constants)
SPLIT_BOOT_DRAWS = 10_000        # Methods, v28d line 267 (iteration3/disjoint/analyze.py B)
SPLIT_BOOT_SEED = 20261003       # not stated in the text; iteration3/disjoint/analyze.py SEED
SPLIT_CI_PERCENTILES = (2.5, 97.5)   # percentile interval, analyze.py
CLIN_BOOT_DRAWS = 50_000         # Methods, v28d line 275 (clinical_refit.py derive)
CLIN_SIGNFLIP_BATCHES = 100      # clinical_refit.py: 100 batches of 2,000 sign assignments
CLIN_SIGNFLIP_BATCH_SIZE = 2_000
CLIN_SIGNFLIPS = CLIN_SIGNFLIP_BATCHES * CLIN_SIGNFLIP_BATCH_SIZE   # Methods, v28d line 275: 200,000
CLIN_SEED = 20260907             # not stated in the text; clinical_refit.py (bootstrap seed; sign-flip seed + endpoint index)
CLIN_TIE_TOL = 1e-12             # clinical_refit.py: |null| >= |observed| - 1e-12
CLIN_ENDPOINTS = [("Human", 100, "percentage points"), ("pg_parsimony", 1, "quantified simplified groups")]
CAMPI_LAKE_STAGE = "lake"        # SFig15 compute_task_records stage of the reservoir build
CAMPI_INFERENCE_ARRAY = "5450700"   # Slurm array of the CAMPI AlphaNovo inference (line 299)
CAMPI_STAGES = {"extraction": "05_persample", "selection": "06_parsimony", "search": "07_sage"}  # line 299
POOLED_EXTRACTION_STAGE = "02_evidence"   # line 299, "Pooled extraction took 11 min"
NUMBER_WORDS = {3: "three", 5: "five"}


def sd(unit, name, **kw):
    """Read ``source_data/<unit>/<name>`` as a tab-separated table."""
    return pd.read_csv(SD / unit / name, sep="\t", **kw)


# ---------------------------------------------------------------- line 215
def stat_reservoirs():
    """Reservoir sizes and assembly FASTA counts (v28d line 215).

    Each reservoir size is the sum of the unique sequences that each source added in
    priority order, as written by the lake builder. Cross-checked against the
    ``reservoir_entries`` column of ``source_data/Fig1/acquisition_sizes_726.tsv``.
    Statistical unit: sequence; n = one reservoir each (descriptive).

    Returns
    -------
    dict
        ``size`` per reservoir, ``assembly_fastas`` per cohort, ``fig1_agrees`` (bool).
    """
    s = sd("methods", "reservoir_sources.tsv")
    size = s.groupby("reservoir").unique_sequences.sum().to_dict()
    fastas = {"PREDICT": int(s[s.source.eq("PREDICT_MG")].files.iloc[0]),
              "CapScan": int(s[s.source.eq("CAPSCAN_MG")].files.iloc[0])}
    a = sd("Fig1", "acquisition_sizes_726.tsv")
    fig1 = a.groupby("cohort").reservoir_entries.unique()
    agree = (list(fig1["HSTOOL"]) == [size["public"]] and list(fig1["PREDICT"]) == [size["PREDICT_cohort"]]
             and list(fig1["CapScan"]) == [size["CapScan_cohort"]] and list(fig1["CAMPI"]) == [size["CAMPI"]])
    return dict(size={k: int(v) for k, v in size.items()}, assembly_fastas=fastas, fig1_agrees=bool(agree))


# ---------------------------------------------------------------- line 227
def stat_comparator_fastas():
    """Records and contaminants of the supplied CAMPI reference and MG FASTAs (line 227).

    Contaminants are the records whose accession occurs in both supplied FASTAs: the
    MG FASTA otherwise holds only assembly proteins (``k99_`` contig accessions), so the
    shared accessions are the appended contaminant block. Statistical unit: sequence
    record; n = 29,673 and 19,435 records (descriptive).

    Returns
    -------
    dict
        ``records`` and ``contaminants`` per database, and consistency checks.
    """
    r = sd("methods", "campi_comparator_records.tsv")
    acc = {db: set(g.accession) for db, g in r.groupby("database")}
    shared = acc["reference"] & acc["MG"]
    out = dict(records={db: int((r.database == db).sum()) for db in ["reference", "MG"]},
               contaminants={db: int((r.database.eq(db) & r.accession.isin(shared)).sum()) for db in ["reference", "MG"]})
    mg_other = r[r.database.eq("MG") & ~r.accession.isin(shared)]
    out["mg_non_contaminants_all_k99"] = bool(mg_other.accession.str.startswith("k99_").all())
    seq = {db: set(r[r.database.eq(db) & r.accession.isin(shared)].sequence_sha256) for db in ["reference", "MG"]}
    out["contaminant_sequences_identical"] = seq["reference"] == seq["MG"]
    out["accessions_unique_within_database"] = bool(not r.duplicated(["database", "accession"]).any())
    return out


# ---------------------------------------------------------------- lines 233 and 303
def stat_reconversion():
    """PREDICT acquisitions searched with reconverted mzML files (lines 233 and 303).

    Statistical unit: acquisition; n = 306 PREDICT public-sequence searches.

    Returns
    -------
    dict
        ``n`` reconverted, their acquisition labels (``MicrobPred_<k>``), whether all of
        them were truncated in the original folder and all others complete.
    """
    p = sd("methods", "predict_reconversion.tsv", dtype=str)
    rec = p[p.searched_reconverted_mzml.eq("True")]
    labels = sorted((a.split("_SA_")[1].split("_")[1] for a in rec.acquisition), key=int)
    return dict(n=len(rec), n_word=NUMBER_WORDS[len(rec)], labels=labels,
                text="MicrobPred_" + ", ".join(labels[:-1]) + " and " + labels[-1],
                reconverted_were_truncated=bool(rec.original_mzml_ending.eq("TRUNCATED").all()),
                others_complete=bool(p[~p.index.isin(rec.index)].original_mzml_ending.eq("COMPLETE").all()),
                all_searched_mzml_closed=bool(p.searched_mzml_closing_tag.eq("True").all()),
                n_searches=len(p))


# ---------------------------------------------------------------- line 267
def split_bootstrap(per_acq=None, draws=SPLIT_BOOT_DRAWS, seed=SPLIT_BOOT_SEED, order_by_rank=True):
    """Participant cluster bootstrap of acquisition-median gains, spectrum-split comparison.

    Ported from ``iteration3/disjoint/analyze.py`` with the logic unchanged: gain =
    100 (FL/MG - 1) per acquisition; median over acquisitions; participants resampled
    with replacement within the stratum keeping all their acquisitions; ``draws`` draws
    shared by all comparisons of a stratum; percentile 95% interval. Participant ids are
    ordered as in the original run: the original formed clusters in ``np.unique`` order
    of the raw participant strings; here each acquisition's cluster key is
    ``participant_sort_rank`` of ``source_data/Fig4/crossfit_per_acquisition.tsv`` (same
    roster index ``i``), which keeps that order without the raw identifiers. With
    ``order_by_rank=False`` the published pseudonyms are used as keys instead.

    Statistical unit: participant (resampling) over acquisitions (estimate); n = 296
    PREDICT acquisitions from 97 patients, 282 CapScan acquisitions from 15 participants.

    Returns
    -------
    pandas.DataFrame
        One row per stratum, length view and comparison, with the unrounded median and
        interval bounds.
    """
    a = sd("methods", "split_spectra_mg_per_acquisition.tsv", dtype={"subject_id": str},
           keep_default_na=False) if per_acq is None else per_acq.copy()
    if order_by_rank:
        rk = sd("Fig4", "crossfit_per_acquisition.tsv", dtype=str, keep_default_na=False)
        rk = rk.set_index(rk.crossfit_index.astype(int))
        a = a.assign(subject_id=[f"{int(rk.participant_sort_rank[int(i)]):05d}" for i in a.i])
        assert (rk.cohort[a.i.astype(int)].to_numpy() == a.cohort.to_numpy()).all()
    rng = np.random.default_rng(seed)
    strata = [("PREDICT", a.cohort == "PREDICT"), ("CAPSCAN_all", a.cohort == "CAPSCAN"),
              ("CAPSCAN_stool", (a.cohort == "CAPSCAN") & (a.group == "stool")),
              ("CAPSCAN_intestinal", (a.cohort == "CAPSCAN") & (a.group != "stool"))]
    out = []
    for name, m in strata:
        g = a[m].reset_index(drop=True)
        subj = g.subject_id.values
        us = np.unique(subj)
        members = [np.flatnonzero(subj == s) for s in us]
        ix_all = [np.concatenate([members[k] for k in rng.integers(0, len(us), len(us))]) for _ in range(draws)]
        for L in ["5", "9"]:
            cols = {c: g[f"{c}_{L}"].astype(float).values for c in
                    ["gain_frozen_fullspectra", "gain_odd_full", "gain_odd_self", "gain_odd_disjoint"]}
            cols["disjoint_minus_full_selection_gain_points"] = cols["gain_odd_disjoint"] - cols["gain_odd_full"]
            for col, v in cols.items():
                bs = np.array([np.median(v[ix]) for ix in ix_all])
                lo, hi = np.percentile(bs, SPLIT_CI_PERCENTILES)
                out.append(dict(stratum=name, min_length=int(L), comparison=col, n_acquisitions=len(g),
                                n_participants=len(us), median_gain_pct=float(np.median(v)),
                                ci95_low=float(lo), ci95_high=float(hi), n_fl_gt_mg=int((v > 0).sum()),
                                n_fl_lt_mg=int((v < 0).sum()), n_equal=int((v == 0).sum())))
    return pd.DataFrame(out)


def compare_split_bootstrap(recomputed):
    """Compare :func:`split_bootstrap` with the original summary (rounded to 2 decimals).

    Returns
    -------
    pandas.DataFrame
        Original and recomputed medians and bounds, and whether each agrees at 2 decimals.
    """
    o = sd("methods", "split_spectra_mg_summary_original.tsv")
    m = o.merge(recomputed, on=["stratum", "min_length", "comparison"], suffixes=("_original", ""))
    for c in ["median_gain_pct", "ci95_low", "ci95_high"]:
        m[f"{c}_agrees"] = np.isclose(m[f"{c}_original"], m[c].round(2), atol=1e-9)
    return m


def stat_crossfit_roster():
    """Matched cross-fitted acquisitions and singular pooled fits (line 267).

    Reads the production cross-fit table (``source_data/Fig4/crossfit_per_acquisition.tsv``)
    and the rerun on the 347-acquisition CapScan roster. Statistical unit: acquisition.

    Returns
    -------
    dict
        ``n`` matched acquisitions, ``n_fallback`` acquisitions with any LDA fallback, for
        both versions, and the fallback arms.
    """
    x = sd("Fig4", "crossfit_per_acquisition.tsv", dtype=str, keep_default_na=False)
    y = sd("methods", "crossfit_lda_fallback_347.tsv", dtype=str, keep_default_na=False)
    fx, fy = x[x.lda_fallback.ne("")], y[y.lda_fallback.ne("")]
    return dict(n=len(x), n_by_cohort=x.cohort.value_counts().to_dict(), n_fallback=len(fx),
                n_fallback_word=NUMBER_WORDS[len(fx)],
                fallback=sorted(zip(fx.cohort, fx.lda_fallback)),
                n_347=len(y), n_fallback_347=len(fy), fallback_347=sorted(zip(fy.cohort, fy.lda_fallback)),
                fallback_index_equal=sorted(fx.crossfit_index.astype(int)) == sorted(fy.i.astype(int)))


def stat_crossfit_interpreter():
    """Python and NumPy versions of the cross-fit pooled-rescoring interpreter (line 267)."""
    r = sd("methods", "crossfit_interpreter.tsv", dtype=str).iloc[0]
    return dict(python=r.python, numpy=r.numpy)


# ---------------------------------------------------------------- line 269
def stat_unions():
    """Exact-sequence unions of the selected FASTAs (line 269).

    The union size is the sum over FASTAs of the sequences each adds to the union
    (exact full sequence, I/L distinct). Statistical unit: protein sequence.

    Returns
    -------
    dict
        Per roster: ``fastas``, ``records``, ``unique``.
    """
    u = sd("methods", "union_contributions.tsv")
    g = u.groupby("roster")
    return {r: dict(fastas=int(len(d)), records=int(d.records.sum()), unique=int(d.new_unique_sequences.sum()))
            for r, d in g}


# ---------------------------------------------------------------- line 275
def _bh(p):
    """Benjamini-Hochberg q values, as in clinical_refit.py ``bh``."""
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    out = np.empty_like(p)
    out[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return out


def clinical_paired_tests(paired):
    """Patient bootstrap and Monte Carlo sign-flip test, antibiotic on minus off (line 275).

    Ported from ``clinical_refit.py`` ``derive`` (paired part) with the logic, seeds and
    tails unchanged. Two-sided test of the mean on-minus-off difference: P = (count + 1) /
    (200,000 + 1), count = sign assignments with |mean| >= |observed mean| - 1e-12;
    percentile 95% interval from 50,000 patient-bootstrap means; BH across the two
    endpoints. Statistical unit: patient; n = 42 patients observed in both states.

    Parameters
    ----------
    paired : pandas.DataFrame
        Controlled ``paired_antibiotic_profiles.tsv`` (84 rows, 42 patients x 2 states).

    Returns
    -------
    pandas.DataFrame
        One row per endpoint (columns as ``clinical_paired_statistics.tsv``).
    """
    assert len(paired) == 84 and paired.patient_id.nunique() == 42
    stats = []
    for idx, (metric, scale, unit) in enumerate(CLIN_ENDPOINTS):
        wide = paired.pivot(index="patient_id", columns="antibiotic_state", values=metric)[["Off", "On"]] * scale
        delta = (wide.On - wide.Off).to_numpy()
        rng = np.random.default_rng(CLIN_SEED)
        boots = delta[rng.integers(0, len(delta), (CLIN_BOOT_DRAWS, len(delta)))].mean(axis=1)
        lo, hi = np.quantile(boots, [.025, .975])
        rng = np.random.default_rng(CLIN_SEED + idx)
        count = 0
        for _ in range(CLIN_SIGNFLIP_BATCHES):
            null = ((rng.integers(0, 2, (CLIN_SIGNFLIP_BATCH_SIZE, len(delta))) * 2 - 1) * delta).mean(axis=1)
            count += int((abs(null) >= abs(delta.mean()) - CLIN_TIE_TOL).sum())
        stats.append(dict(metric=metric, n_patients=len(delta), mean_off=wide.Off.mean(), mean_on=wide.On.mean(),
                          mean_on_minus_off=delta.mean(), CI95_low=lo, CI95_high=hi, units=unit,
                          p_signflip_MC=(count + 1) / (CLIN_SIGNFLIPS + 1), exceed_count=count))
    stats = pd.DataFrame(stats)
    stats["q_two_clinical_endpoints"] = _bh(stats.p_signflip_MC)
    return stats


def stat_clinical():
    """Line 275 counts from the public aggregate table, and the recomputation when allowed.

    Returns
    -------
    dict
        ``n_patients`` (public table), ``draws``, ``signflips``, ``public`` table and
        ``recomputed`` table (None without ``FASTALAKE_CONTROLLED_DIR``).
    """
    pub = sd("Fig4", "clinical_paired_statistics.tsv")
    out = dict(n_patients=int(pub.n_patients.unique()[0]), draws=CLIN_BOOT_DRAWS, signflips=CLIN_SIGNFLIPS,
               public=pub, recomputed=None)
    root = os.environ.get("FASTALAKE_CONTROLLED_DIR")
    path = Path(root) / "paired_antibiotic_profiles.tsv" if root else None
    if path is None or not path.is_file():
        print("SKIPPED: needs controlled input paired_antibiotic_profiles.tsv (see controlled/CONTROLLED_INPUTS.md)")
    else:
        out["recomputed"] = clinical_paired_tests(pd.read_csv(path, sep="\t"))
    return out


# ---------------------------------------------------------------- line 299
def stat_compute():
    """CAMPI compute times and hardware (line 299).

    Statistical unit: Slurm job or task; CAMPI n = 1 reservoir job, 1 pooled extraction
    job, 15 tasks per acquisition stage, 15 inference tasks (descriptive).

    Returns
    -------
    dict
        Minutes, medians, builder threads, node memory (GB, 1 GB = 1000 MB as Slurm
        reports RealMemory in MB) and node features.
    """
    t = sd("SFig15", "compute_task_records.tsv", dtype={"jobid": str})
    s = sd("SFig15", "compute_stage_summary.tsv")
    lake = t[t.cohort.eq("CAMPI") & t.stage.eq(CAMPI_LAKE_STAGE)]
    assert len(lake) == 1
    lake = lake.iloc[0]
    camp = t[t.cohort.eq("CAMPI")]
    med = {k: float((camp[camp.stage.eq(v)].elapsed_s / 60).median()) for k, v in CAMPI_STAGES.items()}
    n_tasks = {k: int(camp.stage.eq(v).sum()) for k, v in CAMPI_STAGES.items()}
    summ = {k: float(s[s.cohort.eq("CAMPI") & s.stage.eq(v)].wall_median_min.iloc[0]) for k, v in CAMPI_STAGES.items()}
    ev = camp[camp.stage.eq(POOLED_EXTRACTION_STAGE)]
    assert len(ev) == 1
    sac = sd("SFig15", "campi_alphanovo_inference_sacct.tsv", dtype=str)
    inf = sac[sac.JobID.str.match(rf"^{CAMPI_INFERENCE_ARRAY}_\d+$")]
    nodes = sd("methods", "compute_nodes.tsv", dtype=str)
    lake_node = nodes[nodes.node.eq(lake.node)].iloc[0]
    inf_nodes = nodes[nodes.node.isin(set(inf.NodeList))]
    return dict(lake_minutes=float(lake.elapsed_s) / 60, lake_jobid=str(lake.jobid), lake_node=lake.node,
                builder_threads=int(sd("methods", "reservoir_sources.tsv").query("reservoir == 'CAMPI'").builder_threads.iloc[0]),
                node_memory_GB=int(lake_node.RealMemory_MB) / 1000, node_features=lake_node.AvailableFeatures,
                node_logical_cpus=int(lake_node.CPUTot),
                pooled_extraction_minutes=float(ev.elapsed_s.iloc[0]) / 60,
                stage_medians=med, stage_n=n_tasks, stage_medians_summary_table=summ,
                inference_median_minutes=float((inf.ElapsedRaw.astype(float) / 60).median()),
                inference_n=len(inf), inference_gres=sorted(set(inf_nodes.Gres)))
