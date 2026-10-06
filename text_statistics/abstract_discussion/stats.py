"""Statistics printed only in the Abstract (v28d line 35) and the Discussion (v28d lines 181-193).

No figure panel belongs to this unit. Every value is recomputed from tables that other
units already ship; following ``source_data/README.md`` they are read where they live and
are not duplicated:

``source_data/Fig4/crossfit_per_acquisition.tsv``
    10-fold cross-fitted peptide counts per matched acquisition (578 rows: 296 PREDICT,
    282 CapScan). Abstract "3.5-5.5% fewer ... in two stool cohorts" and "21-36% more in
    intestinal fluid"; Discussion line 183 ("9-50-residue").
``source_data/Fig1/acquisition_sizes_726.tsv`` and ``source_data/Fig1/human_reference.json``
    Entry counts of the 726 selected FASTAs, the study-reservoir sizes and the 20,431-entry
    human reference. Abstract and Discussion "All 726 databases ... smaller than the human
    reference proteome"; Discussion "233-375 million sequences".
``source_data/SFig15/compute_stage_summary.tsv`` and ``source_data/SFig15/compute_task_records.tsv``
    Slurm accounting of the production chain. Discussion "8-12 min per cohort",
    "1-3 min per acquisition", "less than 1 min".
``source_data/Fig2/member_summary.tsv``
    CAMPI members. Abstract "eight-member community".

The bootstrap is moved unchanged from the original producer
(``LAB_ROOT/.../iteration4/crossfit/run/analyze_final.py``, function ``summarize``), see
PROVENANCE.md. The timing statistics follow the original producer's Computing block
(``LAB_ROOT/.../revision_2026-09-24/build_revision.py`` lines 201-228).
"""

import json

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- choices
MIN_LENGTH = 9                 # Discussion v28d line 183, Methods line 777: 9-50-residue peptides
MAX_LENGTH = 50                # upper bound of every peptide count column used here (Methods)
B = 10_000                     # bootstrap draws; analyze_final.py `B` (Methods v28d line 777)
SEED = 20261004                # analyze_final.py `SEED`; not printed in v28d (OPEN_ITEMS TEXT-GAP)
CI_LEVEL = 95                  # percentile interval, 2.5th and 97.5th percentile of the draws
MAIN_STRATA = ("PREDICT", "CAPSCAN_all", "CAPSCAN_stool", "CAPSCAN_intestinal")  # one RNG stream, this order
REGION_STRATA = ("duodenum", "jejunum", "ileum", "colon", "stool")             # second stream, fresh seed, this order
INTESTINAL_REGIONS = ("duodenum", "jejunum", "ileum", "colon")                 # "intestinal fluid"
LOOKUP_STAGE = "02_evidence"   # "study-level lookup" (Discussion line 193)
SELECTION_STAGES = ("05_persample", "06_parsimony")  # "database selection" (Discussion line 193)
SAGE_STAGE = "07_sage"         # "Sage searches" (Discussion line 193)
COHORTS = ("CAMPI", "HSTOOL", "PREDICT", "CAPSCAN")  # the four cohorts in the compute tables


# ---------------------------------------------------------------- loading
def load_crossfit(sd_fig4):
    """Read the cross-fitted per-acquisition table and add the gain columns.

    Parameters
    ----------
    sd_fig4 : pathlib.Path
        ``source_data/Fig4``.

    Returns
    -------
    pandas.DataFrame
        One row per matched acquisition (578), with ``gain_cf_pooled_5`` and
        ``gain_cf_pooled_9`` = 100 * (FL / MG - 1) as in analyze_final.py.
    """
    a = pd.read_csv(sd_fig4 / "crossfit_per_acquisition.tsv", sep="\t", keep_default_na=False,
                    na_values=[""])
    for L in ("5", "9"):
        fl, mg = a[f"FL_pooled_{L}"].astype(float), a[f"MG_pooled_{L}"].astype(float)
        a[f"gain_cf_pooled_{L}"] = np.where(mg > 0, 100 * (fl / mg - 1), np.nan)
    return a


def strata_masks(a):
    """Return the stratum masks of analyze_final.py, keyed by stratum name."""
    cap = a.cohort.eq("CAPSCAN")
    m = {"PREDICT": a.cohort.eq("PREDICT"), "CAPSCAN_all": cap,
         "CAPSCAN_stool": cap & a.sampling_site.eq("stool"),
         "CAPSCAN_intestinal": cap & a.sampling_site.ne("stool")}
    for g in REGION_STRATA:
        m[f"CAPSCAN_{g}"] = cap & a.sampling_site.eq(g)
    return m


# ---------------------------------------------------------------- cross-fitted gain
def _summarize(a, strata, rng, col):
    """Median and participant-cluster bootstrap interval per stratum (analyze_final.summarize).

    Clusters are participants in the order of the original raw participant strings
    (``participant_sort_rank``); one RNG stream serves all strata in the given order, and
    the B index sets of a stratum are drawn before any column is evaluated, exactly as in
    the original code.
    """
    out = {}
    for name, m in strata:
        g = a[m].reset_index(drop=True)
        subj = g.participant_sort_rank.to_numpy()
        us = np.unique(subj)
        members = [np.flatnonzero(subj == s) for s in us]
        draws = [np.concatenate([members[k] for k in rng.integers(0, len(us), len(us))]) for _ in range(B)]
        v = g[col].to_numpy(dtype=float)
        bs = np.array([np.nanmedian(v[ix]) for ix in draws])
        lo, hi = (100 - CI_LEVEL) / 2, 100 - (100 - CI_LEVEL) / 2
        out[name] = dict(n_acquisitions=len(g), n_participants=len(us), median=float(np.nanmedian(v)),
                         ci_low=float(np.percentile(bs, lo)), ci_high=float(np.percentile(bs, hi)),
                         n_fl_gt_mg=int((v > 0).sum()), n_fl_lt_mg=int((v < 0).sum()),
                         n_nan=int(np.isnan(v).sum()))
    return out


def stat_crossfit_gain(a):
    """Median cross-fitted FastaLake-over-metagenome peptide gain, 9-50 residues.

    Statistical unit: acquisition (estimand: median over acquisitions of
    100 * (FL / MG - 1)); interval by participant-cluster percentile bootstrap,
    B = 10,000, seed 20261004. n: PREDICT 296 acquisitions (97 participants, NOID
    acquisitions as own clusters); CapScan 282 (15); CapScan stool 56 (15); intestinal
    fluid 226 (14); duodenum 56, jejunum 65, ileum 54, colon 51 acquisitions.

    Parameters
    ----------
    a : pandas.DataFrame
        Output of ``load_crossfit``.

    Returns
    -------
    dict
        ``main`` (four strata of the first RNG stream) and ``region`` (five strata of
        the second stream), each stratum -> n, median, ci_low, ci_high and sign counts;
        ``region_range`` (min, max) of the four intestinal-region medians.
    """
    col = f"gain_cf_pooled_{MIN_LENGTH}"
    m = strata_masks(a)
    main = _summarize(a, [(s, m[s]) for s in MAIN_STRATA], np.random.default_rng(SEED), col)
    region = _summarize(a, [(f"CAPSCAN_{g}", m[f"CAPSCAN_{g}"]) for g in REGION_STRATA],
                        np.random.default_rng(SEED), col)
    meds = [region[f"CAPSCAN_{g}"]["median"] for g in INTESTINAL_REGIONS]
    n_stool_cohorts = int(a.loc[a.sampling_site.eq("stool"), "cohort"].nunique())
    return dict(main=main, region=region, region_range=(min(meds), max(meds)),
                n_stool_cohorts=n_stool_cohorts)


# ---------------------------------------------------------------- database sizes
def stat_database_sizes(sd_fig1):
    """Selected-database sizes against the human reference proteome.

    Statistical unit: acquisition (one FASTA each); n = 726 (15 CAMPI, 58 healthy
    stool, 306 PREDICT, 347 CapScan). Descriptive census, no test.

    Parameters
    ----------
    sd_fig1 : pathlib.Path
        ``source_data/Fig1``.

    Returns
    -------
    dict
        ``n_databases``, ``n_below_human`` (entries < human entries), ``max_entries``,
        ``human_entries``, ``n_by_cohort``, ``reservoir_min`` and ``reservoir_max``
        (entries of the four study reservoirs) and the same in millions.
    """
    s = pd.read_csv(sd_fig1 / "acquisition_sizes_726.tsv", sep="\t")
    human = json.loads((sd_fig1 / "human_reference.json").read_text())["entries"]
    res = s.groupby("cohort").reservoir_entries.unique()
    assert all(len(v) == 1 for v in res), "one reservoir per cohort"
    res = {k: int(v[0]) for k, v in res.items()}
    return dict(n_databases=int(len(s)), n_below_human=int((s.entries < human).sum()),
                max_entries=int(s.entries.max()), human_entries=int(human),
                n_by_cohort=s.cohort.value_counts().to_dict(),
                reservoir_entries=res, reservoir_min=min(res.values()), reservoir_max=max(res.values()),
                reservoir_min_millions=min(res.values()) / 1e6, reservoir_max_millions=max(res.values()) / 1e6)


# ---------------------------------------------------------------- computing time
def stat_compute_times(sd_sfig15):
    """Wall times of the production chain (Discussion line 193).

    Statistical unit: Slurm task. Study-level lookup: one task per cohort (4 cohorts);
    the median wall time per cohort is that task's time. Database selection: per
    acquisition, the wall minutes of its ``SELECTION_STAGES`` tasks summed, then the
    median over acquisitions per cohort (n = 15, 58, 306, 347); the sum of the two stage
    medians is reported as the alternative reading. Sage: median wall time of the
    ``SAGE_STAGE`` tasks per cohort.

    Parameters
    ----------
    sd_sfig15 : pathlib.Path
        ``source_data/SFig15``.

    Returns
    -------
    dict
        ``lookup`` (cohort -> min), ``lookup_range``; ``selection_per_acquisition``
        (cohort -> median min), ``selection_n``, ``selection_range``;
        ``selection_sum_of_medians`` and its range; ``sage`` (cohort -> median min) and
        ``sage_max``.
    """
    s = pd.read_csv(sd_sfig15 / "compute_stage_summary.tsv", sep="\t")
    j = pd.read_csv(sd_sfig15 / "compute_task_records.tsv", sep="\t", dtype={"jobid": str})
    lookup = s[s.stage.eq(LOOKUP_STAGE)].set_index("cohort").wall_median_min.reindex(COHORTS)
    assert (s[s.stage.eq(LOOKUP_STAGE)].n_tasks == 1).all()
    sage = s[s.stage.eq(SAGE_STAGE)].set_index("cohort").wall_median_min.reindex(COHORTS)
    g = j[j.stage.isin(SELECTION_STAGES) & j.cohort.isin(COHORTS)].copy()
    g["key"] = g.acquisition.where(g.acquisition.notna(), g.jobid.str.split("_").str[1])
    med, n, alt = {}, {}, {}
    for co in COHORTS:
        x = g[g.cohort.eq(co)]
        per = x.groupby("key").agg(total=("wall_minutes", "sum"), k=("stage", "nunique"))
        assert (per.k == len(SELECTION_STAGES)).all(), co
        med[co], n[co] = float(per.total.median()), int(len(per))
        alt[co] = float(sum(x[x.stage.eq(st)].wall_minutes.median() for st in SELECTION_STAGES))
    return dict(lookup=lookup.to_dict(), lookup_range=(float(lookup.min()), float(lookup.max())),
                selection_per_acquisition=med, selection_n=n,
                selection_range=(min(med.values()), max(med.values())),
                selection_sum_of_medians=alt,
                selection_sum_of_medians_range=(min(alt.values()), max(alt.values())),
                sage=sage.to_dict(), sage_max=float(sage.max()))


# ---------------------------------------------------------------- CAMPI members
def stat_campi_members(sd_fig2):
    """Number of members of the defined CAMPI community (Abstract "eight-member").

    Statistical unit: species; n = distinct ``member`` values in the CAMPI member table.
    """
    m = pd.read_csv(sd_fig2 / "member_summary.tsv", sep="\t")
    return int(m.member.nunique())


WORDS = {2: "two", 8: "eight"}


def number_word(n):
    """Return the English word the manuscript prints for a small count."""
    return WORDS[n]


# ---------------------------------------------------------------- all numbers, NUMBERS.tsv order
def numbers(gain, sizes, times, n_members):
    """Every recomputed value of this unit, in ``NUMBERS.tsv`` row order.

    Parameters
    ----------
    gain, sizes, times : dict
        Outputs of ``stat_crossfit_gain``, ``stat_database_sizes`` and ``stat_compute_times``.
    n_members : int
        Output of ``stat_campi_members``.

    Returns
    -------
    list of (str, object)
        ``(key, value)``; ``None`` marks a row this unit cannot recompute (OPEN).
    """
    m = gain["main"]
    return [
        ("L35 eight-member community", number_word(n_members)),
        ("L35 1% false discovery rate", None),
        ("L35 fewer peptides, PREDICT stool (-median gain, %)", -m["PREDICT"]["median"]),
        ("L35 fewer peptides, CapScan stool (-median gain, %)", -m["CAPSCAN_stool"]["median"]),
        ("L35 two stool cohorts", number_word(gain["n_stool_cohorts"])),
        ("L35 intestinal fluid, smallest region median gain (%)", gain["region_range"][0]),
        ("L35 intestinal fluid, largest region median gain (%)", gain["region_range"][1]),
        ("L35 726 databases (all below human reference)", sizes["n_below_human"]),
        ("L181 726 databases (all below human reference)", sizes["n_below_human"]),
        ("L183 minimum peptide length", MIN_LENGTH),
        ("L183 maximum peptide length", MAX_LENGTH),
        ("L193 smallest reservoir (million sequences)", sizes["reservoir_min_millions"]),
        ("L193 largest reservoir (million sequences)", sizes["reservoir_max_millions"]),
        ("L193 shortest study-level lookup (min)", times["lookup_range"][0]),
        ("L193 longest study-level lookup (min)", times["lookup_range"][1]),
        ("L193 smallest cohort median selection time per acquisition (min)", times["selection_range"][0]),
        ("L193 largest cohort median selection time per acquisition (min)", times["selection_range"][1]),
        ("L193 largest cohort median Sage time (min)", times["sage_max"]),
    ]
