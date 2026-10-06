"""Supplementary Methods 6-7 (validation, error controls, transfer and pooling): statistics.

Recomputes every result number printed in v28d lines 661-777 (and line 633, the
additional-confidence rosters) from ``source_data/supplementary_methods_6_7/``. Design
parameters typed in these lines (seeds, thresholds, length windows, software versions)
are the named constants below; they are not results and have no NUMBERS.tsv row.

Sections and their source_data tables:

* Sequence controls (lines 665-667): scramble_end_to_end.tsv, scramble_lookup_comparison.tsv
* Database sizes (673-675): acquisition_sizes_726.tsv, human_reference.tsv
* Selection sensitivity (679): selection_grid_acquisitions.tsv
* CAMPI recovery (683-687): campi_spectral_inventory.tsv, campi_reference_peptides.tsv,
  campi_member_summary.tsv, campi_reference_fasta_records.tsv
* Candidate clustering (695): clustering_per_acquisition_identifications.tsv
* CAMPI fusion control (713-715): fusion_campi_acquisition_components.tsv
* Reservoir occurrence (721): label_validity_presence.tsv, label_validity_result.tsv
* Donor-selected entrapment (727-729): donor_recount.tsv
* Restricted simulation (733): simulation_all_outcomes.tsv, simulation_design.tsv
* Transfer and pooling (739-745): transfer_pool_comparisons.tsv, hstool_pooling_per_acquisition.tsv
* Common search space (751): common_fasta_union_summary.tsv, representative_presence.tsv,
  union_contributions.tsv
* Spectrum split (755-761): split_half_acquisitions_710.tsv, split_half_model_fallbacks.tsv,
  split_half_zero_yield.tsv
* Joint search (765-769): joint_counts_per_acquisition.tsv, joint_fileq_coverage.tsv,
  microbial_group_pairs.tsv, technical_agreement_summary.tsv
* 10-fold cross-fitting (773-777): crossfit_per_acquisition.tsv, crossfit_lda_fallback_acquisitions.tsv,
  crossfit_q_validation.tsv, crossfit_k1_control.tsv, crossfit_pilot_folds.tsv, crossfit_fold_caches.tsv,
  crossfit_independent_lda_check.tsv,
  confidence_ordinary_all_pairs.tsv, minlen7_benchmark_per_acquisition.tsv
* Additional confidence rosters (633): confidence_ordinary_all_pairs.tsv, confidence_fusion_all.tsv,
  confidence_pep01_zero_yield_cases.tsv
"""

from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- design constants
SCRAMBLE_SEEDS = (20260901, 20260913, 20260914, 20260915, 20260916)  # v28d line 665
Q_CUTOFF = 0.01                    # peptide q, v28d lines 667, 715, 721, 727, 733, 739, 757, 765
FULL_LENGTH = (5, 50)              # accepted canonical peptides, v28d lines 715, 721, 733, 739, 757
NESTED_LENGTH = (9, 50)            # nested subset, v28d lines 715, 721, 757, 775
SHORT_LENGTH = (5, 8)              # v28d lines 715, 721, 777
LOOKUP_TOP_FRACTION = 0.30         # top-30% production lookup queries, v28d lines 677, 679, 755, 773
SELECTION_GRID_LENGTHS = (6, 12)   # minimum evidence lengths, v28d line 679; asserted in stat_selection_grid
CLUSTER_IDENTITIES = (0.99, 0.97, 0.95, 0.90, 0.75, 0.50)  # Linclust --min-seq-id, v28d line 695; asserted in stat_clustering
FUSION_SEEDS = (20260907, 20260908)   # not printed in v28d; asserted against fusion_campi_acquisition_components.tsv in stat_fusion
DONOR_ENTRAPMENT_SEEDS = (20260910, 20260911, 20260912)  # not printed in v28d; asserted against donor_recount.tsv in stat_donor
FDRBENCH_VERSION = "0.0.4"         # v28d line 727
SIM_BOOTSTRAP_DRAWS = 20_000       # v28d line 733; LAB_ROOT/analyses/fdr_reference_benchmark_2026-09-10/collect.py line 22
SIM_BOOTSTRAP_SEED = 2026091051    # not printed in v28d; collect.py line 70
SIM_CI = (0.025, 0.975)            # 95% percentile interval, v28d line 733
JOINT_PERMUTATION_SEED = 0         # not printed in v28d; iteration2 joint_search extract_fileq.py evidence(seed=0)
CROSSFIT_FOLDS = 10                # v28d line 773; asserted against crossfit_pilot_folds.tsv and crossfit_fold_caches.tsv
CROSSFIT_FOLD_KEY = "FastaLake-D11-crossfit-20261004:<native scan>"  # v28d line 773; asserted against crossfit_fold_caches.tsv
CROSSFIT_BOOTSTRAP_DRAWS = 10_000  # v28d line 777
CROSSFIT_BOOTSTRAP_SEED = 20261004  # v28d line 777
CROSSFIT_CI = (2.5, 97.5)          # percentile interval, v28d line 777
PERTURBATION_MOVE = 0.5            # percentage points, v28d line 777
CROSSFIT_REGIONS = ("duodenum", "jejunum", "ileum", "colon", "stool")  # order of the region bootstrap stream

# Design values printed in lines 661-777 that no source_data table of this unit records. They are
# listed for the reader only; no code uses them and nothing here checks them (OPEN_ITEMS.md).
#   coverage -c 0.8 --cov-mode 0 --cluster-mode 2 (line 695); MSFragger 4.1 + Percolator 3.06.5 (line 699);
#   Sage v0.14.7 with five candidates per spectrum (line 703); shuffle seed 20260925 (line 721).


def load(sd, name, **kw):
    """Read one source_data table (tab-separated).

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.
    name : str
        File name of the table.
    **kw
        Passed to ``pandas.read_csv``.

    Returns
    -------
    DataFrame
    """
    return pd.read_csv(Path(sd) / name, sep="\t", **kw)


# ---------------------------------------------------------------- lines 665-667
def stat_scramble(sd):
    """Terminal-preserving and all-position scramble workflows and lookup fractions.

    Unit: workflow (one representative healthy-stool acquisition; one authentic and ten
    shuffled workflows, five seeds per shuffle type) and distinct query (lookup fractions).

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.

    Returns
    -------
    dict
        selected-FASTA protein ranges, accepted peptide counts, terminal-preserving query-match
        fraction ranges (percent) at lengths 9 and 10 and the largest number of matching
        terminal-preserving queries at lengths >= 13.
    """
    e = load(sd, "scramble_end_to_end.tsv")
    p = e[(e.gate == "q_only") & (e.q == Q_CUTOFF) & (e.min_length == 5) & (e.max_length == 50)]
    allp = p[p.arm == "all_position"]
    term = p[p.arm == "terminal"]
    auth = p[p.arm == "authentic"]
    c = load(sd, "scramble_lookup_comparison.tsv")
    c = c[c.arm == "C_TERMINAL_PRESERVED"]
    f = c.assign(frac=100 * c.matched_queries / c.distinct_queries_terminal)
    return {
        "scr_all_min": int(allp.parsimony_proteins.min()), "scr_all_max": int(allp.parsimony_proteins.max()),
        "scr_term_min": int(term.parsimony_proteins.min()), "scr_term_max": int(term.parsimony_proteins.max()),
        "scr_auth_proteins": int(auth.parsimony_proteins.iloc[0]),
        "scr_shuffled_peptides_max": int(p[p.arm != "authentic"].n_peptides.max()),
        "scr_auth_peptides": int(auth.n_peptides.iloc[0]),
        "scr_term9_min": float(f[f.length == 9].frac.min()), "scr_term9_max": float(f[f.length == 9].frac.max()),
        "scr_term10_min": float(f[f.length == 10].frac.min()), "scr_term10_max": float(f[f.length == 10].frac.max()),
        "scr_term_ge13_matched_max": int(c[c.length >= 13].matched_queries.max()),
        "scr_n_seeds": int(c.seed.nunique()),
    }


# ---------------------------------------------------------------- lines 665, 673-675, 713, 721
def stat_database_sizes(sd):
    """Selected-FASTA counts and reservoir sizes.

    Unit: acquisition (726: 15 CAMPI, 58 healthy stool, 306 PREDICT, 347 CapScan).

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.

    Returns
    -------
    dict
        number of selected FASTAs, human Swiss-Prot entries, reservoir entry counts.
    """
    a = load(sd, "acquisition_sizes_726.tsv")
    h = load(sd, "human_reference.tsv")
    res = a.groupby("cohort").reservoir_entries.unique()
    assert all(len(v) == 1 for v in res)
    return {
        "sizes_n_fastas": int(len(a)),
        "sizes_human_entries": int(h.entries.iloc[0]),
        "reservoir_public": int(res["HSTOOL"][0]),
        "reservoir_campi": int(res["CAMPI"][0]),
        "sizes_n_by_cohort": a.cohort.value_counts().to_dict(),
    }


# ---------------------------------------------------------------- line 679
def stat_selection_grid(sd):
    """Selection-sensitivity grid.

    Unit: acquisition-setting observation (58 healthy-stool acquisitions x 19 settings).

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.

    Returns
    -------
    dict
        acquisitions, settings, length range, score-fraction range (percent), observations.
    """
    g = load(sd, "selection_grid_acquisitions.tsv")
    ln = g[g.kind == "length"]
    assert (int(ln.min_length.min()), int(ln.min_length.max())) == SELECTION_GRID_LENGTHS
    assert (g.loc[g.kind.isin(["length", "reference"]), "top_percent"] == LOOKUP_TOP_FRACTION).all()
    lengths = g.loc[g.kind.isin(["length", "reference"]), "min_length"]
    scores = 100 * g.loc[g.kind == "score", "top_percent"]
    return {
        "grid_acquisitions": int(g["sample"].nunique()), "grid_settings": int(g.arm.nunique()),
        "grid_len_min": int(lengths.min()), "grid_len_max": int(lengths.max()),
        "grid_score_min": float(scores.min()), "grid_score_max": float(scores.max()),
        "grid_observations": int(len(g)),
        "grid_reference_top_percent": float(g.loc[g.kind == "reference", "top_percent"].iloc[0]),
    }


# ---------------------------------------------------------------- lines 683-687
def stat_campi(sd):
    """CAMPI acquisitions, community members, pooled reference peptides and depth draw.

    Unit: acquisition (15) and distinct peptide (pooled reference set).

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.

    Returns
    -------
    dict
        acquisitions, members, pooled distinct reference peptides (5-50 residues, union over
        acquisitions, sequences as listed in the per-acquisition frequency tables), smallest
        MS2 count.
    """
    inv = load(sd, "campi_spectral_inventory.tsv")
    ref = load(sd, "campi_reference_peptides.tsv", dtype={"peptide": str})
    mem = load(sd, "campi_member_summary.tsv")
    rec = load(sd, "campi_reference_fasta_records.tsv", dtype=str, keep_default_na=False)
    return {
        "campi_contaminants": int(rec.loc[rec.header_kind == "other", "header_token"].nunique()),
        "campi_reference_records": int(len(rec)),
        "campi_acquisitions": int(len(inv)),
        "campi_members": int(mem.member.nunique()),
        "campi_ref_pooled": int(ref.peptide.nunique()),
        "campi_ref_pooled_IL": int(ref.peptide.str.replace("I", "L").nunique()),
        "campi_min_ms2": int(inv.ms2.min()),
    }


# ---------------------------------------------------------------- line 695
def stat_clustering(sd):
    """Candidate-clustering configurations.

    Unit: acquisition-arm configuration (15 CAMPI acquisitions x 7 arms).

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.

    Returns
    -------
    dict
        named values (keys as used in NUMBER_SPEC and stats_rows).
    """
    t = load(sd, "clustering_per_acquisition_identifications.tsv")
    assert sorted(t.loc[t.arm != "unclustered", "identity"].unique(), reverse=True) == list(CLUSTER_IDENTITIES)
    return {"cluster_configs": int(len(t)), "cluster_arms": int(t.arm.nunique()),
            "cluster_acquisitions": int(t["sample"].nunique())}


# ---------------------------------------------------------------- line 715
def stat_fusion(sd):
    """CAMPI target-entrapment fusion diagnostic at peptide q <= 0.01.

    Unit: acquisition-level counts summed over the 15 CAMPI acquisitions, per seed (2 seeds).
    Diagnostic = (E + E0 + 2B + U) / N; paired = (E + E0 + 2B) / N; unresolved = U / N.

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.

    Returns
    -------
    dict
        per seed and length window, the three components in percent and the unresolved share
        of the total in percent; plus their ranges across seeds.
    """
    t = load(sd, "fusion_campi_acquisition_components.tsv")
    assert tuple(sorted(int(x) for x in t.seed.unique())) == FUSION_SEEDS
    t = t[t.q_threshold == Q_CUTOFF]
    out, rows = {}, []
    for (seed, lo, hi), g in t.groupby(["seed", "min_residues", "max_residues"]):
        n, pn, u = g.n_canonical.sum(), g.paired_numerator.sum(), g.n_unpaired_accepted.sum()
        assert g["sample"].nunique() == 15
        rows.append(dict(seed=seed, lo=lo, hi=hi, total=100 * (pn + u) / n, paired=100 * pn / n,
                         unresolved=100 * u / n, share=100 * u / (pn + u)))
    r = pd.DataFrame(rows)
    for (lo, hi), g in r.groupby(["lo", "hi"]):
        for k in ("total", "paired", "unresolved", "share"):
            out[f"fus_{lo}_{hi}_{k}_min"] = float(g[k].min())
            out[f"fus_{lo}_{hi}_{k}_max"] = float(g[k].max())
    out["fus_table"] = r
    out["fus_seeds"] = sorted(int(s) for s in t.seed.unique())
    return out


# ---------------------------------------------------------------- line 721
def stat_label_validity(sd):
    """Reservoir occurrence of accepted control sequences.

    Unit: unique sequence (I/L-equivalent), deduplicated across acquisitions and seeds.

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.

    Returns
    -------
    dict
        named values (keys as used in NUMBER_SPEC and stats_rows).
    """
    p = load(sd, "label_validity_presence.tsv", dtype={"peptide": str})
    p["length"] = p.peptide.str.len()
    out = {}
    for arm, key in (("ENTRAP", "ctl"), ("SHUFCTL", "shuf"), ("REALPOS", "pos")):
        for (lo, hi) in (SHORT_LENGTH, NESTED_LENGTH, FULL_LENGTH):
            g = p[(p.arm == arm) & p.length.between(lo, hi)]
            out[f"lab_{key}_{lo}_{hi}_present"] = int(g.present.sum())
            out[f"lab_{key}_{lo}_{hi}_total"] = int(len(g))
            out[f"lab_{key}_{lo}_{hi}_pct"] = float(100 * g.present.mean())
    r = load(sd, "label_validity_result.tsv")
    out["lab_reservoir"] = int(r.reservoir_sequences.iloc[0])
    return out


# ---------------------------------------------------------------- lines 727-729
def stat_donor(sd):
    """Donor-selected explicit-peptide entrapment: roster, seeds and score ties.

    Unit: target acquisition (58 donor-selected, 5 self-selected), 3 entrapment seeds.

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.

    Returns
    -------
    dict
        targets per tier, seeds, exact score ties at q <= 0.01, and the median over target
        acquisitions of the seed-averaged paired estimate (E + E0 + 2B)/(T + E) per length.
    """
    d = load(sd, "donor_recount.tsv")
    assert tuple(sorted(int(x) for x in d.seed.unique())) == DONOR_ENTRAPMENT_SEEDS
    out = {
        "donor_targets": int(d.loc[d.tier == "T1_within_group", "target"].nunique()),
        "self_targets": int(d.loc[d.tier == "T0_self", "target"].nunique()),
        "donor_seeds": int(d.seed.nunique()),
        "donor_seed_values": sorted(int(s) for s in d.seed.unique()),
        "donor_ties_q01": int(d.loc[d.q_threshold == Q_CUTOFF, "exact_score_ties"].sum()),
        "donor_ties_all": int(d.exact_score_ties.sum()),
    }
    q = d[d.q_threshold == Q_CUTOFF]
    for tier, key in (("T1_within_group", "donor"), ("T0_self", "self")):
        for lo in (5, 9):
            m = q[(q.tier == tier) & (q.min_length == lo)].groupby("target").paired_estimate.mean()
            out[f"{key}_paired_median_{lo}"] = float(100 * m.median())
            out[f"{key}_paired_n_{lo}"] = int(len(m))
    return out


# ---------------------------------------------------------------- line 733
def _sim_summary(values, rng):
    """Mean and percentile bootstrap of the mean, as in the original collect.py summarize().

    Parameters
    ----------
    values : array-like of float
        One value per synthetic mixture.
    rng : numpy.random.Generator
        The shared stream (consumed: SIM_BOOTSTRAP_DRAWS x len(values) integers).

    Returns
    -------
    tuple of float
        mean, 2.5% quantile and 97.5% quantile of the bootstrap means.
    """
    values = np.asarray(values, dtype=float)
    n = len(values)
    draws = values[rng.integers(n, size=(SIM_BOOTSTRAP_DRAWS, n))].mean(axis=1)
    return float(values.mean()), float(np.quantile(draws, SIM_CI[0])), float(np.quantile(draws, SIM_CI[1]))


def stat_simulation(sd):
    """Restricted-reference simulation: design counts, mean FDP and bootstrap intervals.

    Unit: synthetic mixture (12 baseline, 12 stress). The bootstrap resamples whole mixtures
    (20,000 draws, numpy default_rng seed 2026091051, one stream consumed in the original
    order: baseline fixed, baseline independent, baseline self, stress fixed, stress
    independent, stress self, then the paired self-minus-independent differences for
    baseline and stress).

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.

    Returns
    -------
    dict
        design values, outcome counts, means (percent), 95% bootstrap intervals (percent)
        and paired differences (percentage points).
    """
    o = load(sd, "simulation_all_outcomes.tsv")
    des = load(sd, "simulation_design.tsv").set_index(["parameter", "regime"]).value
    assert des[("bootstrap_draws", "all")] == SIM_BOOTSTRAP_DRAWS and des[("bootstrap_seed", "all")] == SIM_BOOTSTRAP_SEED
    ordinary = o[o.arm.isin(["fixed", "self", "independent"])]
    rng = np.random.default_rng(SIM_BOOTSTRAP_SEED)
    out = {}
    for (regime, arm), f in ordinary.groupby(["regime", "arm"], sort=True):
        m, lo, hi = _sim_summary(f.realised_FDP, rng)
        out[f"sim_{regime}_{arm}_mean"] = 100 * m
        out[f"sim_{regime}_{arm}_lo"] = 100 * lo
        out[f"sim_{regime}_{arm}_hi"] = 100 * hi
        out[f"sim_{regime}_{arm}_n"] = int(len(f))
    paired = ordinary.pivot(index=["regime", "unit_id"], columns="arm", values="realised_FDP")
    for regime, f in paired.groupby(level="regime", sort=True):
        m, lo, hi = _sim_summary(f.self - f.independent, rng)
        out[f"sim_{regime}_self_minus_indep"] = 100 * m
        out[f"sim_{regime}_self_minus_indep_lo"] = 100 * lo
        out[f"sim_{regime}_self_minus_indep_hi"] = 100 * hi
        out[f"sim_{regime}_self_minus_fixed"] = float(100 * (f.self - f.fixed).mean())
    n = des[("spectra_per_acquisition", "all")]
    out.update({
        "sim_spectra": int(n),
        "sim_bootstrap_draws": int(des[("bootstrap_draws", "all")]),
        "sim_units_baseline": int(o.loc[o.regime == "baseline", "unit_id"].nunique()),
        "sim_units_stress": int(o.loc[o.regime == "stress", "unit_id"].nunique()),
        "sim_len_min": int(des[("generating_min_length", "all")]),
        "sim_len_max": int(des[("generating_max_length", "all")]),
        "sim_reference_proteins": int(des[("reference_proteins", "all")]),
        # the text says "approximately 30,000"; the exact record count is sim_reference_proteins
        "sim_reference_proteins_approx": int(round(des[("reference_proteins", "all")], -3)),
        "sim_ordinary_outcomes": int(len(ordinary)),
        "sim_all_outcomes": int(len(o)),
        "sim_absent_stress_pct": float(100 * (n - des[("present_spectra", "stress")]) / n),
        "sim_absent_baseline_pct": float(100 * (n - des[("present_spectra", "baseline")]) / n),
    })
    return out


# ---------------------------------------------------------------- lines 739-745
def stat_transfer_pool(sd):
    """Transfer donors and shuffle-only acceptance with pooled databases.

    Unit: target acquisition (58 transfer targets; 10 inter-individual) and
    peptide-acquisition observation (12 technical replicates, q <= 0.01, 5-50 residues).

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.

    Returns
    -------
    dict
        named values (keys as used in NUMBER_SPEC and stats_rows).
    """
    t = load(sd, "transfer_pool_comparisons.tsv")
    tr = t[(t.kind == "transfer") & (t.min_length == 5)]
    inter = tr[tr.group == "Interindividual"].pivot(index="sample", columns="arm", values="donor")
    p = load(sd, "hstool_pooling_per_acquisition.tsv")
    p = p[(p.q_threshold == Q_CUTOFF) & (p.min_residues == 5)]
    out = {
        "transfer_targets": int(tr["sample"].nunique()),
        "inter_targets": int(len(inter)),
        "inter_donors_differ": int((inter.T1_within_group != inter.T4_different_donor).sum()),
        "inter_donors_same": int((inter.T1_within_group == inter.T4_different_donor).sum()),
        "pool_replicates": int(p["sample"].nunique()),
    }
    for arm in ("own", "pool_other11"):
        g = p[p.arm == arm]
        out[f"pool_{arm}_shuffle_only"] = int(g.shuffle_only_not_in_original_reservoir.sum())
        out[f"pool_{arm}_observations"] = int(g.n_canonical.sum())
    return out


# ---------------------------------------------------------------- line 751
def stat_common_space(sd):
    """Common healthy-stool FASTA union and representative presence.

    Unit: protein sequence (union) and quantified peptide feature counted per acquisition
    (58 acquisitions, pooled numerators and denominators).

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.

    Returns
    -------
    dict
        named values (keys as used in NUMBER_SPEC and stats_rows).
    """
    u = load(sd, "common_fasta_union_summary.tsv")
    r = load(sd, "representative_presence.tsv")
    c = load(sd, "union_contributions.tsv").groupby("roster")
    unions = c.new_unique_sequences.sum()
    assert unions["HSTOOL58"] == u.unique_sequences.iloc[0]
    return {
        "union_PREDICT306": int(unions["PREDICT306"]), "union_CAPSCAN331": int(unions["CAPSCAN331"]),
        "union_CAPSCAN347": int(unions["CAPSCAN347"]), "union_HSTOOL58_recount": int(unions["HSTOOL58"]),
        "union_sources": int(u.source_fastas.iloc[0]),
        "union_sequences": int(u.unique_sequences.iloc[0]),
        "outside_pct": float(100 * r.outside_local_fasta_features.sum() / r.retained_lfq_features.sum()),
        "outside_acquisitions": int(len(r)),
    }


# ---------------------------------------------------------------- lines 755-761
def stat_spectrum_split(sd, table="split_half_acquisitions_710.tsv"):
    """Spectrum-split roster, participant mapping, reconverted inputs, fallbacks, zero yields.

    Unit: acquisition (710 = 58 + 15 + 306 + 331; the 726 table has 347 CapScan) and
    participant (93 PREDICT patients, 15 CapScan participants).

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.
    table : str
        ``split_half_acquisitions_710.tsv`` (v28d roster, default) or ``split_half_acquisitions_726.tsv``.

    Returns
    -------
    dict
        named values (keys as used in NUMBER_SPEC and stats_rows).
    """
    a = load(sd, table, dtype={"participant": str}, keep_default_na=False)
    a5 = a[a.min_length == 5]
    n = a5.cohort.value_counts().to_dict()
    out = {f"split_n_{k}": int(v) for k, v in n.items()}
    pr = a5[a5.cohort == "PREDICT"]
    cs = a5[a5.cohort == "CAPSCAN"]
    out.update({
        "split_predict_participants": int(pr.loc[pr.participant != "", "participant"].nunique()),
        "split_predict_mapped": int((pr.participant != "").sum()),
        "split_predict_missing": int((pr.participant == "").sum()),
        "split_capscan_participants": int(cs.loc[cs.participant != "", "participant"].nunique()),
        "split_capscan_mapped": int((cs.participant != "").sum()),
        "split_reconverted": int((pr.conversion_recovered.astype(str) == "True").sum()),
        "split_capscan_zero_disjoint": int((cs.disjoint.astype(int) == 0).sum()),
        "split_n_total": int(len(a5)),
    })
    if table == "split_half_acquisitions_710.tsv":
        f = load(sd, "split_half_model_fallbacks.tsv")
        z = load(sd, "split_half_zero_yield.tsv")
        out.update({
            "split_fallback_searches": int(len(f)),
            "split_fallback_acquisitions": int(f["index"].nunique()),
            "split_zero_min_q": float(z.loc[z.arm == "HALF_DISJOINT", "minimum_target_peptide_q"].iloc[0]),
            "split_zero_accepted": int(z.loc[z.arm == "HALF_DISJOINT", "accepted_peptides"].iloc[0]),
        })
    return out


# ---------------------------------------------------------------- lines 765-769
def stat_joint(sd):
    """Joint healthy-stool search: per-acquisition confidence and quantitative agreement.

    Unit: acquisition (58; peptide gains) and acquisition pair (66 technical-replicate pairs
    from 12 acquisitions, 105 independent-preparation pairs from 15 acquisitions).

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.

    Returns
    -------
    dict
        median gains in percent (joint study-level D/T estimator and Sage q, each over the
        acquisition-specific Sage count), permutation-control maxima, and medians across pairs
        of the per-pair median CV (percent) and log2 Pearson r on matched microbial groups.
    """
    c = load(sd, "joint_counts_per_acquisition.tsv")
    own = c[c.arm == "own"].set_index("filename")
    jo = c[c.arm == "joint"].set_index("filename")
    assert len(own) == len(jo) == 58 and set(own.index) == set(jo.index)
    jo = jo.loc[own.index]
    cov = load(sd, "joint_fileq_coverage.tsv")
    out = {
        "joint_acquisitions": int(len(own)),
        "joint_glob_gain": float(100 * (np.median(jo.glob_plain / own.sage) - 1)),
        "joint_sage_gain": float(100 * (np.median(jo.sage / own.sage) - 1)),
        "joint_perm_max": int(cov.negctl_permuted_label_accepted.max()),
        "joint_plus1_max": int(jo.negctl_file_plus1.max()),
    }
    out["joint_plus1_word"] = "none" if out["joint_plus1_max"] == 0 else str(out["joint_plus1_max"])
    p = load(sd, "microbial_group_pairs.tsv")
    p = p[p.definition == "matched"]
    for comp, key in (("original_vs_joint_fileq", "fileq"), ("original_vs_joint_corrected", "study")):
        for tier, tk in (("TechnicalReplicates", "tech"), ("WorkflowReplicates", "prep")):
            g = p[(p.comparison == comp) & (p.tier == tier)]
            for method in g.method.unique():
                m = "orig" if method.startswith("original") else "joint"
                h = g[g.method == method]
                out[f"jq_{key}_{tk}_{m}_cv"] = float(h.median_pair_cv_pct.median())
                out[f"jq_{key}_{tk}_{m}_r"] = float(h.log2_pearson.median())
                out[f"jq_{key}_{tk}_pairs"] = int(len(h))
                out[f"jq_{key}_{tk}_acq"] = int(len(set(h.sample_a) | set(h.sample_b)))
    a = load(sd, "technical_agreement_summary.tsv")
    a = a[(a.method == "fixed_joint_corrected") & (a.threshold == "both_resolved")].set_index("tier")
    out.update({
        "lead_tech_orig": float(100 * a.loc["TechnicalReplicates", "original_agreement"]),
        "lead_tech_joint": float(100 * a.loc["TechnicalReplicates", "test_agreement"]),
        "lead_prep_orig": float(100 * a.loc["WorkflowReplicates", "original_agreement"]),
        "lead_prep_joint": float(100 * a.loc["WorkflowReplicates", "test_agreement"]),
    })
    return out


# ---------------------------------------------------------------- lines 773-777
def _strata(a):
    """Cross-fit summary strata in the order of the original bootstrap stream.

    Parameters
    ----------
    a : DataFrame
        crossfit_per_acquisition.tsv.

    Returns
    -------
    list of (str, boolean Series)
        PREDICT, CapScan all, CapScan stool, CapScan intestinal fluid (non-stool regions).
    """
    cap = a.cohort == "CAPSCAN"
    return [("PREDICT", a.cohort == "PREDICT"), ("CAPSCAN_all", cap),
            ("CAPSCAN_stool", cap & (a.group == "stool")), ("CAPSCAN_intestinal", cap & (a.group != "stool"))]


def crossfit_bootstrap(a, strata, columns, seed=CROSSFIT_BOOTSTRAP_SEED):
    """Participant cluster bootstrap of the acquisition median, as in the original analyze_final.py.

    One numpy default_rng stream runs through the strata in the given order; within a
    stratum the same 10,000 resampled participant sets are used for every column.
    Acquisitions without a participant identifier are their own cluster. Clusters are taken
    in the order of ``cluster_order`` (rank of the original participant identifier within its
    cohort), which is the order the original analysis used; with pseudonyms alone the draws
    differ. The column discloses only that rank, not the identifiers.

    Parameters
    ----------
    a : DataFrame
        crossfit_per_acquisition.tsv.
    strata : list of (name, boolean mask)
    columns : list of str
        Columns to summarise; the stream does not depend on this list.

    Returns
    -------
    DataFrame
        stratum, column, n acquisitions, n participants, median, 2.5 and 97.5 percentiles.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for name, m in strata:
        g = a[m].reset_index(drop=True)
        subj = g.cluster_order.values if "cluster_order" in g else g.subject_id.values
        us = np.unique(subj)
        members = [np.flatnonzero(subj == s) for s in us]
        draws = [np.concatenate([members[k] for k in rng.integers(0, len(us), len(us))])
                 for _ in range(CROSSFIT_BOOTSTRAP_DRAWS)]
        for col in columns:
            v = g[col].values.astype(float)
            bs = np.array([np.nanmedian(v[ix]) for ix in draws])
            rows.append(dict(stratum=name, column=col, n_acquisitions=len(g), n_participants=len(us),
                             median=float(np.nanmedian(v)), ci_low=float(np.percentile(bs, CROSSFIT_CI[0])),
                             ci_high=float(np.percentile(bs, CROSSFIT_CI[1]))))
    return pd.DataFrame(rows)


def stat_crossfit(sd, bootstrap=True):
    """10-fold spectrum cross-fitting: checks, sensitivities and gains.

    Unit: acquisition (578 matched: 296 PREDICT, 282 CapScan = 56 stool + 226 intestinal
    fluid; out-of-fold tables 306 PREDICT and 331 CapScan; 8 acquisitions in the independent
    closed-form check); intervals resample participants (10,000 draws, seed 20261004).

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.
    bootstrap : bool
        Also run the participant cluster bootstraps (``cf_boot``, ``cf_boot_regions``).

    Returns
    -------
    dict
        named values (keys as used in NUMBER_SPEC and stats_rows).
    """
    a = load(sd, "crossfit_per_acquisition.tsv", dtype={"subject_id": str})
    a["lda_fallback"] = a.lda_fallback.fillna("")
    st = _strata(a)
    med = lambda m, col: float(np.nanmedian(a.loc[m, col]))
    out = {"cf_acquisitions": int(len(a)), "cf_n_predict": int((a.cohort == "PREDICT").sum()),
           "cf_n_capscan": int((a.cohort == "CAPSCAN").sum())}
    # pilot folds, q re-implementation, K = 1 control
    f = load(sd, "crossfit_pilot_folds.tsv")
    assert sorted(f.fold.unique()) == list(range(CROSSFIT_FOLDS))
    f["frac"] = 100 * f.fold_ms2 / f.groupby("sample").fold_ms2.transform("sum")
    # out-of-fold prediction tables written for each cohort and fold (line 773)
    oc = load(sd, "crossfit_fold_caches.tsv")
    assert (oc.K == CROSSFIT_FOLDS).all() and (oc.fold_seed + ":<native scan>" == CROSSFIT_FOLD_KEY).all()
    for c, g in oc.groupby("cohort"):
        assert sorted(g.fold.unique()) == list(range(CROSSFIT_FOLDS))
        per_fold = g.groupby("fold")["sample"].apply(frozenset)
        assert per_fold.nunique() == 1 and (g.groupby("fold").n_acquisitions.first() == g.groupby("fold").size()).all()
        out[f"cf_oof_{c}"] = int(g["sample"].nunique())
    # independent closed-form discriminant (line 777): acquisitions without fallback in any arm whose
    # closed-form accepted sets equal the producer's at both lengths in all three arms
    ic = load(sd, "crossfit_independent_lda_check.tsv")
    ok = ic.assign(same=(ic.closed_form_accepted_5 == ic.producer_accepted_5) & (ic.closed_form_accepted_9 == ic.producer_accepted_9)
                   & (ic.closed_form_jaccard == 1.0), lda=ic.producer_rescoring == "lda").groupby("i")[["same", "lda"]].all()
    out["cf_closed_form_checked"] = int(len(ok))
    out["cf_closed_form_no_fallback"] = int(ok.lda.sum())
    out["cf_closed_form_exact"] = int((ok.lda & ok.same).sum())
    v = load(sd, "crossfit_q_validation.tsv")
    k = load(sd, "crossfit_k1_control.tsv")
    p252 = k[(k.arm == "MG") & k["sample"].str.endswith("MicrobPred_252")].set_index("min_length")
    out.update({
        "cf_pilots": int(f["sample"].nunique()), "cf_fold_min": float(f.frac.min()), "cf_fold_max": float(f.frac.max()),
        "cf_q_searches": int(len(v)), "cf_q_equal": int(v.q_reimplementation_equal.sum()),
        "cf_rescoring_equal": int(v.full_rescoring_equal.sum()),
        "cf_k1_selection_equal": int(k.drop_duplicates("sample").selection_equals_archived_fasta.sum()),
        "cf_k1_selection_n": int(k["sample"].nunique()),
        "cf_k1_sage_equal": int(k.sage_equals_frozen.sum()), "cf_k1_comparisons": int(len(k)),
        "cf_k1_pooled_equal": int(k.pooled_rescoring_equals_frozen.sum()),
        "cf_k1_pooled_differ": int((~k.pooled_rescoring_equals_frozen.astype(bool)).sum()),
        "cf_rescoring_differ": int((~v.full_rescoring_equal.astype(bool)).sum()),
        "cf_k1_252_pooled_5": int(p252.loc[5, "pooled_rescoring_k1"]), "cf_k1_252_frozen_5": int(p252.loc[5, "frozen"]),
        "cf_k1_252_pooled_9": int(p252.loc[9, "pooled_rescoring_k1"]), "cf_k1_252_frozen_9": int(p252.loc[9, "frozen"]),
    })
    # calibration, perturbation, union sensitivity
    calib, pert, uni, reg = [], [], [], []
    for name, m in st:
        for L in ("5", "9"):
            calib.append(abs(med(m, f"calib_full_minus_frozen_{L}")))
            base = med(m, f"gain_cf_pooled_{L}")
            pert += [abs(med(m, f"gain_cf_pert{s}_{L}") - base) for s in (1, 2)]
            uni.append(abs(med(m, f"gain_cf_union_{L}") - base))
    cap = a.cohort == "CAPSCAN"
    for g in CROSSFIT_REGIONS[:4]:
        m = cap & (a.group == g)
        for L in ("5", "9"):
            reg.append(abs(med(m, f"gain_cf_union_{L}") - med(m, f"gain_cf_pooled_{L}")))
    moved = [(abs(a[f"gain_cf_pert{s}_5"] - a.gain_cf_pooled_5) > PERTURBATION_MOVE).sum() for s in (1, 2)]
    out.update({
        "cf_calib_max": max(calib), "cf_pert_shift_max": max(pert),
        "cf_moved_pert1": int(moved[0]), "cf_moved_pert2": int(moved[1]),
        "cf_moved_pert1_pct": float(100 * moved[0] / len(a)), "cf_moved_pert2_pct": float(100 * moved[1] / len(a)),
        "cf_union_shift_max": max(uni), "cf_union_region_shift_max": max(reg),
        "cf_mg_union_zero_capscan": int(((a.cohort == "CAPSCAN") & (a.MG_union_5 == 0)).sum()),
    })
    # fallback acquisitions
    fb = a[a.lda_fallback != ""]
    out.update({"cf_fallback_n": int(len(fb)), "cf_fallback_predict": int((fb.cohort == "PREDICT").sum()),
                "cf_fallback_capscan": int((fb.cohort == "CAPSCAN").sum()),
                "cf_fallback_capscan_group": ";".join(fb.loc[fb.cohort == "CAPSCAN", "group"])})
    lf = load(sd, "crossfit_lda_fallback_acquisitions.tsv").set_index("i")
    for i in lf.index:
        out[f"cf_fb_{i}_producer"] = float(lf.loc[i, "gain_cf_pooled_5"])
        out[f"cf_fb_{i}_independent"] = float(lf.loc[i, "gain_independent_refit_5"])
    pr = a.cohort == "PREDICT"
    for L in ("5", "9"):
        s = a.set_index("i")[f"gain_cf_pooled_{L}"].copy()
        out[f"cf_predict_median_{L}"] = float(np.nanmedian(s[pr.values]))
        s.loc[lf.index] = lf[f"gain_independent_refit_{L}"]
        out[f"cf_predict_median_{L}_replaced"] = float(np.nanmedian(s[pr.values]))
        out[f"cf_capscan_median_{L}_replaced"] = float(np.nanmedian(s[cap.values]))
        out[f"cf_capscan_median_{L}"] = med(cap, f"gain_cf_pooled_{L}")
    # gains and counts printed in line 777
    for name, m in st:
        out[f"cf_{name}_n"] = int(m.sum())
        out[f"cf_{name}_pos9"] = int((a.loc[m, "gain_cf_pooled_9"] > 0).sum())
        out[f"cf_{name}_cf_vs_full_5"] = med(m, "cf_vs_full_pooled_5")
        out[f"cf_{name}_odd_disjoint_5"] = med(m, "gain_odd_disjoint_5")
        out[f"cf_{name}_odd_self_5"] = med(m, "gain_odd_self_5")
        out[f"cf_{name}_fewer_5"] = -med(m, "cf_vs_full_pooled_5")
        for L in ("5", "9"):
            out[f"cf_{name}_gain_{L}"] = med(m, f"gain_cf_pooled_{L}")
    if bootstrap:
        cols = ["gain_cf_pooled_5", "gain_cf_pooled_9", "gain_odd_self_5", "gain_odd_disjoint_5"]
        out["cf_boot"] = crossfit_bootstrap(a, st, cols)
        regs = [(f"CAPSCAN_{g}", cap & (a.group == g)) for g in CROSSFIT_REGIONS]
        out["cf_boot_regions"] = crossfit_bootstrap(a, regs, ["gain_cf_pooled_9"])
        b = out["cf_boot"].set_index(["stratum", "column"])
        r = out["cf_boot_regions"].set_index(["stratum", "column"])
        out["cf_stool9_lo"] = b.loc[("CAPSCAN_stool", "gain_cf_pooled_9"), "ci_low"]
        out["cf_stool9_hi"] = b.loc[("CAPSCAN_stool", "gain_cf_pooled_9"), "ci_high"]
        out["cf_stool9_region_lo"] = r.loc[("CAPSCAN_stool", "gain_cf_pooled_9"), "ci_low"]
        out["cf_stool9_region_hi"] = r.loc[("CAPSCAN_stool", "gain_cf_pooled_9"), "ci_high"]
    return out


def stat_same_spectrum(sd):
    """Same-spectrum (all-prediction selection) benchmark at 9-50 and 5-8 residues and minimum length 7.

    Unit: acquisition (296 PREDICT, 282 CapScan matched); gain = 100 x (FL / MG - 1), median
    over acquisitions; q <= 0.01, q-only gate.

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.

    Returns
    -------
    dict
        named values (keys as used in NUMBER_SPEC and stats_rows).
    """
    o = load(sd, "confidence_ordinary_all_pairs.tsv")
    o = o[(o.gate == "q_only") & (o.q == Q_CUTOFF)]
    out = {}
    for c in ("PREDICT", "CAPSCAN"):
        for lo, hi in (NESTED_LENGTH, SHORT_LENGTH, FULL_LENGTH):
            g = o[(o.cohort == c) & (o.min_length == lo) & (o.max_length == hi)]
            out[f"ss_{c}_{lo}_{hi}_median"] = float(g.gain_pct.median())
            out[f"ss_{c}_{lo}_{hi}_pos"] = int((g.n_fl > g.n_mg).sum())
            out[f"ss_{c}_{lo}_{hi}_n"] = int(len(g))
    m = load(sd, "minlen7_benchmark_per_acquisition.tsv")
    m = m[(m.population == "fresh7") & (m.arm == "FL_public") & (m.comparator == "MG")]
    for c in ("PREDICT", "CapScan", "CAPSCAN"):
        g = m[m.cohort == c]
        if len(g):
            key = c.upper()
            out[f"ml7_{key}_gain"] = float(g.gain_pct.median())
            out[f"ml7_{key}_recovery"] = float(g.comparator_identity_recovery_pct.median())
            out[f"ml7_{key}_n"] = int(len(g))
    return out


# ---------------------------------------------------------------- line 633
def stat_confidence_rosters(sd):
    """Additional-confidence rosters: matched acquisitions, CAMPI fusion searches, zero yields.

    Unit: acquisition (296 PREDICT, 282 CapScan) and CAMPI fusion search (15 acquisitions x 2 seeds).

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.

    Returns
    -------
    dict
        named values (keys as used in NUMBER_SPEC and stats_rows).
    """
    o = load(sd, "confidence_ordinary_all_pairs.tsv")
    fz = load(sd, "confidence_fusion_all.tsv")
    z = load(sd, "confidence_pep01_zero_yield_cases.tsv")
    return {
        "conf_predict": int(o.loc[o.cohort == "PREDICT", "sample"].nunique()),
        "conf_capscan": int(o.loc[o.cohort == "CAPSCAN", "sample"].nunique()),
        "conf_fusion_searches": int(fz[["sample", "seed"]].drop_duplicates().shape[0]),
        "conf_zero_yield_predict": int((z.cohort == "PREDICT").sum()),
    }


# ---------------------------------------------------------------- all
def compute_all(sd, bootstrap=True):
    """Run every statistic of the unit.

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.
    bootstrap : bool
        Run the cross-fit participant bootstraps.

    Returns
    -------
    dict
        One flat dict of named values (tables kept under their keys); values of the 726-acquisition
        spectrum split carry the suffix ``_726``.
    """
    out = {}
    for f in (stat_scramble, stat_database_sizes, stat_selection_grid, stat_campi, stat_clustering, stat_fusion,
              stat_label_validity, stat_donor, stat_simulation, stat_transfer_pool, stat_common_space,
              stat_spectrum_split, stat_joint, stat_same_spectrum, stat_confidence_rosters):
        out.update(f(sd))
    out.update(stat_crossfit(sd, bootstrap=bootstrap))
    s726 = stat_spectrum_split(sd, "split_half_acquisitions_726.tsv")
    out.update({f"{k}_726": v for k, v in s726.items()})
    return out


def show(d, keys=None):
    """Print named values, one per line, for the notebook.

    Parameters
    ----------
    d : dict
        Named values, for example the output of a stat_* function.
    keys : iterable of str, optional
        Print only these keys.

    Returns
    -------
    None
    """
    for k, x in d.items():
        if keys is None or k in keys:
            print(f"{k:40s} {x}")


# ---------------------------------------------------------------- printed numbers
WORDS = {0: "zero", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight",
         9: "nine", 10: "ten", 11: "eleven", 12: "twelve", 13: "thirteen", 15: "fifteen"}


def as_word(n, printed):
    """Write count ``n`` as the English word used in the text, keeping the printed capitalisation.

    Parameters
    ----------
    n : int
        Recomputed count.
    printed : str
        The word as printed; "no", "none" and "No" stand for zero.

    Returns
    -------
    str
    """
    if printed.lower() in ("no", "none"):
        word = printed.lower() if n == 0 else str(n)
    else:
        word = WORDS.get(int(n), str(n))
    return word.capitalize() if printed[:1].isupper() else word


SD_ = "source_data/supplementary_methods_6_7/"
C_ = "text_statistics/supplementary_methods_6_7/stats.py:"
# (v28d_line, quote, printed_value, key, function, source_data files, status)
NUMBER_SPEC = [
    (633, "fixed 296-PREDICT and 282-CapScan matched acquisition rosters", "296", "conf_predict", "stat_confidence_rosters", "confidence_ordinary_all_pairs.tsv", "MATCH"),
    (633, "fixed 296-PREDICT and 282-CapScan matched acquisition rosters", "282", "conf_capscan", "stat_confidence_rosters", "confidence_ordinary_all_pairs.tsv", "MATCH"),
    (633, "all 30 broad-reservoir CAMPI fusion searches", "30", "conf_fusion_searches", "stat_confidence_rosters", "confidence_fusion_all.tsv", "MATCH"),
    (633, "All three affected PREDICT FastaLake acquisitions", "three", "conf_zero_yield_predict", "stat_confidence_rosters", "confidence_pep01_zero_yield_cases.tsv", "MATCH"),
    (665, "against the same 233,497,702-sequence reservoir", "233,497,702", "reservoir_public", "stat_database_sizes", "acquisition_sizes_726.tsv", "MATCH"),
    (667, "Selected FASTAs contained 499–557 proteins for all-position shuffles", "499", "scr_all_min", "stat_scramble", "scramble_end_to_end.tsv", "MATCH"),
    (667, "Selected FASTAs contained 499–557 proteins for all-position shuffles", "557", "scr_all_max", "stat_scramble", "scramble_end_to_end.tsv", "MATCH"),
    (667, "489–536 for terminal-preserving shuffles", "489", "scr_term_min", "stat_scramble", "scramble_end_to_end.tsv", "MATCH"),
    (667, "489–536 for terminal-preserving shuffles", "536", "scr_term_max", "stat_scramble", "scramble_end_to_end.tsv", "MATCH"),
    (667, "compared with 12,606 authentic proteins", "12,606", "scr_auth_proteins", "stat_scramble", "scramble_end_to_end.tsv", "MATCH"),
    (667, "Each shuffled arm yielded zero canonical target peptides", "zero", "scr_shuffled_peptides_max", "stat_scramble", "scramble_end_to_end.tsv", "MATCH"),
    (667, "The authentic arm yielded 29,663", "29,663", "scr_auth_peptides", "stat_scramble", "scramble_end_to_end.tsv", "MATCH"),
    (667, "were 8.11–8.65% at nine residues", "8.11", "scr_term9_min", "stat_scramble", "scramble_lookup_comparison.tsv", "MATCH"),
    (667, "were 8.11–8.65% at nine residues", "8.65%", "scr_term9_max", "stat_scramble", "scramble_lookup_comparison.tsv", "MATCH"),
    (667, "0.72–1.03% at ten residues", "0.72", "scr_term10_min", "stat_scramble", "scramble_lookup_comparison.tsv", "MATCH"),
    (667, "0.72–1.03% at ten residues", "1.03%", "scr_term10_max", "stat_scramble", "scramble_lookup_comparison.tsv", "MATCH"),
    (667, "no terminal-preserving query matched at lengths ≥13", "no", "scr_term_ge13_matched_max", "stat_scramble", "scramble_lookup_comparison.tsv", "MATCH"),
    (673, "for all 726 primary acquisitions", "726", "sizes_n_fastas", "stat_database_sizes", "acquisition_sizes_726.tsv", "MATCH"),
    (673, "the 20,431 canonical entries", "20,431", "sizes_human_entries", "stat_database_sizes", "human_reference.tsv", "MATCH"),
    (675, "we used the same 710 selected FASTAs", "710", "sizes_n_fastas", "stat_database_sizes", "acquisition_sizes_726.tsv", "v28f_change:726"),
    (675, "benchmark used 233,497,702 entries", "233,497,702", "reservoir_public", "stat_database_sizes", "acquisition_sizes_726.tsv", "MATCH"),
    (679, "the same 58 healthy-stool acquisitions for 19 settings", "58", "grid_acquisitions", "stat_selection_grid", "selection_grid_acquisitions.tsv", "MATCH"),
    (679, "the same 58 healthy-stool acquisitions for 19 settings", "19", "grid_settings", "stat_selection_grid", "selection_grid_acquisitions.tsv", "MATCH"),
    (679, "minimum evidence lengths of 6–12 at top 30%", "6", "grid_len_min", "stat_selection_grid", "selection_grid_acquisitions.tsv", "MATCH"),
    (679, "minimum evidence lengths of 6–12 at top 30%", "12", "grid_len_max", "stat_selection_grid", "selection_grid_acquisitions.tsv", "MATCH"),
    (679, "score fractions of 5–100% at minimum length nine", "5", "grid_score_min", "stat_selection_grid", "selection_grid_acquisitions.tsv", "MATCH"),
    (679, "score fractions of 5–100% at minimum length nine", "100%", "grid_score_max", "stat_selection_grid", "selection_grid_acquisitions.tsv", "MATCH"),
    (679, "All 1,102 acquisition–setting observations", "1,102", "grid_observations", "stat_selection_grid", "selection_grid_acquisitions.tsv", "MATCH"),
    (683, "set unions across the 15 acquisitions", "15", "campi_acquisitions", "stat_campi", "campi_spectral_inventory.tsv", "MATCH"),
    (683, "in exactly one of the eight supplied member proteomes", "eight", "campi_members", "stat_campi", "campi_member_summary.tsv", "MATCH"),
    (683, "The 116 contaminant entries were excluded", "116", "campi_contaminants", "stat_campi", "campi_reference_fasta_records.tsv", "MATCH"),
    (687, "overlap with the 76,127 distinct pooled reference peptides", "76,127", "campi_ref_pooled", "stat_campi", "campi_reference_peptides.tsv", "MATCH"),
    (687, "when drawing 38,254 MS2 spectra without replacement", "38,254", "campi_min_ms2", "stat_campi", "campi_spectral_inventory.tsv", "MATCH"),
    (687, "matching the smallest of the 15 acquisitions", "15", "campi_acquisitions", "stat_campi", "campi_spectral_inventory.tsv", "MATCH"),
    (695, "all 105 acquisition–arm configurations", "105", "cluster_configs", "stat_clustering", "clustering_per_acquisition_identifications.tsv", "MATCH"),
    (699, "We also searched the 15 acquisition-specific FastaLake", "15", "campi_acquisitions", "stat_campi", "campi_spectral_inventory.tsv", "MATCH"),
    (713, "in the 375,333,654-sequence reservoir", "375,333,654", "reservoir_campi", "stat_database_sizes", "acquisition_sizes_726.tsv", "MATCH"),
    (715, "the pooled diagnostic was 2.59–2.79% for 5–50 residues", "2.59", "fus_5_50_total_min", "stat_fusion", "fusion_campi_acquisition_components.tsv", "MATCH"),
    (715, "the pooled diagnostic was 2.59–2.79% for 5–50 residues", "2.79%", "fus_5_50_total_max", "stat_fusion", "fusion_campi_acquisition_components.tsv", "MATCH"),
    (715, "1.55–1.56% for the nested 9–50-residue subset", "1.55", "fus_9_50_total_min", "stat_fusion", "fusion_campi_acquisition_components.tsv", "MATCH"),
    (715, "1.55–1.56% for the nested 9–50-residue subset", "1.56%", "fus_9_50_total_max", "stat_fusion", "fusion_campi_acquisition_components.tsv", "MATCH"),
    (715, "the paired contribution was 3.09–3.38%", "3.09", "fus_5_8_paired_min", "stat_fusion", "fusion_campi_acquisition_components.tsv", "MATCH"),
    (715, "the paired contribution was 3.09–3.38%", "3.38%", "fus_5_8_paired_max", "stat_fusion", "fusion_campi_acquisition_components.tsv", "MATCH"),
    (715, "the unresolved contribution 6.35–8.30%", "6.35", "fus_5_8_unresolved_min", "stat_fusion", "fusion_campi_acquisition_components.tsv", "MATCH"),
    (715, "the unresolved contribution 6.35–8.30%", "8.30%", "fus_5_8_unresolved_max", "stat_fusion", "fusion_campi_acquisition_components.tsv", "MATCH"),
    (715, "giving 9.72–11.39% in total", "9.72", "fus_5_8_total_min", "stat_fusion", "fusion_campi_acquisition_components.tsv", "MATCH"),
    (715, "giving 9.72–11.39% in total", "11.39%", "fus_5_8_total_max", "stat_fusion", "fusion_campi_acquisition_components.tsv", "MATCH"),
    (715, "Unresolved pairing accounted for 65–73% of that total", "65", "fus_5_8_share_min", "stat_fusion", "fusion_campi_acquisition_components.tsv", "MATCH"),
    (715, "Unresolved pairing accounted for 65–73% of that total", "73%", "fus_5_8_share_max", "stat_fusion", "fusion_campi_acquisition_components.tsv", "MATCH"),
    (715, "the paired contribution was 1.10–1.20%", "1.10", "fus_9_50_paired_min", "stat_fusion", "fusion_campi_acquisition_components.tsv", "MATCH"),
    (715, "the paired contribution was 1.10–1.20%", "1.20%", "fus_9_50_paired_max", "stat_fusion", "fusion_campi_acquisition_components.tsv", "MATCH"),
    (715, "the unresolved contribution 0.36–0.45%", "0.36", "fus_9_50_unresolved_min", "stat_fusion", "fusion_campi_acquisition_components.tsv", "MATCH"),
    (715, "the unresolved contribution 0.36–0.45%", "0.45%", "fus_9_50_unresolved_max", "stat_fusion", "fusion_campi_acquisition_components.tsv", "MATCH"),
    (721, "across all 375,333,654 original reservoir proteins", "375,333,654", "lab_reservoir", "stat_label_validity", "label_validity_result.tsv", "MATCH"),
    (721, "All 4,000 positive sequences sampled", "4,000", "lab_pos_5_50_total", "stat_label_validity", "label_validity_presence.tsv", "MATCH"),
    (721, "All 4,000 positive sequences sampled from original reservoir proteins were recovered", "4,000", "lab_pos_5_50_present", "stat_label_validity", "label_validity_presence.tsv", "MATCH"),
    (721, "1,517/1,544 accepted controls (98.25%)", "1,517", "lab_ctl_5_8_present", "stat_label_validity", "label_validity_presence.tsv", "MATCH"),
    (721, "1,517/1,544 accepted controls (98.25%)", "1,544", "lab_ctl_5_8_total", "stat_label_validity", "label_validity_presence.tsv", "MATCH"),
    (721, "1,517/1,544 accepted controls (98.25%)", "98.25%", "lab_ctl_5_8_pct", "stat_label_validity", "label_validity_presence.tsv", "MATCH"),
    (721, "1,511/1,543 shuffled sequences (97.93%)", "1,511", "lab_shuf_5_8_present", "stat_label_validity", "label_validity_presence.tsv", "MATCH"),
    (721, "1,511/1,543 shuffled sequences (97.93%)", "1,543", "lab_shuf_5_8_total", "stat_label_validity", "label_validity_presence.tsv", "MATCH"),
    (721, "1,511/1,543 shuffled sequences (97.93%)", "97.93%", "lab_shuf_5_8_pct", "stat_label_validity", "label_validity_presence.tsv", "MATCH"),
    (721, "the corresponding counts were 480/3,177 (15.11%)", "480", "lab_ctl_9_50_present", "stat_label_validity", "label_validity_presence.tsv", "MATCH"),
    (721, "the corresponding counts were 480/3,177 (15.11%)", "3,177", "lab_ctl_9_50_total", "stat_label_validity", "label_validity_presence.tsv", "MATCH"),
    (721, "the corresponding counts were 480/3,177 (15.11%)", "15.11%", "lab_ctl_9_50_pct", "stat_label_validity", "label_validity_presence.tsv", "MATCH"),
    (721, "and 350/3,177 (11.02%)", "350", "lab_shuf_9_50_present", "stat_label_validity", "label_validity_presence.tsv", "MATCH"),
    (721, "and 350/3,177 (11.02%)", "3,177", "lab_shuf_9_50_total", "stat_label_validity", "label_validity_presence.tsv", "MATCH"),
    (721, "and 350/3,177 (11.02%)", "11.02%", "lab_shuf_9_50_pct", "stat_label_validity", "label_validity_presence.tsv", "MATCH"),
    (727, "Each of the 58 healthy-stool acquisitions was searched", "58", "donor_targets", "stat_donor", "donor_recount.tsv", "MATCH"),
    (727, "No entrapment peptide received exactly the same score as its paired target", "No", "donor_ties_q01", "stat_donor", "donor_recount.tsv", "MATCH"),
    (729, "Three seeds were averaged within target acquisition", "Three", "donor_seeds", "stat_donor", "donor_recount.tsv", "MATCH"),
    (729, "The five self-selected controls remain a separate stratum", "five", "self_targets", "stat_donor", "donor_recount.tsv", "MATCH"),
    (733, "spectra for 2,048 generating peptides per acquisition", "2,048", "sim_spectra", "stat_simulation", "simulation_design.tsv", "MATCH"),
    (733, "to create twelve baseline and twelve stress mixtures", "twelve", "sim_units_baseline", "stat_simulation", "simulation_all_outcomes.tsv", "MATCH"),
    (733, "and twelve stress mixtures", "twelve", "sim_units_stress", "stat_simulation", "simulation_all_outcomes.tsv", "MATCH"),
    (733, "Generating peptides were 7–30 residues", "7", "sim_len_min", "stat_simulation", "simulation_design.tsv", "MATCH"),
    (733, "Generating peptides were 7–30 residues", "30", "sim_len_max", "stat_simulation", "simulation_design.tsv", "MATCH"),
    (733, "an approximately 30,000-protein CAMPI reference", "30,000", "sim_reference_proteins_approx", "stat_simulation", "simulation_design.tsv", "MATCH"),
    (733, "All 72 ordinary outcomes were retained", "72", "sim_ordinary_outcomes", "stat_simulation", "simulation_all_outcomes.tsv", "MATCH"),
    (733, "separately for the twelve units in each regime", "twelve", "sim_stress_self_n", "stat_simulation", "simulation_all_outcomes.tsv", "MATCH"),
    (733, "(20,000 bootstrap resamples)", "20,000", "sim_bootstrap_draws", "stat_simulation", "simulation_design.tsv", "MATCH"),
    (733, "The 72 ordinary-search outcomes are a defined subset of 288 total simulation outcomes", "72", "sim_ordinary_outcomes", "stat_simulation", "simulation_all_outcomes.tsv", "MATCH"),
    (733, "The 72 ordinary-search outcomes are a defined subset of 288 total simulation outcomes", "288", "sim_all_outcomes", "stat_simulation", "simulation_all_outcomes.tsv", "MATCH"),
    (733, "were 1.25% for the fixed-reference search", "1.25%", "sim_baseline_fixed_mean", "stat_simulation", "simulation_all_outcomes.tsv", "MATCH"),
    (733, "1.06% with selection from the replicate", "1.06%", "sim_baseline_independent_mean", "stat_simulation", "simulation_all_outcomes.tsv", "MATCH"),
    (733, "1.01% with selection from the measurement being searched", "1.01%", "sim_baseline_self_mean", "stat_simulation", "simulation_all_outcomes.tsv", "MATCH"),
    (733, "and 1.50%, 1.20% and 1.91%, respectively", "1.50%", "sim_stress_fixed_mean", "stat_simulation", "simulation_all_outcomes.tsv", "MATCH"),
    (733, "and 1.50%, 1.20% and 1.91%, respectively", "1.20%", "sim_stress_independent_mean", "stat_simulation", "simulation_all_outcomes.tsv", "MATCH"),
    (733, "and 1.50%, 1.20% and 1.91%, respectively", "1.91%", "sim_stress_self_mean", "stat_simulation", "simulation_all_outcomes.tsv", "MATCH"),
    (733, "50% of generating peptides were deliberately absent", "50%", "sim_absent_stress_pct", "stat_simulation", "simulation_design.tsv", "MATCH"),
    (733, "compared with 25% at baseline", "25%", "sim_absent_baseline_pct", "stat_simulation", "simulation_design.tsv", "MATCH"),
    (733, "increased the proportion by 0.71 percentage points", "0.71", "sim_stress_self_minus_indep", "stat_simulation", "simulation_all_outcomes.tsv", "MATCH"),
    (733, "(95% bootstrap interval, 0.35–1.12)", "0.35", "sim_stress_self_minus_indep_lo", "stat_simulation", "simulation_all_outcomes.tsv", "MATCH"),
    (733, "(95% bootstrap interval, 0.35–1.12)", "1.12", "sim_stress_self_minus_indep_hi", "stat_simulation", "simulation_all_outcomes.tsv", "MATCH"),
    (733, "by 0.41 percentage points relative to the fixed-reference search", "0.41", "sim_stress_self_minus_fixed", "stat_simulation", "simulation_all_outcomes.tsv", "MATCH"),
    (739, "For each of the 58 healthy-stool acquisitions", "58", "transfer_targets", "stat_transfer_pool", "transfer_pool_comparisons.tsv", "MATCH"),
    (739, "selected from the ten inter-individual acquisitions", "ten", "inter_targets", "stat_transfer_pool", "transfer_pool_comparisons.tsv", "MATCH"),
    (739, "For the ten inter-individual targets", "ten", "inter_targets", "stat_transfer_pool", "transfer_pool_comparisons.tsv", "MATCH"),
    (739, "differed for nine targets and coincided for one", "nine", "inter_donors_differ", "stat_transfer_pool", "transfer_pool_comparisons.tsv", "MATCH"),
    (739, "differed for nine targets and coincided for one", "one", "inter_donors_same", "stat_transfer_pool", "transfer_pool_comparisons.tsv", "MATCH"),
    (743, "accounted for 2 of 376,412 accepted", "2", "pool_own_shuffle_only", "stat_transfer_pool", "hstool_pooling_per_acquisition.tsv", "MATCH"),
    (743, "accounted for 2 of 376,412 accepted", "376,412", "pool_own_observations", "stat_transfer_pool", "hstool_pooling_per_acquisition.tsv", "MATCH"),
    (743, "12 of 411,354 with the other-eleven unions", "12", "pool_pool_other11_shuffle_only", "stat_transfer_pool", "hstool_pooling_per_acquisition.tsv", "MATCH"),
    (743, "12 of 411,354 with the other-eleven unions", "411,354", "pool_pool_other11_observations", "stat_transfer_pool", "hstool_pooling_per_acquisition.tsv", "MATCH"),
    (745, "For the 58 healthy-stool acquisitions, we built KO profiles", "58", "outside_acquisitions", "stat_common_space", "representative_presence.tsv", "MATCH"),
    (751, "we deduplicated the 58 acquisition-selected FASTAs", "58", "union_sources", "stat_common_space", "common_fasta_union_summary.tsv", "MATCH"),
    (751, "to give a union of 153,267 protein sequences", "153,267", "union_sequences", "stat_common_space", "common_fasta_union_summary.tsv", "MATCH"),
    (751, "contained 555,987 sequences for PREDICT", "555,987", "union_PREDICT306", "stat_common_space", "union_contributions.tsv", "MATCH"),
    (751, "and 452,255 for CapScan", "452,255", "union_CAPSCAN347", "stat_common_space", "union_contributions.tsv", "v28f_change:460,977"),
    (751, "43.2% of quantified peptide features", "43.2%", "outside_pct", "stat_common_space", "representative_presence.tsv", "MATCH"),
    (755, "The fixed roster comprised 58 healthy-stool, 15 CAMPI, 306 PREDICT and 331 CapScan acquisitions", "58", "split_n_HSTOOL", "stat_spectrum_split", "split_half_acquisitions_710.tsv", "MATCH"),
    (755, "The fixed roster comprised 58 healthy-stool, 15 CAMPI, 306 PREDICT and 331 CapScan acquisitions", "15", "split_n_CAMPI_MOCK", "stat_spectrum_split", "split_half_acquisitions_710.tsv", "MATCH"),
    (755, "The fixed roster comprised 58 healthy-stool, 15 CAMPI, 306 PREDICT and 331 CapScan acquisitions", "306", "split_n_PREDICT", "stat_spectrum_split", "split_half_acquisitions_710.tsv", "MATCH"),
    (755, "The fixed roster comprised 58 healthy-stool, 15 CAMPI, 306 PREDICT and 331 CapScan acquisitions", "331", "split_n_CAPSCAN", "stat_spectrum_split", "split_half_acquisitions_710.tsv", "MATCH"),
    (757, "93 PREDICT patients (302 mapped acquisitions, four missing IDs)", "93", "split_predict_participants", "stat_spectrum_split", "split_half_acquisitions_710.tsv", "MATCH"),
    (757, "93 PREDICT patients (302 mapped acquisitions, four missing IDs)", "302", "split_predict_mapped", "stat_spectrum_split", "split_half_acquisitions_710.tsv", "MATCH"),
    (757, "93 PREDICT patients (302 mapped acquisitions, four missing IDs)", "four", "split_predict_missing", "stat_spectrum_split", "split_half_acquisitions_710.tsv", "MATCH"),
    (757, "15 CapScan participants (all 331 acquisitions)", "15", "split_capscan_participants", "stat_spectrum_split", "split_half_acquisitions_710.tsv", "MATCH"),
    (757, "15 CapScan participants (all 331 acquisitions)", "331", "split_capscan_mapped", "stat_spectrum_split", "split_half_acquisitions_710.tsv", "MATCH"),
    (759, "five truncated PREDICT inputs were replaced", "five", "split_reconverted", "stat_spectrum_split", "split_half_acquisitions_710.tsv", "MATCH"),
    (759, "A sensitivity analysis excluding these five acquisitions", "five", "split_reconverted", "stat_spectrum_split", "split_half_acquisitions_710.tsv", "MATCH"),
    (761, "Thirteen Sage arm searches across eight acquisitions", "Thirteen", "split_fallback_searches", "stat_spectrum_split", "split_half_model_fallbacks.tsv", "MATCH"),
    (761, "Thirteen Sage arm searches across eight acquisitions", "eight", "split_fallback_acquisitions", "stat_spectrum_split", "split_half_model_fallbacks.tsv", "MATCH"),
    (761, "One low-yield CapScan disjoint arm had zero accepted peptides", "One", "split_capscan_zero_disjoint", "stat_spectrum_split", "split_half_acquisitions_710.tsv", "MATCH"),
    (761, "One low-yield CapScan disjoint arm had zero accepted peptides", "zero", "split_zero_accepted", "stat_spectrum_split", "split_half_zero_yield.tsv", "MATCH"),
    (761, "(minimum target peptide q 0.0473)", "0.0473", "split_zero_min_q", "stat_spectrum_split", "split_half_zero_yield.tsv", "MATCH"),
    (765, "Applied once across all 58 acquisitions together", "58", "joint_acquisitions", "stat_joint", "joint_counts_per_acquisition.tsv", "MATCH"),
    (765, "gave a median gain of 32.0%", "32.0%", "joint_glob_gain", "stat_joint", "joint_counts_per_acquisition.tsv", "MATCH"),
    (765, "compared with 32.9% from Sage", "32.9%", "joint_sage_gain", "stat_joint", "joint_counts_per_acquisition.tsv", "MATCH"),
    (765, "it accepted at most 17 peptides in any acquisition", "17", "joint_perm_max", "stat_joint", "joint_fileq_coverage.tsv", "MATCH"),
    (765, "With one added to the decoy count, it accepted none", "none", "joint_plus1_max", "stat_joint", "joint_counts_per_acquisition.tsv", "MATCH"),
    (769, "accepted same-acquisition peptide evidence in all twelve technical replicates", "twelve", "jq_fileq_tech_acq", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "across all twelve technical or fifteen preparation replicates", "twelve", "jq_fileq_tech_acq", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "across all twelve technical or fifteen preparation replicates", "fifteen", "jq_fileq_prep_acq", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "between technical replicates from 20.1% to 15.2%", "20.1%", "jq_fileq_tech_orig_cv", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "between technical replicates from 20.1% to 15.2%", "15.2%", "jq_fileq_tech_joint_cv", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "increased Pearson correlation from 0.83 to 0.87", "0.83", "jq_fileq_tech_orig_r", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "increased Pearson correlation from 0.83 to 0.87", "0.87", "jq_fileq_tech_joint_r", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "CV increased from 32.9% to 35.1%", "32.9%", "jq_fileq_prep_orig_cv", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "CV increased from 32.9% to 35.1%", "35.1%", "jq_fileq_prep_joint_cv", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "correlation decreased from 0.76 to 0.72", "0.76", "jq_fileq_prep_orig_r", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "correlation decreased from 0.76 to 0.72", "0.72", "jq_fileq_prep_joint_r", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "the corresponding changes were from 20.3% to 14.5% and 0.83 to 0.87 for technical replicates", "20.3%", "jq_study_tech_orig_cv", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "the corresponding changes were from 20.3% to 14.5% and 0.83 to 0.87 for technical replicates", "14.5%", "jq_study_tech_joint_cv", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "the corresponding changes were from 20.3% to 14.5% and 0.83 to 0.87 for technical replicates", "0.83", "jq_study_tech_orig_r", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "the corresponding changes were from 20.3% to 14.5% and 0.83 to 0.87 for technical replicates", "0.87", "jq_study_tech_joint_r", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "from 33.1% to 33.6% and 0.75 to 0.74 for independent preparations", "33.1%", "jq_study_prep_orig_cv", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "from 33.1% to 33.6% and 0.75 to 0.74 for independent preparations", "33.6%", "jq_study_prep_joint_cv", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "from 33.1% to 33.6% and 0.75 to 0.74 for independent preparations", "0.75", "jq_study_prep_orig_r", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "from 33.1% to 33.6% and 0.75 to 0.74 for independent preparations", "0.74", "jq_study_prep_joint_r", "stat_joint", "microbial_group_pairs.tsv", "MATCH"),
    (769, "remained 74.9% among technical replicates under both strategies", "74.9%", "lead_tech_orig", "stat_joint", "technical_agreement_summary.tsv", "MATCH"),
    (769, "remained 74.9% among technical replicates under both strategies", "74.9%", "lead_tech_joint", "stat_joint", "technical_agreement_summary.tsv", "MATCH"),
    (769, "decreased from 70.8% to 66.5% among independent preparations", "70.8%", "lead_prep_orig", "stat_joint", "technical_agreement_summary.tsv", "MATCH"),
    (769, "decreased from 70.8% to 66.5% among independent preparations", "66.5%", "lead_prep_joint", "stat_joint", "technical_agreement_summary.tsv", "MATCH"),
    (773, "Folds held 9.9–10.1% of MS2 spectra in the six pilot acquisitions", "9.9", "cf_fold_min", "stat_crossfit", "crossfit_pilot_folds.tsv", "MATCH"),
    (773, "Folds held 9.9–10.1% of MS2 spectra in the six pilot acquisitions", "10.1%", "cf_fold_max", "stat_crossfit", "crossfit_pilot_folds.tsv", "v28f_change:10.2%"),
    (773, "Folds held 9.9–10.1% of MS2 spectra in the six pilot acquisitions", "six", "cf_pilots", "stat_crossfit", "crossfit_pilot_folds.tsv", "MATCH"),
    (773, "out-of-fold prediction tables of all 306 PREDICT and 331 CapScan acquisitions", "306", "cf_oof_PREDICT", "stat_crossfit", "crossfit_fold_caches.tsv", "MATCH"),
    (773, "out-of-fold prediction tables of all 306 PREDICT and 331 CapScan acquisitions", "331", "cf_oof_CAPSCAN", "stat_crossfit", "crossfit_fold_caches.tsv", "MATCH"),
    (775, "On 36 archived unsplit searches", "36", "cf_q_searches", "stat_crossfit", "crossfit_q_validation.tsv", "MATCH"),
    (775, "reproduced Sage's accepted sets in 36 of 36", "36", "cf_q_equal", "stat_crossfit", "crossfit_q_validation.tsv", "MATCH"),
    (775, "the full rescoring in 33 of 36", "33", "cf_rescoring_equal", "stat_crossfit", "crossfit_q_validation.tsv", "v28f_change:31"),
    (775, "the full rescoring in 33 of 36", "36", "cf_q_searches", "stat_crossfit", "crossfit_q_validation.tsv", "MATCH"),
    (775, "selection reproduced the archived databases in six of six pilot acquisitions", "six", "cf_k1_selection_equal", "stat_crossfit", "crossfit_k1_control.tsv", "MATCH"),
    (775, "Sage reproduced the original accepted peptide sets in 24 of 24 comparisons", "24", "cf_k1_sage_equal", "stat_crossfit", "crossfit_k1_control.tsv", "MATCH"),
    (775, "Through the pooled rescoring, 22 of 24 were identical", "22", "cf_k1_pooled_equal", "stat_crossfit", "crossfit_k1_control.tsv", "MATCH"),
    (775, "Through the pooled rescoring, 22 of 24 were identical", "24", "cf_k1_comparisons", "stat_crossfit", "crossfit_k1_control.tsv", "MATCH"),
    (775, "with 3,047 against 3,100 peptides at 5–50 residues", "3,047", "cf_k1_252_pooled_5", "stat_crossfit", "crossfit_k1_control.tsv", "MATCH"),
    (775, "with 3,047 against 3,100 peptides at 5–50 residues", "3,100", "cf_k1_252_frozen_5", "stat_crossfit", "crossfit_k1_control.tsv", "MATCH"),
    (775, "2,361 against 2,429 at 9–50 residues", "2,361", "cf_k1_252_pooled_9", "stat_crossfit", "crossfit_k1_control.tsv", "MATCH"),
    (775, "2,361 against 2,429 at 9–50 residues", "2,429", "cf_k1_252_frozen_9", "stat_crossfit", "crossfit_k1_control.tsv", "MATCH"),
    (775, "Across all 578 acquisitions", "578", "cf_acquisitions", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (775, "with median differences of at most 0.02 points in every cohort stratum", "0.02", "cf_calib_max", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "changed stratum medians by at most 0.16 points", "0.16", "cf_pert_shift_max", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "but changed 25 and 28 of 578 acquisitions (4–5%)", "25", "cf_moved_pert1", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "but changed 25 and 28 of 578 acquisitions (4–5%)", "28", "cf_moved_pert2", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "but changed 25 and 28 of 578 acquisitions (4–5%)", "578", "cf_acquisitions", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "but changed 25 and 28 of 578 acquisitions (4–5%)", "4", "cf_moved_pert1_pct", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "but changed 25 and 28 of 578 acquisitions (4–5%)", "5%", "cf_moved_pert2_pct", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "The three differences in the rescoring check and the two in the one-fold control", "three", "cf_rescoring_differ", "stat_crossfit", "crossfit_q_validation.tsv", "v28f_change:five"),
    (777, "The three differences in the rescoring check and the two in the one-fold control", "two", "cf_k1_pooled_differ", "stat_crossfit", "crossfit_k1_control.tsv", "MATCH"),
    (777, "In three acquisitions (two PREDICT and one CapScan region 3)", "three", "cf_fallback_n", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "In three acquisitions (two PREDICT and one CapScan region 3)", "two", "cf_fallback_predict", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "In three acquisitions (two PREDICT and one CapScan region 3)", "one", "cf_fallback_capscan", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "which reproduced the accepted sets exactly in five checked acquisitions without fallback", "five", "cf_closed_form_exact", "stat_crossfit", "crossfit_independent_lda_check.tsv", "MATCH"),
    (777, "gave gains of 75.1%, 13.3% and 11.7% at 5–50 residues", "75.1%", "cf_fb_168_independent", "stat_crossfit", "crossfit_lda_fallback_acquisitions.tsv", "MATCH"),
    (777, "gave gains of 75.1%, 13.3% and 11.7% at 5–50 residues", "13.3%", "cf_fb_185_independent", "stat_crossfit", "crossfit_lda_fallback_acquisitions.tsv", "MATCH"),
    (777, "gave gains of 75.1%, 13.3% and 11.7% at 5–50 residues", "11.7%", "cf_fb_496_independent", "stat_crossfit", "crossfit_lda_fallback_acquisitions.tsv", "MATCH"),
    (777, "instead of 31.5%, −10.7% and −74.3%", "31.5%", "cf_fb_168_producer", "stat_crossfit", "crossfit_lda_fallback_acquisitions.tsv", "MATCH"),
    (777, "instead of 31.5%, −10.7% and −74.3%", "−10.7%", "cf_fb_185_producer", "stat_crossfit", "crossfit_lda_fallback_acquisitions.tsv", "MATCH"),
    (777, "instead of 31.5%, −10.7% and −74.3%", "−74.3%", "cf_fb_496_producer", "stat_crossfit", "crossfit_lda_fallback_acquisitions.tsv", "MATCH"),
    (777, "changed the PREDICT medians from 6.50% to 6.55% at 5–50 residues", "6.50%", "cf_predict_median_5", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "changed the PREDICT medians from 6.50% to 6.55% at 5–50 residues", "6.55%", "cf_predict_median_5_replaced", "stat_crossfit", "crossfit_per_acquisition.tsv; crossfit_lda_fallback_acquisitions.tsv", "MATCH"),
    (777, "from −3.54% to −3.49% at 9–50 residues", "−3.54%", "cf_predict_median_9", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "from −3.54% to −3.49% at 9–50 residues", "−3.49%", "cf_predict_median_9_replaced", "stat_crossfit", "crossfit_per_acquisition.tsv; crossfit_lda_fallback_acquisitions.tsv", "MATCH"),
    (777, "(−8.4 to −2.5% instead of −8.4 to −2.3%)", "−8.4", "cf_stool9_region_lo", "crossfit_bootstrap", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "(−8.4 to −2.5% instead of −8.4 to −2.3%)", "−2.5%", "cf_stool9_region_hi", "crossfit_bootstrap", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "(−8.4 to −2.5% instead of −8.4 to −2.3%)", "−8.4", "cf_stool9_lo", "crossfit_bootstrap", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "(−8.4 to −2.5% instead of −8.4 to −2.3%)", "−2.3%", "cf_stool9_hi", "crossfit_bootstrap", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "changed the cohort, stool and intestinal-fluid medians by at most 2 points", "2", "cf_union_shift_max", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "single sampling-region medians by up to 6.4 points", "6.4", "cf_union_region_shift_max", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "In two CapScan acquisitions every fold's metagenome search accepted no peptides", "two", "cf_mg_union_zero_capscan", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "in 106 of 296 PREDICT", "106", "cf_PREDICT_pos9", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "in 106 of 296 PREDICT", "296", "cf_PREDICT_n", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "16 of 56 CapScan stool", "16", "cf_CAPSCAN_stool_pos9", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "16 of 56 CapScan stool", "56", "cf_CAPSCAN_stool_n", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "206 of 226 intestinal-fluid acquisitions", "206", "cf_CAPSCAN_intestinal_pos9", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "206 of 226 intestinal-fluid acquisitions", "226", "cf_CAPSCAN_intestinal_n", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "identified a median 9.8% (PREDICT) and 14.2% (CapScan) fewer peptides", "9.8%", "cf_PREDICT_fewer_5", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "identified a median 9.8% (PREDICT) and 14.2% (CapScan) fewer peptides", "14.2%", "cf_CAPSCAN_all_fewer_5", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "gave lower gains of −1.0%, −10.5% and 22.9% at 5–50 residues", "−1.0%", "cf_PREDICT_odd_disjoint_5", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "gave lower gains of −1.0%, −10.5% and 22.9% at 5–50 residues", "−10.5%", "cf_CAPSCAN_stool_odd_disjoint_5", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "gave lower gains of −1.0%, −10.5% and 22.9% at 5–50 residues", "22.9%", "cf_CAPSCAN_intestinal_odd_disjoint_5", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "gave median gains of 12.4% in PREDICT, 11.6% in CapScan stool and 52.9% in intestinal fluid", "12.4%", "cf_PREDICT_odd_self_5", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "gave median gains of 12.4% in PREDICT, 11.6% in CapScan stool and 52.9% in intestinal fluid", "11.6%", "cf_CAPSCAN_stool_odd_self_5", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "gave median gains of 12.4% in PREDICT, 11.6% in CapScan stool and 52.9% in intestinal fluid", "52.9%", "cf_CAPSCAN_intestinal_odd_self_5", "stat_crossfit", "crossfit_per_acquisition.tsv", "MATCH"),
    (777, "reduced the median gains to 14.5% (PREDICT) and 45.0% (CapScan)", "14.5%", "ml7_PREDICT_gain", "stat_same_spectrum", "minlen7_benchmark_per_acquisition.tsv", "MATCH"),
    (777, "reduced the median gains to 14.5% (PREDICT) and 45.0% (CapScan)", "45.0%", "ml7_CAPSCAN_gain", "stat_same_spectrum", "minlen7_benchmark_per_acquisition.tsv", "MATCH"),
    (777, "identity recovery remained 85.4% and 84.0%", "85.4%", "ml7_PREDICT_recovery", "stat_same_spectrum", "minlen7_benchmark_per_acquisition.tsv", "MATCH"),
    (777, "identity recovery remained 85.4% and 84.0%", "84.0%", "ml7_CAPSCAN_recovery", "stat_same_spectrum", "minlen7_benchmark_per_acquisition.tsv", "MATCH"),
    (777, "reduced the gains to 9.4% and 42.4%", "9.4%", "ss_PREDICT_9_50_median", "stat_same_spectrum", "confidence_ordinary_all_pairs.tsv", "MATCH"),
    (777, "reduced the gains to 9.4% and 42.4%", "42.4%", "ss_CAPSCAN_9_50_median", "stat_same_spectrum", "confidence_ordinary_all_pairs.tsv", "MATCH"),
    (777, "in 272 of 296 and 281 of 282 acquisitions", "272", "ss_PREDICT_9_50_pos", "stat_same_spectrum", "confidence_ordinary_all_pairs.tsv", "MATCH"),
    (777, "in 272 of 296 and 281 of 282 acquisitions", "296", "ss_PREDICT_9_50_n", "stat_same_spectrum", "confidence_ordinary_all_pairs.tsv", "MATCH"),
    (777, "in 272 of 296 and 281 of 282 acquisitions", "281", "ss_CAPSCAN_9_50_pos", "stat_same_spectrum", "confidence_ordinary_all_pairs.tsv", "MATCH"),
    (777, "in 272 of 296 and 281 of 282 acquisitions", "282", "ss_CAPSCAN_9_50_n", "stat_same_spectrum", "confidence_ordinary_all_pairs.tsv", "MATCH"),
    (777, "showed gains of 67% and 96%", "67%", "ss_PREDICT_5_8_median", "stat_same_spectrum", "confidence_ordinary_all_pairs.tsv", "MATCH"),
    (777, "showed gains of 67% and 96%", "96%", "ss_CAPSCAN_5_8_median", "stat_same_spectrum", "confidence_ordinary_all_pairs.tsv", "MATCH"),
]


def recomputed(v, key, printed):
    """Return the recomputed value of ``key`` as a full-precision string (word if printed as a word).

    Parameters
    ----------
    v : dict
        Output of compute_all().
    key : str
        Key of the value.
    printed : str
        Value as printed in v28d.

    Returns
    -------
    str
    """
    x = v[key]
    if printed.strip().isalpha():
        return as_word(int(x), printed.strip())
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    return repr(float(x))


def numbers_rows(v):
    """Rows of NUMBERS.tsv, in the order of NUMBER_SPEC.

    Parameters
    ----------
    v : dict
        Output of compute_all().

    Returns
    -------
    list of dict
        One dict per NUMBERS.tsv row.
    """
    rows = []
    for line, quote, printed, key, func, src, status in NUMBER_SPEC:
        rows.append(dict(
            v28d_line=str(line), quote=quote, printed_value=printed,
            recomputed_value="" if key is None else recomputed(v, key, printed),
            code="" if func is None else C_ + func,
            source_data="; ".join(SD_ + s.strip() for s in src.split(";")) if src else "",
            status=status))
    return rows


# ---------------------------------------------------------------- statistical tests and intervals
STATS_COLUMNS = ["claim", "v28d_line", "test", "sided", "unit", "n", "distinct_or_repeated", "covariates", "correction",
                 "family", "family_size", "estimate", "uncertainty", "uncertainty_construction", "statistic", "df",
                 "exact_P", "effect_size", "effect_size_how", "draws", "seed", "software_version", "code"]


def _row(**kw):
    """One STATS.tsv row with the defaults for a descriptive estimate (no test).

    Parameters
    ----------
    **kw
        Column values; any of STATS_COLUMNS.

    Returns
    -------
    dict
        All STATS_COLUMNS as strings.
    """
    base = dict(sided="NA", covariates="none", correction="none", family="NA", family_size="NA", statistic="NA",
                df="NA", exact_P="NA (no test)", effect_size="NA", effect_size_how="NA", draws="NA", seed="NA")
    base.update(kw)
    return {c: str(base.get(c, "NA")) for c in STATS_COLUMNS}


def stats_rows(v, numpy_version):
    """Rows of STATS.tsv.

    No line of Supplementary Methods 6-7 reports a hypothesis test; every row is an interval
    or a descriptive estimator whose construction the text states only in part.

    Parameters
    ----------
    v : dict
        Output of compute_all().
    numpy_version : str
        NumPy version of the recomputation, written into software_version.

    Returns
    -------
    list of dict
        One dict per STATS.tsv row.
    """
    rows = []
    pdv = f"pandas {pd.__version__}"
    sw_sim = f"numpy {numpy_version} default_rng (PCG64); original run CPython 3.12.4, numpy 2.2.6, scipy 1.16.3"
    sw_cf = f"numpy {numpy_version} default_rng (PCG64); original run Python 3.12.4, numpy 2.2.6, pandas 2.3.3"
    for regime in ("baseline", "stress"):
        for arm, label in (("fixed", "fixed reference"), ("independent", "selection from independently noised replicate"),
                           ("self", "selection from the searched measurement")):
            k = f"sim_{regime}_{arm}"
            rows.append(_row(
                claim=f"Mean simulated false-discovery proportion, {regime}, {label}", v28d_line="733",
                test="none (mean over mixtures with percentile bootstrap interval)", unit="acquisition",
                n=f"{v[k + '_n']} synthetic mixtures", distinct_or_repeated="distinct mixtures; the three arms are repeated searches of each mixture",
                estimate=f"mean FDP {v[k + '_mean']:.4f}%", uncertainty=f"95% CI {v[k + '_lo']:.4f} to {v[k + '_hi']:.4f}%",
                uncertainty_construction=f"percentile bootstrap (2.5% and 97.5% quantiles), {SIM_BOOTSTRAP_DRAWS:,} draws, resampling whole mixtures; one random stream shared by the eight summaries in fixed order",
                draws=SIM_BOOTSTRAP_DRAWS, seed=SIM_BOOTSTRAP_SEED, software_version=sw_sim, code=C_ + "stat_simulation"))
        k = f"sim_{regime}_self_minus_indep"
        rows.append(_row(
            claim=f"Paired FDP difference, selection from searched measurement minus replicate, {regime}", v28d_line="733",
            test="none (mean paired difference with percentile bootstrap interval)", unit="acquisition",
            n="12 synthetic mixtures", distinct_or_repeated="repeated: two search arms within each mixture (paired difference)",
            estimate=f"mean difference {v[k]:.4f} percentage points", uncertainty=f"95% CI {v[k + '_lo']:.4f} to {v[k + '_hi']:.4f} percentage points",
            uncertainty_construction=f"percentile bootstrap, {SIM_BOOTSTRAP_DRAWS:,} draws, resampling whole mixtures (pairs kept together)",
            effect_size=f"{v[k]:.4f}", effect_size_how="mean of per-mixture (self FDP - independent FDP), percentage points",
            draws=SIM_BOOTSTRAP_DRAWS, seed=SIM_BOOTSTRAP_SEED, software_version=sw_sim, code=C_ + "stat_simulation"))
    rows.append(_row(
        claim="Paired FDP difference, selection from searched measurement minus fixed reference, stress", v28d_line="733",
        test="none (mean paired difference)", unit="acquisition", n="12 synthetic mixtures",
        distinct_or_repeated="repeated: two search arms within each mixture (paired difference)",
        estimate=f"mean difference {v['sim_stress_self_minus_fixed']:.4f} percentage points",
        uncertainty="none computed by the original code", uncertainty_construction="NA",
        effect_size=f"{v['sim_stress_self_minus_fixed']:.4f}", effect_size_how="mean of per-mixture (self FDP - fixed FDP), percentage points",
        software_version=sw_sim, code=C_ + "stat_simulation"))
    b = v["cf_boot"].set_index(["stratum", "column"])
    names = {"PREDICT": "PREDICT", "CAPSCAN_stool": "CapScan stool", "CAPSCAN_intestinal": "CapScan intestinal fluid"}
    labels = {"gain_cf_pooled_5": ("10-fold cross-fitted gain over metagenome, 5-50 residues", "777"),
              "gain_cf_pooled_9": ("10-fold cross-fitted gain over metagenome, 9-50 residues", "777"),
              "gain_odd_self_5": ("Gain over metagenome, selection from odd-scan predictions, odd-scan spectra, 5-50 residues", "777"),
              "gain_odd_disjoint_5": ("Gain over metagenome, selection from even-scan predictions, odd-scan spectra, 5-50 residues", "777")}
    for col, (lab, line) in labels.items():
        for st, nm in names.items():
            r = b.loc[(st, col)]
            rows.append(_row(
                claim=f"{lab}; {nm}", v28d_line=line, test="none (median with participant cluster bootstrap interval)",
                unit="acquisition", n=f"{int(r.n_acquisitions)} acquisitions ({int(r.n_participants)} participants)",
                distinct_or_repeated="repeated: acquisitions within participants",
                estimate=f"median gain {r['median']:.4f}%", uncertainty=f"95% CI {r.ci_low:.4f} to {r.ci_high:.4f}%",
                uncertainty_construction=(f"percentile bootstrap (2.5th and 97.5th percentiles), {CROSSFIT_BOOTSTRAP_DRAWS:,} draws, "
                                          "resampling participants with replacement (acquisitions without a participant are their own cluster); "
                                          "one random stream through the strata PREDICT, CapScan all, CapScan stool, CapScan intestinal fluid"),
                effect_size=f"{r['median']:.4f}", effect_size_how="median over acquisitions of 100 x (FL / MG - 1)",
                draws=CROSSFIT_BOOTSTRAP_DRAWS, seed=CROSSFIT_BOOTSTRAP_SEED, software_version=sw_cf, code=C_ + "crossfit_bootstrap"))
    rg = v["cf_boot_regions"].set_index(["stratum", "column"])
    for g in CROSSFIT_REGIONS:
        r = rg.loc[(f"CAPSCAN_{g}", "gain_cf_pooled_9")]
        rows.append(_row(
            claim=f"10-fold cross-fitted gain over metagenome, 9-50 residues; CapScan {g} (region stream)", v28d_line="777",
            test="none (median with participant cluster bootstrap interval)", unit="acquisition",
            n=f"{int(r.n_acquisitions)} acquisitions ({int(r.n_participants)} participants)", distinct_or_repeated="repeated: acquisitions within participants",
            estimate=f"median gain {r['median']:.4f}%", uncertainty=f"95% CI {r.ci_low:.4f} to {r.ci_high:.4f}%",
            uncertainty_construction=(f"percentile bootstrap, {CROSSFIT_BOOTSTRAP_DRAWS:,} draws, resampling participants; separate random stream "
                                      "restarted with the same seed through duodenum, jejunum, ileum, colon, stool"),
            effect_size=f"{r['median']:.4f}", effect_size_how="median over acquisitions of 100 x (FL / MG - 1)",
            draws=CROSSFIT_BOOTSTRAP_DRAWS, seed=CROSSFIT_BOOTSTRAP_SEED, software_version=sw_cf, code=C_ + "crossfit_bootstrap"))
    t = v["fus_table"]
    for (lo, hi), g in t.groupby(["lo", "hi"]):
        rows.append(_row(
            claim=f"CAMPI fusion diagnostic, peptide q <= 0.01, {lo}-{hi} residues", v28d_line="715",
            test="none (descriptive estimator (E + E0 + 2B + U) / N, acquisition counts summed before division)", unit="peptide",
            n="15 acquisitions x 2 control seeds",
            distinct_or_repeated="repeated: two control seeds of the same 15 acquisitions",
            estimate="diagnostic " + " and ".join(f"{x:.4f}%" for x in g.total) + f" (seeds {', '.join(str(s) for s in g.seed)})",
            uncertainty=f"range across seeds {g.total.min():.4f} to {g.total.max():.4f}% (not a confidence interval)",
            uncertainty_construction="two seeded positional-permutation controls; no resampling",
            effect_size="paired " + " and ".join(f"{x:.4f}%" for x in g.paired) + "; unresolved " + " and ".join(f"{x:.4f}%" for x in g.unresolved),
            effect_size_how="paired = (E + E0 + 2B) / N; unresolved = U / N",
            seed="; ".join(str(s) for s in g.seed), software_version=f"{pdv} (summation only); numpy {numpy_version}", code=C_ + "stat_fusion"))
    for L in (9, 10):
        rows.append(_row(
            claim=f"Terminal-preserving scrambled-query match fraction at {L} residues", v28d_line="667",
            test="none (observed fraction in a randomised control)", unit="sequence", n="distinct queries per seed (scramble_lookup_comparison.tsv)",
            distinct_or_repeated="repeated: five scramble seeds of the same acquisition",
            estimate=f"{v[f'scr_term{L}_min']:.4f} to {v[f'scr_term{L}_max']:.4f}%", uncertainty="range across five seeds (not a confidence interval)",
            uncertainty_construction="five seeded scrambles; no resampling", seed="; ".join(str(s) for s in SCRAMBLE_SEEDS),
            software_version=pdv, code=C_ + "stat_scramble"))
    rows.append(_row(
        claim="Terminal-preserving scrambled queries matching at lengths >= 13", v28d_line="667",
        test="none (observed count in a randomised control)", unit="sequence", n="distinct queries per seed and length",
        distinct_or_repeated="repeated: five scramble seeds of the same acquisition",
        estimate=f"maximum matched queries {v['scr_term_ge13_matched_max']}", uncertainty="none", uncertainty_construction="NA",
        seed="; ".join(str(s) for s in SCRAMBLE_SEEDS), software_version=pdv, code=C_ + "stat_scramble"))
    rows.append(_row(
        claim="Joint search, target/decoy label permutation within acquisition, per-acquisition D/T estimator", v28d_line="765",
        test="none (negative control: one label permutation)", unit="acquisition", n="58 acquisitions",
        distinct_or_repeated="distinct", estimate=f"maximum accepted peptides {v['joint_perm_max']}; with (D+1)/T {v['joint_plus1_max']}",
        uncertainty="none (single permutation)", uncertainty_construction="NA", draws="1 permutation per acquisition",
        seed=f"{JOINT_PERMUTATION_SEED} (numpy default_rng, one stream through acquisitions in sorted file-name order)",
        software_version="original run: iteration2 joint_search extract_fileq.py (numpy default_rng)", code=C_ + "stat_joint"))
    for key, conf in (("fileq", "per-acquisition confidence"), ("study", "study-level confidence")):
        for tk, tier in (("tech", "technical replicates"), ("prep", "independent preparations")):
            rows.append(_row(
                claim=f"Microbial protein-group pair CV and log2 Pearson r, acquisition-specific vs joint search, {conf}, {tier}",
                v28d_line="769", test="none (median across acquisition pairs of per-pair median CV; median of per-pair Pearson r)",
                unit="acquisition", n=f"{v[f'jq_{key}_{tk}_pairs']} pairs of {v[f'jq_{key}_{tk}_acq']} acquisitions",
                distinct_or_repeated="repeated: each acquisition occurs in several pairs; the two strategies are repeated analyses of the same pairs",
                estimate=(f"CV {v[f'jq_{key}_{tk}_orig_cv']:.4f}% to {v[f'jq_{key}_{tk}_joint_cv']:.4f}%; "
                          f"r {v[f'jq_{key}_{tk}_orig_r']:.4f} to {v[f'jq_{key}_{tk}_joint_r']:.4f}"),
                uncertainty="none reported", uncertainty_construction="NA",
                effect_size_how="CV = population SD / mean x 100 of two raw intensities, median over rows within a pair; Pearson r on log2 intensities of matched rows",
                software_version=pdv, code=C_ + "stat_joint"))
    rows.append(_row(
        claim="Leading annotated species agreement, technical replicates and independent preparations", v28d_line="769",
        test="none (function-averaged pairwise agreement)", unit="acquisition", n="12 technical and 15 preparation acquisitions",
        distinct_or_repeated="repeated: acquisitions in several pairs",
        estimate=(f"technical {v['lead_tech_orig']:.4f}% to {v['lead_tech_joint']:.4f}%; "
                  f"preparations {v['lead_prep_orig']:.4f}% to {v['lead_prep_joint']:.4f}%"),
        uncertainty="none reported", uncertainty_construction="NA",
        software_version="read from technical_agreement_summary.tsv (producer A21/A33 leaders step; not recomputed here)",
        code=C_ + "stat_joint"))
    for key, lab in (("donor", "donor-selected (58 targets)"), ("self", "self-selected (5 targets)")):
        for lo in (5, 9):
            rows.append(_row(
                claim=f"FDRBench paired entrapment estimate, {lab}, {lo}-50 residues, q <= 0.01", v28d_line="727-729",
                test="none (descriptive estimator (E + E0 + 2B) / (T + E))", unit="acquisition",
                n=f"{v[f'{key}_paired_n_{lo}']} target acquisitions x 3 seeds",
                distinct_or_repeated="repeated: three entrapment seeds averaged within target acquisition",
                estimate=f"median over targets {v[f'{key}_paired_median_{lo}']:.4f}%", uncertainty="none reported",
                uncertainty_construction="NA", seed="; ".join(str(s) for s in v["donor_seed_values"]),
                software_version=f"FDRBench v{FDRBENCH_VERSION} (control generation); estimator recomputed with {pdv}",
                code=C_ + "stat_donor"))
    rows.append(_row(
        claim="CAMPI acquisition-depth standardisation to 38,254 MS2 spectra", v28d_line="687",
        test="none (analytic hypergeometric expectation of distinct accepted peptides)", unit="acquisition", n="15 acquisitions",
        distinct_or_repeated="distinct", estimate="no value printed in lines 661-777", uncertainty="NA",
        uncertainty_construction="analytic expectation; no random draws", draws="NA (analytic)", seed="NA (analytic)",
        software_version="NA", code=C_ + "stat_campi"))
    return rows


# ---------------------------------------------------------------- comparison with the Source Data package summaries
def compare_with_package(sd, v):
    """Compare recomputed values with the summary sheets copied from the Source Data package.

    Parameters
    ----------
    sd : str or Path
        Directory ``source_data/supplementary_methods_6_7``.
    v : dict
        Output of compute_all().

    Returns
    -------
    DataFrame
        one row per comparison: package sheet, quantity, recomputed, package value, absolute difference.
    """
    rows = []

    def add(sheet, what, mine, theirs):
        """Append one comparison row (package sheet, quantity, recomputed, package value)."""
        rows.append(dict(package_sheet=sheet, quantity=what, recomputed=float(mine), package=float(theirs),
                         abs_diff=abs(float(mine) - float(theirs))))

    u = load(sd, "simulation_uncertainty.tsv").set_index(["regime", "endpoint"])
    for regime in ("baseline", "stress"):
        for arm in ("fixed", "independent", "self"):
            k = f"sim_{regime}_{arm}"
            for mine, col in (("_mean", "mean"), ("_lo", "bootstrap_low"), ("_hi", "bootstrap_high")):
                add("Supplementary_Figure_05/simulation_uncertainty.tsv", f"{regime} {arm} {col} (percent)",
                    v[k + mine], 100 * u.loc[(regime, arm), col])
        k = f"sim_{regime}_self_minus_indep"
        for mine, col in (("", "mean"), ("_lo", "bootstrap_low"), ("_hi", "bootstrap_high")):
            add("Supplementary_Figure_05/simulation_uncertainty.tsv", f"{regime} self_minus_independent {col} (pp)",
                v[k + mine], 100 * u.loc[(regime, "self_minus_independent"), col])
    s = load(sd, "crossfit_summary.tsv")
    s["col"] = s.comparison + "_" + s.min_length.astype(str)
    s = s.set_index(["stratum", "col"])
    for r in v["cf_boot"].itertuples():
        for mine, col in (("median", "median_exact"), ("ci_low", "ci95_low_exact"), ("ci_high", "ci95_high_exact")):
            add("Supplementary_Figure_06/metagenome_crossfit/summary.tsv", f"{r.stratum} {r.column} {col}",
                getattr(r, mine), s.loc[(r.stratum, r.column), col])
    s = load(sd, "crossfit_summary_by_region.tsv")
    s["col"] = s.comparison + "_" + s.min_length.astype(str)
    s = s.set_index(["stratum", "col"])
    for r in v["cf_boot_regions"].itertuples():
        for mine, col in (("median", "median_exact"), ("ci_low", "ci95_low_exact"), ("ci_high", "ci95_high_exact")):
            add("Supplementary_Figure_06/metagenome_crossfit/summary_capscan_by_region.tsv", f"{r.stratum} {r.column} {col}",
                getattr(r, mine), s.loc[(r.stratum, r.column), col])
    d = load(sd, "fusion_decomposition_current375.tsv").set_index(["seed", "min", "max"])
    for r in v["fus_table"].itertuples():
        for mine, col in (("total", "total_pct"), ("paired", "paired_pct"), ("unresolved", "unresolved_pct"),
                          ("share", "unresolved_share_pct")):
            add("Supplementary_Figure_05/fusion_decomposition_current375.tsv", f"seed {r.seed} {r.lo}-{r.hi} {col}",
                getattr(r, mine), d.loc[(r.seed, r.lo, r.hi), col])
    ls = load(sd, "scramble_lookup_summary.tsv").set_index("length")
    for L in (9, 10):
        add("Review_v16/Terminal_scramble/scramble_lookup_summary.tsv", f"length {L} terminal_min (percent)",
            v[f"scr_term{L}_min"], 100 * ls.loc[L, "terminal_min"])
        add("Review_v16/Terminal_scramble/scramble_lookup_summary.tsv", f"length {L} terminal_max (percent)",
            v[f"scr_term{L}_max"], 100 * ls.loc[L, "terminal_max"])
    ms = load(sd, "microbial_group_pair_summary.tsv")
    ms = ms[ms.definition == "matched"].set_index(["comparison", "method", "tier"])
    for comp, key in (("original_vs_joint_fileq", "fileq"), ("original_vs_joint_corrected", "study")):
        for tier, tk in (("TechnicalReplicates", "tech"), ("WorkflowReplicates", "prep")):
            for method in ms.loc[comp].index.get_level_values(0).unique():
                m = "orig" if method.startswith("original") else "joint"
                add("Supplementary_Figure_07/per_acquisition_confidence/microbial_group_pair_summary.tsv",
                    f"{comp} {method} {tier} median_pair_cv_pct", v[f"jq_{key}_{tk}_{m}_cv"],
                    ms.loc[(comp, method, tier), "median_pair_cv_pct"])
                add("Supplementary_Figure_07/per_acquisition_confidence/microbial_group_pair_summary.tsv",
                    f"{comp} {method} {tier} log2_pearson", v[f"jq_{key}_{tk}_{m}_r"],
                    ms.loc[(comp, method, tier), "log2_pearson"])
    os_ = load(sd, "confidence_ordinary_summary.tsv")
    os_ = os_[(os_.gate == "q_only") & (os_.q == Q_CUTOFF)].set_index(["cohort", "min_length", "max_length"])
    for c in ("PREDICT", "CAPSCAN"):
        for lo, hi in (NESTED_LENGTH, SHORT_LENGTH, FULL_LENGTH):
            add("Review_v16/Confidence/ordinary_summary.tsv", f"{c} {lo}-{hi} median_gain_pct",
                v[f"ss_{c}_{lo}_{hi}_median"], os_.loc[(c, lo, hi), "median_gain_pct"])
    rs = load(sd, "campi_recovery_summary.tsv")
    add("Figure2/recovery_summary.tsv", "pooled_reference (REF, 5-50)", v["campi_ref_pooled"],
        rs.loc[(rs.arm == "REF") & (rs.min_length == 5), "pooled_reference"].iloc[0])
    return pd.DataFrame(rows)
