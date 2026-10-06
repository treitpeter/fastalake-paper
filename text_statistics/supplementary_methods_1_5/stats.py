"""Statistics printed in Supplementary Methods 1-5 (v28d lines 493-659).

Every function reads tables from ``source_data/`` only and returns the numbers the text
prints, at full precision. The groups follow the order of the text:

========================================  ===========  ==============================================
statistic group                           v28d line    source_data tables
========================================  ===========  ==============================================
AlphaNovo permuted-prediction control     553          Fig2/alphanovo_accuracy_per_acquisition.tsv
bag-aware and mass-equivalent lookup      561, 579     SFig02/bags_bounded48_per_acquisition.tsv,
                                                       SFig02/bags_bounded48_orderings.tsv
healthy-stool study-supported space       567          Fig1/stage1_stats.json, hstool97_*.tsv,
                                                       search_input_completeness.tsv
parser audit                              573          parser_audit_acquisitions_726.tsv
CAMPI target databases                    609          campi_target_fasta_records.tsv
CAMPI member presence in the reservoir    611          Fig2/reservoir_protein_presence_by_proteome.tsv
search input completeness                 629          search_input_completeness.tsv
human-assigned share with missing origin  651          origin_rows_*.tsv
annotation sources and seed-ortholog      657          mgnify_seed_ortholog_check_by_genome.tsv,
check                                                  annotation_mgnify_coverage_by_arm.tsv,
                                                       campi_supplemental_annotation_queries.tsv
functional-allocation sensitivity         657          allocation_*.tsv
========================================  ===========  ==============================================

Tables owned by other units (``Fig1``, ``Fig2``, ``SFig02``) are read from their folders and
are not duplicated. ``python stats.py --write`` regenerates ``NUMBERS.tsv`` and ``STATS.tsv``.
"""
import csv
import json
import platform
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.spatial.distance import squareform
from scipy.stats import rankdata

UNIT = Path(__file__).resolve().parent
REPO = UNIT.parents[1]
SD = REPO / "source_data" / UNIT.name
SDR = REPO / "source_data"

# ---------------------------------------------------------------- choices (named constants)
ACCURACY_WINDOW = "9-50"            # reference peptide length window, v28d line 553
BAG_ARM = "bags_all"                # recorded bags without an ordering bound, line 561
BAG_BOUNDED_ARM = "bags_le8"        # bags with at most 8 orderings, line 561
BAG_TIGHT_ARM = "bags_le2"          # bags with at most 2 orderings, line 561
MASS_ARM = "mass_all"               # full mass-equivalent swap table, line 561
ORDERING_BOUND = 8                  # ordering bound of BAG_BOUNDED_ARM, line 561
PARSER_ORIGINAL_FILES = 710         # files of the original audit (indices 0-709), line 573
COMPLETENESS_THRESHOLD_TEXT = 0.95  # scan-range fraction stated in the text, line 629
COMPLETENESS_THRESHOLD_CODE = 0.90  # fraction used by the original check (sage_input_integrity.py)
SEED_ORTHOLOG_MIN_RATIO = 0.5       # DIAMOND bitscore / reported bitscore, line 657
RESTRICTED_RULES = ("unique", "complete", "unique_complete")   # allocation rules, line 657
EQUAL_RULE = "equal"                # production allocation rule, line 657
OWN_STOOL_SEED = 20260910           # alloc_capscan347.py, permutation draws
OWN_STOOL_DRAWS = 99_999            # alloc_capscan347.py
OWN_STOOL_TOL = 1e-14               # tie tolerance of the producer's P computation
PREVALENCE_FILTERS = (0.0, 0.1, 0.2)  # OG prevalence filters of the own-stool test
MANTEL_SEED = 20260907              # check_prevalence.py SEED
MANTEL_PERMUTATIONS = 1999          # check_prevalence.py NPERM
PREDICT_LAYERS = ("metaG", "metaT", "microbial_metaP", "host_metaP", "metabolites")  # evaluate() order
MICROBIAL = "microbial_metaP"
HOST = "host_metaP"
LEVELS = ("KO", "OG", "Pathway")
HSTOOL_TECHNICAL = "TechnicalReplicates"
HSTOOL_INTER = "Interindividual"

WORDS = {1: "one", 3: "three", 7: "seven", 8: "eight", 14: "fourteen", 15: "fifteen"}


def R(x):
    """Full-precision text of a number (shortest repr that round-trips)."""
    return repr(float(x))


def word(n, capital=False):
    """Spell a small count as the text does (for example ``8`` -> ``eight``)."""
    w = WORDS[int(n)]
    return w.capitalize() if capital else w


def superscript_power(p):
    """Format a power of ten as printed, for example ``1e-05`` -> ``1 × 10⁻⁵``."""
    e = int(round(np.log10(p)))
    assert np.isclose(p, 10.0 ** e, rtol=1e-12), p
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"1 × 10{str(e).translate(sup)}"


def tsv(path, **kw):
    """Read a tab-separated source_data table."""
    return pd.read_csv(path, sep="\t", **kw)


# ---------------------------------------------------------------- line 553

def stat_alphanovo_permuted():
    """Exact-match rate of permuted AlphaNovo predictions against reference identifications.

    Statistical unit: peptide-spectrum identification (reference Sage, target, rank 1,
    peptide q <= 0.01, 9-50 residues), n = 801,613 in 15 CAMPI acquisitions; one
    within-acquisition permutation of the predicted sequences per acquisition.

    Returns
    -------
    dict
        ``acquisitions``, ``reference_spectra``, ``exact_shuffled``, ``pct_shuffled`` (pooled,
        %), ``pct_shuffled_min``/``_max`` (per acquisition, %), ``pct_exact_all`` (pooled, %).
    """
    d = tsv(SDR / "Fig2/alphanovo_accuracy_per_acquisition.tsv")
    d = d[d.window == ACCURACY_WINDOW]
    pa = 100 * d.n_exact_shuffled / d.n_reference_spectra
    return dict(acquisitions=int(d["sample"].nunique()), reference_spectra=int(d.n_reference_spectra.sum()),
                exact_shuffled=int(d.n_exact_shuffled.sum()),
                pct_shuffled=100 * d.n_exact_shuffled.sum() / d.n_reference_spectra.sum(),
                pct_shuffled_min=float(pa.min()), pct_shuffled_max=float(pa.max()),
                pct_exact_all=100 * d.n_exact_all.sum() / d.n_reference_spectra.sum())


# ---------------------------------------------------------------- lines 561, 579

def stat_bags():
    """Gains and losses of peptide-acquisition identifications against exact lookup.

    Statistical unit: acquisition (n = 48: 15 CAMPI, 32 healthy stool, 1 CapScan); counts are
    summed peptide-acquisition identifications at nominal peptide q <= 0.01.

    Returns
    -------
    dict
        Per arm ``<arm>_gained``, ``_lost``, ``_net``, ``_positive`` and healthy-stool mass
        values, ordering coverage and roster counts.
    """
    d = tsv(SDR / "SFig02/bags_bounded48_per_acquisition.tsv")
    o = tsv(SDR / "SFig02/bags_bounded48_orderings.tsv")
    out = dict(acquisitions=int(d[["cohort", "sample"]].drop_duplicates().shape[0]), arm_rows=len(d),
               arms=int(d.arm.nunique()))
    for c in ["CAMPI", "HSTOOL", "CAPSCAN"]:
        out["n_" + c] = int(d.loc[d.cohort == c, "sample"].nunique())
    out["n_PREDICT"] = int(d.loc[d.cohort == "PREDICT", "sample"].nunique())
    for arm in [BAG_ARM, BAG_BOUNDED_ARM, BAG_TIGHT_ARM, MASS_ARM]:
        a = d[d.arm == arm]
        out.update({f"{arm}_gained": int(a.gained.sum()), f"{arm}_lost": int(a.lost.sum()),
                    f"{arm}_net": int(a.net.sum()), f"{arm}_positive": int((a.net > 0).sum())})
        assert (a.gained - a.lost == a.net).all()
    h = d[(d.arm == MASS_ARM) & (d.cohort == "HSTOOL")]
    out.update(hstool_mass_net=int(h.net.sum()), hstool_mass_positive=int((h.net > 0).sum()))
    el = o.drop_duplicates(["cohort", "sample"]).lookup_eligible_predictions.sum()
    assert el == o.predictions.sum()
    out.update(eligible_rows=int(el), rows_within_bound=int(o.loc[o.orderings <= ORDERING_BOUND, "predictions"].sum()))
    out["pct_within_bound"] = 100 * out["rows_within_bound"] / out["eligible_rows"]
    return out


# ---------------------------------------------------------------- line 567

def stat_hstool_spaces():
    """Study-supported protein space of the 58-acquisition and the original 97-acquisition lookup.

    Statistical unit: sequence (reservoir protein); acquisitions only define the pooled queries.

    Returns
    -------
    dict
        ``space_58``, ``space_97``, ``pooled_files_97``, ``excluded_dc`` (pooled files not in
        the 58-acquisition roster), ``reservoir``.
    """
    s58 = json.loads((SDR / "Fig1/stage1_stats.json").read_text())
    s97 = tsv(SD / "hstool97_stage1_stats.tsv", dtype=str).set_index("key").value
    files = tsv(SD / "hstool97_pooled_prediction_files.tsv").prediction_file.str.replace("_predictions.csv", "", regex=False)
    roster58 = set(tsv(SD / "search_input_completeness.tsv", dtype=str).query("lineage == 'healthy_stool_58'")["sample"])
    inside = files.isin(roster58)
    assert int(inside.sum()) == len(roster58) == 58
    assert s58["proteins_searched"] == int(s97["proteins_searched"])
    return dict(space_58=int(s58["proteins_matched"]), space_97=int(s97["proteins_matched"]),
                pooled_files_97=int(len(files)), roster_58=int(inside.sum()), excluded_dc=int((~inside).sum()),
                excluded_named_dc=int(files[~inside].str.contains("DC|Centrifug", regex=True).sum()),
                reservoir=int(s58["proteins_searched"]))


# ---------------------------------------------------------------- line 573

def stat_parser_audit():
    """Archived versus current Rust parser on every retained prediction file.

    Statistical unit: acquisition (prediction file); n = 710 in v28d, 726 with the 16
    participant-15 CapScan files added for v28f.

    Returns
    -------
    dict
        files and rows for the 710 and 726 sets, all-files-identical flags, and the number of
        lookup queries absent from selection evidence.
    """
    d = tsv(SD / "parser_audit_acquisitions_726.tsv")
    o = d[d["index"] < PARSER_ORIGINAL_FILES]
    return dict(files_710=len(o), rows_710=int(o.rows.sum()), files_726=len(d), rows_726=int(d.rows.sum()),
                lookup_identical_726=bool(d.lookup_rows_exact.all()), selection_identical_726=bool(d.selection_exact.all()),
                lookup_only_queries_710=int(o.lookup_only_queries.sum()), lookup_only_queries_726=int(d.lookup_only_queries.sum()),
                lookup_only_files_726=int((d.lookup_only_queries > 0).sum()))


# ---------------------------------------------------------------- line 609

def stat_campi_targets():
    """Records and contaminants in the supplied CAMPI reference and MG target FASTAs.

    Statistical unit: sequence (FASTA record).

    Returns
    -------
    dict
        ``reference_records``, ``mg_records`` (records in the file), ``mg_distinct_records``
        (distinct header token and sequence), contaminant counts.
    """
    d = tsv(SD / "campi_target_fasta_records.tsv", dtype=str)
    r, m = d[d.database == "reference"], d[d.database == "MG"]
    cr = r[r.header_kind == "other"]
    cm = m[m.header_kind == "other"]
    return dict(reference_records=len(r), reference_distinct=int(r[["header_token", "sequence_sha256"]].drop_duplicates().shape[0]),
                reference_contaminants=int(cr.header_token.nunique()), reference_uniprot=int((r.header_kind == "uniprot").sum()),
                mg_records=len(m), mg_distinct_records=int(m[["header_token", "sequence_sha256"]].drop_duplicates().shape[0]),
                mg_contaminants=int(cm.header_token.nunique()), mg_contaminant_records=len(cm),
                shared_contaminants=len(set(zip(cr.header_token, cr.sequence_sha256)) & set(zip(cm.header_token, cm.sequence_sha256))))


# ---------------------------------------------------------------- line 611

def stat_member_presence():
    """Exact presence of the supplied member reference proteomes in the 375 M reservoir.

    Statistical unit: protein (member reference record); n = 29,557 in eight proteomes.

    Returns
    -------
    dict
        ``members``, ``member_proteins``, ``pct_min``, ``pct_max``, ``shuffled_present``.
    """
    d = tsv(SDR / "Fig2/reservoir_protein_presence_by_proteome.tsv")
    real = d[~d.label.str.endswith("_SHUFFLED")]
    pct = 100 * real.n_present_exact / real.n_proteins
    return dict(members=len(real), member_proteins=int(real.n_proteins.sum()), pct_min=float(pct.min()),
                pct_max=float(pct.max()), shuffled_present=int(d.loc[d.label.str.endswith("_SHUFFLED"), "n_present_exact"].sum()))


# ---------------------------------------------------------------- line 629

def stat_completeness():
    """Sage input completeness: closing tag and highest scored scan / declared spectra.

    Statistical unit: acquisition (search); n = 653 public-sequence (306 PREDICT, 347 CapScan)
    and 58 healthy-stool searches.

    Returns
    -------
    dict
        counts per group, lowest coverage (%), and pass counts at the text (95%) and code (90%)
        thresholds.
    """
    d = tsv(SD / "search_input_completeness.tsv")
    d["coverage"] = d.max_scan / d.declared_spectra
    pub = d[d.lineage == "primary_public_current"]
    hs = d[d.lineage == "healthy_stool_58"]
    ok_tag = d.closing_tag.astype(str).eq("True")
    return dict(public=len(pub), public_PREDICT=int((pub.cohort == "PREDICT").sum()), public_CAPSCAN=int((pub.cohort == "CAPSCAN").sum()),
                hstool=len(hs), all_closing_tag=bool(ok_tag.all()), min_coverage_pct=100 * float(d.coverage.min()),
                pass_text_threshold=int((ok_tag & (d.coverage >= COMPLETENESS_THRESHOLD_TEXT)).sum()),
                pass_code_threshold=int((ok_tag & (d.coverage >= COMPLETENESS_THRESHOLD_CODE)).sum()), searches=len(d),
                producer_flags=int(d.flag.astype(str).eq("True").sum()))


# ---------------------------------------------------------------- line 651

def stat_partial_human():
    """Share of pooled total LFQ intensity that is human-assigned although some compatible
    proteins lacked an origin annotation.

    Statistical unit: LFQ feature intensity pooled over all acquisitions (58 healthy stool,
    306 PREDICT, 331 or 347 CapScan); origin rows are the annotated-origin classes.

    Returns
    -------
    dict
        ``pct_<cohort>`` and ``pct_human_<cohort>`` for HSTOOL, PREDICT, CAPSCAN331, CAPSCAN347.
    """
    out = {}
    for c in ["HSTOOL", "PREDICT", "CAPSCAN331", "CAPSCAN347"]:
        d = tsv(SD / f"origin_rows_{c}.tsv")
        tot = d.signal_pooled.sum()
        h = d.human_species.astype(bool)
        out["pct_" + c] = 100 * float(d.loc[h & d.state.eq("annotated_parent_lca"), "signal_pooled"].sum() / tot)
        out["pct_human_" + c] = 100 * float(d.loc[h, "signal_pooled"].sum() / tot)
        out["acquisitions_" + c] = int(d.acquisitions.iloc[0])
    return out


# ---------------------------------------------------------------- line 657, annotations

def stat_annotations():
    """MGnify seed-ortholog check, MGnify coverage and the supplementary CAMPI annotation run.

    Statistical unit: protein (MGnify gene annotation row; compatible protein).

    Returns
    -------
    dict
        ``rows``, ``failed``, ``no_alignment`` and their percentages, ``genomes``, the shares of
        failures with a short query or an adjacent-gene seed, ``campi_queries`` and the
        lowest/highest MGnify coverage of compatible proteins (v28f roster).
    """
    g = tsv(SD / "mgnify_seed_ortholog_check_by_genome.tsv")
    rows, fail, na = int(g.rows.sum()), int(g.fail_ratio_below_half.sum()), int(g.no_alignment.sum())
    cov = tsv(SD / "annotation_mgnify_coverage_by_arm.tsv")
    v28f = pd.concat([cov[(cov.roster == "v28d") & (cov.arm != "CAPSCAN_noclust")], cov[cov.roster == "347"]])
    v28d = cov[cov.roster == "v28d"]
    pf = 100 * v28f.native_gene_sequence_verified / v28f.compatible_parents
    pd_ = 100 * v28d.native_gene_sequence_verified / v28d.compatible_parents
    q = tsv(SD / "campi_supplemental_annotation_queries.tsv")
    return dict(rows=rows, failed=fail, no_alignment=na, pct_failed=100 * fail / rows, pct_no_alignment=100 * na / rows,
                genomes=len(g), failed_short_query_pct=100 * g.fail_query_shorter_than_half_seed.sum() / fail,
                failed_adjacent_seed_pct=100 * g.fail_seed_equals_adjacent.sum() / fail,
                passed_short_query_pct=100 * g.pass_query_shorter_than_half_seed.sum() / (rows - fail),
                passed_adjacent_seed_pct=100 * g.pass_seed_equals_adjacent.sum() / (rows - fail),
                campi_queries=len(q), campi_annotated=int(q.has_eggnog_record.sum()),
                coverage_min_pct=float(pf.min()), coverage_max_pct=float(pf.max()),
                coverage_min_pct_v28d=float(pd_.min()), coverage_max_pct_v28d=float(pd_.max()))


# ---------------------------------------------------------------- line 657, allocation sensitivity

def stat_allocated_fraction():
    """Median share of total LFQ intensity allocated to functions under each allocation rule.

    Statistical unit: acquisition (58 healthy stool, 306 PREDICT, 346 non-empty CapScan of
    347); medians per cohort, level and rule, profile excluding human-assigned signal.

    Returns
    -------
    tuple of (pandas.DataFrame, dict)
        Medians (``arm, level, variant, allocated_fraction``) and the printed ranges in %.
    """
    rows = []
    for name in ["HSTOOL", "PREDICT", "CAPSCAN347"]:
        c = tsv(SD / f"allocation_coverage_{name}.tsv")
        c = c[c.profile == "other_unresolved"].copy()
        c["allocated_fraction"] = c.allocated / c.total_LFQ.where(c.total_LFQ > 0)
        rows.append(c.groupby(["arm", "level", "variant"]).allocated_fraction.median().reset_index())
    med = pd.concat(rows, ignore_index=True)
    eq, re_ = med[med.variant == EQUAL_RULE].allocated_fraction, med[med.variant.isin(RESTRICTED_RULES)].allocated_fraction
    return med, dict(restricted_min_pct=100 * re_.min(), restricted_max_pct=100 * re_.max(),
                     equal_min_pct=100 * eq.min(), equal_max_pct=100 * eq.max())


def own_stool_test(dd, seed=OWN_STOOL_SEED, draws=OWN_STOOL_DRAWS):
    """Permutation test that intestinal samples are closer to their own participant's stool.

    Port of ``stats()`` in ``alloc_capscan347.py`` starting from its distance matrix.

    Parameters
    ----------
    dd : numpy.ndarray
        n x n matrix; ``dd[i, j]`` is the mean Hellinger distance between the four
        intestinal-region profiles of participant i and the stool profile of participant j.
    seed, draws : int
        Random draws of stool-label permutations (``numpy.random.default_rng(seed)``).

    Returns
    -------
    dict
        ``own_closer`` (participants), ``mean_own``, ``mean_other``, ``mean_advantage`` and the
        two-sided permutation ``p`` = (1 + #|null| >= |observed|) / (draws + 1).
        Statistical unit: participant, n = dd.shape[0].
    """
    n = dd.shape[0]
    perms = np.argsort(np.random.default_rng(seed).random((draws, n)), axis=1)
    own = np.diag(dd); other = (dd.sum(1) - own) / (n - 1)
    obs = float((other - own).mean())
    pdg = dd[np.arange(n)[None, :], perms].mean(1)
    null = (dd.sum(1).mean() - n * pdg) / (n - 1)
    p = (1 + int((abs(null) >= abs(obs) - OWN_STOOL_TOL).sum())) / (draws + 1)
    return dict(own_closer=int((own < other).sum()), mean_own=float(own.mean()), mean_other=float(other.mean()),
                mean_advantage=obs, p=p, participants=n)


def stat_own_stool():
    """CapScan own-stool test under each allocation rule and OG prevalence filter (347 roster).

    Statistical unit: participant, n = 15.

    Returns
    -------
    pandas.DataFrame
        One row per rule and prevalence filter (plus the rotated-stool negative control).
    """
    d = tsv(SD / "allocation_capscan_own_stool_distances_347.tsv")
    rows = []
    for (v, prev), g in d.groupby(["variant", "prevalence"], sort=False):
        parts = sorted(g.intestinal_participant.unique())
        dd = g.pivot(index="intestinal_participant", columns="stool_participant", values="mean_hellinger_distance").loc[parts, parts].to_numpy()
        rows.append(dict(variant=v, prevalence=prev, features=int(g.features.iloc[0]), **own_stool_test(dd)))
    return pd.DataFrame(rows)


def mantel(dist_a, dist_b, perms):
    """Spearman Mantel correlation of two condensed distance vectors, two-sided permutation P.

    Port of ``evaluate()`` in ``check_prevalence.py``: distances are ranked, centred and scaled
    to unit length; rho is their dot product; the null permutes the units of ``dist_b``.

    Parameters
    ----------
    dist_a, dist_b : numpy.ndarray
        Condensed distance vectors (``scipy.spatial.distance.pdist`` order) over the same n units.
    perms : numpy.ndarray
        Permutations of range(n), one per row.

    Returns
    -------
    tuple of (float, float)
        rho and P = (1 + #|null| >= |rho|) / (1 + number of permutations).
    """
    def rk(x):
        r = rankdata(x); r -= r.mean(); r /= np.linalg.norm(r); return r
    ra, rb = rk(dist_a), rk(dist_b)
    rho = float(ra @ rb)
    n = squareform(rb).shape[0]
    ix = np.triu_indices(n, 1); bm = squareform(rb); ext = 0
    for start in range(0, len(perms), 200):
        pp = perms[start:start + 200]
        ext += np.count_nonzero(np.abs(bm[pp[:, ix[0]], pp[:, ix[1]]] @ ra) >= abs(rho))
    return rho, (1 + ext) / (1 + len(perms))


def stat_predict_mantel():
    """PREDICT protein-profile distance correlations under each allocation rule.

    Statistical unit: patient (earliest visit with complete molecular measurements), n = 77.

    Returns
    -------
    pandas.DataFrame
        ``variant, level, layer_a, layer_b, rho, p`` for every pair that involves the microbial
        protein profile, in the producer's order.
    """
    d = tsv(SD / "allocation_predict_distances.tsv")
    n = int(d.unit_j.str[1:].astype(int).max())
    rng = np.random.default_rng(MANTEL_SEED)
    perms = np.array([rng.permutation(n) for _ in range(MANTEL_PERMUTATIONS)])
    get = {k: g.distance.to_numpy() for k, g in d.groupby(["level", "layer", "variant"], sort=False)}
    rows = []
    for v in [EQUAL_RULE, *RESTRICTED_RULES, "NEGATIVE_CONTROL_protein_patients_permuted"]:
        for lev in LEVELS:
            vec = {k: get[(lev, k, "all")] for k in PREDICT_LAYERS if k != MICROBIAL}
            vec[MICROBIAL] = get[(lev, MICROBIAL, v)]
            for a, b in combinations(PREDICT_LAYERS, 2):
                if lev != "OG" and {a, b} == {HOST, "metabolites"}:
                    continue
                if MICROBIAL not in (a, b):
                    continue
                rho, p = mantel(vec[a], vec[b], perms)
                rows.append(dict(variant=v, level=lev, layer_a=a, layer_b=b, n_units=n, rho=rho, p=p))
    return pd.DataFrame(rows)


def stat_predict_ranges(m=None):
    """Printed rho and P ranges of the PREDICT sensitivity (from ``stat_predict_mantel``)."""
    m = stat_predict_mantel() if m is None else m
    host = (m.layer_a == HOST) | (m.layer_b == HOST)
    eq, rs = m.variant == EQUAL_RULE, m.variant.isin(RESTRICTED_RULES)
    return dict(other_restricted_min=m[rs & ~host].rho.min(), other_restricted_max=m[rs & ~host].rho.max(),
                other_equal_min=m[eq & ~host].rho.min(), other_equal_max=m[eq & ~host].rho.max(),
                other_p_max=m[(rs | eq) & ~host].p.max(),
                host_equal_min=m[eq & host].rho.min(), host_equal_max=m[eq & host].rho.max(),
                host_restricted_min=m[rs & host].rho.min(), host_restricted_max=m[rs & host].rho.max(),
                host_restricted_p_max=m[rs & host].p.max(),
                control_p_min=m[m.variant.str.startswith("NEGATIVE")].p.min())


def stat_hstool_r():
    """Healthy-stool KO correlation by design tier under each allocation rule.

    Statistical unit: acquisition pair (66 technical-replicate pairs from 12 acquisitions, 45
    inter-individual pairs from 10); value = median over pairs of the Pearson r of log2 KO
    intensity over KOs positive in both acquisitions. Descriptive; no test.

    Returns
    -------
    tuple of (pandas.DataFrame, dict)
        Medians per profile, rule and tier, and the printed ranges (both profiles pooled).
    """
    d = tsv(SD / "allocation_hstool_ko_pairs.tsv")
    med = d.groupby(["profile", "variant", "design"], sort=False).r_per_pair.median().reset_index()
    real = med[~med.variant.str.startswith("NEGATIVE")]
    t = real[real.design == HSTOOL_TECHNICAL]; i = real[real.design == HSTOOL_INTER]
    ie = i[i.variant == EQUAL_RULE].set_index("profile").r_per_pair
    ir = i[i.variant.isin(RESTRICTED_RULES)].r_per_pair
    return med, dict(technical_min=t.r_per_pair.min(), technical_max=t.r_per_pair.max(),
                     inter_equal_all=float(ie["all"]), inter_equal_other_unresolved=float(ie["other_unresolved"]),
                     inter_restricted_min=ir.min(), inter_restricted_max=ir.max(),
                     pairs_technical=int((d.design == HSTOOL_TECHNICAL).sum() / d.groupby(["profile", "variant"]).ngroups),
                     pairs_inter=int((d.design == HSTOOL_INTER).sum() / d.groupby(["profile", "variant"]).ngroups))


# ---------------------------------------------------------------- comparison with producer outputs

def check_against_producers():
    """Compare every recomputed allocation-sensitivity table with the producer's own output.

    Returns
    -------
    dict
        Maximum absolute differences; raises AssertionError on any disagreement.
    """
    out = {}
    med, _ = stat_allocated_fraction()
    ref = tsv(SD / "allocation_fraction_medians_producer.tsv")
    m = med.merge(ref, on=["arm", "level", "variant"], suffixes=("", "_ref"))
    assert len(m) == len(ref) == 36
    out["allocated_fraction_max_abs_diff"] = float((m.allocated_fraction - m.allocated_fraction_ref).abs().max())
    own = stat_own_stool()
    ref = tsv(SD / "allocation_capscan_own_stool_producer.tsv")
    m = own.merge(ref, on=["variant", "prevalence"])
    assert len(m) == 13 and (m.own_closer == m.own_closer_347).all() and (m.features == m.features_347).all()
    out["own_stool_p_max_abs_diff"] = float((m.p - m.p_347).abs().max())
    out["own_stool_advantage_max_abs_diff"] = float((m.mean_advantage - m.mean_advantage_347).abs().max())
    pm = stat_predict_mantel()
    ref = tsv(SD / "allocation_predict_crossomics_producer.tsv")
    m = pm.merge(ref, on=["variant", "level", "layer_a", "layer_b"], suffixes=("", "_ref"))
    assert len(m) == len(ref) == 60
    out["predict_rho_max_abs_diff"] = float((m.rho - m.rho_ref).abs().max())
    out["predict_p_max_abs_diff"] = float((m.p - m.p_ref).abs().max())
    hm, _ = stat_hstool_r()
    ref = tsv(SD / "allocation_hstool_ko_precision_producer.tsv")
    m = hm.merge(ref, on=["profile", "variant", "design"], suffixes=("", "_ref"))
    assert len(m) == len(ref) == 45
    out["hstool_r_max_abs_diff"] = float((m.r_per_pair - m.r_per_pair_ref).abs().max())
    for k, v in out.items():
        assert v <= 1e-12, (k, v)
    return out


# ---------------------------------------------------------------- NUMBERS.tsv and STATS.tsv

def all_numbers():
    """Return every NUMBERS.tsv row as a list of dicts, in text order."""
    a = stat_alphanovo_permuted(); bg = stat_bags(); hs = stat_hstool_spaces(); pa = stat_parser_audit()
    ct = stat_campi_targets(); mp = stat_member_presence(); co = stat_completeness(); ph = stat_partial_human()
    an = stat_annotations(); _, al = stat_allocated_fraction(); own = stat_own_stool(); pr = stat_predict_ranges()
    _, hr = stat_hstool_r()
    real_own = own[own.variant.isin([EQUAL_RULE, *RESTRICTED_RULES])]
    S = lambda f: f"text_statistics/{UNIT.name}/stats.py:{f}"
    D = lambda t: f"source_data/{UNIT.name}/{t}"
    rows = []

    def row(line, quote, printed, value, code, sd, status="MATCH"):
        rows.append(dict(v28d_line=str(line), quote=quote, printed_value=printed, recomputed_value=R(value) if isinstance(value, (float, np.floating)) else str(value),
                         code=code, source_data=sd, status=status))

    acc = "source_data/Fig2/alphanovo_accuracy_per_acquisition.tsv"
    row(553, "each reference-database Sage identification in the 15 CAMPI acquisitions", "15", a["acquisitions"], S("stat_alphanovo_permuted"), acc)
    row(553, "reduced the match rate to 0.006%", "0.006%", a["pct_shuffled"], S("stat_alphanovo_permuted"), acc)
    bsd = "source_data/SFig02/bags_bounded48_per_acquisition.tsv"
    osd = "source_data/SFig02/bags_bounded48_orderings.tsv"
    q1 = "Over 48 acquisitions (15 CAMPI, 32 healthy stool and one CapScan, with no PREDICT)"
    row(561, q1, "48", bg["acquisitions"], S("stat_bags"), bsd)
    row(561, q1, "15", bg["n_CAMPI"], S("stat_bags"), bsd)
    row(561, q1, "32", bg["n_HSTOOL"], S("stat_bags"), bsd)
    row(561, q1, "one", word(bg["n_CAPSCAN"]), S("stat_bags"), bsd)
    q2 = "admitting recorded bags added 6,961 peptide–acquisition events and removed 5,941 (net +1,020, positive in 38 of 48)"
    row(561, q2, "6,961", bg["bags_all_gained"], S("stat_bags"), bsd)
    row(561, q2, "5,941", bg["bags_all_lost"], S("stat_bags"), bsd)
    row(561, q2, "+1,020", bg["bags_all_net"], S("stat_bags"), bsd)
    row(561, q2, "38", bg["bags_all_positive"], S("stat_bags"), bsd)
    q3 = "which covers 908,037 of 976,861 eligible prediction rows (93.0%)"
    row(561, q3, "908,037", bg["rows_within_bound"], S("stat_bags"), osd)
    row(561, q3, "976,861", bg["eligible_rows"], S("stat_bags"), osd)
    row(561, q3, "93.0%", bg["pct_within_bound"], S("stat_bags"), osd)
    q4 = "(+5,051, −3,318, net +1,733, positive in 37 of 48)"
    row(561, q4, "+5,051", bg["bags_le8_gained"], S("stat_bags"), bsd)
    row(561, q4, "−3,318", -bg["bags_le8_lost"], S("stat_bags"), bsd)
    row(561, q4, "+1,733", bg["bags_le8_net"], S("stat_bags"), bsd)
    row(561, q4, "37", bg["bags_le8_positive"], S("stat_bags"), bsd)
    row(561, "Restricting it to at most 2 orderings was net negative (−793)", "−793", bg["bags_le2_net"], S("stat_bags"), bsd)
    q5 = "added 31,111 and removed 33,306 events (net −2,195, positive in 13 of 48)"
    row(561, q5, "31,111", bg["mass_all_gained"], S("stat_bags"), bsd)
    row(561, q5, "33,306", bg["mass_all_lost"], S("stat_bags"), bsd)
    row(561, q5, "−2,195", bg["mass_all_net"], S("stat_bags"), bsd)
    row(561, q5, "13", bg["mass_all_positive"], S("stat_bags"), bsd)
    q6 = "In the 32 healthy-stool acquisitions the mass arm was net −5,242 and positive in 3"
    row(561, q6, "−5,242", bg["hstool_mass_net"], S("stat_bags"), bsd)
    row(561, q6, "3", bg["hstool_mass_positive"], S("stat_bags"), bsd)
    row(579, "The 48-acquisition comparison contained seven search arms per acquisition", "seven", word(bg["arms"]), S("stat_bags"), bsd)
    row(579, "All 336 arm–acquisition rows", "336", bg["arm_rows"], S("stat_bags"), bsd)
    row(567, "produced a study-supported protein space of 18,007,876 proteins", "18,007,876", hs["space_58"], S("stat_hstool_spaces"), "source_data/Fig1/stage1_stats.json")
    q7 = "The original production run pooled 97 acquisitions, including 39 differential-centrifugation acquisitions"
    row(567, q7, "97", hs["pooled_files_97"], S("stat_hstool_spaces"), D("hstool97_pooled_prediction_files.tsv"))
    row(567, q7, "39", hs["excluded_dc"], S("stat_hstool_spaces"), D("hstool97_pooled_prediction_files.tsv"))
    row(567, "and produced 26,881,038 study-supported proteins", "26,881,038", hs["space_97"], S("stat_hstool_spaces"), D("hstool97_stage1_stats.tsv"))
    q8 = "compared all 710 retained prediction files (112,835,982 rows)"
    pad = D("parser_audit_acquisitions_726.tsv")
    row(573, q8, "710", pa["files_726"], S("stat_parser_audit"), pad, "v28f_change:726")
    row(573, q8, "112,835,982", pa["rows_726"], S("stat_parser_audit"), pad, "v28f_change:114,971,345")
    row(573, "Three lookup sequences absent from selection evidence", "Three", word(pa["lookup_only_queries_726"], True), S("stat_parser_audit"), pad)
    q9 = "The searched target FASTAs contained 29,673 reference and 19,435 MG records, each including 116 contaminants."
    cad = D("campi_target_fasta_records.tsv")
    row(609, q9, "29,673", ct["reference_records"], S("stat_campi_targets"), cad)
    row(609, q9, "19,435", ct["mg_distinct_records"], S("stat_campi_targets"), cad)
    row(609, q9, "116", ct["reference_contaminants"], S("stat_campi_targets"), cad)
    q10 = "For each of the eight supplied member reference proteomes (29,557 member proteins)"
    mpd = "source_data/Fig2/reservoir_protein_presence_by_proteome.tsv"
    row(611, q10, "eight", word(mp["members"]), S("stat_member_presence"), mpd)
    row(611, q10, "29,557", mp["member_proteins"], S("stat_member_presence"), mpd)
    q11 = "for the 637 public-sequence CapScan and PREDICT searches and the 58 acquisition-specific healthy-stool searches"
    cod = D("search_input_completeness.tsv")
    row(629, q11, "637", co["public"], S("stat_completeness"), cod, "v28f_change:653")
    row(629, q11, "58", co["hstool"], S("stat_completeness"), cod)
    row(629, "All passed (lowest coverage 99.6%)", "99.6%", co["min_coverage_pct"], S("stat_completeness"), cod)
    q12 = "0.12% (healthy stool), 0.60% (CapScan) and 1.96% (PREDICT)"
    row(651, q12, "0.12%", ph["pct_HSTOOL"], S("stat_partial_human"), D("origin_rows_HSTOOL.tsv"))
    row(651, q12, "0.60%", ph["pct_CAPSCAN347"], S("stat_partial_human"), D("origin_rows_CAPSCAN347.tsv"), "v28f_change:0.61%")
    row(651, q12, "1.96%", ph["pct_PREDICT"], S("stat_partial_human"), D("origin_rows_PREDICT.tsv"))
    mg = D("mgnify_seed_ortholog_check_by_genome.tsv")
    row(657, "(4,744 species representatives", "4,744", an["genomes"], S("stat_annotations"), mg)
    row(657, "run on 17,935 CAMPI reference proteins", "17,935", an["campi_queries"], S("stat_annotations"), D("campi_supplemental_annotation_queries.tsv"))
    q13 = "Of 9,368,708 MGnify annotations, 372,277 (4.0%) failed this check, including 191,565 (2.0%) with no alignment."
    row(657, q13, "9,368,708", an["rows"], S("stat_annotations"), mg)
    row(657, q13, "372,277", an["failed"], S("stat_annotations"), mg)
    row(657, q13, "4.0%", an["pct_failed"], S("stat_annotations"), mg)
    row(657, q13, "191,565", an["no_alignment"], S("stat_annotations"), mg)
    row(657, q13, "2.0%", an["pct_no_alignment"], S("stat_annotations"), mg)
    q14 = "MGnify tables supplied annotations for 88–95% of compatible proteins"
    row(657, q14, "88", an["coverage_min_pct"], S("stat_annotations"), D("annotation_mgnify_coverage_by_arm.tsv"))
    row(657, q14, "95%", an["coverage_max_pct"], S("stat_annotations"), D("annotation_mgnify_coverage_by_arm.tsv"))
    q15 = "These rules allocated a median of 7–61% instead of 62–81% of total intensity"
    acd = D("allocation_coverage_*.tsv")
    row(657, q15, "7", al["restricted_min_pct"], S("stat_allocated_fraction"), acd)
    row(657, q15, "61%", al["restricted_max_pct"], S("stat_allocated_fraction"), acd)
    row(657, q15, "62", al["equal_min_pct"], S("stat_allocated_fraction"), acd)
    row(657, q15, "81%", al["equal_max_pct"], S("stat_allocated_fraction"), acd)
    q16 = "All fourteen CapScan participants remained closer to their own stool (P≤1 × 10⁻⁵)"
    osd2 = D("allocation_capscan_own_stool_distances_347.tsv")
    assert (real_own.own_closer == real_own.participants).all()
    row(657, q16, "fourteen", word(int(real_own.own_closer.min())), S("stat_own_stool"), osd2, "v28f_change:fifteen")
    row(657, q16, "1 × 10⁻⁵", superscript_power(real_own.p.max()), S("stat_own_stool"), osd2)
    q17 = "ρ=0.16–0.64 against 0.16–0.51 with equal allocation, all P≤0.026"
    pdd = D("allocation_predict_distances.tsv")
    row(657, q17, "0.16", pr["other_restricted_min"], S("stat_predict_ranges"), pdd)
    row(657, q17, "0.64", pr["other_restricted_max"], S("stat_predict_ranges"), pdd)
    row(657, q17, "0.16", pr["other_equal_min"], S("stat_predict_ranges"), pdd)
    row(657, q17, "0.51", pr["other_equal_max"], S("stat_predict_ranges"), pdd)
    row(657, q17, "0.026", pr["other_p_max"], S("stat_predict_ranges"), pdd)
    q18 = "fell from 0.43–0.51 to 0.28–0.35 (P≤0.002)"
    row(657, q18, "0.43", pr["host_equal_min"], S("stat_predict_ranges"), pdd)
    row(657, q18, "0.51", pr["host_equal_max"], S("stat_predict_ranges"), pdd)
    row(657, q18, "0.28", pr["host_restricted_min"], S("stat_predict_ranges"), pdd)
    row(657, q18, "0.35", pr["host_restricted_max"], S("stat_predict_ranges"), pdd)
    row(657, q18, "0.002", pr["host_restricted_p_max"], S("stat_predict_ranges"), pdd)
    hpd = D("allocation_hstool_ko_pairs.tsv")
    row(657, "(r=0.89–0.91)", "0.89", hr["technical_min"], S("stat_hstool_r"), hpd)
    row(657, "(r=0.89–0.91)", "0.91", hr["technical_max"], S("stat_hstool_r"), hpd)
    q19 = "the inter-individual value fell from 0.81 to 0.69–0.73"
    row(657, q19, "0.81", hr["inter_equal_all"], S("stat_hstool_r"), hpd)
    row(657, q19, "0.69", hr["inter_restricted_min"], S("stat_hstool_r"), hpd)
    row(657, q19, "0.73", hr["inter_restricted_max"], S("stat_hstool_r"), hpd)
    return rows


def software():
    """Interpreter and package versions used by this module."""
    return f"CPython {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, pandas {pd.__version__}"


def all_stats():
    """Return every STATS.tsv row as a list of dicts."""
    S = lambda f: f"text_statistics/{UNIT.name}/stats.py:{f}"
    sw = software()
    rows = []
    blank = dict(covariates="none", correction="none", family="NA", family_size="NA", uncertainty="NA",
                 uncertainty_construction="NA", statistic="NA", df="NA", exact_P="NA", effect_size="NA",
                 effect_size_how="NA", draws="NA", seed="NA", software_version=sw)
    a = stat_alphanovo_permuted()
    rows.append({**blank, "claim": "Permuting AlphaNovo predictions across spectra reduced the exact-match rate to 0.006%",
                 "v28d_line": "553", "test": "permutation negative control (predicted sequences shuffled across spectra within acquisition; no test)",
                 "sided": "NA", "unit": "peptide", "n": f"{a['reference_spectra']} reference identifications in {a['acquisitions']} acquisitions",
                 "distinct_or_repeated": "repeated: spectra within acquisition",
                 "estimate": f"pooled exact-match rate after permutation {R(a['pct_shuffled'])}% ({a['exact_shuffled']} of {a['reference_spectra']})",
                 "uncertainty": f"per-acquisition range {R(a['pct_shuffled_min'])}% to {R(a['pct_shuffled_max'])}%",
                 "uncertainty_construction": "range over the 15 acquisitions",
                 "draws": "1 permutation per acquisition",
                 "seed": "20261005 + sum of character codes of the acquisition name (Python random.Random.shuffle)",
                 "software_version": "producer: CPython 3.10.9 (anaconda 2023.03); recomputation here: " + sw,
                 "code": S("stat_alphanovo_permuted") + " (permutation done upstream: campi_alphanovo_accuracy accuracy.py:run)"})
    own = stat_own_stool()
    for r in own.itertuples():
        neg = r.variant.startswith("NEGATIVE")
        rows.append({**blank, "claim": ("Negative control: stool profiles rotated across participants" if neg else
                                        f"CapScan intestinal samples closer to own stool, allocation rule '{r.variant}', OG prevalence >= {r.prevalence}"),
                     "v28d_line": "657", "test": "permutation test of the own-stool advantage (stool labels permuted across participants)",
                     "sided": "two-sided", "unit": "participant", "n": str(r.participants),
                     "distinct_or_repeated": "repeated: four intestinal regions and stool within participant; acquisitions averaged within specimen and specimens within participant-region before the test",
                     "family": "NA" if neg else "own-stool sensitivity tests (4 allocation rules x 3 OG prevalence filters)",
                     "family_size": "NA" if neg else "12",
                     "estimate": f"mean advantage {R(r.mean_advantage)} (Hellinger distance); own stool closer for {r.own_closer} of {r.participants}",
                     "statistic": f"mean advantage = {R(r.mean_advantage)}", "exact_P": R(r.p),
                     "effect_size": R(r.mean_advantage),
                     "effect_size_how": "mean over participants of (mean Hellinger distance from the participant's four intestinal-region OG profiles to other participants' stool minus to own stool)",
                     "draws": str(OWN_STOOL_DRAWS), "seed": str(OWN_STOOL_SEED), "code": S("own_stool_test")})
    pm = stat_predict_mantel()
    for r in pm.itertuples():
        neg = r.variant.startswith("NEGATIVE")
        rows.append({**blank, "claim": (f"Negative control: protein profiles permuted across patients, {r.level} {r.layer_a} vs {r.layer_b}" if neg else
                                        f"PREDICT {r.level} {r.layer_a} vs {r.layer_b} profile distances, allocation rule '{r.variant}'"),
                     "v28d_line": "657", "test": "Mantel test, Spearman correlation of Euclidean distances (Hellinger-transformed profiles; metabolites log10, >=70% present, median-imputed, z-scored)",
                     "sided": "two-sided", "unit": "patient", "n": str(r.n_units), "distinct_or_repeated": "distinct",
                     "family": "NA" if neg else f"PREDICT protein-profile rows under allocation rule '{r.variant}' (3 levels x 4 layer pairs)",
                     "family_size": "NA" if neg else "12",
                     "estimate": f"Mantel rho {R(r.rho)}", "statistic": f"rho = {R(r.rho)}", "exact_P": R(r.p),
                     "effect_size": R(r.rho), "effect_size_how": "Spearman Mantel rho",
                     "draws": str(MANTEL_PERMUTATIONS), "seed": f"{MANTEL_SEED} (numpy default_rng, permutations of patients)",
                     "code": S("mantel")})
    d = tsv(SD / "allocation_hstool_ko_pairs.tsv")
    for (prof, var, des), g in d[d.design.isin([HSTOOL_TECHNICAL, HSTOOL_INTER])].groupby(["profile", "variant", "design"], sort=False):
        q1, q3 = np.percentile(g.r_per_pair, [25, 75])
        rows.append({**blank, "claim": f"Healthy-stool KO correlation, {des}, profile '{prof}', allocation rule '{var}'",
                     "v28d_line": "657", "test": "none (descriptive)", "sided": "NA", "unit": "acquisition",
                     "n": f"{len(g)} acquisition pairs", "distinct_or_repeated": "repeated: each acquisition enters several pairs",
                     "estimate": f"median Pearson r of log2 KO intensity {R(g.r_per_pair.median())}",
                     "uncertainty": f"IQR {R(q1)} to {R(q3)}", "uncertainty_construction": "numpy percentile (linear) over pairs",
                     "code": S("stat_hstool_r")})
    return rows


NUMBER_COLUMNS = ["v28d_line", "quote", "printed_value", "recomputed_value", "code", "source_data", "status"]
STAT_COLUMNS = ["claim", "v28d_line", "test", "sided", "unit", "n", "distinct_or_repeated", "covariates", "correction",
                "family", "family_size", "estimate", "uncertainty", "uncertainty_construction", "statistic", "df", "exact_P",
                "effect_size", "effect_size_how", "draws", "seed", "software_version", "code"]


def write_tables():
    """Write NUMBERS.tsv and STATS.tsv beside this module."""
    for name, rows, cols in [("NUMBERS.tsv", all_numbers(), NUMBER_COLUMNS), ("STATS.tsv", all_stats(), STAT_COLUMNS)]:
        with open(UNIT / name, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(rows)
        print("wrote", name, len(rows))


if __name__ == "__main__":
    if "--write" in sys.argv:
        print(check_against_producers())
        write_tables()
