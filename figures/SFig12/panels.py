"""Supplementary Figure 12 (v30): functional allocation and leading-species attribution.

Panels
------
a  fraction of accepted feature intensity allocated to a KEGG ortholog (KO), per
   acquisition, for healthy stool (HSTOOL), PREDICT and CapScan
b  healthy stool: expected number of distinct leading species per function after
   rarefaction to eight resolved acquisitions, per acquisition design
c  CapScan: per-participant cumulative distribution of the leader share within
   species-resolved signal, over eligible functions
d  CapScan: the same for the species-resolved share of non-host functional signal
e  healthy stool: leader agreement under two resolution filters and two designs, direct
   peptide assignment only
f  CapScan: donor-balanced leader agreement (intestinal region versus own stool) under
   four eligibility rules, direct peptide assignment only

All panels use direct peptide assignment; the original protein-grouping series of the v28f
panels e and f was removed in v30 (its rows stay in the two sensitivity tables, which are
Source Data, but are not drawn).

Source data
-----------
source_data/SFig12/<COHORT>_KO_all.tsv                           panel a (KO x acquisition matrices)
source_data/SFig12/KO_allocation_coverage.tsv                    panel a (derived per-acquisition table)
source_data/SFig12/healthy_stool_leader_diversity.tsv            panel b
source_data/Fig5/contributor_cells.tsv                           panels c, d, f recomputation (owned by Fig5)
source_data/SFig12/healthy_stool_grouping_leader_sensitivity.tsv panel e
source_data/SFig12/CapScan_grouping_leader_sensitivity.tsv       panel f
source_data/SFig12/capscan_contributor_region_cells_direct.tsv   panel f recomputation

``render`` ports the v30 build script (build_sfig12.py): one 180 x 172 mm canvas in the
repository fl_style at 6 pt, with the drawing code and row selection unchanged; paths are
arguments. ``coverage_from_matrices`` is the producer step that derives the panel a table
(build_revision.py lines 174-178), unchanged.

``four_rules`` is the function of the same name from the producer
check_s10c_eligibility.py, moved unchanged except that the comparison with the saved
agreement cells (an internal assertion of the producer; those 19 MB tables are not in
this repository) is left out.
"""

import hashlib
import io
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---- choices (named constants) ------------------------------------------------------
COHORTS = ["HSTOOL", "PREDICT", "CAPSCAN"]          # producer line 175; panel a order
RAREFACTION_DEPTH = 8          # legend, v30 paragraph 243; upstream derive_hstool_pc.py (math.comb(n, 8))
ELIGIBILITY_FRACTION = 0.8     # upstream derive_hstool_pc.py: leader in >= ceil(0.8 n) acquisitions of a design
MIN_PARTICIPANTS_PER_KO = 12   # check_s10c_eligibility.py four_rules: KO eligible if >= 12 participants
STRICT_RESOLVED_MIN = 0.2      # four_rules: species-resolved fraction >= 0.2 (both sides of a pair)
STRICT_GAP_MIN = 0.05          # four_rules: leader gap over runner-up >= 0.05 (both sides of a pair)
FIG_MM = (180.0, 172.0)        # v30 build: authored at the placed size
FONT_PT = 6.0                  # v30 build: 6 pt text
PNG_DPI = 300                  # v30 build: PNG resolution

# sha256 of every input, checked before use.
SHA256 = {
    "HSTOOL_KO_all.tsv": "e443ff28bb1ee485de65c55489f5f98d9b4c86037d3f392e8385edf7dc8248af",
    "PREDICT_KO_all.tsv": "663e1213de6c7fcdec2304bef076d420a4f5bf7320908c4d369136e16f454892",
    "CAPSCAN_KO_all.tsv": "1fcc2eb3fa9a634e9ef2f4f6f7630b8d835a3743c7afcf027312373af507eb23",
    "KO_allocation_coverage.tsv": "05a9a9d09028d3835706be3fad3d4f71903c59288cba70056fd92744b1648b01",
    "healthy_stool_leader_diversity.tsv": "93ed4ead5ea85f8d4fd8bf09750bb07dc8a66c0fe8565129f20c85988f496147",
    "contributor_cells.tsv": "6a704174763ce561aa2189fba7b137f3dc480edeebf7f52bdd16cabbb492dd84",
    "healthy_stool_grouping_leader_sensitivity.tsv": "fb4e1e3c34a233f627aad13fd78e33af77c74d9b785959ce7f934cbd8d589716",
    "CapScan_grouping_leader_sensitivity.tsv": "593d0251c4db0947172fd0437b342908546ea1910cac6cdfa18dd801fe5af9b2",
    "capscan_contributor_region_cells_direct.tsv": "e7cb87a78397858a6487ab4c81475f3dfdc7b465cf9fe8643963af4f6e6ee78c",
    "participant_leaders_resolved.tsv": "12c3fd7510af01a8e530aa56620b77477c4e78f8c50116d127c23e9e11962ef4",
    "participant_leaders_strict_both_methods.tsv": "adc373ee4aaf11076f9a53be8d5be60864c30534e14793effd721fae3bbaef38",
}

_NUMBER_WORDS = {8: "eight", 15: "fifteen"}


def check_sha256(path):
    """Raise if the sha256 of ``path`` differs from ``SHA256[path.name]``.

    Parameters
    ----------
    path : pathlib.Path
        Input file.

    Returns
    -------
    str
        The sha256.
    """
    s = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    if s != SHA256[Path(path).name]:
        raise ValueError(f"{Path(path).name}: sha256 {s} differs from expected {SHA256[Path(path).name]}")
    return s


def load_tables(sd, sd_fig5):
    """Read the five plotted tables exactly as the producer did (``pd.read_csv(sep='\\t')``).

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig12``.
    sd_fig5 : pathlib.Path
        ``source_data/Fig5`` (holds ``contributor_cells.tsv``).

    Returns
    -------
    dict
        ``coverage`` (a), ``leaders`` (b), ``cells`` (c, d), ``agree`` (e), ``capagree`` (f).
    """
    paths = dict(coverage=sd / "KO_allocation_coverage.tsv", leaders=sd / "healthy_stool_leader_diversity.tsv",
                 cells=sd_fig5 / "contributor_cells.tsv", agree=sd / "healthy_stool_grouping_leader_sensitivity.tsv",
                 capagree=sd / "CapScan_grouping_leader_sensitivity.tsv")
    out = {}
    for k, p in paths.items():
        check_sha256(p)
        out[k] = pd.read_csv(p, sep='\t')
    # producer line 180: only the current-analysis rows are plotted in f
    out["capagree"] = out["capagree"][out["capagree"].source.str.startswith('A33_')]
    return out


# ---- v30 figure (build_sfig12.py), drawing unchanged ------------------------------------

def render(t, out_stem):
    """Draw Supplementary Fig. 12 (v30) on one 180 x 172 mm canvas, as build_sfig12.py did.

    Ported from the v30 build script with paths replaced; drawing and row selection are
    unchanged. Panel e keeps only the direct-peptide-assignment rows (``method ==
    'simplified'``, subset ``all_eligible``) and panel f only ``A33_CAPSCAN__simplified``;
    a-d are unchanged in content.

    Parameters
    ----------
    t : dict
        Output of ``load_tables`` (``coverage``, ``leaders``, ``cells``, ``agree``,
        ``capagree``).
    out_stem : pathlib.Path
        Output path without extension; ``.pdf`` and ``.png`` (``PNG_DPI``) are written.

    Returns
    -------
    pandas.DataFrame
        The plotted summary, as the v30 build wrote plotted_summary.tsv.
    """
    from style import fl_style as fl

    coverage, leaders, cells, agree, capagree = t["coverage"], t["leaders"], t["cells"], t["agree"], t["capagree"]
    FS = FONT_PT
    fl.init()
    W, H = FIG_MM
    fig = plt.figure(figsize=fl.mm(W, H))

    def ax_mm(x, y, w, h):
        a = fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H]); fl.style(a); a.tick_params(labelsize=FS); return a

    def letter(x, y, s):
        fig.text(x / W, 1 - y / H, s, fontsize=7, fontweight='bold', color=fl.INK, va='top', ha='left')

    def labels(a, x=None, y=None):
        if x: a.set_xlabel(x, fontsize=FS)
        if y: a.set_ylabel(y, fontsize=FS)
    COL = [fl.BLUE, fl.ORANGE, fl.GREEN, fl.VERM, fl.SKY]
    plotted = []

    def dots(ax, groups, percent=False, panel=''):
        for i, (label, values) in enumerate(groups):
            v = np.asarray(values, dtype=float); v = v[np.isfinite(v)] * (100 if percent else 1)
            x = i + np.linspace(-.16, .16, len(v)); ax.scatter(x, v, s=4, alpha=.48, color=COL[i % 5], lw=0, rasterized=len(v) > 100)
            if len(v): ax.hlines(np.median(v), i - .25, i + .25, color='black', lw=1.5, zorder=4)
            plotted.append(dict(panel=panel, group=label, n=len(v), median=float(np.median(v)) if len(v) else np.nan))
        ax.set_xticks(range(len(groups)), [g[0] for g in groups])

    # a
    ax = ax_mm(16, 4, 66, 40); letter(1, 3, 'a')
    dots(ax, [(co, coverage[coverage.cohort.eq(co)].allocated_fraction) for co in COHORTS], True, 'a')
    ax.set_xticklabels(['Healthy stool', 'PREDICT', 'CapScan'])
    labels(ax, y='Accepted feature intensity\nallocated to KO (%)')
    # b
    ax = ax_mm(108, 4, 66, 40); letter(93, 3, 'b')
    dots(ax, [(k.replace('Replicates', ''), v.expected_distinct_leaders_in_8_resolved_acquisitions) for k, v in leaders.groupby('tier')], False, 'b')
    ax.tick_params(axis='x', rotation=30); [x.set_ha('right') for x in ax.get_xticklabels()]
    labels(ax, y='Expected distinct leaders\nin 8 acquisitions')
    # c, d
    for (x, lx, l), (col, lab) in zip([(16, 1, 'c'), (108, 93, 'd')], [('leader_share', 'Leader / species-resolved signal (%)'), ('species_resolved_fraction', 'Resolved / non-host functional signal (%)')]):
        ax = ax_mm(x, 64, 66, 40); letter(lx, 63, l)
        for participant, g in cells.groupby('participant'):
            v = np.sort(g[col].dropna() * 100); ax.plot(v, np.arange(1, len(v) + 1) / len(v), alpha=.45, lw=1, color=fl.BLUE)
        ax.set(xlim=(0, 100)); labels(ax, lab, 'Cumulative fraction of functions')
    plotted.append(dict(panel='c,d', group='participants', n=cells.participant.nunique(), median=np.nan))
    # e: direct peptide assignment only
    g = agree[agree.subset.eq('all_eligible') & agree.method.eq('simplified')]
    assert len(g) == 4
    ax = ax_mm(16, 120, 66, 36); letter(1, 119, 'e')
    ax.scatter(np.arange(len(g)), g.equal_KO_agreement * 100, s=14, color=fl.BLUE, lw=0, zorder=3)
    ax.set_xticks(np.arange(len(g)), [x.replace('Replicates', '') + '\n' + f.replace('_', ' ') for x, f in zip(g.tier, g['filter'])])
    ax.set_xlim(-.5, len(g) - .5); labels(ax, y='Equal-function leader\nagreement (%)')
    plotted += [dict(panel='e', group=f'{x}/{f}', n=int(k), value_pct=float(v * 100)) for x, f, k, v in zip(g.tier, g['filter'], g.n_kos, g.equal_KO_agreement)]
    # f: direct peptide assignment only
    c = capagree[capagree.source.eq('A33_CAPSCAN__simplified')]
    assert len(c) == 4
    ax = ax_mm(108, 120, 66, 36); letter(93, 119, 'f')
    ax.plot(np.arange(len(c)), c.donor_balanced_agreement * 100, 'o-', color=fl.BLUE, ms=3.5)
    ax.set_xticks(np.arange(len(c)), ['primary', 'strict pair\nfilter', 'strict donor\ncell eligibility', 'strict paired\ndonor eligibility'])
    ax.set_xlim(-.4, len(c) - .6); labels(ax, y='Donor-balanced leader\nagreement (%)')
    plotted += [dict(panel='f', group=r, n=int(k), value_pct=float(v * 100)) for r, k, v in zip(c.rule, c.n_kos, c.donor_balanced_agreement)]
    assert list(c.rule) == ['primary', 'strict_pair_filter', 'strict_donor_cell_eligibility', 'strict_paired_donor_eligibility']
    fl.save(fig, out_stem, dpi=PNG_DPI, formats=('pdf', 'png'))
    plt.close(fig)
    return pd.DataFrame(plotted)


def write_coverage(coverage, path):
    """Write the derived coverage table as the producer's ``export`` did (pandas, tab, no index).

    Parameters
    ----------
    coverage : pandas.DataFrame
        Output of ``coverage_from_matrices``.
    path : pathlib.Path

    Returns
    -------
    bool
        Whether the written file is byte-identical to ``source_data/SFig12/KO_allocation_coverage.tsv``.
    """
    coverage.to_csv(path, sep='\t', index=False)
    return hashlib.sha256(Path(path).read_bytes()).hexdigest() == SHA256["KO_allocation_coverage.tsv"]


# ---- statistics ------------------------------------------------------------------------
def stat_a_coverage(coverage):
    """Per-cohort summary of panel a (descriptive).

    Statistical unit: acquisition; n = 58 HSTOOL, 306 PREDICT, 346 CapScan (one of the
    347 CapScan acquisitions has no accepted features and is not in the allocation).

    Parameters
    ----------
    coverage : pandas.DataFrame

    Returns
    -------
    dict
        ``{cohort: {n, median, min, max}}`` of ``allocated_fraction``.
    """
    return {co: dict(n=int(g.shape[0]), median=float(np.median(g.allocated_fraction)),
                     min=float(g.allocated_fraction.min()), max=float(g.allocated_fraction.max()))
            for co, g in coverage.groupby('cohort')}


def coverage_from_matrices(paths):
    """Producer step that derives ``KO_allocation_coverage.tsv`` from the KO matrices.

    Moved unchanged from build_revision.py lines 174-178. Each matrix is KO x acquisition
    intensity with a ``FunctionallyUnassigned`` row (the remainder not allocated to a KO);
    the denominator of an acquisition is its column sum. Statistical unit: acquisition.

    Parameters
    ----------
    paths : dict
        ``{cohort: path to source_data/SFig12/<cohort>_KO_all.tsv}``.

    Returns
    -------
    pandas.DataFrame
        Columns cohort, sample, total_intensity, unassigned_intensity, allocated_fraction.
    """
    coverage = []
    for co in COHORTS:
        check_sha256(paths[co])
        t = pd.read_csv(paths[co], sep='\t', index_col=0)
        total = t.sum(axis=0); missing = t.loc['FunctionallyUnassigned']
        assert np.isfinite(t.values).all() and (t.values >= 0).all() and (total > 0).all(), co + ' allocation matrix'
        coverage.extend(dict(cohort=co, sample=s, total_intensity=float(total[s]), unassigned_intensity=float(missing[s]), allocated_fraction=float(1 - missing[s] / total[s])) for s in t.columns)
    return pd.DataFrame(coverage)


def stat_b_depth(leaders):
    """Rarefaction depth of panel b, read from the source-data column name.

    The expectation itself was computed upstream (A33 derive_hstool_pc.py) analytically,
    E = sum over leading species of 1 - C(n - v, 8) / C(n, 8), with n resolved
    acquisitions of the function in that design and v those led by the species; no random
    draws. Statistical unit: function within acquisition design.

    Parameters
    ----------
    leaders : pandas.DataFrame

    Returns
    -------
    dict
        ``depth`` (int), ``word`` (as printed), ``n_rows`` and ``n_per_design``.
    """
    cols = [c for c in leaders.columns if c.startswith('expected_distinct_leaders_in_')]
    assert len(cols) == 1
    depth = int(cols[0].split('_in_')[1].split('_')[0])
    return dict(depth=depth, word=_NUMBER_WORDS.get(depth, str(depth)), n_rows=int(len(leaders)),
                n_per_design=leaders.groupby('tier').size().to_dict())


def stat_cd_counts(cells):
    """Counts printed in the legend for panels c,d.

    Statistical unit: function x participant cell; n = 21,168 cells.

    Parameters
    ----------
    cells : pandas.DataFrame
        ``source_data/Fig5/contributor_cells.tsv``.

    Returns
    -------
    dict
        ``cells``, ``functions``, ``participants``, ``participants_word``.
    """
    p = int(cells.participant.nunique())
    return dict(cells=int(len(cells)), functions=int(cells.KO.nunique()), participants=p,
                participants_word=_NUMBER_WORDS.get(p, str(p)))


def stat_cd_median_resolved(cells):
    """Median species-resolved share of non-host functional intensity over eligible cells.

    Statistical unit: function x participant cell, n = 21,168 cells (1,467 functions,
    15 participants); plain median over all cells, no weighting, no interval.

    Parameters
    ----------
    cells : pandas.DataFrame

    Returns
    -------
    float
        Median of ``species_resolved_fraction`` (a fraction, not %).
    """
    return float(np.median(cells.species_resolved_fraction))


def read_cells(p):
    """Producer ``read`` of check_s10c_eligibility.py (participant as string, no NA parsing)."""
    return pd.read_csv(p, sep='\t', dtype={'participant': str}, keep_default_na=False)


def mean(pairs):
    """Producer ``mean``: agreement averaged within participant-region, then region, then participant."""
    return float(pairs.groupby(['participant', 'region']).same_leader.mean().groupby('participant').mean().mean())


def four_rules(cells_path, region_cells_path):
    """Leader agreement of intestinal regions with own stool under four eligibility rules.

    Moved from check_s10c_eligibility.py ``four_rules`` (see module docstring).
    Statistical unit: function x participant x intestinal region pair, summarised as the
    donor-balanced mean (pairs -> participant-region -> participant -> mean over 15
    participants) and the pooled fraction of pairs. Descriptive; no test, no interval.

    Parameters
    ----------
    cells_path : pathlib.Path
        Function x participant cells (eligibility: KO in >= 12 participants).
    region_cells_path : pathlib.Path
        Function x participant x region cells (four intestinal regions and stool).

    Returns
    -------
    tuple
        ``(table, n_eligible_kos, n_cells)``; ``table`` has one row per rule.
    """
    c = read_cells(cells_path)
    s = read_cells(region_cells_path)
    assert not c.duplicated(['KO', 'participant']).any()
    assert not s.duplicated(['KO', 'participant', 'region']).any()
    eligible = set(c.groupby('KO').participant.nunique().loc[lambda x: x >= MIN_PARTICIPANTS_PER_KO].index)
    stool = s[s.region.eq('stool')].drop(columns='region')
    p = s[s.region.ne('stool') & s.KO.isin(eligible)].merge(
        stool, on=['KO', 'participant'], suffixes=('_capsule', '_stool'), validate='m:1')
    p['same_leader'] = p.leading_species_capsule.eq(p.leading_species_stool)
    p['well_resolved'] = (p[['species_resolved_fraction_capsule', 'species_resolved_fraction_stool']].min(axis=1).ge(STRICT_RESOLVED_MIN)
                          & p[['gap_capsule', 'gap_stool']].min(axis=1).ge(STRICT_GAP_MIN))
    strict = p[p.well_resolved]
    strict_c = c[c.species_resolved_fraction.ge(STRICT_RESOLVED_MIN) & c.gap.ge(STRICT_GAP_MIN)]
    keepA = set(strict_c.groupby('KO').participant.nunique().loc[lambda x: x >= MIN_PARTICIPANTS_PER_KO].index)
    keepB = set(strict.groupby('KO').participant.nunique().loc[lambda x: x >= MIN_PARTICIPANTS_PER_KO].index)
    cohorts = {'primary': p, 'strict_pair_filter': strict,
               'strict_donor_cell_eligibility': strict[strict.KO.isin(keepA)],
               'strict_paired_donor_eligibility': strict[strict.KO.isin(keepB)]}
    rows = []
    for rule, f in cohorts.items():
        d = f.groupby(['participant', 'region']).size()
        rows.append(dict(rule=rule, donor_balanced_agreement=mean(f), pooled_cell_agreement=float(f.same_leader.mean()),
                         n_pairs=len(f), n_kos=f.KO.nunique(), n_donors=f.participant.nunique(), n_donor_regions=len(d)))
    return pd.DataFrame(rows), len(eligible), len(c)


def stat_f_recompute(sd, sd_fig5):
    """Recompute the plotted panel f rows (direct peptide assignment) from the cell tables.

    Parameters
    ----------
    sd, sd_fig5 : pathlib.Path
        ``source_data/SFig12`` and ``source_data/Fig5``.

    Returns
    -------
    tuple
        ``(table, max_abs_diff, counts_equal)``: the recomputed four-row table with a
        ``source`` column as in the plotted table, the largest absolute difference of the
        two agreement columns to the ``A33_CAPSCAN__simplified`` rows of
        ``CapScan_grouping_leader_sensitivity.tsv``, and whether all count columns are equal.
        Statistical unit: participant (donor-balanced mean), n = 15.
    """
    check_sha256(sd / "capscan_contributor_region_cells_direct.tsv")
    check_sha256(sd_fig5 / "contributor_cells.tsv")
    new, nk, nc = four_rules(sd_fig5 / "contributor_cells.tsv", sd / "capscan_contributor_region_cells_direct.tsv")
    new.insert(0, 'source', 'A33_CAPSCAN__simplified')
    ref = pd.read_csv(sd / "CapScan_grouping_leader_sensitivity.tsv", sep='\t')
    ref = ref[ref.source.eq('A33_CAPSCAN__simplified')].reset_index(drop=True)
    assert list(new.source) == list(ref.source) and list(new.rule) == list(ref.rule)
    diff = float(max(np.abs(new[c] - ref[c]).max() for c in ['donor_balanced_agreement', 'pooled_cell_agreement']))
    counts = all(list(new[c]) == list(ref[c]) for c in ['n_pairs', 'n_kos', 'n_donors', 'n_donor_regions'])
    return new, diff, counts


def stat_participant_leaders(sd_fig5):
    """Participant-averaged leader agreement printed on v28d line 157 (Fig. 5f tables).

    Agreement per participant is already averaged over functions and intestinal regions
    (upstream fig5ef_run.py); this function averages over the 15 participants.
    Statistical unit: participant, n = 15. Descriptive; no interval.

    Parameters
    ----------
    sd_fig5 : pathlib.Path

    Returns
    -------
    dict
        ``{threshold: {'original': mean, 'simplified': mean, 'participants': n}}``.
    """
    out = {}
    for thr, name in [('resolved', 'participant_leaders_resolved.tsv'), ('strict_both_methods', 'participant_leaders_strict_both_methods.tsv')]:
        check_sha256(sd_fig5 / name)
        t = pd.read_csv(sd_fig5 / name, sep='\t')
        out[thr] = dict(original=float(t.agreement_original.mean()), simplified=float(t.agreement_simplified.mean()),
                        participants=int(t.participant.nunique()))
    return out
