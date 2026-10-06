"""Figure 2: recovery of the defined CAMPI SIHUMIx community (v30, 9-50-residue panels).

Panels a-h are drawn on one 180 x 250 mm canvas by the original producer
(``LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/figure2_redesign/build_campi_figure.py``,
sha256 dc20ed0b...). This module moves that producer's main-figure code here with paths
made repository-relative. The plotting logic, draw order, sizes and colours are unchanged;
the producer's companion figure, SVG export and side-written copies of its inputs are not
ported because they are not part of the manuscript figure.

v30 producer: the len950 rebuild of the module above (``panels_len950.py``, see PROVENANCE.md). Panels b, d, e, f and g read the
nested 9-50-residue population (``PLOT_MIN_LENGTH``): the min_length == 9 rows of the A14
tables, i.e. the accepted 5-50 lists filtered to length >= 9, and an UpSet recount from the 15
membership tables (``upset_at``). Panel c (protein groups) and panel h (lengths by class) stay
5-50. The panel d and e labels quoting a count are computed from the plotted table; the panel b
y-axis and panel e x-axis read "accepted" (v29 relabel). The composite PDF edits that follow
the canvas are in ``assembly.py``.

The producer labelled its review layout a, b, c, d-f, g, h, i. The manuscript letters
(v28d) were assigned afterwards by PDF edits (see ASSEMBLY.md). Mapping, producer -> v28d:
a -> a, b -> b, c -> c, d-f left (loss stages) -> d, d-f right (pooled composition and
recovery) -> e, g (member heatmap) -> f, h (UpSet) -> g, i (lengths) -> h.

Source data (all under ``source_data/``):

- ``Fig2/paired_reference_recovery.tsv``, ``Fig2/recovery_summary.tsv``,
  ``Fig2/member_summary.tsv``, ``Fig2/stage_summary.tsv`` (A14 CAMPI reference gap)
- ``Fig2/arm_<cohort>__reaggregated_per_acquisition.tsv`` and
  ``Fig2/arm_<cohort>__direct_acquisition_summary.tsv`` (simplified grouping, 5 arms)
- ``Fig2/metrics_per_acquisition.tsv`` (database entries per acquisition and arm)
- ``Fig2/length_by_member_class_a14.tsv``, ``Fig2/upset_intersections_a14.tsv`` (5-50),
  ``Fig2/upset_intersections_FL_9_50.tsv`` (9-50; checked against the recount, panel g)
- ``SFig03/CAMPI_membership_<sample>.tsv`` (15 A14 peptide membership tables; shared
  with Supplementary Figure 3 and stored under that unit)
- ``Fig1/acquisition_sizes_726.tsv`` (reservoir size; shared, stored under Figure 1)
- ``Fig2/reservoir_protein_presence_by_proteome.tsv`` (member proteome containment)
- ``Fig2/campi_gain_availability.tsv``, ``Fig2/campi_new_quantified_peptides.tsv``,
  ``Fig2/campi_depth_standardization.tsv``, ``Fig2/campi_FL_5_frequencies_<sample>.tsv``
  (CAMPI variability and add-back mechanism)
- ``Fig2/campi_discrepant_prediction_pooled_summary.tsv`` (same-spectrum audit, A03)
- ``Fig2/alphanovo_accuracy_pooled.tsv``, ``Fig2/alphanovo_accuracy_per_acquisition.tsv``
"""

import hashlib
import math
from collections import Counter
from pathlib import Path

import fitz
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.special import gammaln

from style import fl_style as fl

# ---- choices (each value is used below; README Choices table lists them) ------------
MIN_LENGTH = 5                 # legend, v30 paragraph 42 ("5-50 residues in c and h")
MIN_LENGTH_SENSITIVITY = 9     # Results, v30 paragraph 40 ("only six at 9-50 residues")
MAX_LENGTH = 50                # legend, v28d line 87
JITTER_SEED = 20260922         # producer constant; horizontal point jitter in b and c only
JITTER_HALF_WIDTH = 0.1        # producer constant; display only
GROUPING_METHOD = "simplified"  # reaggregated rows used for c (Methods, protein grouping)
PAGE_MM = (180, 250)           # producer canvas
BOTTLE_RECT_MM = (3, 7, 50, 40.5)  # producer placement of the community artwork in a
PRODUCER_LINE_PT = 0.6        # producer-era house line weight (axes, ticks, lines); see _producer_rc
RENDER_DPI = 300               # producer PNG refresh (FINAL_REVIEW.json png_refresh)
N_ACQUISITIONS = 15            # Results, v28d line 75
N_MEMBERS = 8                  # Results, v28d line 75 ("eight-member community")
PLOT_MIN_LENGTH = 9            # legend, v30 paragraph 42 ("9-50 residues in b and d-g"); 5 reproduces v28f
N_LEAST_REPRESENTED = 3        # Results, v30 paragraph 40: members with the fewest reference-search
                               # member-specific peptides (9-50 residues)
N_MOST_REPRESENTED = 5         # Results, v30 paragraph 40 ("four of the five most represented members")
MG_MISSED_MEMBER = 'Anaerostipes caccae'  # Results, v30 paragraph 40 ("which MG largely missed")

ARMS = ['FL', 'REF', 'MG', 'FL_R', 'FL_R_M']
LABEL = dict(zip(ARMS, ['FL', 'REF', 'MG', 'FL + REF', 'FL + REF + MG']))
COHORT = dict(zip(ARMS, ['CAMPI_FL_UNIPROT', 'CAMPI_REF', 'CAMPI_MG', 'CAMPI_FL_REF', 'CAMPI_FL_REF_MG']))
METRICS_KEY = {'FL_R': 'FLR', 'FL_R_M': 'FLRM'}
MEMBERS = ['Bacteroides thetaiotaomicron', 'Blautia producta', 'Erysipelatoclostridium ramosum',
           'Escherichia coli', 'Anaerostipes caccae', 'Bifidobacterium longum',
           'Clostridium butyricum', 'Lactobacillus plantarum']
SAMPLES = ['S01', 'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10',
           'S11_Fraction1', 'S11_Fraction2', 'S11_Fraction3', 'S11_Fraction4', 'S12']
STAGES = ['absent_from_acquisition_candidates', 'candidate_digest_compatible_removed_by_selection',
          'candidate_sequence_present_digest_incompatible_absent_final',
          'final_sequence_present_digest_incompatible', 'final_digest_compatible_not_accepted']
#: Loss stages at which the reference sequence was not in the selected (final) FASTA.
STAGES_SEQUENCE_ABSENT = STAGES[:3]
#: Same-spectrum audit categories in which the retained prediction disagrees with the reference.
DISCREPANT_CATEGORIES = ['adjacent_swap', 'length_difference', 'other_disagreement',
                         'rearrangement', 'substitution']
SHORT = lambda m: m.split()[0][0] + '. ' + ' '.join(m.split()[1:])  # noqa: E731  (producer)

#: Panel crop boxes on the producer canvas, millimetres (x0, y0, x1, y1), v28d letters.
PANEL_BOXES_MM = {
    'a': (0, 0, 180, 51),
    'b': (0, 51, 90, 101),
    'c': (90, 51, 180, 101),
    'd': (0, 101, 92, 150),
    'e': (92, 101, 180, 150),
    'f': (0, 150, 180, 196),
    'g': (0, 196, 134, 250),
    'h': (134, 196, 180, 250),
}

#: sha256 of the shared size table stored under Figure 1 (checked before use).
SIZES_SHA256 = 'fd9ade3bb4ebd2ee6c129f69819b9697cf20abd7c10c77a12644bef25396e386'
MEMBERSHIP_SHA256 = {
    'S01': 'fe3f27ac42f772792572e8400c30e2665eb7ca8b3e9d075f0b5f6e0735035fb2',
    'S02': '4a64e2930e4d27bd74b913b6863e742a8663195a5364f01c86422ce39b55c78c',
    'S03': 'c065d0f2714bfd2e7999666185c1f299ca49ad65a7acb7beaa576fec9f865205',
    'S04': '2831b3686f6ad449fba4e94657827f719993b2fb91e7eb0e76e5d5770fe1ef17',
    'S05': 'cce302fcbfbd71bc4bb9677863bfc5aa68ad5b5d730cfeab6e822e8ee4f1beb5',
    'S06': '4b81dbb3d5c4ea81e250e274971fd3501f41811d366b6fde8e96845ef18b32ee',
    'S07': 'ec883caeb75013b3ca778ee8d4255b633fde058d99b5cf11e5bc2e38355979a7',
    'S08': '5f1f92c20b3746aec6927465b56f17755545abf88ac8de64afce1a4be9ddd26f',
    'S09': 'dc2aef1c5d7ac110fc202d6b65d8287f0aaa6d28f0542c289d3491ca356feea5',
    'S10': '7520e41bca7aecc7468f7606359c1f80816c51f95396a67e44d377c95031aee0',
    'S11_Fraction1': '55cbde9f72cff8bb5aa72836abcee99047ba581792203d32c5aeb79a7052ee65',
    'S11_Fraction2': '54be84a4f676d0d4d81a902c0c739c20c91338eeeec997405fb13a1c22af7d1b',
    'S11_Fraction3': '0620cde139e1acf8eb4cd47c997619dbdbccdeba6917a8d507302bb073df3b32',
    'S11_Fraction4': 'e93e188f04caf0aa359bb5158b2ced259286c12bddd59caa872b32fa9e807d05',
    'S12': 'd68c1485743dfc9283263e3c32f8392bd4712a556181c37acede4b00cd2f16a6',
}
NUMBER_WORDS = {3: 'three', 4: 'four', 5: 'five', 6: 'six', 7: 'seven', 8: 'eight', 12: 'twelve', 13: 'thirteen', 15: 'fifteen'}


def sha256(path):
    """Return the hex sha256 of a file."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def as_word(n):
    """Spell a small integer as the manuscript prints it (e.g. 12 -> 'twelve')."""
    if float(n) != int(n) or int(n) not in NUMBER_WORDS:
        raise ValueError(f'no word for {n!r}')
    return NUMBER_WORDS[int(n)]


def read(p):
    """Read a tab-separated table as the producer did."""
    return pd.read_csv(p, sep='\t', comment='#')


# ---------------------------------------------------------------- inputs ----------
def load_inputs(sd_root):
    """Load and validate every table the canvas needs (producer lines 53-83).

    Parameters
    ----------
    sd_root : pathlib.Path
        ``source_data/`` directory of the repository.

    Returns
    -------
    dict
        Tables and derived matrices: ``paired_all``, ``paired``, ``rs_all``, ``rs``, ``ms_all``,
        ``ms``, ``ss_all``, ``ss``, ``gp`` (group rows), ``pp`` (peptides, 15 x 5),
        ``gg`` (quantified groups, 15 x 5), ``rr`` (same-acquisition reference recovery %,
        15 x 5), ``support`` (median % groups with >=2 peptides), ``db`` (median entries),
        ``mt`` (pooled member-specific peptides, 8 x 5), ``lt``, ``up``. Statistical unit:
        acquisition (n = 15) for pp, gg, rr; pooled peptides for rs, ms, mt.
    """
    sd = Path(sd_root) / 'Fig2'
    paired_all = read(sd / 'paired_reference_recovery.tsv')
    paired = paired_all[paired_all.min_length.eq(PLOT_MIN_LENGTH)].copy()
    rs_all = read(sd / 'recovery_summary.tsv')
    rs = rs_all[rs_all.min_length.eq(PLOT_MIN_LENGTH)].set_index('arm').loc[ARMS]
    ms_all = read(sd / 'member_summary.tsv')
    ms = ms_all[ms_all.min_length.eq(PLOT_MIN_LENGTH)]
    ss_all = read(sd / 'stage_summary.tsv')
    ss = ss_all[ss_all.min_length.eq(PLOT_MIN_LENGTH)].set_index('stage').loc[STAGES]
    gp = []
    for arm in ARMS:
        g = read(sd / f'arm_{COHORT[arm]}__reaggregated_per_acquisition.tsv')
        g = g[g.method.eq(GROUPING_METHOD)].copy().set_index('sample').sort_index()
        direct = read(sd / f'arm_{COHORT[arm]}__direct_acquisition_summary.tsv').set_index('sample').sort_index()
        assert len(g) == N_ACQUISITIONS and g.index.equals(direct.index)
        assert g.quantified_groups.astype(int).equals(direct.quantified_study_groups.astype(int))
        g['arm'] = arm
        gp.append(g.reset_index())
    gp = pd.concat(gp, ignore_index=True)
    pp = paired.pivot(index='sample', columns='arm', values='peptides').loc[:, ARMS]
    gg = gp.pivot(index='sample', columns='arm', values='quantified_groups').loc[:, ARMS]
    rr = paired.pivot(index='sample', columns='arm', values='reference_recovery').loc[:, ARMS] * 100
    assert pp.shape == gg.shape == rr.shape == (N_ACQUISITIONS, 5)
    assert pp.index.equals(gg.index) and not pp.isna().any().any()
    support = gp.groupby('arm').groups_with_at_least_two_assigned_peptides_pct.median().loc[ARMS]
    db = read(sd / 'metrics_per_acquisition.tsv').groupby('strategy').db_records.median()
    db = {a: int(db[METRICS_KEY.get(a, a)]) for a in ARMS}
    mt = ms.pivot(index='member', columns='arm', values='pooled_specific_peptides').loc[MEMBERS, ARMS]
    lt = read(sd / 'length_by_member_class_a14.tsv')
    up_table = read(sd / 'upset_intersections_a14.tsv')
    up = upset_at(sd_root, PLOT_MIN_LENGTH)
    if PLOT_MIN_LENGTH == MIN_LENGTH:   # control: the recount must equal the published table
        a = up.set_index('members').n_peptides.sort_index(); b = up_table.set_index('members').n_peptides.sort_index()
        assert a.equals(b)
        assert set(up[up.drawn].members) == set(up_table[up_table.drawn.eq(True)].members)
        up = up_table
    elif PLOT_MIN_LENGTH == MIN_LENGTH_SENSITIVITY:   # v30: the recount must equal the shipped 9-50 table
        t9 = read(sd / 'upset_intersections_FL_9_50.tsv')
        assert up[['members', 'degree', 'n_peptides', 'drawn']].reset_index(drop=True).equals(
            t9[['members', 'degree', 'n_peptides', 'drawn']].reset_index(drop=True))
    return dict(paired_all=paired_all, paired=paired, rs_all=rs_all, rs=rs, ms_all=ms_all, ms=ms,
                ss_all=ss_all, ss=ss, gp=gp, pp=pp, gg=gg, rr=rr, support=support, db=db, mt=mt,
                lt=lt, up=up)


def upset_at(sd_root, min_len):
    """len950: UpSet intersections of pooled FL peptides of length >= ``min_len`` (<= 50).

    Recounted from the 15 accepted-peptide membership tables (``SFig03/CAMPI_membership_<S>.tsv``,
    identical to A14_4_3_campi_gap/results/<S>/peptide_membership.tsv). Drawn: the 8
    single-member sets and the 6 largest shared sets, as in the published panel g.
    """
    pool, membership = set(), {}
    for s in SAMPLES:
        p = Path(sd_root) / 'SFig03' / f'CAMPI_membership_{s}.tsv'
        assert sha256(p) == MEMBERSHIP_SHA256[s], p.name
        q = read(p)
        for row in q.itertuples(index=False):
            names = tuple(sorted(str(row.members).split(';'))) if row.n_members else ()
            membership[row.peptide] = names
            if 'FL' in row.arms.split(';') and min_len <= len(row.peptide) <= MAX_LENGTH:
                pool.add(row.peptide)
    cnt = Counter(membership[p] for p in pool)
    rows = [dict(members='|'.join(k) if k else '(none of the eight)', degree=len(k), n_peptides=v) for k, v in cnt.items()]
    up = pd.DataFrame(rows).sort_values(['n_peptides', 'members'], ascending=[False, True]).reset_index(drop=True)
    shared = up[up.degree > 1]
    assert shared.n_peptides.iloc[5] > shared.n_peptides.iloc[6], 'tie at the sixth shared set'
    up['drawn'] = up.degree.eq(1) | up.index.isin(shared.index[:6])
    assert up.drawn.sum() == 14 and up.degree.eq(1).sum() == N_MEMBERS
    up.insert(0, 'rank', range(1, len(up) + 1))
    return up


def recount_membership(sd_root, d):
    """Recount identities and member assignments from the 15 membership tables.

    This is the producer's independent recount (producer lines 85-146), unchanged. It
    rebuilds the pooled peptide sets per arm from ``SFig03/CAMPI_membership_<S>.tsv`` and
    asserts that every per-acquisition and pooled count in the summary tables matches.

    Parameters
    ----------
    sd_root : pathlib.Path
        ``source_data/`` directory.
    d : dict
        Output of :func:`load_inputs`.

    Returns
    -------
    dict
        ``pool`` (arm -> set of peptides), ``membership`` (peptide -> tuple of members),
        ``acq_members`` (arm -> Counter of acquisitions with an L. plantarum-specific
        peptide), ``length_stats`` (class -> n and median length of FL peptides),
        ``checks`` (counts of verified rows). Statistical unit: peptide (pooled).
    """
    pp, rr, rs, ms_all, ss, paired, lt, up = (d[k] for k in ('pp', 'rr', 'rs', 'ms_all', 'ss', 'paired', 'lt', 'up'))
    pool = {a: set() for a in ARMS}
    pool_all = {a: set() for a in ARMS}   # len950: unfiltered 5-50 pools (member and length checks, panel h)
    membership = {}
    acq_members = {a: Counter() for a in ARMS}
    files = [Path(sd_root) / 'SFig03' / f'CAMPI_membership_{s}.tsv' for s in SAMPLES]
    assert len(files) == N_ACQUISITIONS
    for p, sample in zip(files, SAMPLES):
        assert sha256(p) == MEMBERSHIP_SHA256[sample], p.name
        q = read(p)
        local = {a: set() for a in ARMS}
        for row in q.itertuples(index=False):
            pep = row.peptide
            names = tuple(str(row.members).split(';')) if row.n_members else ()
            assert len(names) == row.n_members
            if pep in membership:
                assert membership[pep] == names
            membership[pep] = names
            for arm in row.arms.split(';'):
                pool_all[arm].add(pep)
                if len(pep) >= PLOT_MIN_LENGTH:   # len950: plotted population
                    pool[arm].add(pep)
                    local[arm].add(pep)
        for a in ARMS:
            assert len(local[a]) == pp.loc[sample, a]
            if any(membership[x] == ('Lactobacillus plantarum',) for x in local[a]):
                acq_members[a]['Lactobacillus plantarum'] += 1
            rec = len(local[a] & local['REF']) / len(local['REF'])
            assert abs(rec * 100 - rr.loc[sample, a]) < 1e-8
    for a in ARMS:
        assert len(pool[a]) == rs.loc[a, 'pooled_peptides']
        assert len(pool[a] & pool['REF']) == rs.loc[a, 'pooled_overlap']
        for min_len in [MIN_LENGTH, MIN_LENGTH_SENSITIVITY]:
            counts = Counter(membership[p][0] for p in pool_all[a] if len(p) >= min_len and len(membership[p]) == 1)
            for m in MEMBERS:
                saved = ms_all[(ms_all.arm == a) & (ms_all.min_length == min_len) & (ms_all.member == m)]
                assert counts[m] == int(saved.pooled_specific_peptides.iloc[0])
        classes = Counter(0 if not membership[p] else 1 if len(membership[p]) == 1 else 2 for p in pool[a])
        assert classes[0] == rs.loc[a, 'absent_eight_member_reference']
        assert classes[1] == rs.loc[a, 'known_member_specific']
        assert classes[2] == rs.loc[a, 'known_member_shared']
    assert ss.peptide_acquisition_observations.sum() == paired[paired.arm.eq('FL')].lost.sum()
    length_stats = {}
    for cls, field in [('One member', 1), ('Shared by members', 2), ('No member', 0)]:
        lengths = [len(p) for p in pool_all['FL'] if (0 if not membership[p] else 1 if len(membership[p]) == 1 else 2) == field]
        length_stats[cls] = {'n': len(lengths), 'median': float(np.median(lengths)),
                             'min': min(lengths), 'max': max(lengths)}
        expected = Counter(lengths)
        for row in lt[lt['class'].eq(cls)].itertuples(index=False):
            assert row.n_peptides == expected[row.length]
            assert abs(row.pct_within_class - 100 * row.n_peptides / len(lengths)) < 1e-8
    intersections = Counter('|'.join(membership[p]) for p in pool['FL'])   # plotted population (panel g)
    for row in up.itertuples(index=False):
        key = '|'.join(sorted(str(row.members).split('|'))) if row.degree else ''
        expected = sum(v for k, v in intersections.items() if '|'.join(sorted(k.split('|'))) == key)
        assert row.n_peptides == expected
    checks = dict(raw_membership_files=N_ACQUISITIONS, pooled_member_counts_verified=80,
                  paired_peptide_rows_verified=75, direct_group_rows_verified=75,
                  length_and_upset_recount_pass=True)
    return dict(pool=pool, membership=membership, acq_members=acq_members,
                length_stats=length_stats, checks=checks)


# ---------------------------------------------------------------- canvas ----------
def _producer_rc():
    """Apply the producer's rcParams on top of ``fl.init()`` (producer lines 24-28).

    The repository style module is a later version of the house style that raised axis,
    tick and line weights to 1 pt (a floor added 2026-08-11 in the code repository) and set
    ``patch.linewidth``/``grid.linewidth``. The producer ran with the earlier module, where
    these were 0.6 pt and the matplotlib defaults. They are restored here so that the panels
    are drawn exactly as the producer drew them; nothing is restyled.
    """
    plt.rcParams.update({'lines.linewidth': PRODUCER_LINE_PT, 'axes.linewidth': PRODUCER_LINE_PT,
                         'xtick.major.width': PRODUCER_LINE_PT, 'ytick.major.width': PRODUCER_LINE_PT,
                         'xtick.minor.width': PRODUCER_LINE_PT, 'ytick.minor.width': PRODUCER_LINE_PT,
                         'patch.linewidth': 1.0, 'grid.linewidth': 0.8})
    plt.rcParams.update({'savefig.bbox': None, 'font.size': 7, 'axes.labelsize': 7,
                         'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
                         'text.color': fl.INK, 'axes.labelcolor': fl.INK,
                         'xtick.color': fl.INK, 'ytick.color': fl.INK})


def producer_style(ax):
    """``fl_style.style`` as the producer imported it (tick width 0.6 pt)."""
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.tick_params(length=2.6, width=PRODUCER_LINE_PT, labelsize=7)
    ax.xaxis.label.set_size(7)
    ax.yaxis.label.set_size(7)
    return ax


def canvas(h):
    """Producer helper: a 180 mm wide figure with mm-placed axes and text."""
    fig = plt.figure(figsize=fl.mm(180, h))

    def ax(x, y, w, height, letter=None):
        a = fig.add_axes([x/180, 1-(y+height)/h, w/180, height/h])
        producer_style(a)
        if letter:
            fl.panel(a, letter)
        return a

    def text(x, y, s, **kw):
        return fig.text(x/180, 1-y/h, s, **kw)
    return fig, ax, text


def paired_plot(ax, table, factor, ylabel, ylim, medians=True):
    """Producer helper for b and c: paired points per acquisition with median bars.

    Statistical unit: acquisition (n = 15 per database). The horizontal jitter is drawn from
    ``numpy.random.default_rng(JITTER_SEED)`` and is for display only.
    """
    col = dict(zip(ARMS, [fl.FL, fl.GREY, fl.ORANGE, fl.GREEN, fl.PURPLE]))
    jitter = np.random.default_rng(JITTER_SEED).uniform(-JITTER_HALF_WIDTH, JITTER_HALF_WIDTH, len(table))
    for j, (_, row) in zip(jitter, table.iterrows()):
        ax.plot(np.arange(5) + j, row.to_numpy()*factor, color='#d9d9d9', lw=.5, zorder=0)
    for i, a in enumerate(ARMS):
        values = table[a].to_numpy()*factor
        ax.scatter(i+jitter, values, s=10, c=col[a], alpha=.85, lw=0, zorder=2)
        ax.plot([i-.23, i+.23], [np.median(values)]*2, color=fl.INK, lw=1.15, zorder=3)
        if medians:
            ax.text(i, 1.045, f'{int(table[a].median()):,}', transform=ax.get_xaxis_transform(),
                    ha='center', va='bottom', fontsize=6.5, color=col[a])
    ax.set_xticks(range(5), [LABEL[a] for a in ARMS])
    ax.set_xlim(-.5, 4.5)
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis='x', length=0, labelsize=6)
    ax.tick_params(axis='y', labelsize=6.5)


def panel_a(fig, ax, tx, d):
    """Panel a: database design table (producer lines 194-213; community artwork added at save).

    Statistical unit: acquisition (median database entries over n = 15 per database).
    """
    col = dict(zip(ARMS, [fl.FL, fl.GREY, fl.ORANGE, fl.GREEN, fl.PURPLE]))
    db = d['db']
    tx(3, 5, 'a', fontsize=8)
    t = ax(54, 4, 123, 35)
    t.axis('off'); t.set_xlim(0, 123); t.set_ylim(35, 0)
    headers = [(0, 'Database'), (48, 'FL-selected\nproteins'), (71, 'Eight-member\nREF proteomes'),
               (94, 'Metagenome\nproteins'), (121, 'Entries\n(median)')]
    for x, s in headers:
        t.text(x, 1, s, ha='left' if x == 0 else 'right' if x == 121 else 'center', va='top', fontsize=6.3)
    flags = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 1, 1]]
    for i, a in enumerate(ARMS):
        y = 11 + i * 4.8 + (1.5 if i >= 3 else 0)
        t.text(0, y, LABEL[a], va='center', color=col[a], fontsize=7)
        for x, flag in zip([48, 71, 94], flags[i]):
            if flag:
                t.scatter(x, y, s=17, color=fl.INK)
            else:
                t.plot([x-.9, x+.9], [y, y], color='#bcbcbc', lw=.7)
        t.text(121, y, f'{db[a]:,}', va='center', ha='right', fontsize=6.8)
    t.plot([0, 123], [23, 23], color=fl.GREY, lw=.6)
    tx(54, 42, 'FL: selection from 375 M gut-catalogue and UniProtKB proteins;', fontsize=6.4)
    tx(54, 45.5, 'no community composition, CAMPI REF or specimen MG supplied.', fontsize=6.4)
    tx(3, 50, 'REF and MG were added to the selected FL database after parsimony.', fontsize=7)


def panel_b(fig, ax, tx, d):
    """Panel b: accepted peptides per acquisition (producer lines 215-217).

    Statistical unit: acquisition, n = 15 per database; bars are medians.
    """
    b = ax(18, 59, 64, 28)
    tx(3, 54, 'b', fontsize=8)
    paired_plot(b, d['pp'], .001, 'Accepted peptides\nper acquisition (thousands)', (0, 40))


def panel_c(fig, ax, tx, d):
    """Panel c: quantified protein groups per acquisition and % groups with >=2 peptides.

    Producer lines 218-224. Statistical unit: acquisition, n = 15 per database; bars and the
    percentage row are medians over acquisitions.
    """
    c = ax(111, 59, 64, 28)
    tx(96, 54, 'c', fontsize=8)
    paired_plot(c, d['gg'], .001, 'Quantified protein groups\nper acquisition (thousands)', (0, 7))
    sp = ax(111, 94, 64, 5)
    sp.axis('off'); sp.set_xlim(-.5, 4.5); sp.set_ylim(0, 1)
    sp.text(-.62, .4, '2 or more peptides\n(median %)', ha='right', va='center', fontsize=6)
    for i, a in enumerate(ARMS):
        sp.text(i, .4, f'{d["support"][a]:.1f}', ha='center', va='center', fontsize=6.5)


def panel_d(fig, ax, tx, d):
    """Panel d: stage at which unrecovered reference observations were lost (producer 228-245).

    Statistical unit: peptide-acquisition observation (n = 53,737 summed over 15 acquisitions).
    """
    ss = d['ss']
    tx(3, 104, 'd–f', fontsize=8)
    loss = ax(43, 112, 34, 27)
    stage_labels = ['Not a candidate', 'Tryptic candidate;\nremoved by parsimony',
                    'Non-tryptic candidate;\nabsent from final FASTA',
                    'Final FASTA;\nnot tryptic', 'Final FASTA;\ntryptic, not accepted']
    nd = ss.no_direct_prediction_observations / 1000
    yes = (ss.peptide_acquisition_observations - ss.no_direct_prediction_observations) / 1000
    loss.barh(range(5), nd, color=fl.FL, height=.70, label='No exported prediction')
    loss.barh(range(5), yes, left=nd, color=fl.ORANGE, height=.70, label='Prediction exists')
    for i, n in enumerate(ss.peptide_acquisition_observations):
        loss.text(n/1000+.8, i, f'{int(n):,}', va='center', fontsize=5.8)
    loss.set_yticks(range(5), stage_labels, fontsize=5.9)
    loss.invert_yaxis(); loss.set_xlim(0, 45); loss.set_xticks([0, 20, 40])
    loss.set_xlabel(f'Unrecovered REF peptide–acquisition\nobservations (thousands), n = {int(ss.peptide_acquisition_observations.sum()):,}', fontsize=6)
    loss.tick_params(axis='y', length=0); loss.tick_params(axis='x', labelsize=6)
    handles, labels = loss.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(16/180, 1-101.5/250),
               ncol=2, frameon=False, fontsize=5.7, handlelength=1, columnspacing=.8, borderaxespad=0)


def panel_e(fig, ax, tx, d):
    """Panel e: pooled peptides by member assignment and pooled reference recovery (247-268).

    Statistical unit: distinct peptide pooled over 15 acquisitions; recovery denominator
    76,127 distinct reference peptides.
    """
    rs = d['rs']
    comp = ax(114, 112, 43, 27)
    left = np.zeros(5)
    class_styles = [('known_member_specific', 'One member', fl.FL),
                    ('known_member_shared', 'Several members', fl.GREY),
                    ('absent_eight_member_reference', 'No member', '#dddddd')]
    for field, label, color in class_styles:
        values = rs[field].to_numpy()/1000
        comp.barh(range(5), values, left=left, color=color, height=.7, label=label)
        left += values
    for i, arm in enumerate(ARMS):
        comp.text(2, i, f'{int(rs.loc[arm,"pooled_peptides"]):,}', color='white', fontsize=5.8, va='center')
        comp.text(1.32, i, f'{100*rs.loc[arm,"pooled_reference_recovery"]:.1f}',
                  transform=comp.get_yaxis_transform(), ha='center', va='center', fontsize=6.5)
    comp.text(1.32, 1.14, f'REF recovered\n(% of {int(rs.loc["REF", "pooled_peptides"]):,})', transform=comp.transAxes,
              ha='center', va='bottom', fontsize=6)
    comp.set_yticks(range(5), [LABEL[a] for a in ARMS], fontsize=6.1)
    comp.invert_yaxis(); comp.set_xlim(0, 92); comp.set_xticks([0, 40, 80])
    comp.set_xlabel('Pooled accepted peptides (thousands)', fontsize=6)
    comp.tick_params(axis='y', length=0); comp.tick_params(axis='x', labelsize=6)
    handles, labels = comp.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(96/180, 1-100.5/250),
               ncol=3, frameon=False, fontsize=5.6, handlelength=1, columnspacing=.7, borderaxespad=0)


def panel_f(fig, ax, tx, d):
    """Panel f: pooled member-specific peptides per member and database (producer 270-290).

    Statistical unit: distinct peptide pooled over 15 acquisitions, 8 members x 5 databases.
    The orange frames, the 'Members with 10 or more peptides' row and its key are drawn here
    as the producer drew them; they were removed from the v28d figure by a later PDF edit
    (ASSEMBLY.md step 2).
    """
    mt = d['mt']
    tx(3, 153, 'g', fontsize=8)
    hm = ax(42, 157, 104, 26)
    hm.imshow(np.log10(mt.to_numpy(float)+1), cmap=LinearSegmentedColormap.from_list('member', ['white', fl.FL]),
              vmin=0, vmax=4.8, aspect='auto')
    for i in range(8):
        for j in range(5):
            v = int(mt.iloc[i, j])
            hm.text(j, i, f'{v:,}', ha='center', va='center', fontsize=6.5, color='white' if v >= 1000 else fl.INK)
            if v < 10:
                hm.add_patch(plt.Rectangle((j-.5, i-.5), 1, 1, fill=False, ec=fl.ORANGE, lw=1.2))
    hm.set_yticks(range(8), [SHORT(m) for m in MEMBERS], fontsize=6.5, fontstyle='italic')
    hm.set_xticks(range(5), [LABEL[a] for a in ARMS], fontsize=6.1)
    hm.xaxis.tick_top(); hm.tick_params(length=0)
    for s in hm.spines.values():
        s.set_visible(False)
    count = ax(42, 185, 104, 5)
    count.axis('off'); count.set_xlim(-.5, 4.5); count.set_ylim(0, 1)
    count.text(-.6, .5, 'Members with\n10 or more peptides', ha='right', va='center', fontsize=6.4)
    for i, a in enumerate(ARMS):
        count.text(i, .5, f'{int((mt[a]>=10).sum())}/8', ha='center', va='center', fontsize=7)
    key = ax(151, 160, 26, 23); key.axis('off')
    key.text(0, .95, 'Pooled member-\nspecific peptides', fontsize=6.3, va='top')
    key.add_patch(plt.Rectangle((0, .22), .11, .12, fill=False, ec=fl.ORANGE, lw=1.2, transform=key.transAxes))
    key.text(.16, .28, 'Fewer than 10', fontsize=6.1, va='center')


def panel_g(fig, ax, tx, d):
    """Panel g: UpSet of FL peptides among members, 8 single-member and 6 largest shared sets.

    Producer lines 293-314. Statistical unit: distinct FL peptide pooled over 15 acquisitions.
    """
    up = d['up']
    tx(3, 199, 'h', fontsize=8)
    up_ax = ax(36, 204, 89, 17)
    u = up[up.drawn.eq(True)].copy(); u['shared'] = (u.degree > 1).astype(int)
    u = u.sort_values(['shared', 'n_peptides'], ascending=[True, False])
    assert len(u) == 14 and (PLOT_MIN_LENGTH != MIN_LENGTH or u.n_peptides.sum() == 65007)
    up_ax.bar(range(14), u.n_peptides/1000, color=fl.FL, width=.64)
    for i, n in enumerate(u.n_peptides):
        up_ax.text(i, n/1000+.8, f'{int(n):,}', ha='center', fontsize=5.6)
    up_ax.set_xlim(-.6, 13.6); up_ax.set_ylim(0, 42); up_ax.set_xticks([])
    up_ax.set_ylabel('FL peptides (thousands)', fontsize=6.3)
    up_ax.tick_params(axis='y', labelsize=6)
    mat_main = ax(36, 223, 89, 22)
    mat_main.set_xlim(-.6, 13.6); mat_main.set_ylim(7.5, -.5)
    for i, names in enumerate(u.members):
        included = names.split('|')
        inds = [j for j, m in enumerate(MEMBERS) if m in included]
        if len(inds) > 1:
            mat_main.plot([i, i], [min(inds), max(inds)], color=fl.INK, lw=.9)
        for j, m in enumerate(MEMBERS):
            mat_main.scatter(i, j, s=11, color=fl.INK if m in included else fl.GRID, lw=0)
    mat_main.set_yticks(range(8), [SHORT(m) for m in MEMBERS], fontstyle='italic', fontsize=6.1)
    mat_main.set_xticks([]); mat_main.tick_params(length=0)
    for s in mat_main.spines.values():
        s.set_visible(False)


def panel_h(fig, ax, tx, d):
    """Panel h: FL peptide length distribution by member-assignment class (producer 316-326).

    Statistical unit: distinct FL peptide pooled over 15 acquisitions (n = 63,640 one member,
    3,842 several members, 10,160 no member).
    """
    lt = d['lt']
    tx(136, 199, 'i', fontsize=8)
    length_ax = ax(145, 206, 31, 29)
    for cls, color in [('One member', fl.FL), ('Shared by members', fl.GREY), ('No member', '#bdbdbd')]:
        q = lt[lt['class'].eq(cls)].set_index('length').reindex(range(MIN_LENGTH, MAX_LENGTH + 1), fill_value=0)
        length_ax.plot(q.index, q.pct_within_class, color=color, lw=1.1,
                       label=cls.replace('Shared by members', 'Several members'))
    length_ax.set_xlim(MIN_LENGTH, MAX_LENGTH); length_ax.set_xticks([10, 30, 50])
    length_ax.set_xlabel('Peptide length (residues)', fontsize=6.1)
    length_ax.set_ylabel('Peptides in class (%)', fontsize=6.2)
    length_ax.tick_params(labelsize=6)
    length_ax.legend(loc='upper right', frameon=False, fontsize=5.7, handlelength=1)


def build_canvas(sd_root, bottle_pdf, out_dir, stem='CAMPI_Figure2_review'):
    """Draw the producer canvas, place the community artwork and render it.

    Mirrors the producer's main figure and its ``save(..., bottle=True)`` (producer lines
    164-373): matplotlib PDF, ``fl.enforce_nbt``, community artwork placed with PyMuPDF at
    ``BOTTLE_RECT_MM``, then a raster at ``RENDER_DPI``.

    Parameters
    ----------
    sd_root : pathlib.Path
        ``source_data/`` directory.
    bottle_pdf : pathlib.Path
        ``figures/Fig2/artwork/bottle_61.pdf``.
    out_dir : pathlib.Path
        Output directory (``build/Fig2``).

    Returns
    -------
    dict
        ``pdf`` and ``png`` paths and the loaded inputs ``d``.
    """
    out_dir = Path(out_dir)
    _producer_rc()
    d = load_inputs(sd_root)
    fig, ax, tx = canvas(PAGE_MM[1])
    for draw in (panel_a, panel_b, panel_c, panel_d, panel_e, panel_f, panel_g, panel_h):
        draw(fig, ax, tx, d)
    fl.enforce_nbt(fig)
    pdf = out_dir / (stem + '.pdf')
    fig.savefig(pdf)
    plt.close(fig)
    base = fitz.open(pdf); bottle = fitz.open(bottle_pdf)
    final = fitz.open(); page = final.new_page(width=base[0].rect.width, height=base[0].rect.height)
    page.show_pdf_page(page.rect, base, 0)
    rect = fitz.Rect(*(v * 72/25.4 for v in BOTTLE_RECT_MM))
    page.show_pdf_page(rect, bottle, 0)
    base.close(); final.save(pdf, garbage=4, deflate=True)
    final.close(); bottle.close()
    png = out_dir / (stem + '.png')
    with fitz.open(pdf) as doc:
        doc[0].get_pixmap(dpi=RENDER_DPI, alpha=False).save(png)
        text = doc[0].get_text()
    assert all(x in text for x in LABEL.values())
    return dict(pdf=pdf, png=png, d=d)


def crop_panel(png, letter, out_png):
    """Crop one panel from a canvas raster rendered at ``RENDER_DPI``.

    The same pixel box is used for the regenerated canvas and for the producer's own PNG,
    so the two crops are directly comparable.
    """
    from PIL import Image
    x0, y0, x1, y1 = (int(round(v / 25.4 * RENDER_DPI)) for v in PANEL_BOXES_MM[letter])
    with Image.open(png) as im:
        im.crop((x0, y0, x1, y1)).save(out_png)
    return out_png


def crop_panel_pdf(pdf, letter, out_pdf):
    """Write one panel of the canvas PDF as its own vector PDF (clip, no rescaling)."""
    x0, y0, x1, y1 = (v * 72 / 25.4 for v in PANEL_BOXES_MM[letter])
    clip = fitz.Rect(x0, y0, x1, y1)
    with fitz.open(pdf) as src, fitz.open() as dst:
        page = dst.new_page(width=clip.width, height=clip.height)
        page.show_pdf_page(page.rect, src, 0, clip=clip)
        dst.save(out_pdf, garbage=4, deflate=True)
    return out_pdf


# ---------------------------------------------------------------- statistics ------
def stat_a_database_entries(d):
    """Median database entries per arm (Fig. 2a; v28d line 77).

    Statistical unit: acquisition, n = 15 per database. Returns dict arm -> median.
    """
    return dict(d['db'])


def stat_b_peptide_medians(d):
    """Median accepted 9-50-residue peptides per acquisition (Fig. 2b; v30 paragraph 37).

    Statistical unit: acquisition, n = 15 per database. Returns dict arm -> median.
    """
    return {a: float(d['pp'][a].median()) for a in ARMS}


def stat_c_group_medians(d):
    """Median quantified protein groups per acquisition and per-acquisition comparisons.

    Fig. 2c; v28d line 77. Statistical unit: acquisition, n = 15 per database.

    Returns
    -------
    dict
        ``median`` (arm -> median), ``augmented_exceed_all`` (number of acquisitions in which
        FL + REF and FL + REF + MG each exceed both REF and MG), ``support`` (median % groups
        with >=2 assigned peptides).
    """
    gg = d['gg']
    exceed = int(((gg.FL_R > gg.REF) & (gg.FL_R > gg.MG) & (gg.FL_R_M > gg.REF) & (gg.FL_R_M > gg.MG)).sum())
    return dict(median={a: float(gg[a].median()) for a in ARMS}, augmented_exceed_all=exceed,
                support={a: float(d['support'][a]) for a in ARMS})


def stat_e_recovery(d):
    """Reference identity recovery, 9-50 residues (Fig. 2e; v30 paragraphs 38 and 40).

    Statistical unit: acquisition (n = 15) for the within-acquisition median; distinct
    peptide pooled over 15 acquisitions for pooled recovery (denominator 63,348).

    Returns
    -------
    dict
        ``median_within_FL`` (%), ``pooled`` (arm -> %), ``n_reference_pooled``,
        ``fl_one_member_fraction`` (% of pooled FL peptides assigned to one member),
        ``fl_one_member``, ``fl_pooled``.
    """
    rs = d['rs']
    return dict(median_within_FL=float(d['rr']['FL'].median()),
                pooled={a: float(100 * rs.loc[a, 'pooled_reference_recovery']) for a in ARMS},
                n_reference_pooled=int(rs.loc['REF', 'pooled_peptides']),
                fl_one_member=int(rs.loc['FL', 'known_member_specific']),
                fl_pooled=int(rs.loc['FL', 'pooled_peptides']),
                fl_one_member_fraction=float(100 * rs.loc['FL', 'known_member_specific'] / rs.loc['FL', 'pooled_peptides']))


def stat_d_loss_stages(d):
    """Unrecovered reference observations and the share lost because the sequence was absent.

    Fig. 2d; v30 paragraph 38. Statistical unit: peptide-acquisition observation, summed over 15
    acquisitions (FL arm, 9-50 residues).

    Returns
    -------
    dict
        ``n_unrecovered``, ``n_sequence_absent``, ``pct_sequence_absent`` and per-stage counts.
    """
    ss = d['ss']
    n = int(ss.peptide_acquisition_observations.sum())
    absent = int(ss.loc[STAGES_SEQUENCE_ABSENT, 'peptide_acquisition_observations'].sum())
    return dict(n_unrecovered=n, n_sequence_absent=absent, pct_sequence_absent=100 * absent / n,
                stages={s: int(ss.loc[s, 'peptide_acquisition_observations']) for s in STAGES})


def stat_f_members(d, rc):
    """Member-specific peptide counts at 9-50 residues (Fig. 2f; v30 paragraph 40).

    Statistical unit: distinct peptide pooled over 15 acquisitions; species (n = 8 members).
    Members are ranked by reference-search (REF) member-specific peptides.

    Returns
    -------
    dict
        ``fl_members_with_specific_9`` and ``mg_members_with_specific_9`` (members with any
        member-specific peptide), ``mg_zero_members_9``; for the four of the
        ``N_MOST_REPRESENTED`` most represented members other than ``MG_MISSED_MEMBER``, the
        minimum and maximum FL/REF and FL/MG percentages; FL, MG and REF counts for
        ``MG_MISSED_MEMBER`` and its FL/REF percentage; for the ``N_LEAST_REPRESENTED``
        members, whether FL >= MG and whether FL + REF >= REF; the acceptance windows
        (checked against every pooled peptide).
    """
    ms = d['ms_all']
    piv9 = ms[ms.min_length.eq(MIN_LENGTH_SENSITIVITY)].pivot(index='member', columns='arm', values='pooled_specific_peptides')
    ranked = piv9['REF'].sort_values(ascending=False).index.tolist()
    most = ranked[:N_MOST_REPRESENTED]
    assert MG_MISSED_MEMBER in most
    four = [m for m in most if m != MG_MISSED_MEMBER]
    least = ranked[-N_LEAST_REPRESENTED:]
    fl_ref = {m: float(100 * piv9.loc[m, 'FL'] / piv9.loc[m, 'REF']) for m in four}
    fl_mg = {m: float(100 * piv9.loc[m, 'FL'] / piv9.loc[m, 'MG']) for m in four}
    # The length windows are acceptance filters, not observed extremes: check that every
    # pooled peptide of every arm lies inside the 5-50 window before reporting it.
    lengths = [len(p) for a in ARMS for p in rc['pool'][a]]
    assert MIN_LENGTH_SENSITIVITY <= min(lengths) and max(lengths) <= MAX_LENGTH
    return dict(fl_members_with_specific_9=int((piv9['FL'] > 0).sum()),
                mg_members_with_specific_9=int((piv9['MG'] > 0).sum()),
                mg_zero_members_9=sorted(piv9.index[piv9['MG'].eq(0)].tolist()),
                n_four=len(four), most_represented=most,
                fl_ref_pct=fl_ref, fl_ref_min=min(fl_ref.values()), fl_ref_max=max(fl_ref.values()),
                fl_mg_pct=fl_mg, fl_mg_min=min(fl_mg.values()), fl_mg_max=max(fl_mg.values()),
                missed_fl=int(piv9.loc[MG_MISSED_MEMBER, 'FL']), missed_mg=int(piv9.loc[MG_MISSED_MEMBER, 'MG']),
                missed_fl_ref_pct=float(100 * piv9.loc[MG_MISSED_MEMBER, 'FL'] / piv9.loc[MG_MISSED_MEMBER, 'REF']),
                least_represented=least,
                fl_ge_mg={m: bool(piv9.loc[m, 'FL'] >= piv9.loc[m, 'MG']) for m in least},
                flr_exceeds_ref={m: bool(piv9.loc[m, 'FL_R'] >= piv9.loc[m, 'REF']) for m in least},
                window_5=f'{MIN_LENGTH}–{MAX_LENGTH}', window_9=f'{MIN_LENGTH_SENSITIVITY}–{MAX_LENGTH}')


def stat_h_lengths(d, rc):
    """Median peptide length by member-assignment class, 5-50 residues (Fig. 2h; v30 paragraph 40).

    Computed twice: from the pooled FL peptide recount, and as the weighted median of the
    length table; both must agree. Statistical unit: distinct FL peptide pooled over 15
    acquisitions.
    """
    lt = d['lt']
    out = {}
    for cls in ['One member', 'Shared by members', 'No member']:
        q = lt[lt['class'].eq(cls)].sort_values('length')
        cum = q.n_peptides.cumsum().to_numpy(); n = int(q.n_peptides.sum())
        lens = q.length.to_numpy()
        lo = lens[np.searchsorted(cum, (n + 1) // 2)]
        hi = lens[np.searchsorted(cum, n // 2 + 1)]
        table_median = (lo + hi) / 2
        assert table_median == rc['length_stats'][cls]['median']
        out[cls] = dict(n=n, median=float(table_median))
    return out


def stat_reservoir(sd_root):
    """Reservoir size for CAMPI (v28d line 75). Unit: sequence; one reservoir."""
    p = Path(sd_root) / 'Fig1' / 'acquisition_sizes_726.tsv'
    assert sha256(p) == SIZES_SHA256, 'shared size table changed'
    t = read(p)
    v = t[t.cohort.eq('CAMPI')].reservoir_entries.unique()
    assert len(v) == 1 and t.cohort.eq('CAMPI').sum() == N_ACQUISITIONS
    return int(v[0])


def stat_containment(sd_root):
    """Share of each supplied member reference proteome present in the reservoir (v28d line 75).

    Exact whole-protein sequence equality after I/L equivalence, counted per reference
    record. Statistical unit: species (n = 8 members); shuffled control rows are excluded
    here and returned separately.

    Returns
    -------
    dict
        ``per_member`` (% present), ``min``, ``max``, ``shuffled_max`` (% present among the
        shuffled control proteomes).
    """
    t = read(Path(sd_root) / 'Fig2' / 'reservoir_protein_presence_by_proteome.tsv')
    shuffled = t.label.str.endswith('_SHUFFLED')
    real = t[~shuffled]
    assert len(real) == N_MEMBERS
    pct = (100 * real.n_present_exact / real.n_proteins).set_axis(real.label)
    return dict(per_member=pct.to_dict(), min=float(pct.min()), max=float(pct.max()),
                shuffled_max=float((100 * t[shuffled].n_present_exact / t[shuffled].n_proteins).max()))


def stat_gain_mechanism(sd_root, destination='FLR'):
    """Where the add-back gains come from (v28d line 77).

    Ported from ``build_mechanistic_review.py`` (CAMPI variability analysis). Statistical
    unit: peptide-acquisition observation, summed over 15 acquisitions, for reference
    proteome add-back (destination ``FLR`` = FL + REF).

    Returns
    -------
    dict
        ``gained``, ``unavailable`` (gains not searchable in the selected FL database),
        ``pct_unavailable``, ``new_quantified``, ``new_identified_and_quantified``,
        ``pct_new_identified``.
    """
    sd = Path(sd_root) / 'Fig2'
    a = read(sd / 'campi_gain_availability.tsv')
    q = read(sd / 'campi_new_quantified_peptides.tsv')
    for _, x in a.groupby(['sample', 'destination']):
        assert x.gained_peptides.sum() == x.total_gained.iloc[0]
    for _, x in q.groupby(['sample', 'destination']):
        assert x.peptides.sum() == x.total_destination_only_quantified.iloc[0]
    x = a[a.destination == destination]; z = q[q.destination == destination]
    assert x['sample'].nunique() == N_ACQUISITIONS
    gained = int(x.gained_peptides.sum())
    unavailable = int(x[x.availability != 'digest_compatible_in_FL_fasta'].gained_peptides.sum())
    newq = int(z.peptides.sum())
    newid = int(z[z.category == 'newly_identified_and_quantified'].peptides.sum())
    return dict(gained=gained, unavailable=unavailable, pct_unavailable=100 * unavailable / gained,
                new_quantified=newq, new_identified_and_quantified=newid,
                pct_new_identified=100 * newid / newq)


def expected_unique(freq, N, k):
    """Expected distinct peptides when k of N spectra are drawn without replacement.

    Ported unchanged from the CAMPI variability ``analyze.py``. Each peptide with f accepted
    PSMs among N spectra is missed with probability C(N-f, k) / C(N, k) (hypergeometric).
    """
    freq = np.asarray(freq, dtype=np.int64); assert 0 <= k <= N and np.all((freq > 0) & (freq <= N))
    if k == N:
        return float(len(freq))
    if k == 0:
        return 0.
    possible = N - freq >= k; no = np.zeros(len(freq), float); f = freq[possible]
    no[possible] = np.exp(gammaln(N-f+1)-gammaln(N-f-k+1)-gammaln(N+1)+gammaln(N-k+1))
    return float(np.sum(1-np.clip(no, 0, 1)))


def stat_depth_standardization(sd_root, arm='FL'):
    """FL peptide-yield range before and after post-search subsampling to equal MS2 counts.

    v28d line 77. Statistical unit: acquisition, n = 15. Draw size k is the smallest MS2
    count over all 150 acquisition x arm x length rows (as in ``analyze.py``). The expected
    yield is recomputed here from the per-acquisition peptide PSM frequencies and checked
    against the stored table. No random draws: the expectation is analytic.

    Returns
    -------
    dict
        ``k``, ``observed_fold`` (max/min observed peptides), ``expected_fold`` (max/min
        expected peptides at k spectra), per-sample expected values.
    """
    sd = Path(sd_root) / 'Fig2'
    t = read(sd / 'campi_depth_standardization.tsv')
    k = int(t.ms2.min())
    assert (t.equal_ms2_draw == k).all()
    r = t[(t.arm == arm) & (t.min_length == MIN_LENGTH)].set_index('sample').loc[SAMPLES]
    exp = {}
    for s in SAMPLES:
        f = read(sd / f'campi_FL_5_frequencies_{s}.tsv').accepted_psms.to_numpy()
        assert len(f) == r.loc[s, 'peptides']
        exp[s] = expected_unique(f, int(r.loc[s, 'ms2']), k)
        assert math.isclose(exp[s], r.loc[s, 'expected_peptides_equal_ms2'], rel_tol=1e-9)
    e = pd.Series(exp)
    return dict(k=k, observed_fold=float(r.peptides.max() / r.peptides.min()),
                expected_fold=float(e.max() / e.min()), expected=exp)


def stat_discrepant_predictions(sd_root):
    """Same-spectrum audit: searchable despite de novo disagreement (v28d line 81).

    Population ``production_retained_reference_paired``: spectra whose retained production
    de novo prediction is paired with a reference identification; discrepant = every category
    except exact I/L match. Statistical unit: spectrum (observation), pooled over 15
    acquisitions.

    Returns
    -------
    dict
        ``n_discrepant``, ``n_reference_in_final`` and its %, ``n_direct_lookup`` and its %.
    """
    t = read(Path(sd_root) / 'Fig2' / 'campi_discrepant_prediction_pooled_summary.tsv')
    x = t[(t.population == 'production_retained_reference_paired') & t.category.isin(DISCREPANT_CATEGORIES)]
    assert len(x) == len(DISCREPANT_CATEGORIES)
    assert (x.reference_present_final_evaluable == x.n).all()
    n = int(x.n.sum()); present = int(x.reference_present_final.sum()); look = int(x.any_lookup.sum())
    return dict(n_discrepant=n, n_reference_in_final=present, pct_reference_in_final=100 * present / n,
                n_direct_lookup=look, pct_direct_lookup=100 * look / n)


def stat_alphanovo_accuracy(sd_root, window='9-50'):
    """Top-ranked AlphaNovo prediction accuracy against the reference search (v28d line 81).

    Statistical unit: spectrum for pooled fractions; acquisition (n = 15) for the median and
    range. The median and range are recomputed from the per-acquisition table and checked
    against the pooled table.
    """
    sd = Path(sd_root) / 'Fig2'
    pooled = read(sd / 'alphanovo_accuracy_pooled.tsv').set_index('window').loc[window]
    per = read(sd / 'alphanovo_accuracy_per_acquisition.tsv')
    per = per[per.window == window]
    assert len(per) == N_ACQUISITIONS
    n_ref = int(per.n_reference_spectra.sum()); n_exact = int(per.n_exact_all.sum())
    n_ret = int(per.n_retained_top30.sum()); n_exact_ret = int(per.n_exact_retained.sum())
    assert (n_ref, n_exact, n_ret, n_exact_ret) == (int(pooled.n_reference_spectra), int(pooled.n_exact_all),
                                                    int(pooled.n_retained_top30), int(pooled.n_exact_retained))
    frac = per.n_exact_all / per.n_reference_spectra
    assert math.isclose(frac.median(), pooled.median_frac_exact_all, rel_tol=1e-12)
    return dict(n_reference=n_ref, n_exact=n_exact, pct_exact=100 * n_exact / n_ref,
                median_pct=float(100 * frac.median()), min_pct=float(100 * frac.min()),
                max_pct=float(100 * frac.max()), n_retained=n_ret, n_exact_retained=n_exact_ret,
                pct_exact_retained=100 * n_exact_ret / n_ret)


def recompute_numbers(sd_root, d, rc):
    """Every number in NUMBERS.tsv, in file order, as (v28d_line, name, value).

    Rows follow the v30 text (paragraphs 36-40 and the legend, paragraph 42; v28d lines 75-87).
    Values are full precision as computed. ``None`` marks a printed number that the
    repository cannot recompute (status OPEN in NUMBERS.tsv; reason in OPEN_ITEMS.md).
    """
    a = stat_a_database_entries(d); b = stat_b_peptide_medians(d); c = stat_c_group_medians(d)
    dd = stat_d_loss_stages(d); e = stat_e_recovery(d); f = stat_f_members(d, rc); h = stat_h_lengths(d, rc)
    cont = stat_containment(sd_root); gain = stat_gain_mechanism(sd_root)
    depth = stat_depth_standardization(sd_root); disc = stat_discrepant_predictions(sd_root)
    acc = stat_alphanovo_accuracy(sd_root)
    n_acq = len(d['pp'])
    return [
        ('75', 'acquisitions', n_acq),
        ('75', 'members (word)', as_word(len(cont['per_member']))),
        ('75', 'reservoir (millions)', stat_reservoir(sd_root) / 1e6),
        ('75', 'containment min (%)', cont['min']),
        ('75', 'containment max (%)', cont['max']),
        ('77', 'median entries FL', a['FL']),
        ('77', 'median entries REF', a['REF']),
        ('77', 'median entries MG', a['MG']),
        ('77', 'median peptides FL (9-50)', b['FL']),
        ('77', 'length window', f['window_9']),
        ('77', 'median peptides REF (9-50)', b['REF']),
        ('77', 'median peptides MG (9-50)', b['MG']),
        ('77', 'median groups FL', c['median']['FL']),
        ('77', 'median groups REF', c['median']['REF']),
        ('77', 'median groups MG', c['median']['MG']),
        ('77', 'median groups FL + REF', c['median']['FL_R']),
        ('77', 'median groups FL + REF + MG', c['median']['FL_R_M']),
        ('77', 'acquisitions where both augmented arms exceed REF and MG', c['augmented_exceed_all']),
        ('77', 'length window (add-back gains)', f['window_5']),
        ('77', 'gains not searchable in FL database (%)', gain['pct_unavailable']),
        ('77', 'length window (MS2 subsampling)', f['window_5']),
        ('77', 'observed FL peptide max/min', depth['observed_fold']),
        ('77', 'expected FL peptide max/min at equal MS2', depth['expected_fold']),
        ('79', 'median within-acquisition recovery FL (%)', e['median_within_FL']),
        ('79', 'length window', f['window_9']),
        ('79', 'pooled recovery FL (%)', e['pooled']['FL']),
        ('79', 'pooled reference identities', e['n_reference_pooled']),
        ('79', 'unrecovered observations', dd['n_unrecovered']),
        ('79', 'sequence absent from selected database (%)', dd['pct_sequence_absent']),
        ('79', 'pooled recovery FL + REF (%)', e['pooled']['FL_R']),
        ('79', 'pooled recovery FL + REF + MG (%)', e['pooled']['FL_R_M']),
        ('81', 'length window', f['window_9']),
        ('81', 'AlphaNovo exact (%)', acc['pct_exact']),
        ('81', 'reference-identified spectra', acc['n_reference']),
        ('81', 'per-acquisition median exact (%)', acc['median_pct']),
        ('81', 'per-acquisition min exact (%)', acc['min_pct']),
        ('81', 'per-acquisition max exact (%)', acc['max_pct']),
        ('81', 'filter-passing spectra', acc['n_retained']),
        ('81', 'score filter fraction (%)', None),
        ('81', 'exact among filter-passing (%)', acc['pct_exact_retained']),
        ('81', 'discrepant observations', disc['n_discrepant']),
        ('81', 'discrepant but searchable (%)', disc['pct_reference_in_final']),
        ('81', 'discrepant and searchable', disc['n_reference_in_final']),
        ('81', 'discrepant predictions with direct lookup (%)', disc['pct_direct_lookup']),
        ('83', 'members with FL-specific peptides at 9-50 (word)', as_word(f['fl_members_with_specific_9'])),
        ('83', 'MG members with specific peptides at 9-50 (word)', as_word(f['mg_members_with_specific_9'])),
        ('83', 'length window', f['window_9']),
        ('83', 'most represented members excluding A. caccae (word)', as_word(f['n_four'])),
        ('83', 'most represented members (word)', as_word(len(f['most_represented']))),
        ('83', 'FL/REF min (%)', f['fl_ref_min']),
        ('83', 'FL/REF max (%)', f['fl_ref_max']),
        ('83', 'FL/MG min (%)', f['fl_mg_min']),
        ('83', 'FL/MG max (%)', f['fl_mg_max']),
        ('83', 'A. caccae FL member-specific', f['missed_fl']),
        ('83', 'A. caccae MG member-specific', f['missed_mg']),
        ('83', 'A. caccae FL/REF (%)', f['missed_fl_ref_pct']),
        ('83', 'least represented members (word)', as_word(len(f['least_represented']))),
        ('83', 'FL peptides assigned to one member (%)', e['fl_one_member_fraction']),
        ('83', 'length window', f['window_9']),
        ('83', 'median length shared (5-50)', h['Shared by members']['median']),
        ('83', 'median length one member (5-50)', h['One member']['median']),
        ('83', 'length window', f['window_5']),
        ('87', 'acquisitions (word)', as_word(n_acq).capitalize()),
        ('87', 'databases (word)', as_word(len(ARMS))),
        ('87', 'minimum assigned peptides per group', None),
        ('87', 'members (word)', as_word(len(cont['per_member']))),
        ('87', 'pooled reference identities', e['n_reference_pooled']),
        ('87', 'length window b, d-g', f['window_9']),
        ('87', 'length window c, h', f['window_5']),
        ('87', 'peptide q threshold', None),
    ]
