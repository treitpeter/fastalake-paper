"""Current 717-acquisition bag/mass and fragment analysis (v42, panels e,f).

Drawing follows the manuscript producer. Totals are independently recomputed from
per-acquisition rows; historical 48-acquisition inputs remain in source_data.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory
from style import fl_style as fl
from figures.SFig02 import panels

ARMS = ('bags_le2', 'bags_le4', 'bags_le8', 'bags_all', 'mass_all', 'comb_all')
FR = [('bags_all', 'Bags: unrestricted'), ('bags_le8', 'Bags: ≤8 orderings')]
W_EF, H_EF = 155.0, 117.1


def summarize(sd):
    d = pd.read_csv(Path(sd) / 'bags_bounded717_per_acquisition.tsv', sep='\t')
    f = pd.read_csv(Path(sd) / 'fragment717_per_acquisition.tsv', sep='\t')
    assert len(d) == 717 and not d.duplicated(['cohort', 'sample']).any()
    assert d.groupby('cohort').size().to_dict() == {'CAMPI': 15, 'HSTOOL': 50, 'CAPSCAN': 346, 'PREDICT': 306}
    assert not f.duplicated(['cohort', 'sample', 'arm']).any()
    assert set(map(tuple, f[['cohort', 'sample']].values)) == set(map(tuple, d[['cohort', 'sample']].values))
    totals = {}
    for arm in ARMS:
        g, l, n = (d[arm + '_' + k] for k in ('gained', 'lost', 'net'))
        assert (g - l == n).all()
        totals[arm] = dict(gained=int(g.sum()), lost=int(l.sum()), net=int(n.sum()), positive=int((n > 0).sum()))
    assert (f.zero_diagnostic + f.positive_diagnostic == f.composition_preserving_pairs).all()
    fragment = {}
    for arm, _ in FR:
        sub = f[f.arm.eq(arm)]
        assert len(sub) == 717
        fragment[arm] = {k:int(sub[k].sum()) for k in ('zero_diagnostic', 'positive_diagnostic', 'composition_preserving_pairs')}
    return pd.DataFrame(totals).T, fragment


def numbers(sd):
    summary, fr = summarize(sd)
    result = {'acquisitions': 717, 'arms': 7}
    roster = pd.read_csv(Path(sd) / 'bags_bounded717_per_acquisition.tsv', sep='\t')
    result.update({'cohort_' + k: int(v) for k, v in roster.groupby('cohort').size().items()})
    for arm, row in summary.iterrows():
        result.update({arm + '_' + k: int(v) for k, v in row.items()})
    for arm, row in fr.items():
        result.update({arm + '_' + k: int(v) for k, v in row.items()})
        result[arm + '_supported_pct'] = 100 * row['positive_diagnostic'] / row['composition_preserving_pairs']
    return result


def render(sd, out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    summary, fr = summarize(sd)
    mpl.rcdefaults()
    fl.init()
    fig2, axes2 = plt.subplots(2, 1, figsize=fl.mm(W_EF, H_EF), gridspec_kw={'height_ratios': [4, 1.3]})
    fig2.subplots_adjust(left=0.205, right=0.915, top=0.955, bottom=0.135, hspace=0.42)
    e = axes2[0]
    order = ['bags_le2', 'bags_le4', 'bags_le8', 'bags_all', 'mass_all', 'comb_all']
    labels = ['Bags: ≤2 orderings', 'Bags: ≤4 orderings', 'Bags: ≤8 orderings', 'Bags: unrestricted', 'Mass substitutions', 'Bags + mass']
    ss = summary.loc[order]; yy = np.arange(len(order))
    e.barh(yy, ss.gained, color=fl.FL, label='Gained')
    e.barh(yy, -ss.lost, color=fl.VERM, label='Lost')
    e.axvline(0, color=fl.INK, lw=1)
    e.set_yticks(yy, [lab if '≤' not in lab else '' for lab in labels]); e.invert_yaxis()
    # production: xlim (-42000, 48000), label offset 600, ticks +-30,000 for totals up to ~33,000;
    # same proportions for totals up to ~569,000
    SCALE = 15
    e.set_xlim(-760000, 640000)   # wider on the lost side so the 568,964 label clears the tick labels
    for i, (_, r) in enumerate(ss.iterrows()):
        e.text(-r.lost - 600 * SCALE, i, f'{int(r.lost):,}', ha='right', va='center', fontsize=6.5, color=fl.INK)
        e.text(r.gained + 600 * SCALE, i, f'{int(r.gained):,}', va='center', fontsize=6.5, color=fl.INK)
    e.set_xticks([-500000, 0, 500000], ['500,000 lost', '0', '500,000 gained'])
    e.set_xlabel('Peptide–acquisition identifications relative to exact lookup')
    e.legend(loc='upper left', frameon=False, fontsize=7)
    fl.style(e)
    net = e.twinx(); net.set_ylim(e.get_ylim()); net.set_yticks(yy, [f'{int(r.net):+,}' for _, r in ss.iterrows()])
    for sp in ('top', 'right', 'left', 'bottom'):
        net.spines[sp].set_visible(False)
    net.tick_params(axis='y', length=0, labelsize=6.5, pad=2); net.set_ylabel('Net', rotation=270, labelpad=9)

    f = axes2[1]
    fy = np.arange(len(FR))
    xmax = max(v['composition_preserving_pairs'] for v in fr.values())
    for j, (k, _) in enumerate(FR):
        z, p = fr[k]['zero_diagnostic'], fr[k]['positive_diagnostic']
        f.barh([j], [z], color=fl.GREY, label='No diagnostic ion' if j == 0 else None)
        f.barh([j], [p], left=[z], color=fl.GREEN, label='≥1 diagnostic ion' if j == 0 else None)
        f.text(z / 2, j, f'{z:,}', ha='center', va='center', fontsize=7, color='white')
        f.text(z + p + 60, j, f'{p:,} ({100 * p / (z + p):.1f}%)', ha='left', va='center', fontsize=6.5, color=fl.INK)
    f.set_xlim(0, 8000); f.set_yticks(fy, ['' for _ in FR]); f.invert_yaxis()
    f.set_xticks([0, 2000, 4000, 6000, 8000], ['0', '2,000', '4,000', '6,000', '8,000'])
    f.set_xlabel('Order-related same-spectrum comparison pairs')
    fl.style(f)
    fl.panel(e, 'e'); fl.panel(f, 'f')
    fl.normalise_type(fig2, 'pre-align', verbose=False)
    fig2.canvas.draw()
    tick_tr = blended_transform_factory(e.transAxes, e.transData)
    pad_ax = (e.yaxis.get_major_ticks()[0].get_pad() + 2.6) / 72 / (e.get_position().width * fig2.get_figwidth())
    for i, lab in enumerate(labels):
        if '≤' in lab:
            pre, suf = lab.split('≤')
            panels.draw_rel_label(e, -pad_ax, i, pre, 'le', suf, tick_tr, fontsize=7, color=fl.INK, ha='right')
    ftr = blended_transform_factory(f.transAxes, f.transData)
    fpad = (f.yaxis.get_major_ticks()[0].get_pad() + 2.6) / 72 / (f.get_position().width * fig2.get_figwidth())
    for j, (_, lab) in enumerate(FR):
        if '≤' in lab:
            pre, suf = lab.split('≤')
            panels.draw_rel_label(f, -fpad, j, pre, 'le', suf, ftr, fontsize=7, color=fl.INK, ha='right')
        else:
            f.text(-fpad, j, lab, transform=ftr, ha='right', va='center', fontsize=7, color=fl.INK)
    lt = f.transAxes
    for x0, col, txt in [(0.27, fl.GREY, 'No diagnostic ion'), (0.56, fl.GREEN, None)]:
        key = Rectangle((x0, -0.73), 0.035, 0.14, transform=lt, facecolor=col, edgecolor='none', clip_on=False)
        f.add_artist(key)
        if txt:
            f.text(x0 + 0.05, -0.66, txt, transform=lt, ha='left', va='center', fontsize=7, color=fl.INK)
        else:
            panels.draw_rel_label(f, x0 + 0.05, -0.66, '', 'ge', '1 diagnostic ion', lt, fontsize=7, color=fl.INK, ha='left')
    panels.align_letters(fig2, [[e, f]])
    assert net.get_ylim() == e.get_ylim()

    # ---- plotted values, read back from the artists -------------------------------------------
    bars = [p for p in e.patches if isinstance(p, Rectangle)]
    gained_drawn = [b.get_width() for b in bars[:6]]; lost_drawn = [-b.get_width() for b in bars[6:12]]
    plotted = []
    for i, arm in enumerate(order):
        r = ss.loc[arm]
        assert gained_drawn[i] == r.gained and lost_drawn[i] == r.lost
        plotted.append(dict(panel='e', arm=arm, label=labels[i], gained=int(gained_drawn[i]), lost=int(lost_drawn[i]),
                            net_label=net.get_yticklabels()[i].get_text(), net=int(r.net), acquisitions_positive=int(r.positive)))
    fb = [p for p in f.patches if isinstance(p, Rectangle) and not p.get_transform() == lt]
    for j, (k, lab) in enumerate(FR):
        z, p = fr[k]['zero_diagnostic'], fr[k]['positive_diagnostic']
        plotted.append(dict(panel='f', arm=k, label=lab, pairs=z + p, no_diagnostic=z, with_diagnostic=p,
                            supported_pct=round(100 * p / (z + p), 1)))
    pd.DataFrame(plotted).to_csv(out / 'PLOTTED_VALUES.tsv', sep='\t', index=False)
    print(pd.DataFrame(plotted).to_string())

    pdf, _ = fl.save(fig2, out / 'SFig02_continued_v42', formats=('pdf',))
    png = out / 'SFig02_continued_v42.png'
    px = panels.render_png(pdf, png, 300)
    plt.close(fig2)
    return {"pdf": pdf, "png": png, "plotted": plotted}
