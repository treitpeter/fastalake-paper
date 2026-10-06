"""NEW_VALUES.tsv for Figure 2 (len950). Reads only repository source_data; 9-50 = min_length 9 rows
(accepted 5-50 lists filtered to length >= 9), UpSet recounted by panels_len950.upset_at."""
import os, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
REPO = Path((os.environ.get('FASTALAKE_PAPER_REPO') or str(Path(__file__).resolve().parents[3])))
os.environ.setdefault('FASTALAKE_FONT_DIR', os.environ.get('FASTALAKE_FONT_DIR', ''))
sys.path[:0] = [str(HERE), str(REPO)]
import panels_len950 as P
SD = REPO / 'source_data'
rows = []
def add(panel, q, old, new, src):
    rows.append((panel, q, old, new, src))
D = {}
for L in (5, 9):
    P.PLOT_MIN_LENGTH = L
    d = P.load_inputs(SD); rc = P.recount_membership(SD, d)
    ms9 = d['ms_all'][d['ms_all'].min_length.eq(L)].pivot(index='member', columns='arm', values='pooled_specific_peptides')
    D[L] = dict(d=d, b=P.stat_b_peptide_medians(d), e=P.stat_e_recovery(d), dd=P.stat_d_loss_stages(d), ms=ms9)
f5 = lambda v: f'{v:.1f}'
PR = 'source_data/Fig2/paired_reference_recovery.tsv'; RS = 'source_data/Fig2/recovery_summary.tsv'
SS = 'source_data/Fig2/stage_summary.tsv'; MS = 'source_data/Fig2/member_summary.tsv'
flt = lambda s: f'{s}, min_length == 9 rows (accepted 5-50 list filtered to length >= 9)'
for a in P.ARMS:
    add('b', f'median accepted peptides per acquisition, {P.LABEL[a]}', int(D[5]['b'][a]), int(D[9]['b'][a]), flt(PR) + ', median of peptides over 15 acquisitions')
for a in P.ARMS:
    add('e', f'pooled peptides, {P.LABEL[a]}', int(D[5]['d']['rs'].loc[a, 'pooled_peptides']), int(D[9]['d']['rs'].loc[a, 'pooled_peptides']), flt(RS))
    for c in ['known_member_specific', 'known_member_shared', 'absent_eight_member_reference']:
        add('e', f'{c}, {P.LABEL[a]}', int(D[5]['d']['rs'].loc[a, c]), int(D[9]['d']['rs'].loc[a, c]), flt(RS))
    add('e', f'pooled reference recovery %, {P.LABEL[a]}', f5(D[5]['e']['pooled'][a]), f5(D[9]['e']['pooled'][a]), flt(RS))
add('e', 'pooled reference identities (label denominator)', D[5]['e']['n_reference_pooled'], D[9]['e']['n_reference_pooled'], flt(RS) + ', REF pooled_peptides')
add('e/text', 'FL pooled peptides assigned to one member %', f5(D[5]['e']['fl_one_member_fraction']), f5(D[9]['e']['fl_one_member_fraction']), flt(RS))
add('text 79', 'median within-acquisition reference recovery FL %', f5(D[5]['e']['median_within_FL']), f5(D[9]['e']['median_within_FL']), flt(PR) + ', median reference_recovery')
for s in P.STAGES:
    add('d', f'unrecovered observations, {s}', D[5]['dd']['stages'][s], D[9]['dd']['stages'][s], flt(SS))
    nd = lambda L: int(D[L]['d']['ss'].loc[s, 'no_direct_prediction_observations'])
    add('d', f'  of which no exported prediction, {s}', nd(5), nd(9), flt(SS))
add('d', 'unrecovered observations total (axis label n)', D[5]['dd']['n_unrecovered'], D[9]['dd']['n_unrecovered'], flt(SS))
add('d/text', 'sequence absent from selected database %', f5(D[5]['dd']['pct_sequence_absent']), f5(D[9]['dd']['pct_sequence_absent']), flt(SS) + ', first three stages / total')
for m in P.MEMBERS:
    for a in P.ARMS:
        add('f', f'pooled member-specific peptides, {m}, {P.LABEL[a]}', int(D[5]['ms'].loc[m, a]), int(D[9]['ms'].loc[m, a]), flt(MS))
for L in (5, 9):
    ms = D[L]['ms']; least = ms['REF'].sort_values().index[:3].tolist()
    D[L]['least'] = least; D[L]['fold'] = {m: ms.loc[m, 'REF'] / ms.loc[m, 'FL'] for m in least}
    D[L]['flr'] = all(ms.loc[m, 'FL_R'] > ms.loc[m, 'REF'] for m in least)
add('f/text 83', 'three least represented members (by REF)', '; '.join(D[5]['least']), '; '.join(D[9]['least']), flt(MS))
add('f/text 83', 'REF/FL fold, min', f'{min(D[5]["fold"].values()):.2f}', f'{min(D[9]["fold"].values()):.2f}', flt(MS))
add('f/text 83', 'REF/FL fold, max', f'{max(D[5]["fold"].values()):.2f}', f'{max(D[9]["fold"].values()):.2f}', flt(MS))
add('f/text 83', 'FL + REF exceeds REF for all three', D[5]['flr'], D[9]['flr'], flt(MS))
add('f/text 83', 'members with >=1 FL-specific peptide', int((D[5]['ms']['FL'] > 0).sum()), int((D[9]['ms']['FL'] > 0).sum()), flt(MS))
up5 = P.upset_at(SD, 5); up9 = P.upset_at(SD, 9)
u = lambda up: up[up.drawn].assign(sh=lambda x: (x.degree > 1).astype(int)).sort_values(['sh', 'n_peptides'], ascending=[True, False])
for up, L in ((up5, 5), (up9, 9)):
    D[L]['up'] = u(up)
add('g', 'drawn FL peptides total (14 sets)', int(D[5]['up'].n_peptides.sum()), int(D[9]['up'].n_peptides.sum()), 'recount of source_data/SFig03/CAMPI_membership_<S>.tsv, FL arm, length 5-50 vs 9-50')
for i in range(14):
    r5, r9 = D[5]['up'].iloc[i], D[9]['up'].iloc[i]
    add('g', f'bar {i+1}: 5-50 set [{r5.members}] | 9-50 set [{r9.members}]', int(r5.n_peptides), int(r9.n_peptides), 'recount of source_data/SFig03/CAMPI_membership_<S>.tsv, FL arm, length filter')
up9.to_csv(HERE / 'upset_intersections_FL_9_50.tsv', sep='\t', index=False)
with open(HERE / 'NEW_VALUES.tsv', 'w') as fh:
    fh.write('panel\tquantity\told_5_50\tnew_9_50\tsource_file_and_filter\n')
    for r in rows:
        fh.write('\t'.join(map(str, r)) + '\n')
print(open(HERE / 'NEW_VALUES.tsv').read())
