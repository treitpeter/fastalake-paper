"""A4 collector: old (production distance_refit, 12 participants, 331 roster) vs new (run347, 13 participants)."""
import json, hashlib
from pathlib import Path
import pandas as pd
A = Path(__file__).resolve().parents[1]
B = (Path(__import__('os').environ['FASTALAKE_DATA_DIR'])/'review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit')
ctl = json.loads((A/'control331/COMPLETE.json').read_text()); new = json.loads((A/'run347/COMPLETE.json').read_text())
assert ctl['status'] == new['status'] == 'PASS' and new['mode'] == 'fixed347' and len(new['donors']) == 13
for d in ('control331', 'run347'):
    c = json.loads((A/d/'COMPLETE.json').read_text())
    for f, h in c['outputs'].items(): assert hashlib.sha256((A/d/f).read_bytes()).hexdigest() == h, (d, f)
old = pd.read_csv(B/'distance_refit/controlled_176_comparisons.tsv', sep='\t'); assert len(old) == 176
oc = old[old.cohort == 'CAPSCAN']; nc = pd.read_csv(A/'run347/controlled_CAPSCAN_96.tsv', sep='\t')
newall = pd.concat([old[old.cohort == 'PREDICT'], nc[old.columns]], ignore_index=True); assert len(newall) == 176
newall.to_csv(A/'controlled_176_comparisons_347.tsv', sep='\t', index=False)
keys = ['features', 'distance', 'adjustment', 'outer_correlation', 'layer_a', 'layer_b']
m = oc.merge(nc, on=['cohort'] + keys, suffixes=('_331', '_347'), validate='one_to_one'); assert len(m) == 96
m.to_csv(A/'CAPSCAN_controlled_331_vs_347.tsv', sep='\t', index=False)
def summary(t):
    r = t[(t.adjustment == 'region_adjusted') & (t.distance == 'rank')]; h = t[(t.adjustment == 'region_adjusted') & (t.distance == 'hellinger')]
    ns = r[r.q >= .05]; pairs = ns.groupby(['layer_a', 'layer_b']).size().to_dict()
    return dict(rank_adj_n=len(r), rank_adj_not_sig=len(ns), rank_adj_not_sig_by_pair={' vs '.join(k): v for k, v in pairs.items()},
                hell_adj_n=len(h), hell_adj_sig=int((h.q < .05).sum()), n_participants=sorted(set(t.n_participants)), n_profiles=sorted(set(t.n_profiles)),
                all_capscan_sig=int((t.q < .05).sum()), all_capscan=len(t))
so, sn = summary(oc), summary(nc)
pair_names = [('Non-host protein', 'Metabolites'), ('Host protein', 'Metabolites'), ('DNA', 'Metabolites'), ('DNA', 'Host protein'), ('DNA', 'Non-host protein'), ('Non-host protein', 'Host protein')]
L = ['# A4 controlled distance sensitivity (Supplementary Methods 8, Supplementary Fig. 11e,f): 331 vs 347', '',
     'Old = production biology_refit/distance_refit (12 DNA-complete participants, 331 roster; reproduced exactly by control331, atol 1e-10).',
     'New = run347 (13 DNA-complete participants incl. participant 15, 347 roster, C2 groups_347). PREDICT rows unchanged (no CapScan input).', '',
     '| printed statement | old (331) | new (347) |', '|---|---|---|',
     f"| controlled estimates retained | 176 | {len(newall)} (unchanged by construction: grid size does not depend on n) |",
     f"| region-adjusted rank CapScan comparisons | {so['rank_adj_n']} | {sn['rank_adj_n']} |",
     f"| of these, not meeting q<0.05 ('11 of 24') | {so['rank_adj_not_sig']} of {so['rank_adj_n']} | {sn['rank_adj_not_sig']} of {sn['rank_adj_n']} |"]
for a, b in pair_names:
    k = a + ' vs ' + b
    L.append(f"| not significant, {a}–{b} (of 4) | {so['rank_adj_not_sig_by_pair'].get(k, 0)} | {sn['rank_adj_not_sig_by_pair'].get(k, 0)} |")
L += [f"| region-adjusted Hellinger comparisons meeting q<0.05 ('All 24') | {so['hell_adj_sig']} of {so['hell_adj_n']} | {sn['hell_adj_sig']} of {sn['hell_adj_n']} |",
      f"| participants / profiles in CapScan rows | {so['n_participants']} / {so['n_profiles']} | {sn['n_participants']} / {sn['n_profiles']} |",
      f"| all 96 CapScan controlled rows with q<0.05 | {so['all_capscan_sig']} | {sn['all_capscan_sig']} |", '',
      'Classification changes (q<0.05 in one population only), all 96 controlled CapScan rows:', '']
ch = m[(m.q_331 < .05) != (m.q_347 < .05)]
L += ['| features | distance | adjustment | outer | layers | rho 331 | q 331 | rho 347 | q 347 |', '|---|---|---|---|---|---|---|---|---|']
for _, r in ch.iterrows():
    L.append(f"| {r.features} | {r.distance} | {r.adjustment} | {r.outer_correlation} | {r.layer_a}–{r.layer_b} | {r.coefficient_331:.3f} | {r.q_331:.4f} | {r.coefficient_347:.3f} | {r.q_347:.4f} |")
L += ['', 'Full row-level comparison: CAPSCAN_controlled_331_vs_347.tsv; 176-row table for Source Data: controlled_176_comparisons_347.tsv.',
      'Distance cross-checks (347): DNA and metabolite region-adjusted Hellinger distances equal cohort_reconciliation current331_plus16_all: ' + json.dumps(new['distance_errors'])]
(A/'RESULTS.md').write_text('\n'.join(L) + '\n'); print('\n'.join(L))
