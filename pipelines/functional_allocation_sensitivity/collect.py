"""A3 collector: allocation sensitivity, old (331, production job 6048208) vs new (347, CapScan rerun; HSTOOL/PREDICT unchanged)."""
import json, hashlib
from pathlib import Path
import pandas as pd
A = Path(__file__).resolve().parents[1]
P = A.parents[2]/'iteration1/proposals_methods/alloc_sensitivity/out'
for d in ('control331', 'run347'):
    c = json.loads((A/d/'COMPLETE.json').read_text()); assert c['status'] == 'PASS'
    for f, h in c['outputs'].items(): assert hashlib.sha256((A/d/f).read_bytes()).hexdigest() == h, (d, f)
new = json.loads((A/'run347/COMPLETE.json').read_text()); assert new['mode'] == 'fixed347' and len(new['participants']) == 15
old_m = pd.read_csv(P/'allocated_fraction_medians.tsv', sep='\t')
cap = pd.read_csv(A/'run347/allocated_fraction_medians_CAPSCAN.tsv', sep='\t')
new_m = pd.concat([old_m[old_m.arm != 'CAPSCAN_noclust'], cap], ignore_index=True); new_m.to_csv(A/'allocated_fraction_medians_347.tsv', sep='\t', index=False)
def rng(m):
    eq = m[m.variant == 'equal'].allocated_fraction; ot = m[m.variant != 'equal'].allocated_fraction
    return f"{100*ot.min():.1f}–{100*ot.max():.1f}% instead of {100*eq.min():.1f}–{100*eq.max():.1f}%"
oo = pd.read_csv(P/'capscan_own_stool.tsv', sep='\t'); no = pd.read_csv(A/'run347/capscan_own_stool.tsv', sep='\t')
m = oo.merge(no, on=['variant', 'prevalence'], suffixes=('_331', '_347'))
m.to_csv(A/'capscan_own_stool_331_vs_347.tsv', sep='\t', index=False)
real = m[~m.variant.str.startswith('NEGATIVE')]
L = ['# A3 functional-allocation sensitivity (Supplementary Methods 5): 331 vs 347', '',
     'Producer: claude_review_2026-10-03/iteration1/proposals_methods/alloc_sensitivity/alloc_sensitivity.py (job 6048208).',
     'Old = production out/ (CapScan 14 participants, 331 roster; reproduced exactly by control331). New = CapScan rebuilt from the C2 347',
     'annotated-origin consensus (15 participants); HSTOOL and PREDICT rows are the production values (no CapScan input).', '',
     '| printed statement | old (331) | new (347) |', '|---|---|---|',
     f"| allocated median of total intensity, restricted rules vs equal allocation (all cohorts and levels) | {rng(old_m)} | {rng(new_m)} |",
     f"| same, CapScan only | {rng(old_m[old_m.arm=='CAPSCAN_noclust'])} | {rng(cap)} |",
     f"| CapScan participants closer to own stool, all variants and prevalence filters | {real.own_closer_331.min()}–{real.own_closer_331.max()} of 14 | {real.own_closer_347.min()}–{real.own_closer_347.max()} of 15 |",
     f"| largest own-stool P across variants/filters | {real.p_331.max():g} | {real.p_347.max():g} |",
     '| PREDICT rho 0.16–0.64 vs 0.16–0.51; human-assigned 0.43–0.51 to 0.28–0.35 | unchanged (PREDICT only) | unchanged |',
     '| healthy stool r 0.89–0.91, 0.81 to 0.69–0.73 | unchanged (HSTOOL only) | unchanged |', '',
     '| variant | prevalence | own closer 331 | advantage 331 | P 331 | own closer 347 | advantage 347 | P 347 |', '|---|---|---|---|---|---|---|---|']
for _, r in m.iterrows():
    L.append(f"| {r.variant} | {r.prevalence} | {r.own_closer_331} | {r.mean_advantage_331:.4f} | {r.p_331:g} | {r.own_closer_347} | {r.mean_advantage_347:.4f} | {r.p_347:g} |")
L += ['', '| arm | level | variant | allocated median 331 | allocated median 347 |', '|---|---|---|---|---|']
cm = old_m[old_m.arm == 'CAPSCAN_noclust'].merge(cap, on=['arm', 'level', 'variant'], suffixes=('_331', '_347'))
for _, r in cm.iterrows(): L.append(f"| {r.arm} | {r.level} | {r.variant} | {100*r.allocated_fraction_331:.2f}% | {100*r.allocated_fraction_347:.2f}% |")
L += ['', 'Files: capscan_own_stool_331_vs_347.tsv, allocated_fraction_medians_347.tsv, run347/.']
(A/'RESULTS.md').write_text('\n'.join(L) + '\n'); print('\n'.join(L))
