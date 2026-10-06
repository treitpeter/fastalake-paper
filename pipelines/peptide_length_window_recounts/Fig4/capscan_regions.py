"""v28d line 129 same-spectrum CapScan values by compartment at 5-50 and 9-50 (FL_public vs MG rows)."""
import os  # release: paths come from environment variables
from pathlib import Path
import pandas as pd, json, sys
REPO = (os.environ.get('FASTALAKE_PAPER_REPO') or str(Path(__file__).resolve().parents[3]))
b = pd.read_csv(f'{REPO}/source_data/SFig06/same_spectrum_benchmark_1_capscan.tsv', sep='\t')
b = b[b.test_arm.eq('FL_public') & b.reference_arm.eq('MG')]
out = {}
for L in (5, 9):
    d = b[b.min_length.eq(L)]; assert len(d) == 282
    for grp, g in d.groupby('group'):
        out[f'{L}_{grp}'] = dict(n=len(g), gain=100 * (g.count_ratio.median() - 1), fl_gt_mg=int((g.n_test > g.n_reference).sum()))
    f = d[d.group.ne('stool')]
    out[f'{L}_fluid'] = dict(n=len(f), gain=100 * (f.count_ratio.median() - 1), fl_gt_mg=int((f.n_test > f.n_reference).sum()))
print(json.dumps(out, indent=1)); open(sys.argv[1], 'w').write(json.dumps(out, indent=1))
