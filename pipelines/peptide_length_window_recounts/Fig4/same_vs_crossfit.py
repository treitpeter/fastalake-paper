"""Fig. 4 text checks at 9-50: same-spectrum vs cross-fitted gain counts (v28d line 133), from
source_data/Fig4/crossfit_per_acquisition.tsv; same rule as panels.stat_bc_same_vs_crossfit, only L changed.
Also confirms Fig. 4b,c plot the _9 columns (L_BC)."""
import os  # release: paths come from environment variables
import sys, json
from pathlib import Path
import pandas as pd
REPO = Path((os.environ.get('FASTALAKE_PAPER_REPO') or str(Path(__file__).resolve().parents[3])))
sys.path.insert(0, str(REPO / 'figures/Fig4'))
import panels as P
x = pd.read_csv(REPO / 'source_data/Fig4/crossfit_per_acquisition.tsv', sep='\t')
g_ = lambda fl, mg: fl / mg - 1
out = {'L_BC': P.L_BC}
for L in ('5', '9'):
    for coh, g in x.groupby('cohort'):
        same = g_(g[f'frozen_fl_{L}'], g[f'frozen_mg_{L}']); cf = g_(g[f'FL_pooled_{L}'], g[f'MG_pooled_{L}'])
        out[f'{coh}_{L}'] = [int((same > cf).sum()), len(g)]
out['repo_function_5'] = P.stat_bc_same_vs_crossfit(x)
print(json.dumps(out)); Path(sys.argv[1]).write_text(json.dumps(out, indent=1))
