"""Fig. 3 text check: transfer and pooling medians at 5-50 and 9-50 residues from
source_data/Fig3/transfer_pool_comparisons.tsv, by setting panels.MIN_LENGTH (no new search)."""
import os  # release: paths come from environment variables
from pathlib import Path
import sys; sys.path.insert(0,((os.environ.get('FASTALAKE_PAPER_REPO') or str(Path(__file__).resolve().parents[3])) + '/figures/Fig3'))
import pandas as pd, panels as P
c=pd.read_csv(((os.environ.get('FASTALAKE_PAPER_REPO') or str(Path(__file__).resolve().parents[3])) + '/source_data/Fig3/transfer_pool_comparisons.tsv'),sep='\t')
for L in (5,9):
    P.MIN_LENGTH=L
    t=P.stat_bc_transfer(c); print(L); print(t[['arm','design','field','median','minimum','maximum']].to_string())
    for k,v in P.stat_de_pool(c).items(): print(k, {kk:(vv if not isinstance(vv,dict) else (vv['median'],vv['minimum'],vv['maximum'])) for kk,vv in v.items()})
    print(P.stat_methods_donor_assignments(c)['differ'])
