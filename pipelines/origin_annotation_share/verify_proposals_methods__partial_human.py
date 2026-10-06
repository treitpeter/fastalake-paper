import json,numpy as np
from scipy import sparse
from pathlib import Path
A=Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/github_workflow_reconciliation_2026-09-30/resumed_audit/annotation_repair/annotated_origin_consensus_2026-10-01/full')
for c in ['HSTOOL_noclust','CAPSCAN_noclust','PREDICT_noclust']:
    tx=json.load(open(A/c/'taxonomy_rows.json')); m=sparse.load_npz(A/c/'taxonomy_evidence_signal.npz').tocsr()
    tot=m.sum(); r=np.asarray(m.sum(1)).ravel()
    h=np.array([t['supported_lineage'][-1:]==['s__Homo sapiens'] for t in tx])
    miss=np.array([t['state']=='annotated_parent_lca' for t in tx])
    print(c, len(tx), m.shape, 'human %.3f%%'%(100*r[h].sum()/tot), 'human&missing %.3f%%'%(100*r[h&miss].sum()/tot))
