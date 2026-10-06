from pathlib import Path
import json,sys,hashlib
from datetime import datetime,timezone
import numpy as np
from scipy import sparse

H=Path(__file__).resolve().parent
mode=sys.argv[1]
assert mode in ('pilot','full')
spec=json.loads((H/'INPUTS.json').read_text())
out=H/('pilot_gate' if mode=='pilot' else 'full_collection')
out.mkdir(exist_ok=False)
rows=[]
for task in spec['tasks']:
    d=H/mode/task['id'];r=json.loads((d/'COMPLETE.json').read_text())
    assert r['status']=='PASS' and r['arm']==task['id'] and r['mode']==mode
    assert r['acquisitions']==task['acquisitions']
    z=np.load(d/'sample_totals.npz')
    for level in ('KO','OG','Pathway'):
        matrix=sparse.load_npz(d/f'{level}_evidence_signal.npz')
        assert matrix.shape[1]==task['acquisitions']
        assert np.allclose(np.asarray(matrix.sum(axis=0)).ravel(),z['intensity'],rtol=1e-9,atol=1e-6)
    assert r['frozen_complete_parent_graph_preserved']
    if mode=='pilot' and task['id']=='CAMPI_REF':
        original=json.loads((H.parent.parent/'evidence_all_cohorts/full/CAMPI_REF/COMPLETE.json').read_text())
        assert r['query_peptides']==original['query_peptides']
        assert np.array_equal(z['features'],z['source_features'])
        assert np.allclose(z['intensity'],z['source_intensity'],rtol=1e-9,atol=1e-6)
    rows.append({k:r[k] for k in ('arm','acquisitions','positive_acquisitions','query_peptides','compatible_parents','parent_edges','quantified_features','total_intensity','job_id')})
result=dict(status='PASS',utc=datetime.now(timezone.utc).isoformat(),mode=mode,arms=len(rows),
    acquisition_arm_records=sum(x['acquisitions'] for x in rows),results=rows,
    interpretation='Quantitative evidence reconstruction passed. Biological claims and figure disposition are not inferred from successful computation.',
    remaining=['Assess function/taxonomy coverage and missing annotation namespaces.','Reconcile biological specimen metadata, then use participant-level aggregation.','Rebuild species/evenness and functional conclusions only at supported resolution; synchronize portable code and manuscript.'])
assert result['acquisition_arm_records']==770 and result['arms']==8
(out/'COMPLETE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
