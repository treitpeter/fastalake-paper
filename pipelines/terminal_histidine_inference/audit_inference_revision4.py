#!/usr/bin/env python3
"""Independent real-donor aggregation and exact-test audit of the balanced gate."""
import collections,itertools,json,os,sys
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
RUN=Path(sys.argv[1]).resolve();sys.path.insert(0,str(RUN/'frozen'))
from lib import dump,meta,verify
manifest=json.loads((RUN/'inputs/independent_audit_manifest_revision4.json').read_text())
for x in manifest['code']+manifest['inputs']:verify(x)
base=RUN/'results_revision4/statistics_smoke_CAPSCAN'
receipt=json.loads((base/'COMPLETE.json').read_text());assert receipt['status']=='PASS' and receipt['acquisitions']==32
for x in receipt['outputs']:verify(x)
output=RUN/'results_revision4/independent_audit';output.mkdir(exist_ok=True)
assert not (output/'COMPLETE.json').exists()
def read(name):return pd.read_csv(base/name,sep='\t',keep_default_na=False,na_values=[''])
def numeric(row,key):return float(row[key])
def aggregate(frame,keys,values):
    # Explicit nested arithmetic, independently of the collector's groupby chain.
    specimens=collections.defaultdict(lambda:collections.defaultdict(list))
    for row in frame.to_dict('records'):
        key=tuple(row[k] for k in keys)+(str(row['participant']),row['region'],row['specimen_id'])
        for name,formula in values.items():
            v=formula(row)
            if np.isfinite(v):specimens[key][name].append(v)
    donors=collections.defaultdict(lambda:collections.defaultdict(list))
    for key,value in specimens.items():
        for name,v in value.items():donors[key[:-1]][name].append(sum(v)/len(v))
    return {key:{name:sum(v)/len(v) for name,v in values.items()} for key,values in donors.items()}
def fraction(a,b):return a/b if np.isfinite(b) and b>0 else np.nan
source=aggregate(read('acquisition_enrichment.tsv.gz'),['route','side','origin','residue'],
    {'O':lambda r:fraction(r['observed'],r['denominator']),'E':lambda r:fraction(r['expected'],r['denominator'])})
for v in source.values():
    v['difference']=v.get('O',np.nan)-v.get('E',np.nan);v['enrichment']=fraction(v.get('O',np.nan),v.get('E',np.nan))
pe=read('participant_region_enrichment.tsv');source_checks=0
for row in pe.to_dict('records'):
    key=tuple(row[k] for k in ['route','side','origin','residue'])+(str(row['participant']),row['region'])
    expected=source[key]
    for k in ['O','E','difference','enrichment']:
        assert np.isclose(float(row[k]),expected.get(k,np.nan),atol=2e-12,rtol=2e-10,equal_nan=True),(key,k,row[k],expected)
        source_checks+=1
motif=aggregate(read('acquisition_motif_counts.tsv.gz'),['route','side','origin','population','position','residue'],
    {'difference':lambda r:fraction(r['observed'],r['observed_total'])-fraction(r['expected'],r['expected_total'])})
pm=read('participant_motif_differences.tsv.gz');motif_checks=0
for row in pm.to_dict('records'):
    key=tuple(row[k] for k in ['route','side','origin','population','position','residue'])+(str(row['participant']),row['region'])
    expected=motif.get(key,{}).get('difference',np.nan)
    assert np.isclose(float(row['difference']),expected,atol=2e-12,rtol=2e-10,equal_nan=True),(key,row['difference'],expected)
    motif_checks+=1
reciprocal=read('participant_reciprocal_controls.tsv')
REGIONS=['duodenum','jejunum','ileum','colon','stool']
subjects=sorted({key[-2] for key in source})
source_lookup=collections.defaultdict(dict)
for key,value in source.items():source_lookup[key[:-2]][key[-2:]]=value
motif_lookup=collections.defaultdict(dict)
for key,value in motif.items():motif_lookup[key[:-2]][key[-2:]]=value.get('difference',np.nan)
reciprocal_lookup=collections.defaultdict(list)
for row in reciprocal.to_dict('records'):reciprocal_lookup[(row['route'],row['origin'],row['region'])].append(float(row['difference']))
def regional(cells,donor,region,metric):
    group={'proximal':REGIONS[:2],'distal':REGIONS[2:4]}.get(region,[region])
    return sum(cells.get((donor,r),{}).get(metric,np.nan) for r in group)/len(group)
sign_vectors={n:np.asarray(list(itertools.product((-1.,1.),repeat=n))) for n in range(8,9)}
def exact(v):
    v=np.asarray(v,float);v=v[np.isfinite(v)];n=len(v)
    assert n<=8,'Balanced gate unexpectedly contains more than8 donors'
    mean=float(sum(v)/n) if n else np.nan
    if n<8:return n,mean,1.
    # Enumerate all 2^n sign assignments in a different order from the worker.
    means=np.sum(sign_vectors[n]*v,axis=1)/n
    return n,mean,float(np.count_nonzero(abs(means)>=abs(mean)-1e-12)/len(means))
tests=read('complete_statistical_families.tsv');counts=tests.groupby('family').size().to_dict()
assert counts=={'source_P1_all_residues':7800,'all_motif_positions_residues':21600,'H_spatial_all_routes':2184,'all_reciprocal_controls':60},counts
for row in tests.to_dict('records'):
    endpoint=row['endpoint']
    if endpoint=='source_enrichment':
        cells=source_lookup[(row['route'],row['side'],row['origin'],row['residue'])]
        values=[cells.get((d,row['contrast']),{}).get('difference',np.nan) for d in subjects]
    elif endpoint.startswith(('spatial_H_','ordinal_fluid_H_')):
        cells=source_lookup[(row['route'],row['side'],row['origin'],'H')]
        metric=endpoint.split('_H_',1)[1]
        if endpoint.startswith('spatial_'):
            a,b=row['contrast'].split('__');values=[regional(cells,d,a,metric)-regional(cells,d,b,metric) for d in subjects]
        else:values=[sum((i-1.5)*regional(cells,d,r,metric) for i,r in enumerate(REGIONS[:4]))/5 for d in subjects]
    elif endpoint=='motif':
        cells=motif_lookup[(row['route'],row['side'],row['origin'],row['population'],int(row['position']),row['residue'])]
        group={'proximal':REGIONS[:2],'distal':REGIONS[2:4],'stool':['stool']}[row['contrast']]
        values=[sum(cells.get((d,r),np.nan) for r in group)/len(group) for d in subjects]
    elif endpoint=='matched_reciprocal':values=reciprocal_lookup[(row['route'],row['origin'],row['contrast'])]
    else:raise AssertionError('Unrecognized endpoint '+endpoint)
    n,mean,p=exact(values)
    assert int(row['n'])==n,(endpoint,row,n)
    assert np.isclose(float(row['mean']),mean,atol=2e-12,rtol=2e-10,equal_nan=True),(endpoint,row,mean)
    assert np.isclose(float(row['p']),p,atol=1e-14,rtol=0),(endpoint,row,p)
qchecks=0
for family,g in tests.groupby('family'):
    q=multipletests(g.p.to_numpy(float),method='fdr_bh')[1]
    assert np.allclose(q,g.q.to_numpy(float),atol=1e-14,rtol=1e-12),family
    qchecks+=len(g)
estimable=tests[tests.estimable.eq(True)]
assert len(estimable)>0
for lo,hi in [('CI95_low','CI95_high'),('enrichment_CI95_low','enrichment_CI95_high')]:
    both=tests[[lo,hi]].dropna();assert len(both)>0 and (both[lo]<=both[hi]).all()
record=dict(status='PASS',source_donor_values_reconstructed=source_checks,motif_donor_values_reconstructed=motif_checks,
    effects_counts_and_exact_p_reconstructed=len(tests),q_values_independently_reconstructed=qchecks,
    family_counts=counts,estimable_by_family=estimable.groupby('family').size().to_dict(),
    limitations='Validation subset only; CIs checked for finite ordered bounds, not independently bootstrap-resampled; full cohort scientific review pending.')
dump(output/'AUDIT.json',record)
dump(output/'COMPLETE.json',dict(status='PASS',job_id=os.environ.get('SLURM_JOB_ID'),source_receipt=meta(base/'COMPLETE.json'),scientific_review='PENDING',outputs=[meta(output/'AUDIT.json')]))
print(json.dumps(record),flush=True)
