#!/usr/bin/env python3
"""Matched HSTOOL identity, feature/group precision and missingness comparisons."""
import gzip,itertools,json,os,sys
from pathlib import Path
import numpy as np
import pandas as pd
RUN=Path(sys.argv[1]).resolve();sys.path.insert(0,str(RUN/'frozen'))
from lib import dump,meta,verify
TIERS=['TechnicalReplicates','WorkflowReplicates','OnePiece','Intra','Interindividual']
def tier(s):return next(t for t in TIERS if t in s)
def corr(x,y):return float(np.corrcoef(x,y)[0,1]) if len(x)>2 and np.std(x)>0 and np.std(y)>0 else np.nan
def write(df,name):df.to_csv(RUN/'results'/name,sep='\t',index=False)
def load(p):
    z=np.load(p);return pd.DataFrame(z['values'],index=z['keys'],columns=z['samples'])
def finite_med(a):return float(np.median(a)) if len(a) else np.nan

m=json.loads((RUN/'inputs/manifest_revision2.json').read_text())
for x in m['frozen_code']+[m['protocol']]:verify(x)
names=[a['name'] for a in m['arms']];samples=m['samples'];receipts=[];ids={}
for arm in names:
    folder=RUN/'results'/arm;c=json.loads((folder/'COMPLETE.json').read_text());assert c['status']=='PASS' and not c['pilot'] and c['acquisitions']==58
    for x in c['outputs']:verify(x)
    receipts.append(meta(folder/'COMPLETE.json'))
    f=pd.read_csv(folder/'accepted_identities.tsv.gz',sep='\t');ids[arm]={s:set(f[f['sample'].eq(s)].peptide) for s in samples}
identity=[]
for arm in names[1:]:
    for s in samples:
        for lo in [5,9]:
            a={p for p in ids['original'][s] if len(p)>=lo};b={p for p in ids[arm][s] if len(p)>=lo}
            identity.append(dict(sample=s,tier=tier(s),test=arm,reference='original',min_length=lo,original=len(a),test_count=len(b),overlap=len(a&b),gained=len(b-a),lost=len(a-b),count_ratio=len(b)/len(a) if a else np.nan,original_identity_recovery=len(a&b)/len(a) if a else np.nan))
write(pd.DataFrame(identity),'identity_comparison.tsv')
write(pd.DataFrame(identity).groupby(['test','tier','min_length']).median(numeric_only=True).reset_index(),'identity_summary.tsv')
coverage=[];pairs=[];same=[];cvs=[];cvdetails=[]
spaces=['features','fixed_groups','fixed_groups.Human','fixed_groups.Microbial','fixed_groups.Mixed','fixed_groups.Unknown']
for space in spaces:
    frames={arm:load(RUN/'results'/arm/('same_file_id.'+space+'.npz')) for arm in names}
    vocab=sorted(set().union(*(set(f.index) for f in frames.values())));a={arm:f.reindex(index=vocab,columns=samples).fillna(0).to_numpy() for arm,f in frames.items()}
    for arm,v in a.items():
        for t in TIERS:
            ix=[i for i,s in enumerate(samples) if tier(s)==t];det=(v[:,ix]>0).sum(axis=1)
            for threshold,fraction in [('any',1/len(ix)),('half',.5),('80pct',.8),('all',1)]:
                required=int(np.ceil(len(ix)*fraction));coverage.append(dict(space=space,arm=arm,tier=t,threshold=threshold,required=required,acquisitions=len(ix),global_union_vocabulary=len(vocab),rows=int((det>=required).sum())))
    for test in names[1:]:
        x,y=a['original'],a[test]
        for i,s in enumerate(samples):
            keep=(x[:,i]>0)&(y[:,i]>0);lx=np.log2(x[keep,i]);ly=np.log2(y[keep,i]);delta=ly-lx
            same.append(dict(space=space,test=test,sample=s,tier=tier(s),shared_positive_rows=int(keep.sum()),log2_pearson=corr(lx,ly),median_log2_ratio=finite_med(delta),centered_log2_mad=finite_med(abs(delta-finite_med(delta)))))
        for t in TIERS:
            ix=[i for i,s in enumerate(samples) if tier(s)==t]
            for i,j in itertools.combinations(ix,2):
                keep=(x[:,i]>0)&(x[:,j]>0)&(y[:,i]>0)&(y[:,j]>0)
                for label,v in [('original',x),(test,y)]:
                    vals=v[keep][:,[i,j]];logs=np.log2(vals);d=logs[:,1]-logs[:,0]
                    paircv=100*np.std(vals,axis=1,ddof=0)/np.mean(vals,axis=1)
                    pairs.append(dict(space=space,comparison=test,method=label,tier=t,sample_a=samples[i],sample_b=samples[j],common_rows=int(keep.sum()),log2_pearson=corr(logs[:,0],logs[:,1]),median_pair_cv_pct=finite_med(paircv),centered_log2_mad=finite_med(abs(d-finite_med(d)))))
            if t not in TIERS[:2]:continue
            keep=(x[:,ix]>0).all(axis=1)&(y[:,ix]>0).all(axis=1)
            for label,v in [('original',x),(test,y)]:
                vals=v[keep][:,ix]
                for normalization in ['none','median_on_same_shared_complete_rows']:
                    data=vals if normalization=='none' or not len(vals) else vals/np.median(vals,axis=0)
                    cv=100*np.std(data,axis=1,ddof=1)/np.mean(data,axis=1) if len(data) else np.array([])
                    cvs.append(dict(space=space,comparison=test,method=label,tier=t,normalization=normalization,shared_complete_rows=len(cv),median_cv_pct=finite_med(cv),p25=float(np.percentile(cv,25)) if len(cv) else np.nan,p75=float(np.percentile(cv,75)) if len(cv) else np.nan))
                    if normalization=='none':cvdetails.extend(dict(space=space,comparison=test,method=label,tier=t,key=vocab[k],cv_pct=float(z)) for k,z in zip(np.flatnonzero(keep),cv))
write(pd.DataFrame(coverage),'completeness.tsv');write(pd.DataFrame(pairs),'replicate_pairs.tsv');write(pd.DataFrame(same),'same_acquisition_agreement.tsv');write(pd.DataFrame(cvs),'replicate_cv.tsv');write(pd.DataFrame(cvdetails),'replicate_cv_per_row.tsv.gz')
write(pd.DataFrame(pairs).groupby(['space','comparison','method','tier']).median(numeric_only=True).reset_index(),'replicate_pair_summary.tsv')
write(pd.concat([pd.read_csv(RUN/'results'/a/'transfer_coverage.tsv',sep='\t') for a in names]),'transfer_coverage_all.tsv')
write(pd.concat([pd.read_csv(RUN/'results'/a/'coverage.tsv',sep='\t') for a in names]),'coverage_all.tsv')
report='''# HSTOOL quantitative comparison

Five completed search strategies are independently recounted on exactly the same
58 acquisitions. Identification counts use target rank-1 peptide-q <=0.01 and
canonical bare sequences. Quantitative primary comparisons use positive finite
modified-peptide/combined-charge features with an accepted canonical peptide in
that same acquisition. Modification-specific same-acquisition evidence is saved
as a stricter alternative. Missing values remain zero/unquantified; no imputation.

The pre-existing study dictionary is a restricted diagnostic projection. Every
unmapped feature and positive cell is counted explicitly; no new representative
cover is fitted for this comparison. Local host/microbial contributions are
classified from actual parent sets and union sequence aliases, not the reporting
representative. The original Human/Microbial matrices must reproduce the paper's
validated local-contribution baseline before this report can pass.

Paired precision comparisons use identical rows positive under both methods in
both acquisitions. All-replicate CV uses identical complete rows, separately for
12 technical and 15 workflow replicates, with ddof=1. Pair CV uses ddof=0 to match
the manuscript's original pair statistic. These are different estimators and are
named separately. Overlapping pairs are descriptive, not independent replicates.
Across-person differences are biological similarity, not technical precision.

Unfiltered old LFQ q values are not interpreted as calibrated confidence.
Corrected joint pooled q<=0.01 and q<=0.05 are separate sensitivity exports.
Their independently verified q mapping is pooled across features; it does not
establish a per-acquisition error rate for transferred-only cells. The native
Sage parent-set/directLFQ branch is recorded separately from fixed-dictionary
diagnostics. Scientific/visual review and manuscript integration remain pending.
'''
(RUN/'results/REPORT.md').write_text(report)
dump(RUN/'results/AGGREGATE_VALIDATION.json',dict(status='PASS',arms=5,acquisitions=58,receipts=receipts,job_id=os.environ.get('SLURM_JOB_ID'),scientific_review='PENDING',native_directlfq='SEPARATE_DEPENDENT_STAGE',outputs=[meta(p) for p in (RUN/'results').iterdir() if p.is_file() and p.name!='AGGREGATE_VALIDATION.json']))
print('HSTOOL MATCHED QUANT COMPARISON PASS',flush=True)
