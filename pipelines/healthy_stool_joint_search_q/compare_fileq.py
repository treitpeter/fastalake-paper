#!/usr/bin/env python3
"""D3 / CR-030 follow-up: completeness, identity and microbial protein-group precision with
the joint search gated by per-acquisition peptide confidence. Statistics copy A21_7
frozen/compare_revision2.py (completeness 'all' rows; pair CV ddof=0 median over rows then over
pairs; log2 Pearson; matched rows = positive under both arms in both acquisitions). The own-row
definition (positive in both acquisitions under the arm being scored) is added for Supp Fig 7c,d.
Usage: compare_fileq.py OUTDIR ARMDIR [ARMDIR...]   (an ARMDIR's basename is its label)"""
import itertools,json,sys
from pathlib import Path
import numpy as np,pandas as pd
A21=Path('LAB_ROOT/MM_review_analysis/A21_7_hstool_quant_comparison')
TIERS=['TechnicalReplicates','WorkflowReplicates','OnePiece','Intra','Interindividual']
def tier(s):return next(t for t in TIERS if t in s)
def corr(x,y):return float(np.corrcoef(x,y)[0,1]) if len(x)>2 and np.std(x)>0 and np.std(y)>0 else np.nan
def med(a):return float(np.median(a)) if len(a) else np.nan
def load(p):z=np.load(p);return pd.DataFrame(z['values'],index=z['keys'],columns=z['samples'])
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
dirs={Path(d).name:Path(d) for d in sys.argv[2:]}
m=json.loads((A21/'inputs/manifest_revision2.json').read_text());samples=m['samples']
for k,d in dirs.items():
    c=json.loads((d/'COMPLETE.json').read_text());assert c['status']=='PASS' and not c.get('pilot') and c['acquisitions']==58,k
ids={k:(lambda f:{s:set(f[f['sample'].eq(s)].peptide) for s in samples})(pd.read_csv(d/'accepted_identities.tsv.gz',sep='\t')) for k,d in dirs.items()}
# identities
rows=[];summ=[]
for k in dirs:
    for s in samples:
        a=ids['original'][s];b=ids[k][s]
        rows.append(dict(arm=k,sample=s,tier=tier(s),original=len(a),test=len(b),overlap=len(a&b),count_ratio=len(b)/len(a),original_identity_recovery=len(a&b)/len(a)))
idt=pd.DataFrame(rows);idt.to_csv(out/'identity_per_acquisition.tsv',sep='\t',index=False)
tech=[s for s in samples if tier(s)=='TechnicalReplicates']
for k in dirs:
    x=idt[idt.arm.eq(k)]
    summ.append(dict(arm=k,median_count_ratio=x.count_ratio.median(),n_higher=int((x.count_ratio>1).sum()),median_identity_recovery=x.original_identity_recovery.median(),
        study_distinct=len(set().union(*ids[k].values())),ids_in_all12_technical=len(set.intersection(*[ids[k][s] for s in tech]))))
# completeness
comp=[]
for space in ['features','fixed_groups','fixed_groups.Microbial']:
    for k,d in dirs.items():
        v=load(d/('same_file_id.'+space+'.npz')).reindex(columns=samples).fillna(0)
        for t in TIERS[:2]:
            ix=[s for s in samples if tier(s)==t];comp.append(dict(space=space,arm=k,tier=t,acquisitions=len(ix),rows_complete_all=int((v[ix]>0).all(axis=1).sum())))
comp=pd.DataFrame(comp);comp.to_csv(out/'completeness_all.tsv',sep='\t',index=False)
# precision, fixed_groups.Microbial
space='fixed_groups.Microbial';fr={k:load(d/('same_file_id.'+space+'.npz')) for k,d in dirs.items()}
vocab=sorted(set().union(*(set(f.index) for f in fr.values())));arr={k:f.reindex(index=vocab,columns=samples).fillna(0).to_numpy() for k,f in fr.items()}
pairs=[]
def pairstats(v,keep,i,j):
    vals=v[keep][:,[i,j]];logs=np.log2(vals);cv=100*np.std(vals,axis=1,ddof=0)/np.mean(vals,axis=1)
    return int(keep.sum()),corr(logs[:,0],logs[:,1]),med(cv)
refs=[('original',k) for k in dirs if k!='original']+([('original_fileq','joint_fileq')] if {'original_fileq','joint_fileq'}<=set(dirs) else [])
for t in TIERS[:2]:
    ix=[i for i,s in enumerate(samples) if tier(s)==t]
    for i,j in itertools.combinations(ix,2):
        for ref,test in refs:
            x,y=arr[ref],arr[test];keep=(x[:,i]>0)&(x[:,j]>0)&(y[:,i]>0)&(y[:,j]>0)
            for lab,v in [(ref,x),(test,y)]:
                n,r,c=pairstats(v,keep,i,j);pairs.append(dict(definition='matched',comparison=ref+'_vs_'+test,method=lab,tier=t,sample_a=samples[i],sample_b=samples[j],common_rows=n,log2_pearson=r,median_pair_cv_pct=c))
        for k,v in arr.items():
            keep=(v[:,i]>0)&(v[:,j]>0);n,r,c=pairstats(v,keep,i,j)
            pairs.append(dict(definition='own_rows',comparison=k,method=k,tier=t,sample_a=samples[i],sample_b=samples[j],common_rows=n,log2_pearson=r,median_pair_cv_pct=c))
pairs=pd.DataFrame(pairs);pairs.to_csv(out/'microbial_group_pairs.tsv',sep='\t',index=False)
ps=pairs.groupby(['definition','comparison','method','tier']).agg(pairs=('common_rows','size'),common_rows=('common_rows','median'),log2_pearson=('log2_pearson','median'),median_pair_cv_pct=('median_pair_cv_pct','median')).reset_index()
ps.to_csv(out/'microbial_group_pair_summary.tsv',sep='\t',index=False)
pd.DataFrame(summ).to_csv(out/'identity_summary.tsv',sep='\t',index=False)
res=dict(identity=summ,completeness=comp.to_dict('records'),precision=ps.to_dict('records'))
json.dump(res,open(out/'RESULTS.json','w'),indent=1,default=float)
print(pd.DataFrame(summ).to_string());print(comp.to_string());print(ps.to_string())
