"""CR-030: per-acquisition re-estimation of peptide-level target/decoy q inside the joint HSTOOL search.

For each acquisition we count canonical peptides (modifications stripped, I->L, 5-50 residues) from rank-1
target PSMs passing peptide q<=0.01 under:
  sage      : Sage's own peptide_q column (study-level in the joint search)
  glob_X    : our re-estimate over the whole invocation (should reproduce 'sage' if the estimator matches)
  file_X    : our re-estimate within each acquisition only (per-acquisition confidence)
X = 'plain' (q=D/T) or 'plus1' ((D+1)/T); best PSM score per (modified peptide, label) as in peptide-level FDR.
Own-database (acquisition-specific) searches are processed with the same code as a positive control.
Negative control: labels permuted within each acquisition -> per-file accepted counts should collapse.
Read-only on all inputs; writes only to --out.
"""
import argparse, json, re, sys, time, hashlib
import numpy as np, pandas as pd
from pathlib import Path
COLS=['peptide','filename','rank','label','sage_discriminant_score','peptide_q']
MOD=re.compile(r'\[[^\]]*\]')
def canon(p):
    s=MOD.sub('',p).replace('-','').replace('I','L'); return s
def qvals(score,label,plus1):
    o=np.argsort(-score,kind='mergesort'); lab=label[o]
    T=np.cumsum(lab==1); D=np.cumsum(lab!=1)
    q=((D+(1 if plus1 else 0))/np.maximum(T,1)).astype(float)
    q=np.minimum.accumulate(q[::-1])[::-1]
    out=np.empty_like(q); out[o]=q; return out
TOK=re.compile(r'[A-Z](?:\[[^\]]*\])?|\[[^\]]*\]-?')
def unrev(p):
    t=TOK.findall(p)
    if len(t)<3: return p
    return ''.join([t[0]]+t[1:-1][::-1]+[t[-1]])
PICKED=False
def pep_q(df,plus1):
    best=df.groupby(['peptide','label'],sort=False)['sage_discriminant_score'].max().reset_index()
    if PICKED:
        best['key']=np.where(best.label==1,best.peptide,best.peptide.map(unrev))
        win=best.sort_values('sage_discriminant_score',ascending=False).drop_duplicates('key')
        win=win.copy(); win['q']=qvals(win.sage_discriminant_score.values,win.label.values,plus1)
        best=best.merge(win[['peptide','label','q']],on=['peptide','label'],how='left')
        best['q']=best.q.fillna(1.0)
    else:
        best['q']=qvals(best.sage_discriminant_score.values,best.label.values,plus1)
    return df.merge(best[['peptide','label','q']],on=['peptide','label'],how='left')['q'].values
def count(df,mask):
    t=df[mask&(df.label==1)&(df['rank']==1)]
    c=t.peptide.map(canon); c=c[(c.str.len()>=5)&(c.str.len()<=50)]
    return c.nunique()
def summarize(df,tag,glob=True,negctl=False,seed=0):
    rows=[]
    if glob:
        df=df.copy(); df['q_glob_plain']=pep_q(df,False); df['q_glob_plus1']=pep_q(df,True)
    rng=np.random.default_rng(seed)
    for f,g in df.groupby('filename',sort=True):
        g=g.copy(); r={'arm':tag,'filename':f,'psms':len(g),'sage':count(g,g.peptide_q<=0.01)}
        if glob:
            r['glob_plain']=count(g,g.q_glob_plain<=0.01); r['glob_plus1']=count(g,g.q_glob_plus1<=0.01)
        g['qp']=pep_q(g,False); g['q1']=pep_q(g,True)
        r['file_plain']=count(g,g.qp<=0.01); r['file_plus1']=count(g,g.q1<=0.01)
        if negctl:
            h=g.copy(); h['label']=rng.permutation(h.label.values); h['q1']=pep_q(h,True)
            r['negctl_file_plus1']=count(h,h.q1<=0.01)
        rows.append(r)
    return rows
def read(path,nrows=None):
    return pd.read_csv(path,sep='\t',usecols=COLS,nrows=nrows,dtype={'peptide':str,'filename':str,'rank':np.int32,'label':np.int8,'sage_discriminant_score':np.float64,'peptide_q':np.float64})
ap=argparse.ArgumentParser(); ap.add_argument('--joint'); ap.add_argument('--manifest'); ap.add_argument('--out'); ap.add_argument('--nrows',type=int); ap.add_argument('--own-limit',type=int); ap.add_argument('--picked',action='store_true')
a=ap.parse_args(); PICKED=a.picked; out=Path(a.out); out.mkdir(parents=True,exist_ok=True); t0=time.time()
man=json.load(open(a.manifest)); samples=man['samples']
own_paths=[]
for arm in man.get('arms',man.get('inputs',[])) if isinstance(man.get('arms'),list) else []: pass
txt=open(a.manifest).read()
own_paths=sorted(set(re.findall(r'(LAB_ROOT/[^"]*/HSTOOL/5_sage/razor_hash-acc/[^"/]+/results\.sage\.tsv)',txt)))
print('own paths',len(own_paths),flush=True)
rows=[]
for p in own_paths[:a.own_limit] if a.own_limit else own_paths:
    d=read(p,a.nrows); rows+=summarize(d,'own',glob=False,negctl=True)
print('own done',time.time()-t0,flush=True)
j=read(a.joint,a.nrows); print('joint rows',len(j),'files',j.filename.nunique(),time.time()-t0,flush=True)
rows+=summarize(j,'joint',glob=True,negctl=True)
r=pd.DataFrame(rows); r.to_csv(out/'per_acquisition_counts.tsv',sep='\t',index=False)
r['stem']=r.filename.str.replace(r'\.mzML(\.gz)?$','',regex=True)
own=r[r.arm=='own'].set_index('stem'); jt=r[r.arm=='joint'].set_index('stem')
common=own.index.intersection(jt.index)
summ={'n_own':len(own),'n_joint':len(jt),'n_common':len(common),'elapsed_s':time.time()-t0}
for col in ['sage','glob_plain','glob_plus1','file_plain','file_plus1']:
    if col in jt and len(common):
        summ[f'median_joint_{col}_over_own_sage']=float((jt.loc[common,col]/own.loc[common,'sage']).median())
for col in ['file_plain','file_plus1']:
    if len(common): summ[f'median_joint_{col}_over_own_{col}']=float((jt.loc[common,col]/own.loc[common,col]).median())
if len(own): summ['own_file_plus1_equals_sage_frac']=float((own.file_plus1==own.sage).mean()); summ['own_file_plain_equals_sage_frac']=float((own.file_plain==own.sage).mean())
summ['own_max_abs_rel_diff_plus1']=float(((own.file_plus1-own.sage).abs()/own.sage.clip(lower=1)).max()) if len(own) else None
if len(jt):
    summ['joint_glob_plus1_equals_sage_frac']=float((jt.glob_plus1==jt.sage).mean()); summ['joint_glob_plain_equals_sage_frac']=float((jt.glob_plain==jt.sage).mean())
    summ['median_joint_file_plus1_over_joint_sage']=float((jt.file_plus1/jt.sage).median())
    summ['negctl_joint_max']=int(jt.negctl_file_plus1.max()); summ['negctl_own_max']=int(own.negctl_file_plus1.max()) if len(own) else None
json.dump(summ,open(out/'SUMMARY.json','w'),indent=2); print(json.dumps(summ,indent=2))
