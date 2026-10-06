#!/usr/bin/env python3
"""D3 / CR-030 follow-up: re-extract HSTOOL LFQ features gated by accepted same-acquisition
peptide identifications, where acceptance is defined either by Sage's own peptide_q
('sageq', positive control: must reproduce A21_7 exactly) or by peptide q re-estimated
within each acquisition ('fileq', per-acquisition confidence; identical estimator to the
iteration-1 CR-030 recount: best sage_discriminant_score per (modified peptide, label)
within the acquisition, q = D/T, unpicked, rank-1 targets, q <= 0.01).

The quantitative part copies A21_7 frozen/search_revision2.py (same LFQ table, same
positive-finite rule, same I/L feature keys, same fixed study dictionary, same local
Human/Microbial/Mixed/Unknown classification). The native Sage parent-set branch is not
needed here and is omitted. All inputs are read-only; writes only under --out.

Usage: extract_fileq.py ARM MODE OUT [--pilot] [--nrows N]
  ARM  : manifest arm name (original | joint_corrected | ...)
  MODE : sageq | fileq
"""
import argparse,csv,gzip,json,os,sys,time
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
A21=Path('LAB_ROOT/MM_review_analysis/A21_7_hstool_quant_comparison')
sys.dont_write_bytecode=True
sys.path.insert(0,str(A21/'frozen'))
from lib import dump,meta,verify,table
from peptide import canonical
from parameters import parameters
from origin import classify_members
from digest import fasta

def key(p,c):return p.replace('I','L')+'|charge='+str(c)
def protkey(p):return ';'.join(sorted(set(p.split(';'))))

def qvals(score,label):
    o=np.argsort(-score,kind='mergesort');lab=label[o]
    T=np.cumsum(lab==1);D=np.cumsum(lab!=1)
    q=(D/np.maximum(T,1)).astype(float);q=np.minimum.accumulate(q[::-1])[::-1]
    out=np.empty_like(q);out[o]=q;return out
def pep_q(df):
    best=df.groupby(['peptide','label'],sort=False)['sage_discriminant_score'].max().reset_index()
    best['q']=qvals(best.sage_discriminant_score.values,best.label.values)
    return df.merge(best[['peptide','label','q']],on=['peptide','label'],how='left')['q'].values

COLS=['peptide','filename','rank','label','sage_discriminant_score','peptide_q']
def evidence(path,names,mode,nrows=None,seed=0):
    ids={n:set() for n in names};mods={n:set() for n in names};stats={}
    d=pd.read_csv(path,sep='\t',usecols=COLS,nrows=nrows,dtype={'peptide':str,'filename':str,'rank':np.int32,'label':np.int8,'sage_discriminant_score':np.float64,'peptide_q':np.float64})
    d['filename']=d.filename.map(lambda x:Path(x).name)
    assert set(d.filename.unique())<=set(names),'unexpected filename'
    rng=np.random.default_rng(seed)
    for n,g in d.groupby('filename',sort=True):
        if mode=='sageq':q=g.peptide_q.values
        else:q=pep_q(g)
        ok=(q<=.01)&(g['rank'].values==1)&(g.label.values==1)
        peps=g.peptide.values[ok]
        can=[canonical(p) for p in peps]
        keep=[5<=len(c)<=50 for c in can]
        ids[n]={c for c,k in zip(can,keep) if k};mods[n]={p for p,k in zip(peps,keep) if k}
        st=dict(psm_rows=len(g),accepted_rank1_target_psms=int(ok.sum()),dropped_length=int(len(keep)-sum(keep)))
        if mode=='fileq':
            h=g.copy();h['label']=rng.permutation(h.label.values);hq=pep_q(h)
            st['negctl_permuted_label_accepted']=len({canonical(p) for p in h.peptide.values[(hq<=.01)&(h['rank'].values==1)&(h.label.values==1)]})
            st['sage_q_accepted']=len({canonical(p) for p in g.peptide.values[(g.peptide_q.values<=.01)&(g['rank'].values==1)&(g.label.values==1)]})
        stats[n]=st
    return ids,mods,stats

def main():
    ap=argparse.ArgumentParser();ap.add_argument('arm');ap.add_argument('mode',choices=['sageq','fileq']);ap.add_argument('out')
    ap.add_argument('--pilot',action='store_true');ap.add_argument('--nrows',type=int)
    a=ap.parse_args();t0=time.time()
    m=json.loads((A21/'inputs/manifest_revision2.json').read_text())
    for x in m['frozen_code']+[m['protocol']]:verify(x)
    t=next(x for x in m['arms'] if x['name']==a.arm)
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    samples=m['samples'];selected=samples[:1] if a.pilot else samples;selected_names=[s+'.mzML' for s in selected]
    mapping=pd.read_csv(m['mapping']['path'],sep='\t',dtype=str).set_index('canonical_peptide').study_group.to_dict();verify(m['mapping'])
    aliases=defaultdict(set)
    with gzip.open(m['union_provenance']['path'],'rt') as fh:
        for r in csv.DictReader(fh,delimiter='\t'):aliases[r['union_id']].add(r['original_accession'])
    verify(m['union_provenance'])
    class_table=pd.read_csv(m['local_classification']['path'],sep='\t',usecols=['protein_group','classification'])
    classes={protkey(p):c for p,c in zip(class_table.protein_group,class_table.classification)};verify(m['local_classification'])
    origin_mats={c:{} for c in ['Human','Microbial','Mixed','Unknown']}
    records=[r for r in t['records'] if any(Path(p).name in selected_names for p in json.loads(Path(r['config']['path']).read_text())['mzml_paths'])]
    idmat=defaultdict(dict);allmat=defaultdict(dict);ids_by_sample={};coverage=[];inputs=[]
    for rec in records:
        for x in rec['inputs']:verify(x)
        config=json.loads(Path(rec['config']['path']).read_text());all_names=[Path(p).name for p in config['mzml_paths']]
        assert parameters(config)==m['search_parameters']
        actual_names=[n for n in all_names if n in selected_names]
        ids,modified,stats=evidence(rec['psm']['path'],all_names,a.mode,a.nrows)
        cols=['peptide','charge','proteins','q_value']+all_names
        f=pd.read_csv(rec['lfq']['path'],sep='\t',usecols=cols);assert f.charge.eq(-1).all()
        f['key']=[key(p,c) for p,c in zip(f.peptide,f.charge)];f['canonical']=f.peptide.map(canonical);f['native_group']=f.proteins.map(protkey)
        fa_audit={};fa={h.split()[0]:seq for h,seq in fasta(rec['fasta']['path'],fa_audit)};assert fa_audit['sha256']==rec['fasta']['sha256']
        group_classes={}
        for peptide,group in set(zip(f['canonical'],f['native_group'])):
            members=group.split(';');assert all(pid in fa and peptide in fa[pid] for pid in members)
            expanded=set().union(*(aliases[pid] if pid.startswith('FLU_') else {pid} for pid in members))
            group_classes[group]=classes.get(protkey(';'.join(expanded)),classify_members(expanded))
        rowclasses=f.native_group.map(group_classes).to_numpy()
        raw=f[actual_names].to_numpy(float);positive=np.isfinite(raw)&(raw>0);v=np.where(positive,raw,0.)
        gate=np.stack([f.canonical.isin(ids[n]).to_numpy() for n in actual_names],axis=1)
        for mat,mask in [(idmat,gate),(allmat,np.ones_like(gate))]:
            frame=pd.DataFrame(np.where(mask,v,0.),index=f['key'],columns=[Path(n).stem for n in actual_names]).groupby(level=0).sum()
            for s in frame:mat[s]=frame[s]
        for cls,store in origin_mats.items():
            vals=np.where(gate&(rowclasses==cls)[:,None],v,0.)
            of=pd.DataFrame(vals,index=f['key'],columns=[Path(n).stem for n in actual_names]).groupby(level=0).sum()
            for s in of:store[s]=of[s]
        for j,n in enumerate(actual_names):
            s=Path(n).stem;ids_by_sample[s]=ids[n]
            coverage.append(dict(sample=s,arm=a.arm,mode=a.mode,identified_canonical_peptides=len(ids[n]),positive_source_cells=int(positive[:,j].sum()),
                same_file_id_source_cells=int((positive[:,j]&gate[:,j]).sum()),positive_without_same_file_id=int((positive[:,j]&~gate[:,j]).sum()),**stats[n]))
        inputs.extend(rec['inputs'])
    assert set(idmat)==set(selected)
    proj={}
    for name,mat in [('same_file_id',idmat),('all_positive',allmat)]:
        frame=pd.DataFrame(mat).reindex(columns=selected).fillna(0.);frame=frame[frame.gt(0).any(axis=1)]
        np.savez_compressed(out/(name+'.features.npz'),keys=frame.index.to_numpy(str),samples=np.array(selected),values=frame.to_numpy())
        groups=pd.Series([mapping.get(canonical(k.rsplit('|charge=',1)[0])) for k in frame.index],index=frame.index)
        known=groups.notna();g=frame.loc[known].groupby(groups[known]).sum()
        np.savez_compressed(out/(name+'.fixed_groups.npz'),keys=g.index.to_numpy(str),samples=np.array(selected),values=g.to_numpy())
        proj[name]=dict(features=len(frame),mapped_features=int(known.sum()),groups=len(g))
    for cls,mat in origin_mats.items():
        frame=pd.DataFrame(mat).reindex(columns=selected).fillna(0.)
        groups=pd.Series([mapping.get(canonical(k.rsplit('|charge=',1)[0])) for k in frame.index],index=frame.index)
        known=groups.notna();g=frame.loc[known].groupby(groups[known]).sum();g=g[g.gt(0).any(axis=1)]
        np.savez_compressed(out/('same_file_id.fixed_groups.'+cls+'.npz'),keys=g.index.to_numpy(str),samples=np.array(selected),values=g.to_numpy())
    with gzip.open(out/'accepted_identities.tsv.gz','wt') as fh:
        w=csv.writer(fh,delimiter='\t',lineterminator='\n');w.writerow(['sample','peptide'])
        for s,ps in ids_by_sample.items():w.writerows((s,p) for p in sorted(ps))
    table(out/'coverage.tsv',coverage)
    dump(out/'COMPLETE.json',dict(status='PASS',arm=a.arm,mode=a.mode,pilot=a.pilot,nrows=a.nrows,acquisitions=len(selected),projection=proj,inputs=inputs,
        job_id=os.environ.get('SLURM_JOB_ID'),elapsed_s=time.time()-t0,outputs=[meta(p) for p in out.iterdir() if p.is_file() and p.name!='COMPLETE.json']))
    print('EXTRACT PASS',a.arm,a.mode,len(selected),proj,round(time.time()-t0),flush=True)
if __name__=='__main__':main()
