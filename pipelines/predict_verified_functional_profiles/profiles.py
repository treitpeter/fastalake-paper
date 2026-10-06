"""Conserved, distinct-term reporting profiles on the frozen expanded reference."""
from pathlib import Path
from datetime import datetime,timezone
import csv,gzip,hashlib,json,os,re,sys
import numpy as np
import pandas as pd
from scipy import sparse
H=Path(__file__).resolve().parent
A=H.parent/'github_workflow_reconciliation_2026-09-30/resumed_audit/annotation_repair'
C=A/'annotated_origin_consensus_2026-10-01'
UNASSIGNED='Unassigned_at_this_resolution'
HUMAN='Human_assigned_signal_withheld'
LEVELS=['KO','OG','Pathway']

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()
def load(p):return json.loads(Path(p).read_text())
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def human(t):return t['rank']=='species' and t['supported_lineage'][-1:]==['s__Homo sapiens']
def eligible(f,level):
    if level=='KO':return sorted({t for t in f['terms'] if re.fullmatch(r'K\d{5}',t)})
    if level=='Pathway':return sorted({t for t in f['terms'] if re.fullmatch(r'map\d{5}',t) and not t.startswith('map011')})
    # Each retained parent annotation set supplies one root OG; repeated term support
    # never gives that term additional weight. Do not silently admit multi-OG sets.
    return sorted({s[0] for s in f['annotation_sets'] if len(s)==1 and re.fullmatch(r'[^@]+@(?:2\|Bacteria|2157\|Archaea)',s[0])})

def projector(recs,level,exclude_human):
    assignments=[]
    for r in recs:
        ts=eligible(r['function'],level)
        assignments.append([HUMAN] if exclude_human and human(r['taxonomy']) else (ts or [UNASSIGNED]))
    labels=sorted({t for ts in assignments for t in ts}-{UNASSIGNED,HUMAN})+[UNASSIGNED]+([HUMAN] if exclude_human else [])
    ix={t:i for i,t in enumerate(labels)};rr=[];cc=[];vv=[]
    for j,ts in enumerate(assignments):
        for t in ts:rr.append(ix[t]);cc.append(j);vv.append(1/len(ts))
    p=sparse.coo_matrix((vv,(rr,cc)),shape=(len(labels),len(recs))).tocsr()
    assert np.allclose(p.sum(0).A.ravel(),1,rtol=0,atol=1e-12)
    return labels,p

def tests():
    tax=dict(rank='unresolved',supported_lineage=[])
    def r(ts,missing=False,ann=None):return dict(function=dict(terms=ts,missing=missing,annotation_sets=ann if ann is not None else [[t] for t in ts]),taxonomy=tax)
    rr=[r(['K00001','K00001','K00002'],True),r([]),r(['K00001'])]
    names,p=projector(rr,'KO',False);v=p @ np.array([6.,5.,7.])
    assert dict(zip(names,v))=={'K00001':10.,'K00002':3.,UNASSIGNED:5.}
    rr[0]['taxonomy']=dict(rank='species',supported_lineage=['s__Homo sapiens'])
    names,p=projector(rr,'KO',True);v=p @ np.array([6.,5.,7.])
    assert dict(zip(names,v))=={'K00001':7.,UNASSIGNED:5.,HUMAN:6.}
    assert eligible(r(['map01100','map00010','map00010','bad'])['function'],'Pathway')==['map00010']
    f=r(['A@2|Bacteria','B@2|Bacteria'],ann=[['A@2|Bacteria','B@2|Bacteria']])['function']
    assert eligible(f,'OG')==[]
    f['annotation_sets'].append(['A@2|Bacteria']);assert eligible(f,'OG')==['A@2|Bacteria']

def verify():
    for name,value in load(H/'PROFILE_CODE.json').items():assert sha(H/name)==value,name
    cfg=load(H/'PROFILE_INPUTS.json')
    for path,value in cfg['receipts'].items():assert sha(path)==value,path
    return cfg

def run_arm(arm,dest):
    d=C/'full'/arm;source=load(d/'COMPLETE.json');assert source['status']=='PASS'
    used=['sample_totals.npz','quantity_definitions.jsonl.gz','taxonomy_rows.json','taxonomy_evidence_signal.npz']
    used += [f'{l}_{s}' for l in LEVELS for s in ['quantity_ids.json','evidence_signal.npz']]
    for name in used:assert sha(d/name)==source['outputs'][name]['sha256'],name
    dest.mkdir(parents=True,exist_ok=False)
    z=np.load(d/'sample_totals.npz');total=z['intensity'];samples=list(z['samples'])
    defs={}
    with gzip.open(d/'quantity_definitions.jsonl.gz','rt') as f:
        for r in map(json.loads,f):defs[r['quantity_id']]=r
    tx=load(d/'taxonomy_rows.json');tm=sparse.load_npz(d/'taxonomy_evidence_signal.npz')
    hmask=np.array([human(t) for t in tx]);hmass=tm[hmask].sum(0).A.ravel()
    assert np.allclose(tm.sum(0).A.ravel(),total,rtol=1e-9,atol=1e-5)
    taxclasses=[]
    for t in tx:
        lin=t['supported_lineage'];mic=lin[:1] in [['d__Bacteria'],['d__Archaea']]
        label='Human-assigned' if human(t) else ('Microbial species' if mic and t['rank']=='species' else ('Microbial above species' if mic else ('Other annotated origins' if lin else 'Unresolved origin')))
        taxclasses.append(label)
    taxout=pd.DataFrame(dict(sample=samples,total_LFQ=total))
    for cls in sorted(set(taxclasses)):taxout[cls]=tm[np.array([t==cls for t in taxclasses])].sum(0).A.ravel()
    assert np.allclose(taxout.drop(columns=['sample','total_LFQ']).sum(1),total,rtol=1e-9,atol=1e-5)
    taxout.to_csv(dest/'taxonomy_signal.tsv',sep='\t',index=False)
    cov=[];independent=[]
    for lev in LEVELS:
        ids=load(d/f'{lev}_quantity_ids.json');recs=[defs[k] for k in ids];mat=sparse.load_npz(d/f'{lev}_evidence_signal.npz')
        assert np.allclose(mat.sum(0).A.ravel(),total,rtol=1e-9,atol=1e-5)
        eligible_sets=[eligible(r['function'],lev) for r in recs]
        rowhuman=np.array([human(r['taxonomy']) for r in recs])
        assert np.allclose(mat[rowhuman].sum(0).A.ravel(),hmass,rtol=1e-9,atol=1e-5)
        for exclude,name in [(False,'all'),(True,'other_unresolved')]:
            labels,p=projector(recs,lev,exclude);result=(p@mat).tocsr()
            assert np.allclose(result.sum(0).A.ravel(),total,rtol=1e-9,atol=1e-5)
            keep=~rowhuman if exclude else np.ones(len(recs),bool)
            single=mat[keep & np.array([len(ts)==1 for ts in eligible_sets])].sum(0).A.ravel()
            multi=mat[keep & np.array([len(ts)>1 for ts in eligible_sets])].sum(0).A.ravel()
            no=mat[keep & np.array([not ts for ts in eligible_sets])].sum(0).A.ravel()
            partial=mat[keep & np.array([bool(ts) and r['function']['missing'] for ts,r in zip(eligible_sets,recs)])].sum(0).A.ravel()
            withheld=hmass if exclude else np.zeros_like(hmass)
            assert np.allclose(single+multi+no+withheld,total,rtol=1e-9,atol=1e-5)
            assert np.allclose(result[labels.index(UNASSIGNED)].toarray().ravel(),no,rtol=1e-9,atol=1e-5)
            mask=np.array([t not in [UNASSIGNED,HUMAN] for t in labels]);named=result[mask]
            assert np.allclose(named.sum(0).A.ravel(),single+multi,rtol=1e-9,atol=1e-5)
            stem=f'{lev}_{name}';sparse.save_npz(dest/(stem+'.npz'),named.T.tocsr())
            dump(dest/(stem+'_axes.json'),dict(samples=samples,features=list(np.array(labels)[mask]),representation='equal_share_distinct_eligible_terms',reference='searched-target union plus local full UHGG',human_exclusion=exclude))
            pd.DataFrame(named.toarray(),index=np.array(labels)[mask],columns=samples).to_csv(dest/(stem+'.tsv.gz'),sep='\t',index_label='term')
            for i,s in enumerate(samples):cov.append(dict(arm=arm,level=lev,profile=name,sample=s,total_LFQ=float(total[i]),single_eligible_term=float(single[i]),multiple_eligible_terms=float(multi[i]),no_eligible_term=float(no[i]),human_withheld=float(withheld[i]),positive_terms_with_missing_annotations=float(partial[i])))
            if lev=='KO' and not exclude:
                old=A/'functional_unassigned_review_2026-10-01/full'/arm
                old_receipt=load(old/'COMPLETE.json');assert old_receipt['status']=='PASS'
                for old_name in ['axes.json','proposed_equal_share_KO_signal.npz']:
                    assert sha(old/old_name)==old_receipt['outputs'][old_name],old_name
                oldax=load(old/'axes.json');oldm=sparse.load_npz(old/'proposed_equal_share_KO_signal.npz')
                assert oldax['samples']==samples
                oldlabels=[UNASSIGNED if t=='Unassigned_at_KO_resolution' else t for t in oldax['terms']]
                assert len(oldlabels)==len(set(oldlabels)) and set(oldlabels)==set(labels)
                order=[oldlabels.index(t) for t in labels];delta=result-oldm[order]
                assert np.allclose(delta.data,0,rtol=0,atol=max(1e-5,float(total.max())*1e-12))
                independent.append('New KO all profile reproduces previously verified broader-reference equal-share diagnostic')
    pd.DataFrame(cov).to_csv(dest/'coverage.tsv',sep='\t',index=False)
    dump(dest/'COMPLETE.json',dict(status='PASS',arm=arm,acquisitions=len(samples),accepted_features=int(z['features'].sum()),
        measured_intensity_conserved=True,eligibility_counterexamples_passed=True,independent_checks=independent,
        source_receipt_sha256=sha(d/'COMPLETE.json'),code_manifest_sha256=sha(H/'PROFILE_CODE.json'),
        independent_KO_receipt_sha256=sha(A/'functional_unassigned_review_2026-10-01/full'/arm/'COMPLETE.json'),
        input_manifest_sha256=sha(H/'PROFILE_INPUTS.json'),job_id=os.environ.get('SLURM_JOB_ID'),
        outputs={p.name:sha(p) for p in dest.iterdir() if p.is_file()}))
    print('PASS',arm,len(samples),flush=True)

def main(mode):
    cfg=verify();tests();out=H/('profile_pilot' if mode=='pilot' else 'profiles')
    if mode=='full':
        pilot=load(H/'profile_pilot/COMPLETE.json');assert pilot['status']=='PASS'
        assert pilot['code_manifest_sha256']==sha(H/'PROFILE_CODE.json')
        assert pilot['input_manifest_sha256']==sha(H/'PROFILE_INPUTS.json')
    out.mkdir(exist_ok=False)
    arms=['CAMPI_REF'] if mode=='pilot' else cfg['arms']
    for arm in arms:run_arm(arm,out/arm)
    dump(out/'COMPLETE.json',dict(status='PASS',mode=mode,arms=arms,utc=datetime.now(timezone.utc).isoformat(),
        job_id=os.environ.get('SLURM_JOB_ID'),code_manifest_sha256=sha(H/'PROFILE_CODE.json'),input_manifest_sha256=sha(H/'PROFILE_INPUTS.json'),
        arm_receipts={arm:sha(out/arm/'COMPLETE.json') for arm in arms},
        scope='Reporting refit on frozen acquisition sets; CapScan biological-label finalization is separate'))

if __name__=='__main__':main(sys.argv[1])
