"""Re-derive peptide function evidence, quantity matrices and reporting profiles for one arm after the
native MGnify UHGG eggNOG rows are filtered by the non-circular seed-alignment check.
variant 'replay'      : keep every native row (the original rule); must reproduce the frozen outputs exactly.
variant 'verifiedNN'  : keep a native row only if its own seed ortholog aligns to the attached UHGG sequence
                        with DIAMOND bitscore >= NN/100 x the bitscore MGnify reported for that row.
Taxonomy, parent graph, peptides, acquisition roster and LFQ ledger are frozen and re-checked."""
from pathlib import Path
from collections import defaultdict,Counter
import csv,gzip,hashlib,json,os,sys,importlib.util,glob
import numpy as np, pandas as pd
from scipy import sparse
csv.field_size_limit(1<<30)
W=Path(__file__).resolve().parents[1]
RA=(Path(__import__('os').environ['FASTALAKE_DATA_DIR'])/'review/v15_integration/integration_2026-09-29/analysis_followup/github_workflow_reconciliation_2026-09-30/resumed_audit')
E=RA/'evidence_all_cohorts';N=RA/'annotation_repair';REP=N/'evidence_repaired';CONS=N/'annotated_origin_consensus_2026-10-01'
CP=(Path(__import__('os').environ['FASTALAKE_DATA_DIR'])/'review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01')
sys.path.insert(0,str(E))
from run_projection import records,annotation_index,digest
from evidence_rules import LEVELS,PeptideEvidence,lineage,parse_terms,signature,valid_intensity
spec=importlib.util.spec_from_file_location('cons_policy',CONS/'policy.py');pol=importlib.util.module_from_spec(spec);spec.loader.exec_module(pol)
spec=importlib.util.spec_from_file_location('frozen_profiles',CP/'profiles.py');prof=importlib.util.module_from_spec(spec);spec.loader.exec_module(prof)
arm,variant=sys.argv[1],sys.argv[2]
INP=json.loads((E/'INPUTS.json').read_text());cfg=next(t for t in INP['tasks'] if t['id']==arm)
out=W/'reproject'/variant/arm;out.mkdir(parents=True,exist_ok=False)
pout=W/'profiles'/variant/arm
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,sort_keys=True,default=str)+'\n')
# frozen receipts
for d in [REP/'full'/arm,CONS/'full'/arm]:
    rc=json.loads((d/'COMPLETE.json').read_text());assert rc['status']=='PASS'
    for name in ['peptide_evidence.tsv.gz','sample_totals.npz']:assert digest(d/name)==rc['outputs'][name]['sha256'],(d,name)
rc=json.loads((REP/'full'/arm/'COMPLETE.json').read_text())
for name in ['parent_evidence.tsv.gz','peptide_parent_edges.tsv.gz']:assert digest(REP/'full'/arm/name)==rc['outputs'][name]['sha256'],name
# ---- verification table
thr=None
if variant!='replay':
    assert variant.startswith('verified');thr=int(variant[8:])/100
    vt={}
    for p in sorted(glob.glob(str(W/'verify/prep_real_out/shard_*.tsv.gz'))):
        for r in records(p):
            ns=float(r['native_score']);b=float(r['pair_bitscore']) if r['pair_bitscore'] else 0.
            vt[r['accession']]=(r['seed_ortholog'],b>=thr*ns)
    assert len(glob.glob(str(W/'verify/prep_real_out/shard_*.COMPLETE.json')))==64
print('VERIFY TABLE',None if thr is None else len(vt),flush=True)
# ---- pass 1: parents that received native rows under the original rule
need={}
for r in records(REP/'full'/arm/'parent_evidence.tsv.gz'):
    if r['native_eggnog_status']=='gene_sequence_verified':need[r['accession']]=r['sequence_sha256_IL']
bygenome=defaultdict(list)
for a in need:bygenome[a.split('_')[0]].append(a)
intern={};native_by_sha=defaultdict(set);ncount=Counter()
for g,accs in bygenome.items():
    want=set(accs);seen=set()
    for row in records(INP['native_eggnog'][g]):
        a=row['#query']
        if a not in want or a in seen:continue
        seen.add(a)
        vals=tuple(parse_terms(row.get({'KO':'KEGG_ko','OG':'eggNOG_OGs','Pathway':'KEGG_Pathway'}[l],''),l) for l in LEVELS)
        vals=intern.setdefault(vals,vals)
        if thr is not None:
            seed,ok=vt[a];assert seed==row['seed_ortholog'],(a,seed,row['seed_ortholog'])
            ncount['native_rows_pass' if ok else 'native_rows_fail']+=1
            if not ok:continue
        native_by_sha[need[a]].add(vals)
    assert seen==want,(g,len(want-seen))
ncount['native_accessions']=len(need)
print('NATIVE',dict(ncount),flush=True)
del need,bygenome
mapper,_=annotation_index(cfg['annotations'])
adapter={}
with gzip.open(N/'adapter/functions.jsonl.gz','rt') as f:
    for r in map(json.loads,f):adapter[r['sequence_sha256_IL']]=r['annotations']
ann_ids={};ann_list=[];parent_info={};pchange=Counter()
def norm(d):return {k:frozenset(v) for k,v in d.items() if v}
for r in records(REP/'full'/arm/'parent_evidence.tsv.gz'):
    sh=r['sequence_sha256_IL'];ann={l:set() for l in LEVELS}
    m=mapper.get(bytes.fromhex(sh))
    if m:
        for l in LEVELS:ann[l]|=m.get(l,set())
    for vals in native_by_sha.get(sh,()):
        for i,l in enumerate(LEVELS):ann[l].add(vals[i])
    for l,vs in adapter.get(sh,{}).items():ann[l]|={tuple(v) for v in vs}
    stored=norm({k:{tuple(v) for v in vs} for k,vs in json.loads(r['annotation_sets_json']).items()})
    new=norm(ann)
    if variant=='replay':assert new==stored,(r['parent_id'],new,stored)
    else:
        for l in LEVELS:
            if new.get(l,frozenset())!=stored.get(l,frozenset()):pchange['parents_changed_'+l]+=1
        if new!=stored:pchange['parents_changed_any']+=1
    key=tuple(new.get(l,frozenset()) for l in LEVELS)
    if key not in ann_ids:ann_ids[key]=len(ann_list);ann_list.append({l:set(key[i]) for i,l in enumerate(LEVELS)})
    parent_info[bytes.fromhex(r['parent_id'])]=(lineage(r['accepted_lineage']),ann_ids[key]);pchange['parents']+=1
print('PARENTS',dict(pchange),'unique annotation sets',len(ann_list),flush=True)
del mapper,adapter,native_by_sha,ann_ids
# ---- peptides (frozen consensus order, taxonomy and stored function evidence)
pep=[];stored_tax=[];stored_fun=[];npar=[]
for r in records(CONS/'full'/arm/'peptide_evidence.tsv.gz'):
    assert int(r['peptide_id'])==len(pep);pep.append(r['canonical_peptide']);stored_tax.append(r['taxonomy_json']);stored_fun.append(r['function_evidence_json']);npar.append(int(r['compatible_parents']))
ids={p:i for i,p in enumerate(pep)};evidence=[PeptideEvidence() for _ in pep]
for e in records(REP/'full'/arm/'peptide_parent_edges.tsv.gz'):
    tax,aid=parent_info[bytes.fromhex(e['parent_id'])];evidence[int(e['peptide_id'])].add(tax,ann_list[aid])
del parent_info
catalog={};lookup={l:[] for l in LEVELS};tax_catalog={};tax_ids=[];changed=np.zeros((len(pep),len(LEVELS)),bool);fchange=Counter()
with gzip.open(out/'peptide_evidence.tsv.gz','wt',compresslevel=3) as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['peptide_id','canonical_peptide','compatible_parents','taxonomy_json','function_evidence_json'])
    for i,e in enumerate(evidence):
        assert e.parents==npar[i]
        tax=pol.taxonomy(e.taxonomy(),e.parents);assert json.dumps(tax,sort_keys=True)==stored_tax[i],i
        fun={l:e.functions[l].result() for l in LEVELS};raw=json.dumps(fun,sort_keys=True)
        if variant=='replay':assert raw==stored_fun[i],i
        else:
            old=json.loads(stored_fun[i])
            for k,l in enumerate(LEVELS):
                # compare in JSON space (tuples vs lists); diagnostic only, never feeds the matrices
                if old[l]!=json.loads(json.dumps(fun[l])):changed[i,k]=True;fchange['peptides_changed_'+l]+=1
                if old[l]['terms']!=fun[l]['terms']:fchange['peptides_terms_changed_'+l]+=1
        w.writerow([i,pep[i],e.parents,json.dumps(tax,sort_keys=True),raw])
        tx=json.dumps({k:tax[k] for k in ('state','rank','supported_lineage')},sort_keys=True);tax_ids.append(tax_catalog.setdefault(tx,len(tax_catalog)))
        for l in LEVELS:
            key,rec=signature(l,fun[l],tax);catalog.setdefault(key,rec);lookup[l].append(key)
del evidence,stored_fun,stored_tax
print('PEPTIDES',len(pep),dict(fchange),flush=True)
summary=list(records(cfg['acquisition_summary']));samples=[r['sample'] for r in summary];sx={s:i for i,s in enumerate(samples)}
totals=np.zeros(len(samples));counts=np.zeros(len(samples),np.int64);chg=np.zeros((len(LEVELS),len(samples)))
matrices={l:defaultdict(float) for l in LEVELS};taxmatrix=defaultdict(float)
for r in records(cfg['feature_ledger']):
    si=sx[r['sample']];v=valid_intensity(r['intensity']);pi=ids[r['canonical_peptide'].upper().replace('I','L')]
    totals[si]+=v;counts[si]+=1;chg[:,si]+=changed[pi]*v
    for l in LEVELS:matrices[l][lookup[l][pi],si]+=v
    taxmatrix[tax_ids[pi],si]+=v
prior=np.load(CONS/'full'/arm/'sample_totals.npz');assert list(prior['samples'])==samples
assert np.array_equal(prior['features'],counts) and np.allclose(prior['intensity'],totals,rtol=1e-9,atol=1e-6)
np.savez_compressed(out/'sample_totals.npz',samples=np.array(samples),intensity=totals,features=counts,source_intensity=totals,source_features=counts)
with gzip.open(out/'quantity_definitions.jsonl.gz','wt',compresslevel=3) as f:
    for key in sorted(catalog):f.write(json.dumps(dict(quantity_id=key,**catalog[key]),sort_keys=True)+'\n')
for l,d in matrices.items():
    names=sorted({k for k,s in d});ix={k:i for i,k in enumerate(names)}
    mat=sparse.coo_matrix(([v for v in d.values()],([ix[k] for k,s in d],[s for k,s in d])),shape=(len(names),len(samples))).tocsr()
    assert np.allclose(mat.sum(0).A.ravel(),totals,rtol=1e-9,atol=1e-6)
    sparse.save_npz(out/f'{l}_evidence_signal.npz',mat);(out/f'{l}_quantity_ids.json').write_text(json.dumps(names)+'\n')
    if variant=='replay':
        assert json.loads((CONS/'full'/arm/f'{l}_quantity_ids.json').read_text())==names,l
        old=sparse.load_npz(CONS/'full'/arm/f'{l}_evidence_signal.npz');assert abs(old-mat).max()<=1e-6*max(1.,totals.max()*1e-9),l
tx=sorted(tax_catalog,key=tax_catalog.get);tm=sparse.coo_matrix(([v for v in taxmatrix.values()],([k for k,s in taxmatrix],[s for k,s in taxmatrix])),shape=(len(tx),len(samples))).tocsr()
sparse.save_npz(out/'taxonomy_evidence_signal.npz',tm);(out/'taxonomy_rows.json').write_text(json.dumps([json.loads(t) for t in tx],indent=2)+'\n')
pd.DataFrame(dict(sample=samples,total_LFQ=totals,**{f'LFQ_peptides_function_evidence_changed_{l}':chg[k] for k,l in enumerate(LEVELS)})).to_csv(out/'changed_signal.tsv',sep='\t',index=False)
dump(out/'COMPLETE.json',dict(status='PASS',arm=arm,variant=variant,threshold=thr,job=os.environ.get('SLURM_JOB_ID'),native=dict(ncount),parents=dict(pchange),peptides=dict(fchange),peptide_count=len(pep),
    replay_exact=variant=='replay',outputs={p.name:dict(sha256=digest(p)) for p in out.iterdir() if p.is_file() and p.name!='COMPLETE.json'}))
# ---- reporting profiles (frozen profiles.py functions, same arithmetic as profiles.run_arm)
prof.tests();pout.mkdir(parents=True,exist_ok=False);d=out
defs={}
with gzip.open(d/'quantity_definitions.jsonl.gz','rt') as f:
    for r in map(json.loads,f):defs[r['quantity_id']]=r
txr=json.loads((d/'taxonomy_rows.json').read_text());hmask=np.array([prof.human(t) for t in txr]);hmass=tm[hmask].sum(0).A.ravel();total=totals
taxclasses=[]
for t in txr:
    lin=t['supported_lineage'];mic=lin[:1] in [['d__Bacteria'],['d__Archaea']]
    taxclasses.append('Human-assigned' if prof.human(t) else ('Microbial species' if mic and t['rank']=='species' else ('Microbial above species' if mic else ('Other annotated origins' if lin else 'Unresolved origin'))))
taxout=pd.DataFrame(dict(sample=samples,total_LFQ=total))
for cls in sorted(set(taxclasses)):taxout[cls]=tm[np.array([t==cls for t in taxclasses])].sum(0).A.ravel()
taxout.to_csv(pout/'taxonomy_signal.tsv',sep='\t',index=False)
cov=[];checks={}
for lev in LEVELS:
    idsq=json.loads((d/f'{lev}_quantity_ids.json').read_text());recs=[defs[k] for k in idsq];mat=sparse.load_npz(d/f'{lev}_evidence_signal.npz')
    eligible_sets=[prof.eligible(r['function'],lev) for r in recs];rowhuman=np.array([prof.human(r['taxonomy']) for r in recs])
    for exclude,name in [(False,'all'),(True,'other_unresolved')]:
        labels,p=prof.projector(recs,lev,exclude);result=(p@mat).tocsr()
        assert np.allclose(result.sum(0).A.ravel(),total,rtol=1e-9,atol=1e-5)
        keep=~rowhuman if exclude else np.ones(len(recs),bool)
        single=mat[keep & np.array([len(ts)==1 for ts in eligible_sets])].sum(0).A.ravel()
        multi=mat[keep & np.array([len(ts)>1 for ts in eligible_sets])].sum(0).A.ravel()
        no=mat[keep & np.array([not ts for ts in eligible_sets])].sum(0).A.ravel()
        partial=mat[keep & np.array([bool(ts) and r['function']['missing'] for ts,r in zip(eligible_sets,recs)])].sum(0).A.ravel()
        withheld=hmass if exclude else np.zeros_like(hmass)
        mask=np.array([t not in [prof.UNASSIGNED,prof.HUMAN] for t in labels]);named=result[mask]
        stem=f'{lev}_{name}';sparse.save_npz(pout/(stem+'.npz'),named.T.tocsr())
        (pout/(stem+'_axes.json')).write_text(json.dumps(dict(samples=samples,features=list(np.array(labels)[mask]),representation='equal_share_distinct_eligible_terms',reference='searched-target union plus local full UHGG',human_exclusion=exclude),indent=2,allow_nan=False)+'\n')
        for i,s in enumerate(samples):cov.append(dict(arm=arm,level=lev,profile=name,sample=s,total_LFQ=float(total[i]),single_eligible_term=float(single[i]),multiple_eligible_terms=float(multi[i]),no_eligible_term=float(no[i]),human_withheld=float(withheld[i]),positive_terms_with_missing_annotations=float(partial[i])))
        if variant=='replay':
            fa=json.loads((CP/'profiles'/arm/(stem+'_axes.json')).read_text());fm=sparse.load_npz(CP/'profiles'/arm/(stem+'.npz'))
            assert fa['samples']==samples and fa['features']==list(np.array(labels)[mask]),stem
            err=float(abs(fm-named.T.tocsr()).max()) if fm.nnz else 0.;assert err<=1e-6,(stem,err);checks[stem]=err
pd.DataFrame(cov).to_csv(pout/'coverage.tsv',sep='\t',index=False)
if variant=='replay':
    a=pd.read_csv(pout/'taxonomy_signal.tsv',sep='\t');b=pd.read_csv(CP/'profiles'/arm/'taxonomy_signal.tsv',sep='\t');assert list(a.columns)==list(b.columns);assert np.allclose(a.iloc[:,1:],b.iloc[:,1:],rtol=1e-12)
    a=pd.read_csv(pout/'coverage.tsv',sep='\t');b=pd.read_csv(CP/'profiles'/arm/'coverage.tsv',sep='\t');assert np.allclose(a.select_dtypes('number'),b.select_dtypes('number'),rtol=1e-9,atol=1e-5)
dump(pout/'COMPLETE.json',dict(status='PASS',arm=arm,variant=variant,acquisitions=len(samples),replay_profile_max_abs_error=checks,job=os.environ.get('SLURM_JOB_ID'),
    outputs={p.name:digest(p) for p in pout.iterdir() if p.is_file() and p.name!='COMPLETE.json'}))
print('PASS',arm,variant,flush=True)
