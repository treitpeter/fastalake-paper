#!/usr/bin/env python3
"""Reproduce baseline leaders, then audit fixed-dictionary and native-set alternatives."""
import csv,gzip,hashlib,itertools,json,math,os,re,sys
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
RUN=Path(sys.argv[1]).resolve();sys.path.insert(0,str(RUN/'frozen'))
from lib import dump,meta,verify
from raw_fasta import fasta
TIERS=['TechnicalReplicates','WorkflowReplicates','OnePiece','Intra','Interindividual']
def tier(s):return next(t for t in TIERS if t in s)
def matrix(p):
    # Native identifiers were saved as zero-padded Unicode: 46 GB for the
    # largest keys array, despite much smaller actual strings. Preserve exact
    # strings while decoding bounded chunks; no quantity/inference change.
    import zipfile
    keys=[]
    with zipfile.ZipFile(p) as archive:
        with archive.open('keys.npy') as stream:
            version=np.lib.format.read_magic(stream)
            shape,order,dtype=np.lib.format._read_array_header(stream,version)
            assert len(shape)==1 and not order and dtype.kind=='U'
            chunk=max(1,8*1024*1024//dtype.itemsize)
            for start in range(0,shape[0],chunk):
                count=min(chunk,shape[0]-start)
                block=stream.read(count*dtype.itemsize)
                assert len(block)==count*dtype.itemsize
                keys.extend(np.frombuffer(block,dtype=dtype,count=count).tolist())
            assert stream.read(1)==b''
    with np.load(p) as z:
        values=z['values'];samples=z['samples']
    assert values.shape==(len(keys),len(samples))
    return pd.DataFrame(values,index=keys,columns=samples)
def state(f,annotations):
    profiles=defaultdict(lambda:np.zeros(f.shape[1]));total=f.sum().to_numpy()
    for group,v in zip(f.index,f.to_numpy()):
        a=annotations[group]
        for ko in a['kos']:profiles[ko,a['taxon']]+=v
    rows=[]
    for ko in sorted({k for k,t in profiles}):
        names=sorted(t for k,t in profiles if k==ko);v=np.array([profiles[ko,t] for t in names]);r=np.array([n!='Unassigned' for n in names]);amount=v.sum(axis=0);resolved=v[r].sum(axis=0)
        for i,s in enumerate(f.columns):
            if amount[i]<=0:continue
            order=sorted(np.flatnonzero(r),key=lambda j:(-v[j,i],names[j]));top=order[0] if order and resolved[i]>0 else None
            second=v[order[1],i] if top is not None and len(order)>1 else 0
            rows.append(dict(sample=s,tier=tier(s),ko=ko,normalised_function_intensity=amount[i]/total[i] if total[i]>0 else np.nan,
                resolved_fraction=resolved[i]/amount[i],leading_species=names[top] if top is not None else '',
                leader_margin_over_runnerup=(v[top,i]-second)/resolved[i] if top is not None else 0.,
                state='resolved' if top is not None else 'positive_unresolved'))
    return pd.DataFrame(rows,columns=['sample','tier','ko','normalised_function_intensity','resolved_fraction','leading_species','leader_margin_over_runnerup','state'])

m=json.loads((RUN/'inputs/leader_manifest_revision2.json').read_text())
for x in m['code']+m['inputs']:verify(x)
main=json.loads((RUN/'inputs/manifest_revision2.json').read_text());pilot=int(sys.argv[2])==0;selected=main['samples'][:2] if pilot else main['samples'];out=RUN/'results'/('pilot_leaders' if pilot else 'leaders');out.mkdir(parents=True,exist_ok=True);assert not (out/'COMPLETE.json').exists()
ann=pd.read_csv(m['representative_annotations']['path'],sep='\t',keep_default_na=False)
fixed={r.study_group:dict(kos=set(re.findall(r'K\d{5}',r.KEGG_ko)),taxon=r.representative_taxon) for r in ann.itertuples()}
states={};frames={};receipts=[]
for a in main['arms']:
    name=a['name'];folder=RUN/'results'/name;receipt=json.loads((folder/'COMPLETE.json').read_text());assert receipt['status']=='PASS'
    p=folder/'same_file_id.fixed_groups.npz';verify(next(x for x in receipt['outputs'] if Path(x['path'])==p));f=matrix(p).reindex(columns=selected)
    assert set(f.index)<=set(fixed);label='fixed_'+name;frames[label]=f;states[label]=state(f,fixed);receipts.append(meta(folder/'COMPLETE.json'))
old=pd.read_csv(m['baseline_leaders']['path'],sep='\t',keep_default_na=False);old=old[old['sample'].isin(selected)];now=states['fixed_original']
paired=old.merge(now,on=['sample','ko'],suffixes=('_old','_new'),how='outer',validate='1:1',indicator=True);assert paired['_merge'].eq('both').all()
assert paired.leading_species_old.eq(paired.leading_species_new).all()
for field in ['normalised_function_intensity','resolved_fraction','leader_margin_over_runnerup']:
    np.testing.assert_allclose(paired[field+'_old'].astype(float),paired[field+'_new'].astype(float),rtol=1e-9,atol=1e-12)
print('Historical functional leaders reproduced',len(paired),flush=True)
# Project native union entries through exact sequence, preserving every available
# annotation disagreement. No annotation is attached merely because IDs look alike.
audit={};union={h.split()[0]:s for h,s in fasta(m['union_fasta']['path'],audit)};assert audit['sha256']==m['union_fasta']['sha256']
assert all(pid=='FLU_'+hashlib.sha256(seq.encode()).hexdigest() for pid,seq in union.items())
seq_to_union={s:p for p,s in union.items()};assert len(seq_to_union)==len(union)
audit={};source_seqs={h.split()[0]:s for h,s in fasta(m['annotation_fasta']['path'],audit)};assert audit['sha256']==m['annotation_fasta']['sha256']
source=pd.read_csv(m['annotations']['path'],sep='\t',keep_default_na=False);assert source.accession.is_unique
ko_by_seq=defaultdict(list)
for r in source.itertuples():
    if r.accession in source_seqs and source_seqs[r.accession] in seq_to_union:
        ko_by_seq[seq_to_union[source_seqs[r.accession]]].append(set(re.findall(r'K\d{5}',r.KEGG_ko)))
aliases=defaultdict(set)
with gzip.open(m['union_provenance']['path'],'rt') as f:
    for r in csv.DictReader(f,delimiter='\t'):aliases[r['union_id']].add(r['original_accession'])
assert set(aliases)==set(union)
genomes=pd.read_csv(m['genomes']['path'],sep='\t',keep_default_na=False);fallback={r.Species_rep:r.Lineage for r in genomes.itertuples() if r.Lineage};species={}
for r in genomes.itertuples():
    terms=[x[3:] for x in (r.Lineage or fallback.get(r.Species_rep,'')).split(';') if x.startswith('s__') and len(x)>3]
    if terms:species[r.Genome]=terms[-1]
entry={}
for pid in union:
    ks=ko_by_seq[pid];kos=set.intersection(*ks) if ks else set()
    ts={species[g] for a in aliases[pid] for g in re.findall(r'MGYG\d+',a) if g in species}
    entry[pid]=dict(kos=kos,taxa=ts,taxon=next(iter(ts)) if len(ts)==1 else 'Unassigned',annotation_sources=len(ks),ko_disagreement=bool(ks and any(k!=ks[0] for k in ks)))
native_ann={};auditrows=[]
for mode in ['same_file_id','pooled_q01','pooled_q05']:
    folder=RUN/'results/native'/mode;c=json.loads((folder/'COMPLETE.json').read_text());assert c['status']=='PASS' and not c['pilot']
    p=folder/'protein_matrix.npz';verify(next(x for x in c['outputs'] if Path(x['path'])==p));f=matrix(p).reindex(columns=selected);receipts.append(meta(folder/'COMPLETE.json'))
    for g in f.index:
        if g in native_ann:continue
        parents=g.split(';');assert all(p in entry for p in parents)
        # Require an annotation on every parent and agreement across every parent.
        kos=set.intersection(*(entry[p]['kos'] for p in parents))
        taxons={entry[p]['taxon'] for p in parents};taxon=next(iter(taxons)) if len(taxons)==1 and 'Unassigned' not in taxons else 'Unassigned'
        native_ann[g]=dict(kos=kos,taxon=taxon)
        auditrows.append(dict(protein_set=g,parents=len(parents),parents_with_ko=sum(bool(entry[p]['kos']) for p in parents),parents_with_species=sum(entry[p]['taxon']!='Unassigned' for p in parents),parents_with_disagreeing_ko_aliases=sum(entry[p]['ko_disagreement'] for p in parents),consensus_kos=';'.join(sorted(kos)),consensus_taxon=taxon))
    label='native_'+mode;frames[label]=f;states[label]=state(f,native_ann)
pd.DataFrame(auditrows).to_csv(out/'native_annotation_audit.tsv.gz',sep='\t',index=False)
rows=[];transitions=[];agreements=[];tech=[]
base=states['fixed_original'];samples=selected;kos=sorted(set().union(*(set(f.ko) for f in states.values())))
grid=pd.MultiIndex.from_product([samples,kos],names=['sample','ko']).to_frame(index=False)
for label,f in states.items():
    current=grid.merge(f,on=['sample','ko'],how='left',validate='1:1');current['state']=current.state.fillna('not_detected');current['leading_species']=current.leading_species.fillna('');current['tier']=current['sample'].map(tier);current['method']=label
    current.to_csv(out/(label+'.function_states.tsv.gz'),sep='\t',index=False);rows.append(current)
    if label=='fixed_original':continue
    pp=current.merge(base,on=['sample','ko'],suffixes=('_test','_original'),how='left',validate='1:1');pp['state_original']=pp.state_original.fillna('not_detected')
    pp['both_resolved']=pp.state_test.eq('resolved')&pp.state_original.eq('resolved')
    pp['both_strict']=pp.both_resolved & pp.resolved_fraction_test.ge(.2)&pp.resolved_fraction_original.ge(.2)&pp.leader_margin_over_runnerup_test.ge(.05)&pp.leader_margin_over_runnerup_original.ge(.05)
    pp['leader_agrees']=np.where(pp.both_resolved,pp.leading_species_test.eq(pp.leading_species_original),np.nan)
    pp.to_csv(out/(label+'.matched_baseline.tsv.gz'),sep='\t',index=False)
    tr=pp.groupby(['state_original','state_test']).size().rename('cells').reset_index();tr['method']=label;transitions.append(tr)
    for threshold in ['both_resolved','both_strict']:
        for t in TIERS:
            rr=pp[pp[threshold]&pp.tier_test.eq(t)]
            byacq=rr.groupby('sample').leader_agrees.mean();byko=rr.groupby('ko').leader_agrees.mean()
            agreements.append(dict(method=label,tier=t,threshold=threshold,cells=len(rr),acquisitions=len(byacq),functions=len(byko),equal_acquisition_agreement=byacq.mean(),equal_function_agreement=byko.mean(),pooled_cell_agreement=rr.leader_agrees.mean()))
    # Technical/workflow label stability uses identical eligible functions and
    # actual common resolved pairs, with the same 80% design coverage rule.
    for t in TIERS[:2]:
        names=[s for s in samples if tier(s)==t];need=math.ceil(.8*len(names));rr=pp[pp.tier_test.eq(t)]
        eligible=set(rr[rr.state_test.eq('resolved')].groupby('ko').size().loc[lambda x:x>=need].index)&set(rr[rr.state_original.eq('resolved')].groupby('ko').size().loc[lambda x:x>=need].index)
        for threshold in ['both_resolved','both_strict']:
            perko=[]
            for ko,gg in rr[rr.ko.isin(eligible)&rr[threshold]].groupby('ko'):
                n=len(gg)
                if n<2:continue
                a=gg.leading_species_original.to_numpy();b=gg.leading_species_test.to_numpy();pairs=list(itertools.combinations(range(n),2))
                perko.append(dict(ko=ko,pairs=len(pairs),original_agreement=sum(a[i]==a[j] for i,j in pairs)/len(pairs),test_agreement=sum(b[i]==b[j] for i,j in pairs)/len(pairs)))
            for row in perko:tech.append(dict(method=label,tier=t,threshold=threshold,**row))
pd.concat(transitions).to_csv(out/'resolution_transitions.tsv',sep='\t',index=False);pd.DataFrame(agreements).to_csv(out/'matched_leader_summary.tsv',sep='\t',index=False)
tf=pd.DataFrame(tech);tf.to_csv(out/'technical_function_agreement.tsv',sep='\t',index=False)
if len(tf):tf.groupby(['method','tier','threshold'])[['original_agreement','test_agreement']].mean().reset_index().to_csv(out/'technical_agreement_summary.tsv',sep='\t',index=False)
report='''# HSTOOL leading-species comparison

The original frozen-dictionary function intensities, resolved fractions and leaders
are independently reproduced before comparison. The same annotations then apply
to all five fixed-dictionary quantity matrices. This isolates measured evidence
changes on a restricted common reporting dictionary; unmapped features are outside
that diagnostic and explicitly counted by the parent extraction.

Native Sage protein sets use no representative cover. Annotations are projected
through exact unmodified protein sequence. KOs must agree across every annotated
alias and every parent; missing parent annotation prevents consensus. A species
label requires one available species annotation per parent and agreement across
all parents. Unknown/conflicting sets remain unresolved. Available taxonomy is
sequence-compatible labeling, not unique molecular origin. This conservative
consensus differs deliberately from the reporting-representative convention;
the contrast includes inference/annotation policy as well as measurement.

Each KO receives the full intensity of a compatible group, as in the original
leader analysis; totals across KOs are not additive. All missing, positive-unresolved
and resolved cells are saved. Leader agreement uses common resolved cells; strict
comparisons require 20% resolved signal and 5% leader margin in both methods.
Technical/workflow stability uses the common 80%-eligible function set and exactly
the same observed pairs; equal-function summaries accompany all pair counts.
Overlapping pairs do not supply independent biological replicates. Per-cell transfer
confidence remains unestablished even for the corrected pooled-q branches.

Scientific/visual review and manuscript integration remain pending.
'''
(out/'REPORT.md').write_text(report)
dump(out/'COMPLETE.json',dict(status='PASS',pilot=pilot,acquisitions=len(selected),methods=list(states),original_function_baseline_reproduced=True,source_receipts=receipts,scientific_review='PENDING',job_id=os.environ.get('SLURM_JOB_ID'),outputs=[meta(p) for p in out.iterdir() if p.is_file() and p.name!='COMPLETE.json']))
print('HSTOOL LEADERS PASS',len(states),flush=True)
