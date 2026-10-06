#!/usr/bin/env python3
"""Retain historical high confidence and production top-30% as distinct routes."""
import csv,gzip,json,math,os,sys
from collections import Counter,defaultdict
from pathlib import Path
import ahocorasick
RUN=Path(sys.argv[1]).resolve();sys.path.insert(0,str(RUN/'frozen'))
from lib import COL,clean,low,load_queries,meta,verify,dump,table
from context_helpers import fasta,context
from origin import classify_members,is_classifiable_accession
from historical_denovo import denovo_one
AA='ACDEFGHIKLMNPQRSTVWY'
m=json.loads((RUN/'inputs/manifest.json').read_text());i=int(sys.argv[2]);t=m['tasks'][i]
for x in m['code']+[m['protocol'],t['prediction'],t['fasta']]:verify(x)
out=RUN/'results'/t['cohort']/t['sample'];out.mkdir(parents=True,exist_ok=True);assert not (out/'COMPLETE.json').exists()
high=Counter();eligible=[];all_positive=Counter();qc=Counter();ranks=Counter()
with open(t['prediction']['path']) as f:
    rd=csv.DictReader(f);assert {COL,'score'}<=set(rd.fieldnames)
    for rownum,r in enumerate(rd,1):
        qc['all_rows']+=1;s=clean(r.get(COL) or '')
        ranks[r.get('beam_idx',r.get('beam_index',r.get('beam_rank','FIELD_NOT_FOUND')))]+=1
        if not s:qc['empty_sequence']+=1;continue
        try:score=float(r.get('score') or 'nan')
        except ValueError:qc['invalid_score']+=1;continue
        if not score>0:qc['nonpositive_or_nan_score']+=1;continue
        if not math.isfinite(score):qc['positive_infinite_score']+=1
        all_positive[s]+=1
        if score>=.9:high[s]+=1
        if 9<=len(s)<=50 and not low(s):eligible.append((rownum,s,score))
threshold=sorted((x[2] for x in eligible),reverse=True)[max(1,math.ceil(.3*len(eligible)))-1] if eligible else None
prod=Counter(x[1] for x in eligible if x[2]>=threshold) if eligible else Counter()
if eligible:
    queries,old_filter,kept=load_queries(t['prediction']['path']);canonical=Counter()
    for s,n in prod.items():canonical[s.replace('I','L')]+=n
    assert queries==canonical and len(kept)==sum(prod.values())
historical=denovo_one((t['sample'],t['prediction']['path']))
assert historical['dn_n']==sum(high.values())
assert historical['dn_non']==sum(n for s,n in high.items() if s[-1] not in 'KR')
if t.get('baseline'):
    b=t['baseline'];assert int(b['denovo_n'])==historical['dn_n'] and int(b['denovo_nontryptic'])==historical['dn_non']
forms={'denovo_highconfidence':high,'denovo':prod}
universe={s.replace('I','L') for ss in forms.values() for s in ss if 5<=len(s)<=50}
machine=ahocorasick.Automaton()
for p in sorted(universe):machine.add_word(p,p)
machine.make_automaton();occ=defaultdict(lambda:defaultdict(set));audit={}
for h,s in fasta(t['fasta']['path'],audit):
    pid=h.split()[0]
    if universe:
        for end,p in machine.iter(s.replace('I','L')):occ[p][pid].add(context(s,end+1-len(p),end+1))
assert audit['sha256']==t['fasta']['sha256']
props={};context_rows=[]
for p in sorted({s.replace('I','L') for ss in forms.values() for s in ss}):
    members=set(occ[p]);cc=set().union(*occ[p].values()) if members else set();p1={c[0] for c in cc}
    in_scope=5<=len(p)<=50
    readable=classify_members(members);unknown=sum(not is_classifiable_accession(x) for x in members)
    origin=('Unknown' if unknown else readable) if members else ('NotInSelectedFASTA' if in_scope else 'OutsideContextLengthScope')
    if not members:readable=origin
    v=dict(peptide=p,length=len(p),origin=origin,readable_origin=readable,n_local_parents=len(members),n_unreadable_parents=unknown,
           parents=';'.join(sorted(members)),n_p1_options=';'.join(sorted(p1)),context_length_scope=in_scope)
    props[p]=v;context_rows.append(v)
fields=['peptide','length','origin','readable_origin','n_local_parents','n_unreadable_parents','parents','n_p1_options','context_length_scope']
with gzip.open(out/'peptide_contexts.tsv.gz','wt') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(context_rows)
counts=Counter();denoms=Counter();quick=[]
for route,ss in forms.items():
    for minimum in ([0,5,9] if route=='denovo_highconfidence' else [9]):
        selected=Counter({s:n for s,n in ss.items() if minimum==0 or minimum<=len(s)<=50})
        for unit in ['canonical_peptide','assigned_bare_sequence','accepted_psm_or_prediction_row']:
            items=Counter({s.replace('I','L'):1 for s in selected}) if unit=='canonical_peptide' else (Counter({s:1 for s in selected}) if unit=='assigned_bare_sequence' else selected)
            total=sum(items.values());nk=sum(n for s,n in items.items() if s[-1] not in 'KR');h=sum(n for s,n in items.items() if s[-1]=='H')
            quick.append(dict(route=route,min_length=minimum,max_length=0 if minimum==0 else 50,unit=unit,denominator=total,C_nonKR=nk,C_H=h,C_nonKR_fraction=nk/total if total else None,C_H_fraction_all=h/total if total else None,C_H_fraction_nonKR=h/nk if nk else None))
            for s,n in items.items():
                v=props[s.replace('I','L')]
                for policy in ['conservative','readable_parent_sensitivity']:
                    origin=v['origin' if policy=='conservative' else 'readable_origin']
                    for origin_class in ['All',origin]:
                        for terminus,residues in [('N_first',[s[0]]),('C_last',[s[-1]]),('N_cleavage_P1',v['n_p1_options'].split(';') if v['n_p1_options'] else ['UNMAPPED'])]:
                            key=(route,minimum,policy,unit,origin_class,terminus);denoms[key]+=n
                            for a in residues:
                                if unit=='canonical_peptide' and a in 'IL' and len(a)==1:a='I/L'
                                counts[key+(a,)]+=n/len(residues)
rows=[]
for key,denom in sorted(denoms.items()):
    assert math.isclose(sum(n for k,n in counts.items() if k[:-1]==key),denom,rel_tol=1e-9)
    alphabet=set(AA)|{'PROTEIN_END','UNMAPPED'}|{k[-1] for k in counts if k[:-1]==key}
    if key[3]=='canonical_peptide':alphabet=(alphabet-{'I','L'})|{'I/L'}
    for residue in sorted(alphabet):rows.append(dict(zip(['route','min_length','origin_policy','unit','origin','terminus'],key),residue=residue,count=counts[key+(residue,)],denominator=denom,fraction=counts[key+(residue,)]/denom))
table(out/'terminal_counts.tsv',rows,fields=['route','min_length','origin_policy','unit','origin','terminus','residue','count','denominator','fraction'])
table(out/'threshold_comparison.tsv',quick)
dump(out/'COMPLETE.json',dict(status='PASS',sample=t['sample'],cohort=t['cohort'],design=t['design'],
    qc=dict(qc),beam_ranks=dict(ranks),eligible_production_rows=len(eligible),production_score_threshold=threshold,
    historical_high_confidence=historical,historical_reader_reproduced=True,capscan_baseline_reproduced=bool(t.get('baseline')),
    positive_prediction_rows=sum(all_positive.values()),high_confidence_rows=sum(high.values()),production_rows=sum(prod.values()),
    confidence_definition='Historical score>=0.9 with no length/complexity filter; nested 5-50 and 9-50. Production positive-score non-low-complexity 9-50 top30% including ties. Predictions are exported rank-zero scope, not calibrated peptide q.',
    context_scope='All actual local parents for 5-50; outside-length predictions remain in historical all-length numerator/denominator but are not called absent from FASTA.',
    scientific_review='PENDING',job_id=os.environ.get('SLURM_JOB_ID'),outputs=[meta(p) for p in out.iterdir() if p.is_file() and p.name!='COMPLETE.json']))
print('DENOVO THRESHOLD PASS',t['cohort'],t['sample'],'cutoff',threshold,'high',sum(high.values()),'top30',sum(prod.values()),flush=True)
