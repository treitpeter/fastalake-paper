#!/usr/bin/env python3
"""Matched public/cohort/metagenome identity summaries with declared weighting."""
import csv, json, os, sys
from collections import defaultdict
from pathlib import Path
from statistics import mean, median
RUN=Path(sys.argv[1]).resolve();sys.path.insert(0,str(RUN/'frozen'))
from lib import dump, meta, sha, verify, table

def parameters(c):
    p={k:v for k,v in c.items() if k not in ['version','output_paths','output_directory','mzml_paths','quant']}
    p['database']=dict(p['database']);p['database'].pop('fasta',None)
    e=p['database']['enzyme']
    if e.get('c_terminal') is None: e=dict(e,c_terminal=True)
    if e.get('semi_enzymatic') is None: e=dict(e,semi_enzymatic=False)
    p['database']['enzyme']=e
    return p

def readset(x):
    verify(x);lines=Path(x['path']).read_text().splitlines();s=set(lines)
    assert len(s)==len(lines) and all(5<=len(p)<=50 and 'I' not in p for p in s)
    return s

def pair(a,b):
    return dict(n_test=len(a),n_reference=len(b),overlap=len(a&b),gained=len(a-b),lost=len(b-a),
        count_ratio=len(a)/len(b) if b else None,reference_recovery=len(a&b)/len(b) if b else None,
        test_recovery=len(a&b)/len(a) if a else None)

def run(index):
    m=json.loads((RUN/'inputs/manifest.json').read_text())
    for x in m['frozen_code']+[m['protocol']]: verify(x)
    pilot=index<2;cohort=['PREDICT','CAPSCAN'][index%2]
    tasks=[t for t in m['tasks'] if t['cohort']==cohort]
    if pilot:tasks=tasks[:1]
    out=RUN/'results'/('pilot_'+cohort if pilot else cohort);out.mkdir(parents=True,exist_ok=True)
    assert not (out/'COMPLETE.json').exists()
    rows=[];excluded=[];pooled=defaultdict(lambda:defaultdict(set));comparators=defaultdict(int)
    for t in tasks:
        for x in t['receipts']:verify(x)
        sets={a:readset(x) for a,x in t['sets'].items()}
        for a,x in t['nested_sets'].items():
            assert readset(x)=={p for p in sets[a] if len(p)>=9}
        configs={}
        for a,x in t['configs'].items():
            verify(x);c=json.loads(Path(x['path']).read_text());configs[a]=c
            assert len(c['mzml_paths'])==1 and Path(c['mzml_paths'][0]).resolve()==Path(t['raw']).resolve()
        assert len({json.dumps(parameters(c),sort_keys=True) for c in configs.values()})==1,(cohort,t['sample'],'search parameter mismatch')
        assert {c['version'] for c in configs.values()}=={'0.14.6'}
        if t.get('source_proof'):verify(t['source_proof'])
        d=t['design'];subject=d['subject_id']
        if not subject:excluded.append(dict(sample=t['sample'],reason='missing_subject_id',acquisition_results_retained=True))
        for a,b in [('FL_public','MG'),('FL_public','FL_cohortMG'),('FL_cohortMG','MG')]:
            if b not in sets:continue
            comparators[a,b]+=1
            for lo in (5,9):
                x={p for p in sets[a] if len(p)>=lo};y={p for p in sets[b] if len(p)>=lo}
                pooled[a,b,lo]['test'].update(x);pooled[a,b,lo]['reference'].update(y)
                rows.append(dict(cohort=cohort,sample=t['sample'],subject_id=subject,
                    specimen_id=d['specimen_id'],group=d['group'],aliquot=d['aliquot'],
                    test_arm=a,reference_arm=b,min_length=lo,**pair(x,y)))
    table(out/'per_acquisition.tsv',rows)
    table(out/'participant_exclusions.tsv',excluded,['sample','reason','acquisition_results_retained'])
    metrics=['n_test','n_reference','overlap','gained','lost','count_ratio','reference_recovery','test_recovery']
    # Equal acquisition weight within a recorded specimen; equal specimen weight
    # within participant and stratum. ALL retains the observed specimen mix.
    specimen=defaultdict(list)
    for r in rows:
        if not r['subject_id']:continue
        for stratum in ['ALL','group='+r['group'],'group='+r['group']+';aliquot='+r['aliquot']]:
            specimen[r['test_arm'],r['reference_arm'],r['min_length'],stratum,r['subject_id'],r['specimen_id']].append(r)
    participant=defaultdict(list);srows=[]
    for k,rr in specimen.items():
        sr=dict(test_arm=k[0],reference_arm=k[1],min_length=k[2],stratum=k[3],subject_id=k[4],specimen_id=k[5],n_acquisitions=len(rr))
        for metric in metrics:
            vals=[r[metric] for r in rr if r[metric] is not None];sr[metric]=mean(vals) if vals else None
            sr[metric+'_n_defined']=len(vals)
        srows.append(sr);participant[k[:5]].append(sr)
    table(out/'per_specimen.tsv',srows)
    prows=[]
    for k,rr in participant.items():
        pr=dict(test_arm=k[0],reference_arm=k[1],min_length=k[2],stratum=k[3],subject_id=k[4],n_specimens=len(rr),n_acquisitions=sum(r['n_acquisitions'] for r in rr))
        for metric in metrics:
            vals=[r[metric] for r in rr if r[metric] is not None];pr[metric]=mean(vals) if vals else None
            pr[metric+'_n_defined']=len(vals)
        prows.append(pr)
    table(out/'per_participant.tsv',prows)
    summary=[]
    for a,b,lo in sorted(pooled):
        for population,data in [('acquisition',rows),('participant',prows)]:
            strata=['ALL'] if population=='acquisition' else sorted({r['stratum'] for r in prows})
            for stratum in strata:
                rr=[r for r in data if r['test_arm']==a and r['reference_arm']==b and r['min_length']==lo and (population=='acquisition' or r['stratum']==stratum)]
                if not rr:continue
                for metric in metrics:
                    vals=[r[metric] for r in rr if r[metric] is not None]
                    summary.append(dict(test_arm=a,reference_arm=b,min_length=lo,population=population,stratum=stratum,metric=metric,n_total=len(rr),n_defined=len(vals),mean=mean(vals) if vals else None,median=median(vals) if vals else None,minimum=min(vals) if vals else None,maximum=max(vals) if vals else None))
    table(out/'summary.tsv',summary)
    pooled_rows=[dict(test_arm=a,reference_arm=b,min_length=lo,matched_acquisitions=comparators[a,b],**pair(v['test'],v['reference'])) for (a,b,lo),v in sorted(pooled.items())]
    table(out/'pooled_identity.tsv',pooled_rows)
    outputs=[meta(p) for p in out.iterdir() if p.is_file() and p.name!='COMPLETE.json']
    dump(out/'COMPLETE.json',dict(status='PASS',cohort=cohort,pilot=pilot,acquisitions=len(tasks),
        matched_metagenome=comparators['FL_public','MG'],missing_subjects=len(excluded),
        participant_count=len({t['design']['subject_id'] for t in tasks if t['design']['subject_id']}),
        source_proof='PASS' if cohort=='CAPSCAN' else 'ARCHIVED_NOMETAG_ARM_SOURCE_LEDGER_STILL_REQUIRED',
        search_parameters_matched=True,job_id=os.environ.get('SLURM_JOB_ID'),outputs=outputs))
    print('PUBLIC COMPARISON PASS',cohort,len(tasks),flush=True)

if __name__=='__main__':run(int(sys.argv[2]))
