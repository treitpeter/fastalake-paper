#!/usr/bin/env python3
"""Recount transfer/pooling raw results and verify database availability."""
import csv, json, math, os, sys
from pathlib import Path
RUN=Path(sys.argv[1]).resolve();sys.path.insert(0,str(RUN/'frozen'))
from lib import dump, meta, sha, verify, table
from peptide import canonical
from digest import scan
from parameters import parameters

def fasta(p):
    result={};name=None;seq=[]
    def put():
        if name is not None:
            assert name not in result
            result[name]=''.join(seq)
    with open(p) as f:
        for line in f:
            line=line.strip()
            if line.startswith('>'):put();name=line[1:].split()[0];seq=[]
            elif line:seq.append(line)
    put();return result

def accepted(x):
    verify(x);s=set();modified=set();count=0
    with open(x['path']) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['rank']!='1' or r['label']!='1':continue
            q=float(r['peptide_q']);assert math.isfinite(q) and 0<=q<=1
            if q<=.01:
                p=canonical(r['peptide'])
                if 5<=len(p)<=50:s.add(p);modified.add(r['peptide']);count+=1
    return s,modified,count

def run(index):
    m=json.loads((RUN/'inputs/manifest.json').read_text())
    for x in m['frozen_code']+[m['protocol'],m['historical_transfer'],m['historical_pool']]:verify(x)
    t=m['tasks'][index];out=RUN/'results'/t['id'];out.mkdir(parents=True,exist_ok=True)
    assert not (out/'COMPLETE.json').exists()
    data={};configs={};databases={};sources=[];sourcefastas={}
    for name,x in t['source_fastas'].items():
        verify(x);sourcefastas[name]=fasta(x['path']);sources.append(x)
    for arm,x in t['arms'].items():
        for k in ['config','fasta','search']:
            verify(x[k]);sources.append(meta(x[k]['path']))
        c=json.loads(Path(x['config']['path']).read_text());configs[arm]=c
        assert len(c['mzml_paths'])==1 and Path(c['mzml_paths'][0]).resolve()==Path(t['raw']).resolve()
        assert Path(c['database']['fasta']).resolve()==Path(x['fasta']['path']).resolve()
        db=fasta(x['fasta']['path']);databases[arm]=db
        if t['kind']=='transfer':
            assert db==sourcefastas[x['donor']],(t['id'],arm,'donor FASTA differs')
        else:
            donors=x['donors'];expected={}
            assert len(donors)==len(set(donors))
            if arm=='pool_other11':assert len(donors)==11 and t['sample'] not in donors
            if arm=='pool_all12':assert len(donors)==12 and t['sample'] in donors
            for donor in donors:
                for a,s in sourcefastas[donor].items():
                    assert a not in expected or expected[a]==s
                    expected[a]=s
            assert db==expected,(t['id'],arm,'pool membership differs')
        s,modified,n=accepted(x['search']);data[arm]=s
        if x.get('historical_set'):
            verify(x['historical_set']);old=Path(x['historical_set']['path']).read_text().splitlines()
            assert s=={canonical(p) for p in old if 5<=len(canonical(p))<=50},(t['id'],arm,'historical identity set differs')
        (out/(arm+'_5_50.txt')).write_text(''.join(p+'\n' for p in sorted(s)))
    assert len({json.dumps(parameters(c),sort_keys=True) for c in configs.values()})==1
    assert {c['version'] for c in configs.values()}=={'0.14.6'}
    own=data['own'];rows=[];availability=[]
    for arm,s in data.items():
        present,digest,audit=scan(t['arms'][arm]['fasta']['path'],own|s)
        assert audit['sha256']==t['arms'][arm]['fasta']['sha256']
        assert all(digest[p]>0 for p in s),(t['id'],arm,'accepted peptide digest mismatch')
        for p in sorted(own):
            availability.append(dict(sample=t['sample'],kind=t['kind'],arm=arm,peptide=p,length=len(p),
                final_sequence_present=int(present[p]>0),digest_compatible=int(digest[p]>0),identified=int(p in s)))
        for lo in (5,9):
            a={p for p in s if len(p)>=lo};b={p for p in own if len(p)>=lo}
            rows.append(dict(sample=t['sample'],kind=t['kind'],group=t['group'],arm=arm,min_length=lo,
                own=len(b),peptides=len(a),overlap=len(a&b),lost=len(b-a),gained=len(a-b),
                count_retention_percent=100*len(a)/len(b),own_set_recovery_percent=100*len(a&b)/len(b),
                own_sequence_present=sum(present[p]>0 for p in b),own_digest_compatible=sum(digest[p]>0 for p in b),
                own_digest_compatible_unidentified=sum(digest[p]>0 and p not in a for p in b),
                database_accessions=len(databases[arm]),donor=t['arms'][arm].get('donor','pooled')))
    old=list(csv.DictReader(open(m['historical_transfer' if t['kind']=='transfer' else 'historical_pool']['path']),delimiter='\t'))
    if t['kind']=='transfer':
        r=next(r for r in old if r['run']==t['sample']);a=next(r for r in rows if r['arm']=='T1_within_group' and r['min_length']==5)
        for new,prior in [('own','own_canonical'),('peptides','swap_canonical'),('overlap','intersection'),('lost','lost'),('gained','gained')]:assert a[new]==int(r[prior]),(t['sample'],new,a[new],r[prior])
    else:
        for a in rows:
            if a['min_length']!=5:continue
            r=next(r for r in old if r['sample']==t['sample'] and r['arm']==a['arm'])
            for new,prior in [('own','own_peptides'),('peptides','canonical_peptides'),('overlap','shared_with_own'),('lost','own_peptides_lost'),('gained','peptides_gained_vs_own')]:assert a[new]==int(r[prior]),(t['sample'],new,a[new],r[prior])
    table(out/'comparisons.tsv',rows);table(out/'peptide_availability.tsv',availability)
    dump(out/'COMPLETE.json',dict(status='PASS',kind=t['kind'],sample=t['sample'],arms=list(data),
        historical_counts_and_identities_match=True,databases_and_pools_match=True,search_parameters_match=True,
        exclusion_scope='target final FASTA excluded; shared upstream evidence retained',
        job_id=os.environ.get('SLURM_JOB_ID'),sources=sources,
        outputs=[meta(p) for p in out.iterdir() if p.is_file() and p.name!='COMPLETE.json']))
    print('TRANSFER/POOL PASS',t['id'],flush=True)

if __name__=='__main__':run(int(sys.argv[2]))
