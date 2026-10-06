#!/usr/bin/env python3
"""Frozen single-acquisition searches; validate outputs before releasing dependants."""
import csv,gzip,hashlib,json,math,os,pathlib,re,subprocess,sys,time
from collections import Counter,defaultdict
RUN=pathlib.Path(sys.argv[1]).resolve()
sys.path.insert(0,str(RUN/'frozen'))
from lib import meta,sha,verify,dump,table,command
from peptide import canonical
def read_fasta(path):
    out={};h=None;parts=[]
    for line in open(path):
        if line.startswith('>'):
            if h is not None:
                assert h not in out;out[h]=''.join(parts).upper().replace('I','L')
            h=line[1:].split()[0];parts=[]
        else:parts.append(line.strip())
    if h is not None:assert h not in out;out[h]=''.join(parts).upper().replace('I','L')
    return out
def recount(path):
    z=Counter();sets={n:set() for n in ['peptide_q01','both_q01','length9_50','target_all','decoy_q01']};psms=[]
    with open(path) as f:
        rd=csv.DictReader(f,delimiter='\t')
        assert {'peptide','label','rank','peptide_q','spectrum_q','proteins'}<=set(rd.fieldnames)
        for r in rd:
            z['reported_psms']+=1
            assert r['label'] in ['1','-1']
            for k in ['peptide_q','spectrum_q']:
                v=float(r[k]);assert math.isfinite(v) and 0<=v<=1
            if r['rank']!='1':continue
            p=canonical(r['peptide']);assert p
            if r['label']=='1':sets['target_all'].add(p)
            if float(r['peptide_q'])>.01:continue
            if r['label']=='-1':sets['decoy_q01'].add(p);continue
            sets['peptide_q01'].add(p);psms.append(dict(r,canonical=p))
            z['target_psms_peptide_q01']+=1
            if float(r['spectrum_q'])<=.01:sets['both_q01'].add(p)
            if 9<=len(p)<=50:sets['length9_50'].add(p)
            z['Cterm_'+p[-1]]+=1
    return dict(z),sets,psms
def run(index):
    m=json.loads((RUN/'inputs/manifest.json').read_text())
    for f in m['frozen_code']:verify(f)
    t=m['tasks'][index];case=RUN/'results'/t['id'];case.mkdir(parents=True,exist_ok=True)
    assert not (case/'COMPLETE.json').exists(),'Refuse completed-output overwrite'
    for k in ['fasta','raw','source_config','config']:verify(t[k])
    config=json.loads(pathlib.Path(t['config']['path']).read_text())
    assert config['database']['fasta']==t['fasta']['path']
    assert config['mzml_paths']==[t['raw']['path']]
    assert config['database']['generate_decoys'] is True
    if t['arm']=='unspecific':
        assert config['database']['enzyme']['cleave_at']==''
        assert not config['database']['enzyme'].get('semi_enzymatic')
    verify(m['sage'])
    command([m['sage']['path'],t['config']['path'],'--batch-size','1','--write-pin','--disable-telemetry-i-dont-want-to-improve-sage'],case/'sage')
    path=case/'search/results.sage.tsv';z,sets,psms=recount(path)
    outputs=list((case/'search').glob('*'))
    assert path.exists() and (case/'search/results.json').exists()
    actual=json.loads((case/'search/results.json').read_text())
    assert actual['database']['enzyme']==config['database']['enzyme']
    assert pathlib.Path(actual['database']['fasta']).resolve()==pathlib.Path(t['fasta']['path']).resolve()
    # Independently verify every accepted canonical peptide is in its reported target parent.
    fasta=read_fasta(t['fasta']['path']);relations=defaultdict(set);contexts=defaultdict(set)
    for r in psms:
        p=r['canonical'];parents=[a for a in r['proteins'].split(';') if not a.startswith('rev_')]
        assert parents,(t['id'],r['proteins'])
        matched=False
        for pid in parents:
            assert pid in fasta,('Unresolved protein ID',pid)
            s=fasta[pid];start=s.find(p)
            while start!=-1:
                matched=True;relations[p].add(pid)
                end=start+len(p)
                n=(start==0 or (start==1 and s[0]=='M') or (s[start-1] in 'KR' and s[start]!='P'))
                c=(end==len(s) or (s[end-1] in 'KR' and s[end]!='P'))
                contexts[p].add('fully_tryptic' if n and c else 'semi_tryptic' if n or c else 'fully_nontryptic')
                start=s.find(p,start+1)
        assert matched,('Accepted peptide absent from reported FASTA parents',p)
    rows=[dict(peptide=p,length=len(p),c_terminal=p[-1],non_KR=int(p[-1] not in 'KR'),parents=';'.join(sorted(relations[p])),enzyme_contexts=';'.join(sorted(contexts[p])),ambiguous_context=int(len(contexts[p])>1),both_q01=int(p in sets['both_q01'])) for p in sorted(sets['peptide_q01'])]
    table(case/'accepted_peptide_context.tsv',rows,['peptide','length','c_terminal','non_KR','parents','enzyme_contexts','ambiguous_context','both_q01'])
    for name,peps in sets.items():(case/(name+'.txt')).write_text(''.join(p+'\n' for p in sorted(peps)))
    baseline=t.get('baseline_search')
    comparison={}
    if baseline:
        verify(baseline);_,oldsets,_=recount(baseline['path'])
        for name in ['peptide_q01','both_q01','length9_50']:
            a,b=sets[name],oldsets[name]
            comparison[name]=dict(new=len(a),old=len(b),intersection=len(a&b),gained=len(a-b),lost=len(b-a),old_recovery=len(a&b)/len(b) if b else None)
            for tag,seqs in [('gained',a-b),('lost',b-a)]:(case/(name+'_'+tag+'.txt')).write_text(''.join(p+'\n' for p in sorted(seqs)))
    nontryptic=sum('fully_nontryptic' in x for x in contexts.values())
    summary=dict(status='PASS',sample=t['sample'],arm=t['arm'],task=t['id'],reported=z,peptide_counts={k:len(v) for k,v in sets.items()},baseline_comparison=comparison,fully_nontryptic_compatible_peptides=nontryptic,input_fasta=meta(t['fasta']['path']),input_raw=meta(t['raw']['path']),config=meta(t['config']['path']),engine=m['sage'],job_id=os.environ.get('SLURM_JOB_ID'),scope='Peptide confidence and sequence containment validated; transferred-only LFQ confidence is not inferred',outputs=[meta(p) for p in sorted(case.iterdir()) if p.is_file() and p.name!='COMPLETE.json']+[meta(p) for p in outputs if p.is_file()])
    if t.get('require_nontryptic'):assert nontryptic>0,'Pilot did not demonstrate internal unspecific peptide recovery'
    if config.get('quant',{}).get('lfq'):assert (case/'search/lfq.tsv').exists()
    dump(case/'COMPLETE.json',summary)
    print('VALIDATED',t['id'],summary['peptide_counts'],flush=True)
if __name__=='__main__':run(int(sys.argv[2]))
