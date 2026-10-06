#!/usr/bin/env python3
"""Controlled clustering after lookup; separates peptide availability from search yield."""
import csv,hashlib,json,os,subprocess,sys,time
from pathlib import Path
import ahocorasick

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()
def fasta(p):
    header=None;parts=[]
    with Path(p).open() as f:
        for line in f:
            if line.startswith('>'):
                if header is not None:yield header,''.join(parts)
                header=line[1:].strip();parts=[]
            else:parts.append(line.strip())
    if header is not None:yield header,''.join(parts)
def sequences(p):return {hashlib.sha256(s.encode()).hexdigest() for _,s in fasta(p)}
def run(base,index):
    m=json.loads((base/'manifest.json').read_text());t=m['tasks'][index];out=base/'results'/t['sample'];out.mkdir(parents=True,exist_ok=False)
    sources=[]
    for r in t['sources']+m['binaries']:
        p=Path(r['path']);st=p.stat()
        assert st.st_size==r['bytes'] and st.st_mtime_ns==r['mtime_ns'],str(p)
        h=sha(p)
        if 'sha256' in r:assert h==r['sha256'],str(p)
        assert p.stat().st_mtime_ns==st.st_mtime_ns
        sources.append(dict(r,sha256=h))
    (out/'SOURCE_HASHES.json').write_text(json.dumps(sources,indent=2)+'\n')
    def command(label,args):
        started=time.time()
        with (out/(label+'.stdout')).open('w') as f,(out/(label+'.stderr')).open('w') as e:
            r=subprocess.run(list(map(str,args)),stdout=f,stderr=e)
        with (out/'COMMANDS.jsonl').open('a') as f:f.write(json.dumps(dict(label=label,argv=list(map(str,args)),returncode=r.returncode,seconds=time.time()-started))+'\n')
        assert r.returncode==0,label
    par=m['parsimony'];mm=m['mmseqs'];pred=t['predictions'];candidates=t['candidates']
    def select(label,p):
        dest=out/label;dest.mkdir()
        command(label,[par,'--strategy','razor','--tiebreak','hash-acc','--fasta',p,'--peptides',pred,'--output-dir',dest,'--sample',t['sample'],'--min-length','9','--max-length','50','--normalize-il','-t','2'])
        return dest/(t['sample']+'_razor.fasta')
    baseline=select('unclustered_reproduction',candidates)
    assert sequences(baseline)==sequences(t['final_fasta']),'No-clustering reproduction differs: do not interpret intervention'
    # Stable accession order before MMseqs; candidate set itself is unchanged.
    records=list(fasta(candidates));records.sort(key=lambda x:x[0]);ordered=out/'candidates_sorted.fasta'
    assert len({h.split()[0] for h,_ in records})==len(records),'Duplicate candidate accession'
    with ordered.open('w') as f:
        for h,s in records:f.write('>'+h+'\n'+s+'\n')
    candidate_count=len(records);del records
    db=out/'candidate_db';command('createdb',[mm,'createdb',ordered,db,'--shuffle','0'])
    stages={'candidate_unclustered':Path(candidates),'selected_unclustered':baseline}
    cluster_counts={}
    for threshold in (['0.99']+m['new_thresholds'] if t['sample'] in m['pilot_samples'] else m['new_thresholds']):
        key='cluster_'+threshold.replace('.','');folder=out/key;folder.mkdir()
        clusters=folder/'clusters';tmp=folder/'tmp';representatives=folder/'representatives'
        command(key,[mm,'linclust',db,clusters,tmp,'--min-seq-id',threshold,'-c','0.8','--cov-mode','0','--cluster-mode','2','--threads','2'])
        command(key+'_members',[mm,'createtsv',db,db,clusters,folder/'members.tsv','--threads','2'])
        command(key+'_representatives',[mm,'createsubdb',clusters,db,representatives])
        repfasta=folder/'representatives.fasta'
        command(key+'_fasta',[mm,'convert2fasta',representatives,repfasta])
        members=set();reps=set()
        with (folder/'members.tsv').open() as f:
            for row in csv.reader(f,delimiter='\t'):
                assert len(row)==2 and row[1] not in members
                reps.add(row[0]);members.add(row[1])
        assert len(members)==candidate_count and reps<=members
        cluster_counts[threshold]=dict(candidates=candidate_count,clusters=len(reps),assigned_members=len(members))
        stages[key+'_representatives']=repfasta
        stages[key+'_selected']=select(key+'_parsimony',repfasta)
        if threshold=='0.99':
            assert sequences(repfasta)==sequences(t['prior_h04']+'/cluster_099/representatives.fasta'),'Pilot 99% representatives drifted'
            assert sequences(stages[key+'_selected'])==sequences(t['prior_h04']+'/cluster_099_parsimony/'+t['sample']+'_razor.fasta'),'Pilot 99% selected sequences drifted'
    # Fixed reference and FL peptide identities, independently established by MM A10.
    pep_sets={name:set(Path(path).read_text().splitlines()) for name,path in t['peptide_sets'].items()}
    union=set.union(*pep_sets.values());auto=ahocorasick.Automaton()
    for p in sorted(union):auto.add_word(p.replace('I','L'),p)
    auto.make_automaton();observed={};counts={}
    for name,path in stages.items():
        found=set();n=0
        for _,seq in fasta(path):
            n+=1
            for _,p in auto.iter(seq.replace('I','L')):found.add(p)
        observed[name]=found;counts[name]=n
    rows=[]
    with (out/'peptide_availability.tsv.gz').open('wb') as unused:pass
    import gzip
    with gzip.open(out/'peptide_availability.tsv.gz','wt') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['peptide',*pep_sets,*stages])
        for p in sorted(union):w.writerow([p,*[int(p in v) for v in pep_sets.values()],*[int(p in observed[name]) for name in stages]])
    for universe,peps in pep_sets.items():
        for name,found in observed.items():rows.append(dict(universe=universe,stage=name,n_proteins=counts[name],n_query_peptides=len(peps),available=len(peps&found),lost_relative_to_candidates=len((peps&observed['candidate_unclustered'])-found)))
    with (out/'summary.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    # Never label this as comparative sensitivity/FDR or an upstream-reservoir test.
    result=dict(status='PASS',sample=t['sample'],baseline_sequence_set_reproduced=True,cluster_counts=cluster_counts,summary=rows,scope='Clustering acquisition candidate proteins AFTER identical broad-reservoir lookup, followed by unchanged razor selection. Exact peptide substring availability, not digest-compatible search yield, strain correctness or whole-workflow FDR. Does not test clustering the initial reservoir or compare MetaNovo end to end.',outputs=[dict(path=str(p),sha256=sha(p)) for p in [out/'summary.tsv',out/'peptide_availability.tsv.gz',out/'SOURCE_HASHES.json',out/'COMMANDS.jsonl',*[p for name,p in stages.items() if name!='candidate_unclustered'],*sorted(out.glob('cluster_*/members.tsv'))]])
    (out/'COMPLETE.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':run(Path(sys.argv[1]),int(sys.argv[2]))
