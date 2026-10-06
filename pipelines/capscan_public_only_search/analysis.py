#!/usr/bin/env python3
"""capscan_public_only_search: analysis. Moved from the original analysis as it ran, apart from path scrubbing; see README.md."""
import csv,gzip,hashlib,json,os,pathlib,shutil,sys,time
from collections import Counter
RUN=pathlib.Path(sys.argv[1]).resolve()
sys.path.insert(0,str(RUN/'frozen'))
from lib import meta,sha,verify,dump,table,command
from search_helpers import recount,read_fasta
def raw_records(path):
    h=None;parts=[]
    with open(path) as f:
        for l in f:
            if l.startswith('>'):
                if h is not None:yield h,''.join(parts)
                h=l[1:].split()[0];parts=[]
            else:parts.append(l.strip().upper())
    if h is not None:yield h,''.join(parts)
def extract(m,pool,db,out,threads):
    command([m['extractor']['path'],'--peptides-dir',pool,'--database',db,'--output',out/'selected.fasta','--output-hits',out/'hits.tsv','--top-percent','0.30','--min-length','9','--max-length','50','--normalize-il','true','-t',str(threads)],out/'extract')
    stats=json.loads((out/'selected.stats.json').read_text())
    assert stats['peptide_length_range']==[9,50] and stats['score_filter']=='top_30%_per_sample'
    return stats
def main():
    stage=sys.argv[2];m=json.loads((RUN/'inputs/manifest.json').read_text())
    for x in m['frozen_code']+[m['extractor'],m['parsimony'],m['sage']]:verify(x)
    if stage=='evidence':
        out=RUN/'evidence';out.mkdir(exist_ok=True);pool=out/'predictions';pool.mkdir(exist_ok=True)
        stamps=[]
        for r in m['samples']:
            verify(r['prediction']);stamps.append(meta(r['prediction']['path']))
            link=pool/(r['sample']+'_predictions.csv')
            if not link.exists():link.symlink_to(r['prediction']['path'])
        dump(out/'prediction_hashes.json',stamps)
        verify(m['reservoir']);stats=extract(m,pool,m['reservoir']['path'],out,8)
        assert stats['proteins_searched']==233497702
        # Selected sequence hashes must occur in the public-catalogue source ledger.
        selected={hashlib.sha256(s.encode()).hexdigest() for _,s in raw_records(out/'selected.fasta')}
        remaining=set(selected);counts=Counter();n=0;prev='';h=hashlib.sha256()
        with open(m['source_headers']['path'],'rb') as f:
            line=f.readline();h.update(line)
            assert line.startswith(b'sha256_hash\t')
            for line in f:
                h.update(line);fields=line.rstrip(b'\n').split(b'\t',3)
                assert len(fields)==4,(n,fields)
                key=fields[0].decode();assert len(key)==64 and key>prev;prev=key
                records=fields[3].split(b'|||');assert len(records)==int(fields[1])
                for src in records:
                    label=src.split(b'\t',1)[0].decode()
                    assert label in ['UHGG','GMGC','GMSC','HUMAN_SP'],('Non-public source',label)
                    counts[label]+=1
                remaining.discard(key);n+=1
        assert n==233497702 and not remaining,('Selected sequences missing from public ledger',len(remaining))
        for v in stamps:verify(v)
        verify(m['reservoir']);verify(m['source_headers'])
        dump(out/'COMPLETE.json',dict(status='PASS',reservoir=meta(m['reservoir']['path']),header_source=dict(m['source_headers'],sha256=h.hexdigest()),source_counts=dict(counts),ledger_sequences=n,selected_unique_sequences=len(selected),all_selected_sequences_publicly_provenanced=True,stats=stats,selected=meta(out/'selected.fasta'),predictions=meta(out/'prediction_hashes.json'),outputs=[meta(out/'hits.tsv')]))
        print('PUBLIC EVIDENCE VALIDATED',len(selected),flush=True);return
    index=int(sys.argv[3]);r=m['samples'][index]
    complete=json.loads((RUN/'evidence/COMPLETE.json').read_text());assert complete['status']=='PASS'
    verify(complete['selected']);verify(r['prediction']);verify(r['raw']);verify(r['config'])
    out=RUN/'results'/r['sample'];out.mkdir(parents=True,exist_ok=True)
    assert not (out/'COMPLETE.json').exists()
    selection=out/'selection';selection.mkdir(exist_ok=True);pool=selection/'predictions';pool.mkdir(exist_ok=True)
    link=pool/(r['sample']+'_predictions.csv')
    if not link.exists():link.symlink_to(r['prediction']['path'])
    stat=extract(m,pool,complete['selected']['path'],selection,8)
    par=out/'parsimony';par.mkdir(exist_ok=True)
    command([m['parsimony']['path'],'--strategy','razor','--tiebreak','hash-acc','--fasta',selection/'selected.fasta','--peptides',r['prediction']['path'],'--output-dir',par,'--sample',r['sample'],'--min-length','9','--max-length','50','--normalize-il','-t','8'],par/'selection')
    db=par/(r['sample']+'_razor.fasta');assert db.exists()
    config=json.loads(pathlib.Path(r['config']['path']).read_text())
    for k in ['version','output_paths']:config.pop(k,None)
    config['database']['fasta']=str(db);config['mzml_paths']=[r['raw']['path']];config['output_directory']=str(out/'search')
    cp=out/'sage.json';dump(cp,config)
    command([m['sage']['path'],cp,'--batch-size','1','--write-pin','--disable-telemetry-i-dont-want-to-improve-sage'],out/'sage')
    z,sets,psms=recount(out/'search/results.sage.tsv')
    selected=read_fasta(db)
    for row in psms:
        parents=[p for p in row['proteins'].split(';') if not p.startswith('rev_')]
        assert any(row['canonical'] in selected[p] for p in parents)
    for name,peps in sets.items():(out/(name+'.txt')).write_text(''.join(p+'\n' for p in sorted(peps)))
    baseline={}
    for label,key in [('cohort_metagenomes_available','baseline'),('specimen_metagenome','metagenome')]:
        if r.get(key):
            verify(r[key]);_,old,_=recount(r[key]['path'])
            a,b=sets['peptide_q01'],old['peptide_q01']
            baseline[label]=dict(new=len(a),old=len(b),intersection=len(a&b),gained=len(a-b),lost=len(b-a),comparator_recovery=len(a&b)/len(b) if b else None)
            for tag,v in [('gained',a-b),('lost',b-a)]:(out/(label+'_'+tag+'.txt')).write_text(''.join(p+'\n' for p in sorted(v)))
    dump(out/'COMPLETE.json',dict(status='PASS',sample=r['sample'],public_source_proof=meta(RUN/'evidence/COMPLETE.json'),prediction=meta(r['prediction']['path']),raw=meta(r['raw']['path']),selected_fasta=meta(db),selected_proteins=len(selected),counts={k:len(v) for k,v in sets.items()},comparisons=baseline,job_id=os.environ.get('SLURM_JOB_ID'),outputs=[meta(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='COMPLETE.json']+[meta(p) for p in (out/'search').iterdir() if p.is_file()]))
    print('PUBLIC SAMPLE VALIDATED',r['sample'],len(sets['peptide_q01']),flush=True)
if __name__=='__main__':main()
