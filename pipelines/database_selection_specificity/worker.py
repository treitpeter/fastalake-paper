#!/usr/bin/env python3
import csv,json,os,pathlib,random,shutil,sys,time
from collections import Counter
from lib import *

RUN=pathlib.Path(sys.argv[1]).resolve()
M=json.loads((RUN/'inputs/manifest.json').read_text())
STAGE=sys.argv[2]
def check():
    verify(M['authentic']['frozen']);verify(M['archived_scramble']['frozen'])
    verify(M['extractor']['frozen']);verify(M['composition_binary']);verify(M['reservoir'])
    for r in M['code']:verify(r['frozen'])
def extractor(csvdir,db,out,threads):
    out=pathlib.Path(out);out.mkdir(parents=True,exist_ok=True)
    command([M['extractor']['frozen']['path'],'--peptides-dir',csvdir,'--database',db,'--output',out/'matches.fasta','--output-hits',out/'hits.tsv','--top-percent','0.30','--min-length','9','--max-length','50','--normalize-il','true','-t',str(threads)],out/'extractor')
    return out/'hits.tsv',json.loads((out/'matches.stats.json').read_text())
def shuffle(seed,dst):
    rng=random.Random(seed);same=equivalent=n=0;rowcheck=hashlib.sha256()
    with open(M['authentic']['frozen']['path'],newline='') as fi,open(dst,'w',newline='') as fo:
        reader=csv.DictReader(fi);writer=csv.DictWriter(fo,fieldnames=reader.fieldnames);writer.writeheader()
        for row in reader:
            original=row.get(COL) or '';idx=[i for i,c in enumerate(original) if 'A'<=c<='Z'];res=[original[i] for i in idx];raw=''.join(res)
            if idx:
                for _ in range(100):
                    rng.shuffle(res)
                    if ''.join(res)!=raw:break
                else:same+=1
                assert Counter(res)==Counter(raw)
                s=list(original)
                for i,c in zip(idx,res):s[i]=c
                row[COL]=''.join(s)
                equivalent+=int(''.join(res).replace('I','L')==raw.replace('I','L'))
            writer.writerow(row);n+=1
    # Re-read both tables to check every row and every unchanged field, not a sample.
    with open(M['authentic']['frozen']['path'],newline='') as fi,open(dst,newline='') as fo:
        from itertools import zip_longest
        for a,b in zip_longest(csv.DictReader(fi),csv.DictReader(fo)):
            assert a is not None and b is not None
            assert Counter(clean(a[COL]))==Counter(clean(b[COL]))
            assert {k:v for k,v in a.items() if k!=COL}=={k:v for k,v in b.items() if k!=COL}
    return dict(seed=seed,n_rows=n,unchanged_order=same,unchanged_after_il=equivalent,csv=meta(dst))
def prepare():
    d=RUN/'derived';pool=d/'pool';pool.mkdir(exist_ok=True)
    arms=[];sources=[('REAL',None,pathlib.Path(M['authentic']['frozen']['path']))]
    for seed in M['seeds']:
        arm='SIM_'+str(seed);directory=d/'predictions'/arm;directory.mkdir(parents=True,exist_ok=True)
        p=directory/(M['sample']+'_predictions.csv')
        info=shuffle(seed,p);dump(directory/'SHUFFLE_AUDIT.json',info)
        if seed==20260901:assert sha(p)==M['archived_scramble']['frozen']['sha256'],'Historical shuffle not byte-reproduced'
        sources.append((arm,seed,p))
    realdir=d/'predictions/REAL';realdir.mkdir(parents=True,exist_ok=True)
    shutil.copy2(M['authentic']['frozen']['path'],str(realdir/(M['sample']+'_predictions.csv')))
    sources[0]=('REAL',None,realdir/(M['sample']+'_predictions.csv'))
    combined=d/'all_clean_sequences.fasta'
    all_queries=set()
    with combined.open('wb') as combined_out:
        for arm,seed,p in sources:
            target=d/(arm+'_all_clean.fasta');counts,summary,retained=load_queries(p,target)
            all_queries.update(counts)
            table(d/(arm+'_queries.tsv'),[dict(query=q,length=len(q),retained_rows=n) for q,n in sorted(counts.items())])
            table(d/(arm+'_retained_rows.tsv'),[dict(row=i,query=q,length=len(q),score=s) for i,q,s in retained])
            hits,stats=extractor(p.parent,target,d/'validation'/arm,2)
            observed,_,_=hit_counts(hits)
            assert set(observed)==set(counts),arm+' query reconstruction differs from production'
            assert stats['total_peptides']==len(counts),arm+' production query count differs'
            with target.open('rb') as f:shutil.copyfileobj(f,combined_out,8*1024*1024)
            link=pool/(arm+'_predictions.csv');link.symlink_to(p)
            arms.append(dict(arm=arm,seed=seed,csv=meta(p),queries=meta(d/(arm+'_queries.tsv')),summary=summary,production_fixture_query_set_agrees=True))
    hits,stats=extractor(pool,combined,d/'validation/POOL',2)
    observed,_,_=hit_counts(hits)
    assert set(observed)==all_queries and stats['total_peptides']==len(all_queries),'Pooled per-file filtering differs'
    dump(d/'PREPARED.json',dict(status='PASS',arms=arms,union_queries=len(all_queries),pooled_fixture_agrees=True))
def lookup():
    p=json.loads((RUN/'derived/PREPARED.json').read_text());assert p['status']=='PASS'
    for a in p['arms']:verify(a['csv']);verify(a['queries'])
    hits,stats=extractor(RUN/'derived/pool',M['reservoir']['path'],RUN/'derived/lookup',24)
    assert stats['proteins_searched']==233497702 and stats['total_peptides']==p['union_queries']
    dump(RUN/'derived/LOOKUP.json',dict(status='PASS',hits=meta(hits),statistics=meta(RUN/'derived/lookup/matches.stats.json'),candidate_fasta=meta(RUN/'derived/lookup/matches.fasta'),stats=stats))
def composition():
    prefix=RUN/'derived/reservoir_composition'
    command([M['composition_binary']['path'],M['reservoir']['path'],prefix],RUN/'logs/composition_scan')
    data=json.loads(prefix.with_suffix('.json').read_text())
    assert data['entries']==233497702
    data['source_fasta']=meta(M['reservoir']['path'])
    data['method']='fresh full FASTA residue census; uppercase residues; protein boundaries retained'
    dump(prefix.with_suffix('.json'),data)
check()
if STAGE=='prepare':prepare()
elif STAGE=='lookup':lookup()
elif STAGE=='composition':composition()
elif STAGE=='report':
    from report import build
    build(RUN,M)
else:raise ValueError(STAGE)
check()
dump(RUN/'status'/(STAGE+'.json'),dict(status='PASS',stage=STAGE,job_id=os.environ.get('SLURM_JOB_ID'),utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),manifest_sha256=sha(RUN/'inputs/manifest.json')))
print('VALIDATED',STAGE,flush=True)

