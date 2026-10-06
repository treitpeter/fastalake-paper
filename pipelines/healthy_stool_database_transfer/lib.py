import csv,hashlib,json,math,pathlib,subprocess,time
from collections import Counter

COL='peptide_prediction_detokenized_unmodified'
PYTHON='/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python'
def sha(p):
    h=hashlib.sha256()
    with pathlib.Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def meta(p,digest=True):
    p=pathlib.Path(p).resolve();s=p.stat()
    r=dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
    if digest:r['sha256']=sha(p)
    return r
def verify(r):
    assert meta(r['path'],False)=={k:r[k] for k in ['path','bytes','mtime_ns']},'Changed input '+r['path']
    if 'sha256' in r:assert sha(r['path'])==r['sha256'],'Changed content '+r['path']
def dump(p,obj):
    p=pathlib.Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def table(p,rs,fields=None):
    with pathlib.Path(p).open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
def clean(s):return ''.join(c for c in s if 'A'<=c<='Z')
def low(seq):
    n=len(seq)
    if n==0:return True
    c=Counter(seq)
    return max(c.values())/n>0.6 or (n>=12 and len(c)<=2) or (n>=10 and all(seq[i]==seq[i%2] for i in range(n)))
def load_queries(path,fixture=None):
    rows=[];reasons=Counter();total=0
    fh=open(fixture,'w') if fixture else None
    try:
        with open(path,newline='') as f:
            for i,r in enumerate(csv.DictReader(f),1):
                total+=1;s=clean(r.get(COL) or '')
                if fh and s:fh.write('>row_'+str(i)+'\n'+s+'\n')
                try:score=float(r.get('score') or '')
                except ValueError:reasons['invalid_or_missing_score']+=1;continue
                if not 9<=len(s)<=50:reasons['length']+=1;continue
                if not score>0:reasons['nonpositive_or_nan_score']+=1;continue
                if low(s):reasons['low_complexity']+=1;continue
                rows.append((i,s,score))
    finally:
        if fh:fh.close()
    assert rows,'No eligible queries'
    scores=sorted((r[2] for r in rows),reverse=True)
    threshold=scores[min(len(scores),max(1,math.ceil(len(scores)*.30)))-1]
    kept=[(i,s.replace('I','L'),score) for i,s,score in rows if score>=threshold and not low(s)]
    counts=Counter(s for _,s,_ in kept)
    summary=dict(total_rows=total,eligible_before_top_fraction=len(rows),retained_rows=len(kept),distinct_queries=len(counts),score_threshold=threshold,filter_attrition=dict(reasons))
    return counts,summary,kept
def command(argv,prefix):
    prefix=pathlib.Path(prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
    dump(str(prefix)+'.command.json',dict(argv=list(map(str,argv))))
    with open(str(prefix)+'.out','w') as f,open(str(prefix)+'.err','w') as e:
        subprocess.run(['/usr/bin/time','-v',*map(str,argv)],stdout=f,stderr=e,check=True)
def hit_counts(path):
    # Production hits are contiguous by reservoir entry. Suppress repeated positions
    # within a protein for the supporting entry-multiplicity count.
    occurrences=Counter();entries=Counter();previous=None;seen=set();blocks=0
    with open(path) as f:
        assert f.readline().rstrip('\n')=='peptide\tprotein_id\tposition'
        for line in f:
            q,protein,pos=line.rstrip('\n').split('\t')
            if protein!=previous:
                previous=protein;seen=set();blocks+=1
            occurrences[q]+=1
            if q not in seen:entries[q]+=1;seen.add(q)
    return occurrences,entries,blocks

