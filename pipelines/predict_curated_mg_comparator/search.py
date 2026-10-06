"""predict_curated_mg_comparator: search. Moved from the original analysis as it ran, apart from path scrubbing; see README.md."""
from pathlib import Path
import csv, hashlib, json, os, re, subprocess, sys, time
B=Path(__file__).resolve().parent
csv.field_size_limit(64*1024*1024)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()
def stamp(p):
    p=Path(p).resolve();s=p.stat();return dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns,sha256=sha(p))
def verify(r):
    a=stamp(r['path']);assert all(a[k]==v for k,v in r.items() if k in a),r['path']
def dump(p,j):Path(p).write_text(json.dumps(j,indent=2)+'\n')
def run(i):
    m=json.loads((B/'manifest.json').read_text());assert sha(B/'manifest.json')==(B/'manifest.sha256').read_text().strip()
    assert sha(__file__)==m['worker_sha256']
    t=m['tasks'][i]
    if not t['pilot']:
        for k in range(4):assert json.loads((B/'results'/str(k)/'COMPLETE.json').read_text())['status']=='PASS'
    for rec in [m['sage'],t['fasta'],t['raw'],*t['provenance']]:verify(rec)
    out=B/'results'/str(i);out.mkdir(parents=True,exist_ok=False)
    c=t['config'];c['database']['fasta']=t['fasta']['path'];c['mzml_paths']=[t['raw']['path']];c['output_directory']=str(out)
    dump(out/'config.json',c)
    argv=[m['sage']['path'],str(out/'config.json'),'--batch-size','1','--write-pin','--disable-telemetry-i-dont-want-to-improve-sage']
    dump(out/'COMMAND.json',dict(argv=argv,job_id=os.environ.get('SLURM_JOB_ID'),index=i))
    start=time.monotonic()
    with (out/'stdout.log').open('w') as a,(out/'stderr.log').open('w') as z:
        p=subprocess.run(['/usr/bin/time','-v','-o',str(out/'time.txt'),*argv],stdout=a,stderr=z,env=dict(os.environ,RAYON_NUM_THREADS='8'))
    assert p.returncode==0,p.returncode
    err=(out/'stderr.log').read_text();assert not re.search(r'error parsing|XML.*error|panicked at|premature end',err,re.I),err[-1500:]
    actual=json.loads((out/'results.json').read_text());assert actual['mzml_paths']==c['mzml_paths'] and actual['database']['fasta']==c['database']['fasta']
    assert actual['report_psms']==c['report_psms'] and actual['database']['enzyme']['missed_cleavages']==2
    assert actual['version']=='0.14.6'
    sets=[set(),set()];rows=0;keys=set()
    with (out/'results.sage.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            rows+=1;assert Path(r['filename']).name==Path(t['raw']['path']).name
            k=(r['scannr'],r['rank']);assert k not in keys;keys.add(k)
            assert 1<=int(r['rank'])<=c['report_psms']
            if r['label']=='1' and float(r['peptide_q'])<=.01:
                peptide=re.sub(r'\[[^\]]*\]','',r['peptide']);sets[0].add(peptide);sets[1].add(peptide.replace('I','L'))
    assert rows and sets[0]
    for name,s in zip(['accepted_exact','accepted_canonical'],sets):(out/(name+'.txt')).write_text(''.join(p+'\n' for p in sorted(s)))
    for rec in [t['fasta'],t['raw'],*t['provenance']]:verify(rec)
    audit=dict(status='PASS',index=i,sample=t['sample'],arm=t['arm'],pilot=t['pilot'],seconds=time.monotonic()-start,
        rows=rows,accepted_exact=len(sets[0]),accepted_canonical=len(sets[1]),config=stamp(out/'config.json'),
        task=t,outputs=[stamp(out/n) for n in ['results.json','results.sage.tsv','accepted_exact.txt','accepted_canonical.txt']],
        job_id=os.environ.get('SLURM_JOB_ID'),scientific_review='PENDING')
    if t.get('old_peptides'):
        import gzip
        old={re.sub(r'\[[^\]]*\]','',x.strip()) for x in gzip.open(t['old_peptides'],'rt') if x.strip()}
        audit['historical_identity_difference']=dict(gained=len(sets[0]-old),lost=len(old-sets[0]))
        (out/'historical_gained.txt').write_text('\n'.join(sorted(sets[0]-old)));(out/'historical_lost.txt').write_text('\n'.join(sorted(old-sets[0])))
    dump(out/'COMPLETE.json',audit);print('SUPPLEMENT6 SEARCH PASS',i,t['sample'],t['arm'],len(sets[1]),flush=True)
if __name__=='__main__':run(int(sys.argv[1]))
