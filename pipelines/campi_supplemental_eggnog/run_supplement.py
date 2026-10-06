"""Annotate only never-queried reference proteins with frozen manuscript tools."""
from pathlib import Path
import csv,gzip,hashlib,json,os,subprocess,sys
from collections import Counter
N=Path(__file__).resolve().parent
sys.path.insert(0,str(N))
from extract_native import annotation_rows,fasta,digest,dump

def verify(mode):
    spec=json.loads((N/'SUPPLEMENT_RUN_MANIFEST.json').read_text())
    for p,sha in spec['code'].items():assert digest(p)==sha,p
    hashes={}
    for p,r in spec['database_files'].items():
        st=Path(p).stat();assert st.st_size==r['bytes'] and st.st_mtime_ns==r['mtime_ns'],p
        if mode=='pilot':hashes[p]=digest(p)
    return spec,hashes

def run(mode,index=None):
    spec,hashes=verify(mode);root=N/'supplemental_functions';prep=json.loads((root/'PREPARE_COMPLETE.json').read_text());assert prep['status']=='PASS'
    if mode=='full':
        assert json.loads((root/'pilot/COMPLETE.json').read_text())['status']=='PASS'
        label=f'chunk_{index:03d}'
    else:label='pilot'
    query=root/(label+'.faa');assert digest(query)==prep['outputs'][query.name]
    out=root/(label if mode=='pilot' else 'runs/'+label);out.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ);env['PATH']=spec['path'];env['OMP_NUM_THREADS']='1';env['OPENBLAS_NUM_THREADS']='1'
    command=[spec['python'],spec['mapper'],'-i',str(query),'-o','annotation','--output_dir',str(out),'--data_dir',spec['database_dir'],'-m','diamond','--cpu','4']
    dump(out/'COMMAND.json',dict(argv=command,PATH=env['PATH'],database_hashes_from_pilot=hashes))
    print('COMMAND',json.dumps(command),flush=True)
    subprocess.run(command,env=env,check=True)
    sequences={h.split()[0]:seq for h,seq in fasta(query)}
    assert all(q=='FLSEQ_'+hashlib.sha256(s.encode()).hexdigest() for q,s in sequences.items())
    found=set();positive=Counter()
    for r in annotation_rows(out/'annotation.emapper.annotations'):
        q=r['#query'];assert q in sequences and q not in found;found.add(q)
        for k in ['KEGG_ko','eggNOG_OGs','KEGG_Pathway']:positive[k]+=r[k] not in ('','-')
    if mode=='pilot':assert found, 'No pilot annotation records: inspect runtime/database'
    # Database stat guard after execution; code and data are shared read-only.
    verify('full')
    dump(out/'COMPLETE.json',dict(status='PASS',mode=mode,label=label,queries=len(sequences),annotation_records=len(found),
        positive_fields=dict(positive),query_sha256=digest(query),database_sha256=hashes,
        outputs={p.name:digest(p) for p in out.iterdir() if p.is_file()},job_id=os.environ.get('SLURM_JOB_ID'),
        interpretation='Successful functional annotation execution. No hit or absent KO is preserved; not evidence of biological absence.'))
    print('SUPPLEMENT COMPLETE',label,len(sequences),len(found),dict(positive),flush=True)

def collect():
    spec,_=verify('full');root=N/'supplemental_functions';prep=json.loads((root/'PREPARE_COMPLETE.json').read_text());pilot=json.loads((root/'pilot/COMPLETE.json').read_text());assert pilot['status']=='PASS'
    out=N/'supplemental_native_extract';out.mkdir(exist_ok=False);seen=set();total=annotated=0;receipts={}
    with gzip.open(out/'native_sequences.tsv.gz','wt',compresslevel=3) as sf,gzip.open(out/'native_functions.tsv.gz','wt',compresslevel=3) as af:
        sw=csv.writer(sf,delimiter='\t');aw=csv.writer(af,delimiter='\t')
        sw.writerow(['accession','sequence_sha256_IL','sequence_verified','source_fasta','original_header','terminal_stop_removed'])
        aw.writerow(['accession','sequence_sha256_IL','KEGG_ko','OG','KEGG_Pathway','source_table','source_fasta'])
        for index in range(prep['chunks']):
            label=f'chunk_{index:03d}';d=root/'runs'/label;r=json.loads((d/'COMPLETE.json').read_text());assert r['status']=='PASS' and r['mode']=='full'
            receipts[label]=digest(d/'COMPLETE.json')
            for name,sha in r['outputs'].items():assert digest(d/name)==sha,(label,name)
            query=root/(label+'.faa');assert digest(query)==prep['outputs'][query.name]
            seqs={h.split()[0]:s for h,s in fasta(query)};assert len(seqs)==r['queries'];total+=len(seqs)
            assert not set(seqs)&seen;seen.update(seqs)
            for q,s in seqs.items():
                assert q=='FLSEQ_'+hashlib.sha256(s.encode()).hexdigest()
                sh=hashlib.sha256(s.replace('I','L').encode()).hexdigest();sw.writerow([q,sh,1,str(query),q,0])
            for row in annotation_rows(d/'annotation.emapper.annotations'):
                q=row['#query'];sh=hashlib.sha256(seqs[q].replace('I','L').encode()).hexdigest();annotated+=1
                aw.writerow([q,sh,row['KEGG_ko'],row['eggNOG_OGs'],row['KEGG_Pathway'],str(d/'annotation.emapper.annotations'),str(query)])
    assert total==prep['raw_annotation_queries']
    dump(out/'COMPLETE.json',dict(status='PASS',queries=total,annotation_records=annotated,
        outputs={p.name:digest(p) for p in out.iterdir() if p.is_file()},chunk_receipts=receipts,
        source_manifest_sha256=digest(N/'SUPPLEMENT_RUN_MANIFEST.json'),database_sha256=pilot['database_sha256'],
        job_id=os.environ.get('SLURM_JOB_ID'),scope='Functional annotations for never-queried reference sequences. No taxonomy supplied by eggNOG.'))
    print('SUPPLEMENT COLLECTED',total,annotated,flush=True)

if __name__=='__main__':
    if sys.argv[1]=='collect':collect()
    else:run(sys.argv[1],int(sys.argv[2]) if len(sys.argv)>2 else None)
