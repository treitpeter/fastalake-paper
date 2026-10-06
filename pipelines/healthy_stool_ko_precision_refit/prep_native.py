"""Prepare the non-circular native UHGG eggNOG verification.
Inputs (read-only): the 4,744 MGnify native eggNOG TSVs listed in evidence_all_cohorts/INPUTS.json,
UHGG_all.faa (the reference whose sequences the native rows are attached to), and the eggNOG 5.0.2
DIAMOND database used by the manuscript mapper. Outputs: per-shard query/seed FASTA and native rows."""
import csv,gzip,hashlib,json,sys,subprocess,os,re
from pathlib import Path
csv.field_size_limit(1<<30)
W=Path(__file__).resolve().parents[1];V=W/'verify'
E=Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/github_workflow_reconciliation_2026-09-30/resumed_audit/evidence_all_cohorts')
INP=json.loads((E/'INPUTS.json').read_text())
DMND='LAB_ROOT/eggnog_proteins.dmnd'
DIAMOND='LAB_ROOT/diamond'
NS=int(sys.argv[1]);TEST=len(sys.argv)>2 and sys.argv[2]=='test'
out=V/('prep_test' if TEST else 'prep');out.mkdir(parents=True,exist_ok=False)
genomes=sorted(INP['native_eggnog'])
if TEST:genomes=['MGYG000000543','MGYG000000644']
shard={g:i%NS for i,g in enumerate(genomes)}
seed_shards={};n=0;dup=0;stats={}
with gzip.open(out/'native_rows.tsv.gz','wt') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['genome','shard','accession','seed_ortholog','evalue','score','row_sha256'])
    for g in genomes:
        seen={}
        with open(INP['native_eggnog'][g]) as h:
            for r in csv.DictReader(h,delimiter='\t'):
                q=r['#query'];assert q.startswith(g+'_'),(g,q)
                body='\t'.join(r[k] or '' for k in r if k!='#query')
                rs=hashlib.sha256(body.encode()).hexdigest()
                if q in seen:assert seen[q]==rs,('conflicting duplicate',q);dup+=1;continue
                seen[q]=rs
                w.writerow([g,shard[g],q,r['seed_ortholog'],r['evalue'],r['score'],rs]);n+=1
                seed_shards.setdefault(r['seed_ortholog'],set()).add(shard[g])
stats.update(native_rows=n,identical_duplicate_rows_skipped=dup,genomes=len(genomes),unique_seeds=len(seed_shards))
print(stats,flush=True)
# Query sequences: exactly the UHGG_all.faa records the projection attaches native rows to.
qf={i:gzip.open(out/f'query_{i:03d}.faa.gz','wt',compresslevel=2) for i in set(shard.values())}
nq=0;nqs=0
def fasta(p):
    hdr=None;buf=[]
    with open(p) as f:
        for line in f:
            if line[0]=='>':
                if hdr:yield hdr,''.join(buf)
                hdr=line[1:].split()[0];buf=[]
            else:buf.append(line.strip())
    if hdr:yield hdr,''.join(buf)
for acc,seq in fasta(INP['uhgg']):
    nq+=1;g=acc.split('_')[0]
    if g in shard:
        qf[shard[g]].write(f'>{acc}\n{seq.upper()}\n');nqs+=1
for f in qf.values():f.close()
stats.update(uhgg_records=nq,uhgg_records_in_shards=nqs)
print(stats,flush=True)
# Seed sequences from the eggNOG 5.0.2 DIAMOND database (all records dumped, then filtered).
sf={i:gzip.open(out/f'seeds_{i:03d}.faa.gz','wt',compresslevel=2) for i in set(shard.values())}
p=subprocess.Popen([DIAMOND,'getseq','--db',DMND,'--quiet'],stdout=subprocess.PIPE,text=True,bufsize=1<<20)
found=set();hdr=None;buf=[];ndb=0
def emit(h,b):
    if h in seed_shards:
        s=''.join(b);found.add(h)
        for i in seed_shards[h]:sf[i].write(f'>{h}\n{s}\n')
for line in p.stdout:
    if line[0]=='>':
        if hdr:emit(hdr,buf)
        hdr=line[1:].split()[0];buf=[];ndb+=1
    else:buf.append(line.strip())
if hdr:emit(hdr,buf)
assert p.wait()==0
for f in sf.values():f.close()
missing=sorted(set(seed_shards)-found)
(out/'seeds_missing_from_db.txt').write_text('\n'.join(missing)+'\n')
stats.update(eggnog_db_records=ndb,seeds_found=len(found),seeds_missing=len(missing))
json.dump(dict(status='PASS',shards=NS,test=TEST,job=os.environ.get('SLURM_JOB_ID'),**stats),open(out/'COMPLETE.json','w'),indent=1)
print(stats,flush=True)
