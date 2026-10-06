"""Non-circular check of each native MGnify eggNOG row: does the protein sequence the row is attached to
(UHGG_all.faa record with that gene ID) actually align to the row's own seed ortholog (eggNOG 5.0.2)?
Per genome: DIAMOND blastp of the genome's proteins against the genome's claimed seed set.
Optional 'control' mode plants a negative control by permuting seeds among rows (derangement)."""
import csv,gzip,json,os,subprocess,sys,tempfile,random
from collections import defaultdict
from pathlib import Path
W=Path(__file__).resolve().parents[1];V=W/'verify'
DIAMOND='LAB_ROOT/diamond'
prep=V/sys.argv[1];i=int(sys.argv[2]);mode=sys.argv[3] if len(sys.argv)>3 else 'real'
T=int(os.environ.get('SLURM_CPUS_PER_TASK','4'))
out=V/(sys.argv[1]+'_'+mode+'_out');out.mkdir(exist_ok=True)
dest=out/f'shard_{i:03d}.tsv.gz';assert not dest.exists()
def fasta(p):
    hdr=None;buf=[]
    with gzip.open(p,'rt') as f:
        for line in f:
            if line[0]=='>':
                if hdr:yield hdr,''.join(buf)
                hdr=line[1:].strip();buf=[]
            else:buf.append(line.strip())
    if hdr:yield hdr,''.join(buf)
rows=defaultdict(list)
with gzip.open(prep/'native_rows.tsv.gz','rt') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        if int(r['shard'])==i:rows[r['genome']].append(r)
if mode=='control':
    rng=random.Random(20261003)
    for g,rs in rows.items():
        seeds=[r['seed_ortholog'] for r in rs];perm=seeds[1:]+seeds[:1] if len(seeds)>1 else seeds
        idx=list(range(len(rs)));rng.shuffle(idx)
        sh=[seeds[j] for j in idx];sh=sh[1:]+sh[:1]
        for k,j in enumerate(idx):rs[j]=dict(rs[j],seed_ortholog=sh[k],native_seed=rs[j]['seed_ortholog'])
q={a:s for a,s in fasta(prep/f'query_{i:03d}.faa.gz')}
s={a:x for a,x in fasta(prep/f'seeds_{i:03d}.faa.gz')}
stat=defaultdict(int)
with gzip.open(dest,'wt') as fo, tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR','/tmp')) as td:
    w=csv.writer(fo,delimiter='\t')
    w.writerow(['genome','accession','seed_ortholog','native_score','native_evalue','pair_bitscore','pair_evalue','pair_pident','pair_qcov','pair_scov','best_seed_in_genome_set','best_seed_bitscore','query_length','seed_length','status'])
    for g in sorted(rows):
        rs=rows[g];accs=sorted({r['accession'] for r in rs});seeds=sorted({r['seed_ortholog'] for r in rs if r['seed_ortholog'] in s})
        qa=Path(td)/'q.faa';sa=Path(td)/'s.faa';db=Path(td)/'s';hits=Path(td)/'h.tsv'
        with open(qa,'w') as f:
            for a in accs:
                if a in q:f.write(f'>{a}\n{q[a]}\n')
        pair={};best={}
        if seeds:
            with open(sa,'w') as f:
                for x in seeds:f.write(f'>{x}\n{s[x]}\n')
            subprocess.run([DIAMOND,'makedb','--in',str(sa),'-d',str(db),'--quiet','--threads',str(T)],check=True)
            subprocess.run([DIAMOND,'blastp','-q',str(qa),'-d',str(db),'-o',str(hits),'--ultra-sensitive','--evalue','10','--max-target-seqs','0',
                '--outfmt','6','qseqid','sseqid','bitscore','evalue','pident','qcovhsp','scovhsp','--threads',str(T),'--quiet','--tmpdir',td],check=True)
            with open(hits) as f:
                for line in f:
                    a,x,b,e,p,qc,sc=line.rstrip('\n').split('\t');b=float(b)
                    if (a,x) not in pair or b>pair[(a,x)][0]:pair[(a,x)]=(b,e,p,qc,sc)
                    if a not in best or b>best[a][1]:best[a]=(x,b)
        for r in rs:
            a=r['accession'];x=r['seed_ortholog'];ns=float(r['score'])
            if a not in q:st='query_sequence_absent'
            elif x not in s:st='seed_absent_from_eggnog_db'
            else:
                h=pair.get((a,x))
                st='pair_no_alignment' if h is None else 'pair_aligned'
            h=pair.get((a,x),('', '', '', '', ''));bs=best.get(a,('',''))
            w.writerow([g,a,x,ns,r['evalue'],*h,bs[0],bs[1],len(q.get(a,'')),len(s.get(x,'')),st]);stat[st]+=1
        print(g,len(rs),dict(stat),flush=True)
json.dump(dict(status='PASS',shard=i,mode=mode,counts=stat,job=os.environ.get('SLURM_JOB_ID')),open(out/f'shard_{i:03d}.COMPLETE.json','w'))
