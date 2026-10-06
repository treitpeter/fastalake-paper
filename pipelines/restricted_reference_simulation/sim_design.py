"""RES-47k: design values of the restricted simulation, read from each validation mixture's truth table and GENERATED.json.
usage: python sim_design.py RUNS_VALIDATION_DIR OUT.tsv"""
import csv,json,sys,os
root,out=sys.argv[1:3]
rows=[]
for u in sorted(d for d in os.listdir(root) if d[0]=='u' and d[1:].isdigit()):
    g=json.load(open(f'{root}/{u}/GENERATED.json')); unit=g['unit']
    t=list(csv.DictReader(open(f'{root}/{u}/truth/mixture.tsv'),delimiter='\t'))
    L=[len(r['sequence']) for r in t]; om=sum(r['source_kind']=='omitted' for r in t)
    rows.append(dict(mixture=u,regime=unit.get('regime'),seed=unit.get('seed'),generating_peptides=len(t),
      absent_from_search_reference=om,fraction_absent=round(om/len(t),4),min_length=min(L),max_length=max(L),
      reference_proteins=g['reservoir_proteins']+g['omitted_proteins'],search_reference_proteins=g['reservoir_proteins'],
      omitted_proteins=g['omitted_proteins'],slurm_job=g.get('slurm_job')))
with open(out,'w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
import collections
s=collections.defaultdict(list)
for r in rows: s[r['regime']].append(r)
for k,v in s.items(): print(k,len(v),'absent',sorted({r['fraction_absent'] for r in v}),'len',min(r['min_length'] for r in v),max(r['max_length'] for r in v),'ref',sorted({r['reference_proteins'] for r in v}))
