from pathlib import Path
import sys,json,csv,gzip,statistics,hashlib
from collections import Counter
ROOT=Path.cwd().resolve();O=Path(__file__).resolve().parent.parent;R=O.parent
tasks=list(csv.DictReader((R/'tasks.tsv').open(),delimiter='\t'));data={};summaries=[];missing=[]
for t in tasks:
 p=O/'peptide_tasks'/('%04d'%int(t['index']))
 if not (p/'COMPLETE.json').exists():missing.append(t['index']);continue
 c=json.loads((p/'COMPLETE.json').read_text());assert c['status']=='PASS'
 r=json.loads((p/'SUMMARY.json').read_text());summaries.append(r);data[t['cohort'],t['arm'],t['sample']]=(t,p,r)
D=O/'peptide_summary';D.mkdir(exist_ok=True)
if missing:
 (D/'INCOMPLETE.json').write_text(json.dumps(dict(status='INCOMPLETE',missing=missing),indent=2)+'\n');raise SystemExit('Missing task collectors')
def tab(p,rows):
 with p.open('w',newline='') as f:w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
rows=[]
for r in summaries:
 for pop,n in r['counts'].items():rows.append(dict(index=r['index'],cohort=r['cohort'],arm=r['arm'],sample=r['sample'],population=pop,peptides=n,shared_fresh_nested7=r['fresh_vs_nested7']['shared'],gained_fresh_nested7=r['fresh_vs_nested7']['gained'],lost_fresh_nested7=r['fresh_vs_nested7']['lost']))
tab(D/'per_acquisition_counts.tsv',rows)
def load(k,pop):
 with gzip.open(data[k][1]/(pop+'.txt.gz'),'rt') as f:return set(f.read().split())
a16=json.loads((ROOT/'MM_review_analysis/A16_4_5_public_comparison/inputs/manifest.json').read_text())
pairs=[]
for co in ['PREDICT','CAPSCAN']:
 names={t['sample'] for t in a16['tasks'] if t['cohort']==co and 'MG' in t['configs']}
 assert len(names)==(296 if co=='PREDICT' else 282),(co,len(names))
 for arm in ['FL_public','FL_cohortMG']:
  for name in sorted(names):pairs.append((co,arm,'MG',name))
for arm in ['FL','FL_R','FL_R_M']:
 for comparator in ['REF','MG']:
  for name in sorted(k[2] for k in data if k[:2]==('CAMPI',arm)):pairs.append(('CAMPI',arm,comparator,name))
out=[]
for co,a,b,name in pairs:
 for pop in ['baseline5','nested7','nested9','fresh7','fresh7_nested9']:
  aa,bb=load((co,a,name),pop),load((co,b,name),pop)
  out.append(dict(cohort=co,arm=a,comparator=b,sample=name,population=pop,arm_peptides=len(aa),comparator_peptides=len(bb),shared=len(aa&bb),gain_pct=100*(len(aa)/len(bb)-1) if bb else None,comparator_identity_recovery_pct=100*len(aa&bb)/len(bb) if bb else None))
tab(D/'matched_comparisons.tsv',out)
groups={}
for r in out:groups.setdefault((r['cohort'],r['arm'],r['comparator'],r['population']),[]).append(r)
aggregate=[]
for (co,a,b,pop),g in groups.items():
 aggregate.append(dict(cohort=co,arm=a,comparator=b,population=pop,matched_acquisitions=len(g),median_gain_pct=statistics.median(r['gain_pct'] for r in g if r['gain_pct'] is not None),median_identity_recovery_pct=statistics.median(r['comparator_identity_recovery_pct'] for r in g if r['comparator_identity_recovery_pct'] is not None)))
tab(D/'matched_summary.tsv',aggregate)
(D/'COMPLETE.json').write_text(json.dumps(dict(status='PASS',tasks=len(tasks),matched_rows=len(out),scope='Canonical identification counts/overlaps and acquisition-wise matched summaries. Nested subsets retain baseline q; fresh7 has new search confidence. No final primary-policy decision or workflow-FDR claim.'),indent=2)+'\n')
print('PASS',len(tasks),len(out))
