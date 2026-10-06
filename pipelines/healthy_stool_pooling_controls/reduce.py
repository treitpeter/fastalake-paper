from pathlib import Path
import json,csv,subprocess,hashlib,collections
import pandas as pd
H=Path(__file__).resolve().parent;M=json.loads((H/'MANIFEST.json').read_text());R=Path(M['full_route_parent']);cfg=json.loads((R/'experiment.json').read_text());rows=[];disc=[]
for task in M['tasks']:
 p=Path(task['config']).parent;assert (p/'COMPLETE.json').exists()
 for r in csv.DictReader((p/'summary.tsv').open(),delimiter='\t'):rows.append(dict(sample=task['sample'],arm=task['arm'],**r))
 for r in csv.DictReader((p/'discoveries.tsv').open(),delimiter='\t'):
  if r['assignment']!='target':disc.append(dict(sample=task['sample'],arm=task['arm'],**r))
known={r['peptide']:int(r['present_in_original_reservoir']) for r in csv.DictReader((R/'tables/HSTOOL/original_reservoir_overlap.tsv').open(),delimiter='\t')};missing={r['peptide'] for r in disc}-set(known)
(H/'new_overlap_queries.txt').write_text('\n'.join(sorted(missing))+'\n')
if missing:
 subprocess.run([str(R/'src/scan_overlap'),str(H/'new_overlap_queries.txt'),cfg['cohorts']['HSTOOL']['reservoir'],str(H/'new_original_reservoir_overlap.tsv'),'16'],check=True)
 known.update({r['peptide']:int(r['present_in_original_reservoir']) for r in csv.DictReader((H/'new_original_reservoir_overlap.tsv').open(),delimiter='\t')})
for r in disc:r['present_in_original_reservoir']=known[r['peptide']]
by=collections.defaultdict(list)
for r in disc:by[r['sample'],r['arm'],float(r['q_threshold'])].append(r)
for r in rows:
 rr=by[r['sample'],r['arm'],float(r['q_threshold'])];n=int(r['n_canonical']);floor=int(r['min_residues'])
 for cls in ['shuffle','plant','thermophilic']:
  novel=sum(x['assignment']==cls and not x['present_in_original_reservoir'] and len(x['peptide'])>=floor for x in rr);r[cls+'_only_not_in_original_reservoir']=novel;r[cls+'_only_not_in_original_reservoir_fraction']=novel/n if n else None
new=pd.DataFrame(rows);old=pd.read_csv(R/'tables/HSTOOL/per_acquisition.tsv',sep='\t');samples={r['sample'] for r in M['sources']};old=old[old['sample'].isin(samples)].assign(arm='own');full=pd.concat([old,new],ignore_index=True)
for c in full.columns:
 if c not in ['cohort','sample','arm']:full[c]=pd.to_numeric(full[c])
full.to_csv(H/'per_acquisition.tsv',sep='\t',index=False);pd.DataFrame(disc).to_csv(H/'control_discoveries.tsv',sep='\t',index=False)
primary=full[(full.q_threshold==.01)&(full.min_residues==5)];assert len(primary)==36
summ=[]
for arm,d in primary.groupby('arm'):
 for cls in ['shuffle','plant','thermophilic']:
  n=d[cls+'_only_not_in_original_reservoir'];v=100*d[cls+'_only_not_in_original_reservoir_fraction'];summ.append(dict(arm=arm,control_class=cls,n_acquisitions=len(d),n_positive=int((n>0).sum()),median_count=float(n.median()),total_count=int(n.sum()),median_fraction_pct=float(v.median()),max_fraction_pct=float(v.max())))
pd.DataFrame(summ).to_csv(H/'summary.tsv',sep='\t',index=False)
(H/'VALIDATION.json').write_text(json.dumps({'completed_new_searches':24,'reused_full_route_own_searches':12,'new_sequences_scanned_for_reservoir_overlap':len(missing),'target_final_fasta_excluded_other11':all(r['sample'] not in r['donors'] for r in M['tasks'] if r['arm']=='pool_other11'),'whole_procedure_FDR_calibrated':False,'summary':summ},indent=2));print(json.dumps(summ,indent=2))
