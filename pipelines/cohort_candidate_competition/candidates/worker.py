from pathlib import Path
import json,hashlib,sys,os,csv
import numpy as np
import pandas as pd
H=Path(__file__).resolve().parent;BASE=H.parent
sys.path.insert(0,str(BASE.parent/'confidence'))
from peptide import canonical

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def read(info):
 p=Path(info['results']);st=p.stat();assert (st.st_size,st.st_mtime_ns)==(info['results_bytes'],info['results_mtime_ns'])
 cols=['scannr','rank','label','scored_candidates','peptide_q','peptide','hyperscore','delta_next'];df=pd.read_csv(p,sep='\t',usecols=cols,dtype={'scannr':str,'peptide':str});df=df[df['rank']==1].set_index('scannr');assert df.index.is_unique and (df.scored_candidates>0).all();assert np.isfinite(df.peptide_q).all()
 ids={canonical(s) for s in df.loc[(df.label==1)&(df.peptide_q<=.01),'peptide']};ids={s for s in ids if 5<=len(s)<=50};e=info['expected_sets']['5'];assert sha(e['path'])==e['sha256'];assert ids==set(Path(e['path']).read_text().splitlines())
 h=sha(p);assert (st.st_size,st.st_mtime_ns)==(p.stat().st_size,p.stat().st_mtime_ns)
 return df,len(ids),{'path':str(p),'sha256':h,'baseline_peptide_ids_reproduced':True}
def main():
 mode=sys.argv[1];m=json.loads((BASE/'manifest.json').read_text());rows=m['rows']
 if mode=='pilot':rows=[next(r for r in rows if r['cohort']==c) for c in ['PREDICT','CAPSCAN']]
 assert mode in ['pilot','full'];out=H/mode;out.mkdir(exist_ok=False);summ=[];receipts=[]
 for i,r in enumerate(rows):
  f,nf,rf=read(r['arms']['FL_public']);g,ng,rg=read(r['arms']['MG']);both=f.index.intersection(g.index);assert len(both)>0
  a=f.loc[both];b=g.loc[both];ratios=b.scored_candidates.to_numpy()/a.scored_candidates.to_numpy();reduction=100*(1-a.scored_candidates.to_numpy()/b.scored_candidates.to_numpy())
  item={'cohort':r['cohort'],'sample':r['sample'],**r['design'],'reported_spectra_fl':len(f),'reported_spectra_mg':len(g),'paired_reported_spectra':len(both),'fl_only_reported_spectra':len(f.index.difference(g.index)),'mg_only_reported_spectra':len(g.index.difference(f.index)),'median_candidates_fl_paired':float(np.median(a.scored_candidates)),'median_candidates_mg_paired':float(np.median(b.scored_candidates)),'median_paired_candidate_ratio_mg_over_fl':float(np.median(ratios)),'median_paired_candidate_reduction_pct':float(np.median(reduction)),'peptide_ids_fl':nf,'peptide_ids_mg':ng,'peptide_gain_pct':100*(nf/ng-1) if ng else None}
  summ.append(item);receipts.append({'sample':r['sample'],'FL_public':rf,'MG':rg})
  if mode=='pilot':
   audit=pd.DataFrame({'scannr':both,'scored_candidates_fl':a.scored_candidates.to_numpy(),'scored_candidates_mg':b.scored_candidates.to_numpy(),'mg_over_fl':ratios});audit.to_csv(out/(r['cohort']+'_paired_spectra.tsv.gz'),sep='\t',index=False)
  print(mode,i+1,'/',len(rows),r['sample'],item['median_paired_candidate_ratio_mg_over_fl'],flush=True)
 pd.DataFrame(summ).to_csv(out/'per_acquisition.tsv',sep='\t',index=False)
 df=pd.DataFrame(summ);metrics=['median_candidates_fl_paired','median_candidates_mg_paired','median_paired_candidate_ratio_mg_over_fl','median_paired_candidate_reduction_pct','peptide_gain_pct'];df.groupby('cohort')[metrics].median().to_csv(out/'cohort_median_acquisition.tsv',sep='\t')
 df.groupby(['cohort','subject_id'])[metrics].median().groupby('cohort').median().to_csv(out/'cohort_median_participant.tsv',sep='\t')
 (out/'INPUT_RECEIPTS.json').write_text(json.dumps(receipts,indent=2)+'\n');(out/'COMPLETE.json').write_text(json.dumps({'status':'PASS','mode':mode,'acquisitions':len(rows),'job_id':os.getenv('SLURM_JOB_ID'),'manifest_sha256':sha(BASE/'manifest.json'),'code_sha256':sha(__file__),'metric':'Sage 0.14.6 scored_candidates: preliminary candidates with >=1 fragment match, summed over tested charge/isotope hypotheses','population':'Exact common reported rank-one spectrum IDs, including target and decoy winners, before confidence filtering. Unreported spectra are unavailable, not zero.','identification_endpoint':'Current canonical target peptide IDs, peptide q<=0.01, lengths5–50; every acquisition exactly reproduces current frozen IDs.','interpretation':'Measured search competition associated with yield; not raw theoretical hypotheses, a Bonferroni factor, or isolated causal effect of competition.'},indent=2)+'\n')
if __name__=='__main__':main()
