"""capscan_crossfit_10fold: recovery. Moved from the original analysis as it ran, apart from path scrubbing; see README.md."""
# Identity recovery under the cross-fit: fraction of MG pooled accepted peptides also accepted by FL (cross-fit) and FULL (fold-pooled all-prediction)
import pandas as pd, numpy as np, os
B='LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/'
d=pd.read_csv(B+'iteration4/crossfit/run/PER_ACQUISITION.tsv',sep='\t')
rows=[]
for r in d.itertuples():
    p=f"{B}iteration4/crossfit/outputs/{r.cohort}/{r.sample}/"
    S={a:set(open(p+f"{a}_POOLED_ACCEPTED.txt").read().split()) for a in ['FL','FULL','MG']}
    assert len(S['FL'])==r.FL_pooled_5 and len(S['MG'])==r.MG_pooled_5, r.i
    o={'i':r.i,'cohort':r.cohort,'group':r.group,'subject_id':r.subject_id}
    for L in (5,9):
        T={a:{x for x in s if L<=len(x)<=50} for a,s in S.items()}
        mg=T['MG']
        o[f'rec_cf_{L}']=100*len(T['FL']&mg)/len(mg) if mg else np.nan
        o[f'rec_full_{L}']=100*len(T['FULL']&mg)/len(mg) if mg else np.nan
        o[f'n_mg_{L}']=len(mg)
    rows.append(o)
x=pd.DataFrame(rows)
x['stratum']=np.where(x.cohort=='PREDICT','PREDICT',np.where(x.group=='stool','CAPSCAN_stool','CAPSCAN_intestinal'))
x.to_csv('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/review_0/recovery_per_acquisition.tsv',sep='\t',index=False)
out=[]
for s,g in list(x.groupby('stratum'))+[('CAPSCAN_all',x[x.cohort=='CAPSCAN'])]:
    out.append({'stratum':s,'n':len(g),**{c:round(g[c].median(),2) for c in ['rec_cf_5','rec_full_5','rec_cf_9','rec_full_9']}})
pd.DataFrame(out).to_csv('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/review_0/recovery_summary.tsv',sep='\t',index=False)
print(pd.DataFrame(out).to_string())
