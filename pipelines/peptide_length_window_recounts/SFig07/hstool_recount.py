"""Healthy-stool common-FASTA identities and complete features at 5-50 and 9-50.
Reads the accepted canonical identities and same-acquisition-gated feature matrices written by
A21 (search_revision2.py, extract) and A22 (extract_fileq.py, job 6048637). 9-50 = accepted 5-50
identities/features whose canonical peptide has >= 9 residues; q values unchanged. Statistics copy
compare_fileq.py (job 6048638). Writes per-acquisition and summary tables to OUTDIR."""
import os  # release: paths come from environment variables
import json,re,sys
from pathlib import Path
import numpy as np,pandas as pd
A21=Path((os.environ['FASTALAKE_DATA_DIR'] + '/MM_review_analysis/A21_7_hstool_quant_comparison'))
J=Path((os.environ['FASTALAKE_DATA_DIR'] + '/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/out'))
DIRS={'original':A21/'results/original','union_individual':A21/'results/union_individual','joint_corrected':A21/'results/joint_corrected','joint_fileq':J/'joint_fileq'}
TIERS=['TechnicalReplicates','WorkflowReplicates','OnePiece','Intra','Interindividual']
def tier(s):return next(t for t in TIERS if t in s)
def canonical(s):return re.sub(r'\[[^\]]*\]','',s).replace('I','L')
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
samples=json.loads((A21/'inputs/manifest_revision2.json').read_text())['samples']
ids={}
for k,d in DIRS.items():
    f=pd.read_csv(d/'accepted_identities.tsv.gz',sep='\t',dtype=str)
    ids[k]={s:set(f.loc[f['sample'].eq(s),'peptide']) for s in samples}
rows=[];summ=[];comp=[]
tech=[s for s in samples if tier(s)=='TechnicalReplicates']
for lo in (5,9):
    I={k:{s:{p for p in v if lo<=len(p)<=50} for s,v in ids[k].items()} for k in ids}
    for k in DIRS:
        rr=[]
        for s in samples:
            a=I['original'][s];b=I[k][s]
            rr.append(dict(min_length=lo,arm=k,sample=s,tier=tier(s),original=len(a),test=len(b),overlap=len(a&b),count_ratio=len(b)/len(a),original_identity_recovery=len(a&b)/len(a)))
        x=pd.DataFrame(rr);rows+=rr
        summ.append(dict(min_length=lo,arm=k,median_count_ratio=float(np.median(x.count_ratio)),n_higher=int((x.count_ratio>1).sum()),
            median_identity_recovery=float(np.median(x.original_identity_recovery)),study_distinct=len(set().union(*I[k].values())),
            ids_in_all12_technical=len(set.intersection(*[I[k][s] for s in tech]))))
for k,d in DIRS.items():
    z=np.load(d/'same_file_id.features.npz');v=pd.DataFrame(z['values'],index=z['keys'],columns=z['samples']).reindex(columns=samples).fillna(0)
    L=np.array([len(canonical(key.rsplit('|charge=',1)[0])) for key in v.index])
    for lo in (5,9):
        keep=(L>=lo)&(L<=50)
        for t in TIERS[:2]:
            ix=[s for s in samples if tier(s)==t]
            comp.append(dict(min_length=lo,space='features',arm=k,tier=t,acquisitions=len(ix),rows_complete_all=int((v.loc[keep,ix]>0).all(axis=1).sum())))
pd.DataFrame(rows).to_csv(out/'identity_per_acquisition_len.tsv',sep='\t',index=False)
pd.DataFrame(summ).to_csv(out/'identity_summary_len.tsv',sep='\t',index=False)
pd.DataFrame(comp).to_csv(out/'completeness_len.tsv',sep='\t',index=False)
print(pd.DataFrame(summ).to_string());print(pd.DataFrame(comp).to_string())
