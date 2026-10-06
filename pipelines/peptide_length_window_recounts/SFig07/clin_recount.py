"""Recount accepted peptides per clinical acquisition (own vs common FASTA) at 5-50 and 9-50.
Same acceptance as compare_clinical_common_fasta_20260915.py (rank 1, label 1, peptide_q<=0.01,
canonical = mods stripped, I->L); 9-50 = canonical length >= 9 of the same accepted set. No new search.
Usage: clin_recount.py TASK_INDEX NTASKS OUTDIR"""
import os  # release: paths come from environment variables
import json,re,sys
from pathlib import Path
import pandas as pd
M=(os.environ['FASTALAKE_DATA_DIR'] + '/MM_review_analysis/end_of_day_2026-09-15/common_cache/comparison/inputs/manifest.json')
def canonical(s):return re.sub(r'\[[^\]]*\]','',s).replace('I','L')
def acc(path):
    d=pd.read_csv(path,sep='\t',usecols=['peptide','rank','label','peptide_q'],dtype={'peptide':str})
    d=d[(d['rank']==1)&(d.label==1)&(d.peptide_q<=.01)]
    return {canonical(p) for p in d.peptide}
i,n,out=int(sys.argv[1]),int(sys.argv[2]),Path(sys.argv[3])
tasks=json.load(open(M))['tasks'];rows=[]
for t in tasks[i::n]:
    r=json.load(open(t['common_receipt']));by={Path(x['path']).name:x['path'] for x in r['outputs']}
    old=acc(t['baseline_psm']['path']);new=acc(by['results.sage.tsv'])
    row=dict(cohort=t['cohort'],sample=t['sample'])
    for lo in (5,9):
        a={p for p in old if lo<=len(p)<=50};b={p for p in new if lo<=len(p)<=50}
        row.update({f'original_{lo}':len(a),f'common_{lo}':len(b),f'shared_{lo}':len(a&b)})
    rows.append(row);print(row,flush=True)
pd.DataFrame(rows).to_csv(out/f'part_{i:03d}.tsv',sep='\t',index=False)
