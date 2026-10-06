from pathlib import Path
import ast,csv,gzip,hashlib,json,os,sys
import numpy as np
import pandas as pd
from scipy import sparse
H=Path(__file__).resolve().parent; C=(Path(__import__('os').environ['FASTALAKE_DATA_DIR'])/'review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01'); J=C.parents[1]
R=Path(__import__('os').environ['FASTALAKE_DATA_DIR'])
A33=R/'MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19'
O=C.parent/'github_workflow_reconciliation_2026-09-30/resumed_audit/annotation_repair/annotated_origin_consensus_2026-10-01'
X=J/'source_data/current/Crossomics'
E=C.parent/'eds_reaction_check_2026-09-30'
N=R/'MM_review_analysis/end_of_day_2026-09-15/eds_phenotype'
OLD=R/'analyses/enterotypes_2026-09-07'
REG=['duodenum','jejunum','ileum','colon','stool']; DONORS=[1,2,3,4,6,8,9,10,11,12,13,14]
def load(p):return json.loads(Path(p).read_text())
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,default=str)+'\n')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def rd(p,**kw):return pd.read_csv(p,sep='\t',**kw)
def wr(d,p):d.to_csv(p,sep='\t',index=False)
def human(t):return t['rank']=='species' and t['supported_lineage'][-1:]==['s__Homo sapiens']
def frame(folder,key):
 a=load(folder/(key+'_axes.json'));v=sparse.load_npz(folder/(key+'.npz'))
 return pd.DataFrame(v.toarray(),index=a['samples'],columns=a['features'])
def groupdir(co):
 root=R/('simplified_grouping_results/manuscript_cohorts_2026-09-10_2200' if co=='CAPSCAN' else 'MM_review_analysis/recovery_2026-09-15/biological_grouping_recovery/full')
 return root/'cohorts'/(co+'_noclust')
def functions(path,names,env):
 tree=ast.parse(Path(path).read_text());defs=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
 assert {n.name for n in defs}==set(names)
 exec(compile(ast.Module(body=defs,type_ignores=[]),str(path),'exec'),env);return env

def verify():
 for p,h in load(H/'CODE.json').items():assert sha(H/p)==h,p
 for p,h in load(H/'INPUTS.json').items():assert sha(p)==h,p

def receipt(out,**kw):
 dump(out/'COMPLETE.json',dict(status='PASS',job=os.environ.get('SLURM_JOB_ID'),code=sha(H/'CODE.json'),inputs=sha(H/'INPUTS.json'),**kw,outputs={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='COMPLETE.json'}))
def gate(out):
 d=load(out/'COMPLETE.json');assert d['status']=='PASS' and d['code']==sha(H/'CODE.json') and d['inputs']==sha(H/'INPUTS.json')
 for f,h in d['outputs'].items():assert sha(out/f)==h,f
 return d
