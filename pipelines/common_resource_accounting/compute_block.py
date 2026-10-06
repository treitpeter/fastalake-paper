"""Runs ONLY the 'Computing' block (lines 201-228) of the frozen build_revision.py (byte-identical copy in ../code), with its
own helper definitions (lines 15-21, 32-33 semantics: read/src/export/check), against the frozen inputs in the revision
folder's source_data/ and code/evidence/CAPSCAN.tsv. Modes:
  baseline331: block unchanged.
  fixed347: two substitutions: keepcap from C2 fixed_roster_347.tsv; roster check CAPSCAN 331 -> 347.
  negctl: fixed347 substitutions but keepcap from the 331 roster -> the 347 roster check must fail (exit 1)."""
import sys,hashlib,os
from pathlib import Path
import numpy as np,pandas as pd
mode=sys.argv[1]
O=Path(__file__).resolve().parents[1]
P=Path('LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/supplements/revision_2026-09-24')
code=(O/'code/build_revision.py').read_text().splitlines()
assert code[200].startswith('# Computing: per-acquisition stages') and code[227].startswith('s=pd.DataFrame(summary);ssf=export(s,')
block='\n'.join(code[200:228])
R347='LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/roster/fixed_roster_347.tsv'
if mode in('fixed347','negctl'):
    a="keepcap=set(cp.sample_id)";b="('CAPSCAN',331)"
    assert block.count(a)==1 and block.count(b)==1
    if mode=='fixed347':block=block.replace(a,"keepcap=set(pd.read_csv('"+R347+"',sep='\\t').sample_id);assert len(keepcap)==347")
    block=block.replace(b,"('CAPSCAN',347)")
D=O/mode;DATA=D/'source_data';DATA.mkdir(parents=True,exist_ok=True)
FILES=[];CHECKS=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def src(pattern):
 a=list((P/'source_data').glob(pattern));assert len(a)==1,(pattern,a);return a[0]
def read(pattern):return pd.read_csv(src(pattern),sep='\t')
def export(df,name,origin,meaning):
 f=DATA/(name+'.tsv');df.to_csv(f,sep='\t',index=False);FILES.append(dict(file=str(f),sha256=sha(f),rows=len(df)));return str(f)
def copy(f,folder='source_data/underlying'):pass   # the block's record copying is not part of the summary
def check(name,ok,details):
 CHECKS.append(dict(check=name,passed=bool(ok),details=str(details)));assert ok,(name,details)
(D/'block_executed.py').write_text(block+'\n')
exec(compile(block,str(D/'block_executed.py'),'exec'))
print('CHECKS',CHECKS);print('FILES',FILES)
