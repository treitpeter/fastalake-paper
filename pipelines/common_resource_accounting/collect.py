"""Gate: baseline331 compute_stage_summary.tsv and compute_task_records.tsv byte-identical to the delivered Source Data;
non-CapScan rows of fixed347 identical to baseline. --negative plants one changed cell and must exit 1. Writes RESULTS.md."""
import sys,hashlib
from pathlib import Path
import pandas as pd
O=Path(__file__).resolve().parents[1]
DEL=Path('LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/supplements/revision_2026-09-24/delivery/source_data')
PKG=Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/figure_sd/source_data/package/FastaLake_Source_Data/Supplementary_Figure_15/compute_stage_summary.tsv')
h=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
ok=h(O/'baseline331/source_data/compute_stage_summary.tsv')==h(DEL/'compute_stage_summary.tsv')==h(PKG) and h(O/'baseline331/source_data/compute_task_records.tsv')==h(DEL/'compute_task_records.tsv')
b=pd.read_csv(O/'baseline331/source_data/compute_stage_summary.tsv',sep='\t');f=pd.read_csv(O/'fixed347/source_data/compute_stage_summary.tsv',sep='\t')
if '--negative' in sys.argv:f.loc[0,'n_tasks']+=1
ok=ok and b[b.cohort.ne('CAPSCAN')].equals(f[f.cohort.ne('CAPSCAN')])
ok=ok and int(f[(f.cohort=='CAPSCAN')&(f.stage=='07_sage')].n_tasks.iloc[0])==347
if not ok:print('FAIL');sys.exit(1)
L=['# A7 resource accounting (Supplementary Fig. 15), CapScan 331 vs 347','',
'Producer: build_revision.py (supplements/revision_2026-09-24), Computing block l.201-228, run on its frozen inputs (byte-identical copy, COPY_SHA256.txt). 347 changes two expressions: CapScan keep-set = C2 fixed_roster_347.tsv; roster check 331 -> 347 (see fixed347/block_executed.py vs baseline331/block_executed.py).','',
'Gate PASS: baseline331 compute_stage_summary.tsv and compute_task_records.tsv are byte-identical to the delivered and packaged Source Data; all non-CapScan rows unchanged in 347.','',
'| CapScan stage | n_tasks 331 | n_tasks 347 | memory-sampled 331 | 347 | median wall min 331 | 347 | median memory GiB 331 | 347 | allocated CPU-h 331 | 347 |','|---|---|---|---|---|---|---|---|---|---|---|']
for st in b[b.cohort=='CAPSCAN'].stage:
    x=b[(b.cohort=='CAPSCAN')&(b.stage==st)].iloc[0];y=f[(f.cohort=='CAPSCAN')&(f.stage==st)].iloc[0]
    L.append(f'| {st} | {x.n_tasks} | {y.n_tasks} | {x.n_memory_sampled} | {y.n_memory_sampled} | {x.wall_median_min:.3f} | {y.wall_median_min:.3f} | {x.memory_median_GiB:.3f} | {y.memory_median_GiB:.3f} | {x.allocated_CPU_hours:.2f} | {y.allocated_CPU_hours:.2f} |')
L+=['','Manuscript text: Supplementary Fig. 15 legend and Methods ("306 PREDICT and 331 CapScan (manuscript) acquisitions") -> 347; the shared 368-acquisition CapScan stages (02_evidence, 08_aggregate, 10_fdp_shuffle) are unchanged by construction. The figure x-axis label "CapScan 331" is hard-coded in build_revision.py l.237 (not re-rendered here).']
(O/'RESULTS.md').write_text('\n'.join(L)+'\n');print('PASS')
