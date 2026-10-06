#!/usr/bin/env python3
"""Validate primary resource records and separate known cost from missing scope."""
import csv,json,re,sys
from pathlib import Path
import pandas as pd
import openpyxl
from dateutil.parser import parse as parse_date
from h01_eligibility import sha
base=Path(sys.argv[1]);m=json.loads((base/'manifest.json').read_text());out=base/'results_v2';out.mkdir(exist_ok=False)
for r in m['sources']:assert sha(r['path'])==r['sha256'],r['path']
w=openpyxl.load_workbook(m['workbook'],read_only=True,data_only=True)
sources=[]
it=iter(w['Index'].values);columns=next(it)
for row in it:
    r=dict(zip(columns,row));p=Path(r['path']);assert sha(p)==r['sha256'],str(p)
    sources.append(dict(path=str(p),sha256=sha(p),sheet=r['sheet']))
    if p.suffix=='.tsv':
        d=pd.read_csv(p,sep='\t');d.to_csv(out/(r['sheet']+'.tsv'),sep='\t',index=False)
(out/'PRIMARY_SOURCE_HASHES.json').write_text(json.dumps(sources,indent=2)+'\n')
def sheet(name):
    it=iter(w[name].values);columns=next(it);return pd.DataFrame(it,columns=columns)
historical=sheet('historical_compute_jobs');current=sheet('current_CAMPI_tasks')
for col in ['elapsed_s','cpu_alloc_s','cpu_used_s','ncpus','maxrss_gb']:historical[col]=pd.to_numeric(historical[col],errors='coerce')
summary=[]
for (co,stage),g in historical.groupby(['cohort','stage']):
    summary.append(dict(generation='historical_August',cohort=co,stage=stage,n_jobs=len(g),n_completed=int(g.state.eq('COMPLETED').sum()),sum_elapsed_task_hours=float(g.elapsed_s.sum()/3600),median_elapsed_seconds=float(g.elapsed_s.median()),allocated_CPU_hours=float(g.cpu_alloc_s.sum()/3600),measured_CPU_hours=float(g.cpu_used_s.sum()/3600),n_memory_observations=int(g.maxrss_gb.notna().sum()),maximum_sampled_memory_GB=None if g.maxrss_gb.isna().all() else float(g.maxrss_gb.max())))
pd.DataFrame(summary).to_csv(out/'historical_stage_resource_summary.tsv',sep='\t',index=False)
current.to_csv(out/'current_CAMPI_task_records.tsv',sep='\t',index=False)
# Raw execution timestamps expose prediction wall spans without inventing GPU utilization.
inference=[]
for p0 in m['inference_logs']:
    p=Path(p0);text=p.read_text(errors='replace')
    fields={}
    for key in ['sample','started','finished','checkpoint','rows']:
        match=re.search(r'^'+key+r'\s*:\s*(.+)$',text,re.M)
        fields[key]=match.group(1).strip() if match else None
    fields['population']='CAMPI_MOCK' if (fields['sample'] or '').startswith('S') else 'CAMPI_STOOL' if (fields['sample'] or '').startswith('F') else 'unresolved'
    fields['recorded_execution_span_seconds']=None
    if fields['started'] and fields['finished']:
        fields['recorded_execution_span_seconds']=(parse_date(fields['finished'],tzinfos={'CEST':7200,'CET':3600})-parse_date(fields['started'],tzinfos={'CEST':7200,'CET':3600})).total_seconds()
        assert fields['recorded_execution_span_seconds']>=0
    fields.update(path=str(p),sha256=sha(p));inference.append(fields)
pd.DataFrame(inference).to_csv(out/'CAMPI_prediction_execution_records.tsv',sep='\t',index=False)
# Fresh accounting is retained when the historical scheduler still serves it.
import subprocess
ids=[re.search(r'infer_(\d+)_',Path(p).name).group(1) for p in m['inference_logs'] if re.search(r'infer_(\d+)_',Path(p).name)]
r=subprocess.run(['sacct','-j',','.join(sorted(set(ids))),'--format=JobID,JobName,State,ElapsedRaw,NCPUS,TotalCPU,MaxRSS,Start,End,NodeList','--parsable2','--noheader'],capture_output=True,text=True)
(out/'INFERENCE_SACCT.tsv').write_text(r.stdout);(out/'INFERENCE_SACCT.stderr').write_text(r.stderr)
(out/'REPORT.md').write_text('''# Computational scope and provenance audit

The supplied computational-cost workbook was traced back to every indexed primary
source and its recorded hash. Historical August and current September CAMPI task
records remain separate. Per-stage distributions, actual CPU use, allocated CPU
hours and sampled-memory availability are reported without replacing missing
measurements by zero. Summed task elapsed time is labelled as such; it is not
end-to-end wall time. Acquisition medians are not added into an end-to-end estimate.

CAMPI prediction logs provide actual sample/start/finish/checkpoint/output-row
records. Fresh Slurm accounting is retained if still available; absent scheduler
rows do not imply zero runtime. These July prediction records contain fifteen mock acquisitions and two stool
acquisitions; an explicit population column keeps them separate. They are not
the complete September thirteen-acquisition stool experiment. Prediction GPU utilization, conversion/queueing and other cohorts'
inference require their corresponding records before an all-inclusive speed claim.

The HSTOOL inference scripts are retained in the manifest as primary implementation
evidence: they name the epoch12 checkpoint, ten constrained beams, window precursor
matching and a Hartout worktree. The worktree was inaccessible at preparation,
so its historical commit is not asserted from the path. Matching PREDICT/CapScan
execution scripts were not present alongside the supplied prediction CSVs.

No rerun is needed merely to invent a historical cost. Invoice/service-price,
sequencing-depth and turnaround sources remain separate documentary endpoints
No euro figure or whole-workflow speedup is claimed without a source.
This audit does not close missing records.
''')
(out/'COMPLETE.json').write_text(json.dumps(dict(status='PASS',primary_sources_verified=len(sources),historical_task_rows=len(historical),current_CAMPI_task_rows=len(current),prediction_logs=len(inference),complete_end_to_end_cost=False,outputs=[dict(path=str(p),sha256=sha(p)) for p in sorted(out.iterdir())]),indent=2)+'\n')
