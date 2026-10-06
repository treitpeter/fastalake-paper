"""Submit the frozen, bounded workflow exactly once and retain every job ID."""
from pathlib import Path
import json
import subprocess
from datetime import datetime, timezone
import run_experiment as m

ROOT = Path(__file__).resolve().parent
m.verify_code()
(ROOT/'status/SUBMISSION_LOCK').open('x').close()
record = {'submitted_utc':datetime.now(timezone.utc).isoformat(), 'expected_searches':30,
          'seeds':[20260907,20260908], 'reservoir_sequences':375333654,
          'max_concurrent_tasks':6, 'jobs':{}, 'commands':[]}
def save():
    (ROOT/'status/SUBMITTED.json').write_text(json.dumps(record,indent=2)+'\n')
def submit(stage, name, cpu, memory, limit, dependency=None, array=None):
    cmd=['sbatch','--parsable','--job-name='+name,'--nice=10000',
         '--partition=p.hpcl92,p.hpcl93,p.hpcl8','--cpus-per-task='+str(cpu),
         '--mem='+memory,'--time='+limit,'--chdir='+str(ROOT),
         '--output='+str(ROOT/'logs/%x_%A_%a.out'),'--error='+str(ROOT/'logs/%x_%A_%a.err')]
    if dependency:cmd+=['--dependency=afterok:'+dependency]
    if array:cmd+=['--array='+array]
    cmd += [str(ROOT/'job.sh'),stage]
    job=subprocess.check_output(cmd,text=True).strip().split(';')[0]
    assert job.isdigit(),job
    record['jobs'][name]=job;record['commands'].append(cmd);save()
    print(name,job,flush=True)
    return job

save()
a=submit('preflight','FL_C375_preflight',8,'96G','02:00:00')
b=submit('prepare','FL_C375_prepare',16,'96G','06:00:00',a,'0-1%2')
c=submit('sample','FL_C375_pilots',8,'96G','01:00:00',b,'0-1%2')
d=submit('sample','FL_C375_remaining',8,'96G','01:00:00',c,'2-29%6')
e=submit('collect','FL_C375_collect',2,'8G','00:30:00',d)
record['submission_complete']=True;save()
