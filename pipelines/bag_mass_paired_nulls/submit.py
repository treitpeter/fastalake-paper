from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,subprocess
HERE=Path(__file__).resolve().parent
m=json.loads((HERE/'manifest.json').read_text())
assert not (HERE/'SUBMITTED.json').exists()
for r in json.loads((HERE/'RUN_FREEZE.json').read_text())['files']:
 assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256']
queue=subprocess.check_output(['squeue','-u','treit','-o','%.22i %.38j %.12T %.30R'],text=True)
(HERE/'QUEUE_BEFORE.txt').write_text(queue)
assert 'FL_ANnull_' not in queue, 'Existing null experiment; inspect rather than duplicate'
pilots=set(m['paired_pilot_indices']);full=set(range(len(m['tasks'])))-pilots
assert len(pilots)==4 and len(full)==706 and not pilots & full
# Compact indices without changing membership.
def ranges(indices):
 values=sorted(indices);groups=[];start=last=values[0]
 for v in values[1:]:
  if v==last+1:last=v;continue
  groups.append(str(start) if start==last else f'{start}-{last}');start=last=v
 groups.append(str(start) if start==last else f'{start}-{last}')
 return ','.join(groups)
ids={}
def submit(name,stage,cpus,mem,hours,dependency,indices=None,limit=1):
 cmd=['sbatch','--parsable',f'--job-name=FL_ANnull_{name}',
      '--partition=p.hpcl8,p.hpcl92,p.hpcl93','--nice=10000',f'--cpus-per-task={cpus}',
      f'--mem={mem}G',f'--time={hours}:00:00',f'--chdir={HERE}',
      f'--output={HERE}/logs/%x_%A_%a.out',f'--error={HERE}/logs/%x_%A_%a.err',
      '--kill-on-invalid-dep=yes',f'--dependency={dependency}']
 if indices is not None:cmd.append(f'--array={ranges(indices)}%{limit}')
 cmd.append(str(HERE/'jobs'/f'{stage}.sbatch'))
 result=subprocess.run(cmd,capture_output=True,text=True)
 receipt=dict(name=name,command=cmd,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr,utc=datetime.now(timezone.utc).isoformat())
 with (HERE/'SUBMISSIONS.jsonl').open('a') as h:h.write(json.dumps(receipt)+'\n')
 result.check_returncode();job=result.stdout.strip().split(';')[0];assert job.isdigit()
 ids[name]=job
 (HERE/'logs'/f'scontrol_{name}.txt').write_text(subprocess.check_output(['scontrol','show','job',job],text=True))
 print(name,job,flush=True)
 return job
pc=submit('pilot_candidates','candidate',2,160,24,'aftercorr:5992346',pilots,1)
ps=submit('pilot_searches','search',4,64,24,f'aftercorr:{pc}',pilots,2)
pr=submit('pilot_report','report',2,16,4,f'aftercorr:{ps}:5992347',pilots,1)
pg=submit('pilot_gate','gate',1,4,1,f'afterok:{pr}')
fc=submit('candidates','candidate',2,160,24,f'afterok:{pg},aftercorr:5992349',full,1)
fs=submit('searches','search',4,64,24,f'afterok:{pg},aftercorr:{fc}',full,1)
fr=submit('reports','report',2,16,4,f'afterok:{pg},aftercorr:{fs}:5992350',full,1)
co=submit('collect','collect',1,8,1,f'afterany:{pr}:{fr}')
(HERE/'SUBMITTED.json').write_text(json.dumps(dict(jobs=ids,created_utc=datetime.now(timezone.utc).isoformat(),pilot_indices=sorted(pilots),remaining=len(full),status='SUBMITTED_NOT_COMPLETED',workflow_fdr='UNVERIFIED'),indent=2)+'\n')
(HERE/'QUEUE_AFTER.txt').write_text(subprocess.check_output(['squeue','-j',','.join(ids.values()),'-o','%.22i %.36j %.12T %.35R'],text=True))
