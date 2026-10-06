"""Submit the frozen native AlphaNovo experiment with dependent full-cohort jobs."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import shlex
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REPO = ROOT / 'worktrees/bags'
M = json.loads((HERE / 'manifest.json').read_text())
PY = M['python']

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

assert not (HERE / 'SUBMISSIONS.jsonl').exists(), 'Do not duplicate submissions'
assert json.loads((HERE / 'PREFLIGHT.json').read_text())['tasks'] == 710
for rec in M['code'] + [M['bag_engine'], M['sage']]:
    assert sha(rec['path']) == rec['sha256'], rec['path']
queue = subprocess.check_output(['squeue', '-u', 'treit', '-h', '-o', '%i %j %T %R'], text=True)
assert 'FL_ANmass_' not in queue, 'Native mass experiment is already queued'
(HERE / 'QUEUE_BEFORE.txt').write_text(queue)
manifest_hash = sha(HERE / 'manifest.json')
guard = HERE / 'guard.py'
guard.write_text(
    'from pathlib import Path\nimport hashlib, subprocess, sys\n'
    + 'manifest = Path(' + repr(str(HERE / 'manifest.json')) + ')\n'
    + 'assert hashlib.sha256(manifest.read_bytes()).hexdigest() == ' + repr(manifest_hash) + '\n'
    + 'subprocess.run([sys.executable, "-B", ' + repr(str(REPO / 'tools/run_bags.py'))
    + ', str(manifest), *sys.argv[1:]], check=True)\n'
)
for stage in ['candidate', 'search', 'gate', 'group', 'collect']:
    command = ' '.join(shlex.quote(arg) for arg in [PY, '-B', str(guard), stage])
    if stage in ['candidate', 'search', 'group']:
        command += ' "$SLURM_ARRAY_TASK_ID"'
    script = ('#!/bin/bash\nset -euo pipefail\n'
              'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1\n'
              'exec /usr/bin/time -v ' + command + '\n')
    (HERE / 'jobs' / (stage + '.sbatch')).write_text(script)
freeze = [HERE / 'manifest.json', guard, HERE / 'PREFLIGHT.json', HERE / 'FASTALAKE_TESTS.xml',
          HERE / 'submit.py'] + sorted((HERE / 'jobs').glob('*.sbatch'))
(HERE / 'RUN_FREEZE.json').write_text(json.dumps({str(p): sha(p) for p in freeze}, indent=2) + '\n')

def submit(label, stage, cpus, mem, hours, array=None, dependency=None):
    command = ['sbatch', '--parsable', '--job-name=FL_ANmass_' + label,
               '--partition=p.hpcl8,p.hpcl92,p.hpcl93', '--nice=10000',
               '--cpus-per-task=' + str(cpus), '--mem=' + mem, '--time=' + str(hours) + ':00:00',
               '--chdir=' + str(REPO), '--output=' + str(HERE / 'logs/%x_%A_%a.out'),
               '--error=' + str(HERE / 'logs/%x_%A_%a.err'), '--kill-on-invalid-dep=yes']
    if array:
        command.append('--array=' + array)
    if dependency:
        command.append('--dependency=' + dependency)
    command.append(str(HERE / 'jobs' / (stage + '.sbatch')))
    result = subprocess.run(command, text=True, capture_output=True)
    rec = dict(checked_utc=datetime.now(timezone.utc).isoformat(), label=label, argv=command,
               returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
    if result.returncode == 0:
        rec['job_id'] = result.stdout.strip().split(';')[0]
        assert rec['job_id'].isdigit()
    with (HERE / 'SUBMISSIONS.jsonl').open('a') as f:
        f.write(json.dumps(rec) + '\n')
    print(json.dumps(rec), flush=True)
    result.check_returncode()
    return rec['job_id']

def ranges(values):
    parts = []
    start = previous = values[0]
    for value in values[1:] + [None]:
        if value is not None and value == previous + 1:
            previous = value
            continue
        parts.append(str(start) if start == previous else str(start) + '-' + str(previous))
        start = previous = value
    return ','.join(parts)

pilots = M['paired_pilot_indices']
remaining = [i for i in range(len(M['tasks'])) if i not in pilots]
pilot_array = ','.join(map(str, pilots)) + '%2'
full_array = ranges(remaining) + '%4'
ids = {}
ids['pilot_candidates'] = submit('pilot_candidates', 'candidate', 2, '160G', 24, array=pilot_array)
ids['pilot_searches'] = submit('pilot_searches', 'search', 4, '64G', 24, array=pilot_array,
                               dependency='afterok:' + ids['pilot_candidates'])
ids['pilot_gate'] = submit('pilot_gate', 'gate', 2, '32G', 4,
                           dependency='afterok:' + ids['pilot_searches'])
ids['full_candidates'] = submit('candidates', 'candidate', 2, '160G', 24, array=full_array,
                                dependency='afterok:' + ids['pilot_gate'])
ids['full_searches'] = submit('searches', 'search', 4, '64G', 24, array=full_array,
                              dependency='afterok:' + ids['pilot_gate'] + ',aftercorr:' + ids['full_candidates'])
ids['study_grouping'] = submit('study_groups', 'group', 4, '192G', 24, array='0-15%2',
                              dependency='afterok:' + ids['full_searches'])
ids['collector'] = submit('collect', 'collect', 1, '8G', 2,
                          dependency='afterany:' + ids['full_searches'] + ':' + ids['study_grouping'])
receipt = dict(status='SUBMITTED_NOT_COMPLETED', job_ids=ids, tasks=710, cohorts=M['expected_cohort_counts'],
               pilot_indices=pilots, full_tasks=len(remaining), pilot_concurrency=2, candidate_concurrency=4,
               search_concurrency=4, group_concurrency=2, fasta_lake_revision=M['fasta_lake_revision'],
               alphanovo_revision=M['alphanovo_revision'], manifest_sha256=manifest_hash,
               scientific_review='PENDING', speed_and_identification_results='PENDING',
               checked_utc=datetime.now(timezone.utc).isoformat())
(HERE / 'SUBMITTED.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt, indent=2), flush=True)
