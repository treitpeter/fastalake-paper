"""Print the prepared Slurm chain; --execute submits after successful pilot checks."""
import argparse
import json
from pathlib import Path
import shlex
import subprocess

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--execute', action='store_true')
    p.add_argument('--hours-per-run', type=int, default=4)
    p.add_argument('--parallel', type=int, default=4)
    a = p.parse_args()
    if not 1 <= a.parallel <= 4 or not 1 <= a.hours_per_run <= 8:
        raise ValueError('Use bounded resources: 1–4 concurrent runs and 1–8 hours per run')
    download = BASE.parent / 'capscan_metag_availability_2026-09-09'
    registry = download / 'download_jobs.json'
    jobs = [str(r['job_id']) for r in json.loads(registry.read_text())] if registry.exists() else ['5899719']
    gate = ['sbatch', '--parsable', '--chdir=' + str(ROOT), '--job-name=capscan_download_gate',
            '--partition=p.hpcl92', '--cpus-per-task=1', '--mem=128M', '--time=12:05:00', '--nice=10000',
            '--dependency=afterany:' + ':'.join(jobs), '--output=' + str(BASE / 'logs/download_gate_%j.log'),
            '--wrap=' + shlex.join(['/usr/bin/python3', '-u', str(download / 'verify_download_complete.py'), '--wait-hours', '12'])]
    first = ['sbatch', '--parsable', '--chdir=' + str(ROOT), '--array=0-301%' + str(a.parallel),
             '--dependency=afterok:DOWNLOAD_GATE_ID', '--time=%02d:00:00' % a.hours_per_run, str(BASE / 'profile.sbatch'), 'all']
    if not a.execute:
        print(shlex.join(gate))
        print(shlex.join(first))
        print('Then: sbatch --dependency=afterok:PROFILE_ARRAY_ID ' + str(BASE / 'finish.sbatch'))
        return
    ready = json.loads((BASE / 'readiness.json').read_text())
    if not ready.get('pilots_passed'):
        raise RuntimeError('Both real-data pilots must pass before cohort submission')
    if (BASE / 'production_jobs.json').exists():
        raise RuntimeError('Production jobs already recorded; inspect them before resubmitting')
    gate_job = subprocess.check_output(gate, text=True).strip().split(';')[0]
    record = {'download_gate': gate_job, 'download_gate_command': gate}
    (BASE / 'production_jobs.json').write_text(json.dumps(record, indent=2))
    first = [x.replace('DOWNLOAD_GATE_ID', gate_job) for x in first]
    job = subprocess.check_output(first, text=True).strip().split(';')[0]
    record.update(profile_array=job, profile_command=first)
    # Preserve the submitted job immediately if downstream submission fails.
    (BASE / 'production_jobs.json').write_text(json.dumps(record, indent=2))
    last = ['sbatch', '--parsable', '--chdir=' + str(ROOT), '--dependency=afterok:' + job, str(BASE / 'finish.sbatch')]
    record['comparison_job'] = subprocess.check_output(last, text=True).strip().split(';')[0]
    record['comparison_command'] = last
    (BASE / 'production_jobs.json').write_text(json.dumps(record, indent=2))
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
