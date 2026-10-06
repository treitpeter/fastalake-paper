#!/usr/bin/env python3
"""Exact FASTA unions and a single joint Sage search per retained cohort/arm."""
import argparse
import collections
import copy
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
FL = ROOT.parents[1]
SOURCE = FL / 'simplified_grouping_results/manuscript_cohorts_2026-09-10_2200/manifest.json'
SAGE = Path('LAB_ROOT/sage')
TERMINAL = {'COMPLETED', 'FAILED', 'OUT_OF_MEMORY', 'TIMEOUT', 'CANCELLED',
            'NODE_FAIL', 'BOOT_FAIL', 'DEADLINE', 'PREEMPTED', 'REVOKED'}


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def save(path, obj):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True) + '\n')
    tmp.replace(path)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def stamp(path, hashed=False):
    path = Path(path).resolve(strict=True)
    st = path.stat()
    if not path.is_file() or st.st_size == 0:
        raise ValueError(f'Missing/empty input: {path}')
    out = dict(path=str(path), bytes=st.st_size, mtime_ns=st.st_mtime_ns)
    if hashed:
        out['sha256'] = sha(path)
    return out


def verify(item):
    actual = stamp(item['path'], 'sha256' in item)
    if actual != item:
        raise ValueError(f'Input changed since preparation: {item["path"]}')


def parameters(config):
    out = copy.deepcopy(config)
    for key in ('version', 'output_paths', 'output_directory', 'mzml_paths'):
        out.pop(key, None)
    out['database'].pop('fasta', None)
    return out


def prepare():
    if (ROOT / 'manifest.json').exists():
        raise ValueError('Prepared manifest already exists; do not overwrite a run')
    source = json.loads(SOURCE.read_text())
    previous = json.loads((SOURCE.parent / 'collection/VALIDATION.json').read_text())
    if previous['status'] != 'PASS' or len(previous['tasks']) != 8:
        raise ValueError('Source simplified-grouping run did not pass')
    (ROOT / 'bin').mkdir(exist_ok=True)
    shutil.copy2(SAGE, ROOT / 'bin/sage')
    tasks = []
    for old in source['tasks']:
        folder = ROOT / 'cohorts' / old['id']
        folder.mkdir(parents=True, exist_ok=False)
        (folder / 'source_configs').mkdir()
        names = Path(old['samples_file']).read_text().splitlines()
        if len(names) != old['acquisitions'] or len(set(names)) != len(names):
            raise ValueError(f'Invalid roster: {old["id"]}')
        expected = {str(Path(x['path']).resolve()): x for x in old['search_inputs']}
        records, fastas, configs, base = [], {}, [], None
        for sample in names:
            config_path = Path(old['search_root']) / sample / old['config_name']
            config = json.loads(config_path.read_text())
            source_config = stamp(config_path, True)
            configs.append(source_config)
            shutil.copy2(config_path, folder / 'source_configs' / f'{sample}.json')
            current = parameters(config)
            if base is None:
                base = current
            if current != base:
                raise ValueError(f'Incompatible search parameters: {old["id"]}/{sample}')
            if len(config['mzml_paths']) != 1:
                raise ValueError(f'Expected one original raw per acquisition: {sample}')
            fasta = stamp(config['database']['fasta'])
            for item in (source_config, fasta):
                prev = expected[item['path']]
                if any(item[k] != prev[k] for k in ('bytes', 'mtime_ns')):
                    raise ValueError(f'Changed source simplified-grouping input: {item["path"]}')
            fastas[fasta['path']] = fasta
            raw = stamp(config['mzml_paths'][0])
            records.append(dict(sample=sample, raw=raw, fasta=fasta['path'],
                                source_config=source_config['path']))
        raws = [r['raw']['path'] for r in records]
        if len(set(raws)) != len(raws) or len({Path(p).name for p in raws}) != len(raws):
            raise ValueError('Duplicate raw files or basenames make joint output ambiguous')
        if base['database']['generate_decoys'] is not True:
            raise ValueError('This experiment requires fresh internal Sage decoys')
        config = copy.deepcopy(base)
        config['database']['fasta'] = str(folder / 'union.fasta')
        config['mzml_paths'] = raws
        config['output_directory'] = str(folder / 'search')
        save(folder / 'sage.json', config)
        (folder / 'samples.txt').write_text('\n'.join(names) + '\n')
        with (folder / 'inputs.tsv').open('w') as stream:
            writer = csv.writer(stream, delimiter='\t')
            writer.writerow(['sample', 'mzml', 'source_fasta', 'source_config'])
            writer.writerows((r['sample'], r['raw']['path'], r['fasta'], r['source_config']) for r in records)
        if old['cohort'] == 'HSTOOL':
            mem, cpus, wall = '128G', 24, '04:00:00'
        elif old['cohort'] in ('PREDICT', 'CAPSCAN'):
            mem, cpus, wall = '384G', 32, '08:00:00'
        else:
            mem, cpus, wall = '48G', 16, '02:00:00'
        task = dict(id=old['id'], cohort=old['cohort'], arm=old['arm'],
                    acquisitions=len(records), records=records, fastas=list(fastas.values()),
                    source_configs=configs, folder=str(folder), sage_config=stamp(folder / 'sage.json', True),
                    source_fasta_bytes=sum(f['bytes'] for f in fastas.values()),
                    raw_bytes=sum(r['raw']['bytes'] for r in records),
                    batch_size=4, resources=dict(partition=('p.hpcl93,p.hpcl92' if mem == '384G' else 'p.hpcl93,p.hpcl92,p.hpcl8'), memory=mem, cpus=cpus, walltime=wall))
        tasks.append(task)
    manifest = dict(name='sample_specific_merge_all_raw', prepared_utc=now(),
                    source_manifest=stamp(SOURCE, True), source_validation=stamp(SOURCE.parent / 'collection/VALIDATION.json', True),
                    sage=stamp(ROOT / 'bin/sage', True),
                    sage_version=subprocess.check_output([str(ROOT / 'bin/sage'), '--version'], text=True).strip(),
                    workflow_sha256=sha(__file__), python=sys.executable,
                    max_concurrent_search_jobs=2, begin_local='2026-09-11T22:00:00',
                    deadline_local='2026-09-12T08:00:00', timezone='Europe/Berlin', tasks=tasks)
    save(ROOT / 'manifest.json', manifest)
    print(json.dumps({t['id']: dict(acquisitions=t['acquisitions'], source_fastas=len(t['fastas']),
                                    raw_GB=round(t['raw_bytes']/1e9, 2), resources=t['resources']) for t in tasks}, indent=2))


def fasta_records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    header, parts = None, []
    with opener(path, 'rt') as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts).upper()
                header, parts = line[1:], []
                if not header.strip():
                    raise ValueError(f'Empty FASTA header: {path}')
            elif header is None:
                raise ValueError(f'Sequence before FASTA header: {path}')
            else:
                parts.append(''.join(line.split()))
        if header is not None:
            yield header, ''.join(parts).upper()


def build_union(task):
    folder = Path(task['folder'])
    final = folder / 'union.fasta'
    if final.exists():
        raise ValueError('Union already exists; refusing to overwrite')
    sequences, accession_hashes = {}, {}
    counts = collections.Counter()
    source_hashes, source_counts = [], []
    tag = json.loads((folder / 'sage.json').read_text())['database']['decoy_tag']
    with gzip.open(folder / 'sequence_provenance.tsv.gz', 'wt', compresslevel=3) as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['union_id', 'source_fasta', 'original_accession', 'original_header'])
        for item in sorted(task['fastas'], key=lambda f: f['path']):
            verify(item)
            count = 0
            for header, sequence in fasta_records(item['path']):
                if not sequence or re.fullmatch('[A-Z*]+', sequence) is None:
                    raise ValueError(f'Empty/invalid protein sequence: {item["path"]}: {header}')
                accession = header.split()[0]
                if tag in accession:
                    raise ValueError(f'Unexpected decoy in source target FASTA: {accession}')
                digest = hashlib.sha256(sequence.encode('ascii')).hexdigest()
                existing = sequences.setdefault(digest, sequence)
                if existing != sequence:
                    raise ValueError('SHA256 collision; exact sequence identity is required')
                previous = accession_hashes.setdefault(accession, digest)
                if previous != digest:
                    counts['accession_sequence_conflict_records'] += 1
                writer.writerow(['FLU_' + digest, item['path'], accession, header])
                count += 1
            if count == 0:
                raise ValueError(f'No records in {item["path"]}')
            source_counts.append(dict(path=item['path'], target_records=count))
            source_hashes.append(stamp(item['path'], True))
            verify(item)
            counts['input_records'] += count
            print(f'UNION source={item["path"]} records={count} unique_so_far={len(sequences)}', flush=True)
    tmp = folder / 'union.fasta.tmp'
    with tmp.open('w') as stream:
        for digest, sequence in sorted(sequences.items()):
            stream.write('>FLU_' + digest + '\n')
            for start in range(0, len(sequence), 80):
                stream.write(sequence[start:start+80] + '\n')
    observed = set()
    for header, sequence in fasta_records(tmp):
        digest = hashlib.sha256(sequence.encode('ascii')).hexdigest()
        if header != 'FLU_' + digest or digest in observed or sequences[digest] != sequence:
            raise ValueError('Written union failed exact-sequence roundtrip')
        observed.add(digest)
    if len(observed) != len(sequences):
        raise ValueError('Union lost sequence records')
    tmp.replace(final)
    summary = dict(status='PASS', completed_utc=now(), **counts,
                   unique_sequences=len(sequences), duplicate_records_removed=counts['input_records']-len(sequences),
                   total_residues=sum(map(len, sequences.values())), union=stamp(final, True),
                   source_fastas=source_hashes, source_counts=source_counts,
                   identity='Exact full amino-acid sequence; uppercase and remove wrapping whitespace; retain I/L differences and stop symbols',
                   identifiers='FLU_ + full SHA256(sequence); complete original headers and FASTA membership in sequence_provenance.tsv.gz')
    save(folder / 'union_summary.json', summary)
    print(json.dumps({k: v for k, v in summary.items() if k not in ('source_fastas', 'source_counts')}, indent=2), flush=True)


def load_task(task_id):
    manifest = json.loads((ROOT / 'manifest.json').read_text())
    if sha(__file__) != manifest['workflow_sha256']:
        raise ValueError('Workflow changed after preparation')
    return manifest, next(t for t in manifest['tasks'] if t['id'] == task_id)


def timed(command, prefix, folder):
    command = ['/usr/bin/time', '-v', '-o', str(folder / f'{prefix}.time.txt')] + command
    with (folder / f'{prefix}.stdout.log').open('w') as out, (folder / f'{prefix}.stderr.log').open('w') as err:
        return subprocess.run(command, stdout=out, stderr=err).returncode


def validate_search(task):
    folder = Path(task['folder']) / 'search'
    required = ['results.json', 'results.sage.tsv', 'lfq.tsv', 'results.sage.pin']
    outputs = [stamp(folder / n) for n in required]
    actual = json.loads((folder / 'results.json').read_text())
    expected = json.loads(Path(task['sage_config']['path']).read_text())
    if parameters(actual) != parameters(expected):
        raise ValueError('Sage effective parameters changed')
    if actual['mzml_paths'] != expected['mzml_paths'] or Path(actual['database']['fasta']).resolve() != Path(expected['database']['fasta']).resolve():
        raise ValueError('Sage effective FASTA/raw roster differs from prepared inputs')
    with (folder / 'lfq.tsv').open() as stream:
        fields = next(csv.reader(stream, delimiter='\t'))
    fixed = {'peptide', 'charge', 'proteins', 'q_value', 'score', 'spectral_angle'}
    samples = [f for f in fields if f not in fixed]
    expected_names = {Path(r['raw']['path']).name for r in task['records']}
    if len(samples) != len(expected_names) or set(samples) != expected_names:
        raise ValueError(f'Joint LFQ does not contain every raw column: {samples}')
    with (folder / 'results.sage.tsv').open() as stream:
        fields = next(csv.reader(stream, delimiter='\t'))
    if not {'filename', 'peptide', 'proteins', 'label', 'spectrum_q', 'peptide_q', 'protein_q'} <= set(fields):
        raise ValueError('Unexpected Sage result schema')
    return dict(status='PASS', completed_utc=now(), joint_raw_columns=len(samples), outputs=outputs)


def run(task_id):
    manifest, task = load_task(task_id)
    folder = Path(task['folder'])
    if (folder / 'RUN.json').exists():
        raise ValueError('An attempt exists; do not overwrite it')
    record = dict(task=task_id, state='VERIFYING_INPUTS', started_utc=now(), slurm_job_id=os.environ.get('SLURM_JOB_ID'), host=os.uname().nodename)
    save(folder / 'RUN.json', record)
    start = time.monotonic()
    try:
        verify(manifest['sage'])
        verify(task['sage_config'])
        for item in task['source_configs']:
            verify(item)
        for row in task['records']:
            verify(row['raw'])
        record['state'] = 'BUILDING_UNION'
        save(folder / 'RUN.json', record)
        rc = timed([sys.executable, str(Path(__file__).resolve()), 'build', task_id], 'build', folder)
        if rc:
            raise RuntimeError(f'Union construction exited {rc}; see build.stderr.log')
        (folder / 'search').mkdir(exist_ok=False)
        command = [manifest['sage']['path'], task['sage_config']['path'], '--batch-size', str(task['batch_size']),
                   '--write-pin', '--disable-telemetry-i-dont-want-to-improve-sage']
        record.update(state='SEARCHING', sage_command=command, search_started_utc=now())
        save(folder / 'RUN.json', record)
        rc = timed(command, 'sage', folder)
        record['sage_exit_code'] = rc
        if rc:
            raise RuntimeError(f'Sage exited {rc}; see sage.stderr.log')
        record['state'] = 'VALIDATING'
        save(folder / 'RUN.json', record)
        validation = validate_search(task)
        save(folder / 'VALIDATION.json', validation)
        record['state'] = 'COMPLETED'
    except Exception as exc:
        record.update(state='FAILED', error=str(exc))
        raise
    finally:
        record.update(ended_utc=now(), wall_seconds=time.monotonic()-start)
        save(folder / 'RUN.json', record)


def command_text(args):
    return shlex.join([str(x) for x in args])


def submit():
    manifest = json.loads((ROOT / 'manifest.json').read_text())
    if (ROOT / 'SUBMISSIONS.json').exists():
        raise ValueError('Already submitted; do not duplicate searches')
    submissions = dict(submitted_utc=now(), manifest_sha256=sha(ROOT / 'manifest.json'), jobs=[])
    save(ROOT / 'SUBMISSIONS.json', submissions)
    python = sys.executable
    # Two dependency chains enforce the global cap with per-job resources.
    by_id = {t['id']: t for t in manifest['tasks']}
    order = ['PREDICT_noclust', 'CAPSCAN_noclust', 'HSTOOL_noclust', 'CAMPI_FL_UNIPROT',
             'CAMPI_FL_REF', 'CAMPI_FL_REF_MG', 'CAMPI_REF', 'CAMPI_MG']
    previous_in_lane = [None, None]
    for index, task_id in enumerate(order):
        task = by_id[task_id]
        lane = index % 2
        resource = task['resources']
        path = Path(task['folder']) / 'run.sbatch'
        path.write_text('#!/bin/bash\nset -euo pipefail\n'
                        f'export RAYON_NUM_THREADS={resource["cpus"]}\n'
                        'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\n'
                        f'exec {command_text([python, __file__, "run", task["id"]])}\n')
        command = ['sbatch', '--parsable', '--partition='+resource['partition'], '--nodes=1', '--ntasks=1',
                   '--cpus-per-task='+str(resource['cpus']), '--mem='+resource['memory'], '--time='+resource['walltime'],
                   '--nice=10000', '--begin='+manifest['begin_local'], '--deadline='+manifest['deadline_local'],
                   '--job-name=ssmar_'+task['id'], '--chdir='+str(ROOT),
                   '--output='+str(ROOT / 'logs' / (task['id']+'_%j.out')),
                   '--error='+str(ROOT / 'logs' / (task['id']+'_%j.err'))]
        if previous_in_lane[lane]:
            command.append('--dependency=afterany:'+previous_in_lane[lane])
        command.append(str(path))
        job_id = subprocess.check_output(command, text=True).strip().split(';')[0]
        submissions['jobs'].append(dict(task=task['id'], job_id=job_id, lane=lane,
                                       predecessor=previous_in_lane[lane], command=command))
        previous_in_lane[lane] = job_id
        save(ROOT / 'SUBMISSIONS.json', submissions)
        print(task['id'], job_id, flush=True)
    # Always collect after all searches terminate, including killed/OOM/timeout jobs.
    dependency = 'afterany:' + ':'.join(j['job_id'] for j in submissions['jobs'])
    for mode, extra, wall in [('final', ['--dependency='+dependency], '00:30:00')]:
        script = ROOT / f'{mode}.sbatch'
        script.write_text('#!/bin/bash\nset -euo pipefail\n' +
                          f'exec {command_text([python, __file__, mode])}\n')
        command = ['sbatch', '--parsable', '--partition=p.hpcl92', '--nodes=1', '--ntasks=1',
                   '--cpus-per-task=1', '--mem=256M', '--time='+wall, '--nice=10000',
                   '--job-name=ssmar_'+mode, '--chdir='+str(ROOT),
                   '--output='+str(ROOT / 'logs' / (mode+'_%j.out')),
                   '--error='+str(ROOT / 'logs' / (mode+'_%j.err'))] + extra + [str(script)]
        job_id = subprocess.check_output(command, text=True).strip().split(';')[0]
        submissions[mode] = dict(job_id=job_id, command=command)
        save(ROOT / 'SUBMISSIONS.json', submissions)
        print(mode, job_id, flush=True)


def read_time(path):
    out = {}
    if path.exists():
        content = path.read_text()
        for name, pattern in [('max_rss_kb', r'Maximum resident set size \(kbytes\):\s*(\d+)'),
                              ('elapsed', r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(.+)'),
                              ('exit_status', r'Exit status:\s*(\d+)')]:
            match = re.search(pattern, content)
            if match:
                out[name] = match.group(1).strip()
    return out


def collect(final=False):
    manifest = json.loads((ROOT / 'manifest.json').read_text())
    submissions = json.loads((ROOT / 'SUBMISSIONS.json').read_text())
    ids = ','.join(j['job_id'] for j in submissions['jobs'])
    command = ['sacct', '-j', ids, '--noheader', '--parsable2',
               '--format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS,ReqMem,AllocCPUS,NodeList,Reason']
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError('Slurm accounting unavailable: '+result.stderr)
    accounting = {}
    for line in result.stdout.splitlines():
        values = line.split('|')
        if len(values) >= 10:
            accounting[values[0]] = values
    rows = []
    by_task = {j['task']: j['job_id'] for j in submissions['jobs']}
    for task in manifest['tasks']:
        folder = Path(task['folder'])
        job_id = by_task[task['id']]
        acct = accounting.get(job_id, [''] * 10)
        state = acct[2].split()[0].rstrip('+') if acct[2] else 'ACCOUNTING_PENDING'
        record = json.loads((folder / 'RUN.json').read_text()) if (folder / 'RUN.json').exists() else {}
        union = json.loads((folder / 'union_summary.json').read_text()) if (folder / 'union_summary.json').exists() else {}
        timing = read_time(folder / 'sage.time.txt')
        validated = (folder / 'VALIDATION.json').exists()
        status = 'PASS' if state == 'COMPLETED' and record.get('state') == 'COMPLETED' and validated else ('FAIL' if state in TERMINAL else 'IN_PROGRESS')
        row = dict(task=task['id'], job_id=job_id, scheduler_state=state, phase=record.get('state', 'QUEUED'),
                   status=status, acquisitions=task['acquisitions'], unique_sequences=union.get('unique_sequences', ''),
                   union_MB=round(union['union']['bytes']/1e6, 3) if union else '',
                   sage_elapsed=timing.get('elapsed', ''), sage_peak_RSS_GiB=round(int(timing['max_rss_kb'])/1048576, 3) if 'max_rss_kb' in timing else '',
                   slurm_elapsed=acct[4], slurm_max_RSS=accounting.get(job_id+'.batch', [''] * 10)[5],
                   requested_memory=task['resources']['memory'], exit_code=acct[3],
                   error=record.get('error', ''))
        rows.append(row)
    complete = len(rows) == len(manifest['tasks']) and all(r['scheduler_state'] in TERMINAL for r in rows)
    overall = 'PASS' if complete and all(r['status'] == 'PASS' for r in rows) else ('FAIL' if complete else 'IN_PROGRESS')
    report = dict(updated_utc=now(), status=overall, all_terminal=complete, tasks=rows)
    # Watcher and dependent collector have different output paths, avoiding races.
    dest = ROOT / ('collection' if final else 'monitor')
    dest.mkdir(exist_ok=True)
    save(dest / 'STATUS.json', report)
    (dest / 'sacct.tsv').write_text('JobID|JobName|State|ExitCode|Elapsed|MaxRSS|ReqMem|AllocCPUS|NodeList|Reason\n'+result.stdout)
    with (dest / 'resources.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    lines = ['# sample_specific_merge_all_raw', '', f'Updated {report["updated_utc"]}; status **{overall}**.', '',
             'One full-sequence-deduplicated FASTA and one joint Sage invocation per cohort/arm. Source: the validated simplified-grouping roster.', '',
             '| Cohort/arm | Job | State | Phase | Raw files | Unique sequences | FASTA MB | Sage elapsed | Sage peak GiB |',
             '|---|---|---|---|---:|---:|---:|---|---:|']
    for r in rows:
        lines.append('| '+' | '.join(str(r[k]) for k in ('task','job_id','scheduler_state','phase','acquisitions','unique_sequences','union_MB','sage_elapsed','sage_peak_RSS_GiB'))+' |')
    lines += ['', 'Elapsed time and peak RAM are measured with `/usr/bin/time -v`; Slurm accounting is saved alongside. Requested memory is a ceiling, not an observed result.', '',
              'Sage LFQ output is a peptide-by-raw matrix with matching protein accessions. It is not an inferred protein-group matrix. Protein/species inference is not established by this search benchmark.', '']
    (dest / 'REPORT.md').write_text('\n'.join(lines))
    print(json.dumps(report, indent=2), flush=True)
    return report


def watch():
    end = time.monotonic() + 36 * 3600
    while time.monotonic() < end:
        try:
            report = collect()
            with (ROOT / 'monitor/history.jsonl').open('a') as stream:
                stream.write(json.dumps(report)+'\n')
            if report['all_terminal']:
                return
        except Exception as exc:
            print(now(), 'MONITOR_ERROR', str(exc), flush=True)
        time.sleep(120)
    print('Periodic watch window ended; afterany final collector remains scheduled.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'build', 'run', 'submit', 'status', 'watch', 'final'])
    parser.add_argument('task', nargs='?')
    args = parser.parse_args()
    if args.action == 'prepare':
        prepare()
    elif args.action == 'build':
        _, task = load_task(args.task)
        build_union(task)
    elif args.action == 'run':
        run(args.task)
    elif args.action == 'submit':
        submit()
    elif args.action == 'watch':
        watch()
    else:
        report = collect(final=args.action == 'final')
        if args.action == 'final' and report['status'] != 'PASS':
            raise SystemExit(1)


if __name__ == '__main__':
    main()
