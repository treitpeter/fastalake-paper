"""Submit/run the unchanged reviewed grouping module; validate matched counts."""
from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import importlib.util
import json
import math
import os
import statistics
import traceback
from collections import Counter, defaultdict
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_rows(path):
    with Path(path).open() as stream:
        yield from csv.DictReader(stream, delimiter="\t")


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def preflight(manifest, task):
    for path, expected in [(manifest['runtime_module'], manifest['runtime_module_sha256']),
                           (manifest['samples_file'], manifest['samples_sha256']),
                           (manifest['peptide_reference'], manifest['peptide_reference_sha256'])]:
        if sha(path) != expected:
            raise ValueError(f"Frozen input changed: {path}")
    names = Path(manifest['samples_file']).read_text().splitlines()
    if len(names) != task['acquisitions'] or len(set(names)) != len(names):
        raise ValueError('Acquisition roster changed')
    for record in task['search_inputs']:
        stat = Path(record['path']).stat()
        if (stat.st_size, stat.st_mtime_ns) != (record['bytes'], record['mtime_ns']):
            raise ValueError(f"Search input changed: {record['path']}")
    for name in names:
        actual = Path(task['search_root']) / name
        expected = Path(task['source_search_root']) / name / task['source_subdirectory']
        if actual.resolve() != expected.resolve():
            raise ValueError(f'Input link changed: {actual}')
    return names


def validate_output(manifest, task):
    names = Path(manifest['samples_file']).read_text().splitlines()
    out = Path(task['output'])
    summary = json.loads((out / 'summary.json').read_text())
    if sorted(summary['acquisition_names']) != sorted(names):
        raise ValueError('Output cohort mismatch')
    if summary['q_threshold'] != 0.01 or not summary['representative_sequence_compatibility_checked']:
        raise ValueError('Grouping settings or sequence validation mismatch')
    records = list(read_rows(out / 'acquisition_summary.tsv'))
    if len(records) != len(names) or {r['sample'] for r in records} != set(names):
        raise ValueError('Acquisition summaries missing or duplicated')
    counts, features, intensities = Counter(), Counter(), defaultdict(float)
    previous = None
    for row in read_rows(out / 'study_group_long.tsv'):
        key = row['sample'], row['study_group']
        if previous is not None and key <= previous:
            raise ValueError('Duplicate or unordered group observations')
        previous = key
        value = float(row['intensity'])
        if not math.isfinite(value) or value <= 0 or int(row['n_peptides']) < 1:
            raise ValueError('Invalid quantified group')
        counts[row['sample']] += 1
        features[row['sample']] += int(row['n_features'])
        intensities[row['sample']] += value
    for row in records:
        name = row['sample']
        if counts[name] != int(row['quantified_study_groups']) or features[name] != int(row['lfq_features']):
            raise ValueError(f'Independent group/feature recount failed: {name}')
        if not math.isclose(intensities[name], float(row['intensity']), rel_tol=1e-10, abs_tol=1e-6):
            raise ValueError(f'Independent intensity recount failed: {name}')
    if not math.isclose(summary['intensity_in'], summary['intensity_out'], rel_tol=1e-10, abs_tol=1e-6):
        raise ValueError('Intensity conservation failed')
    module_records = [r for r in json.loads((out/'input_manifest.json').read_text())
                      if Path(r['path']).resolve() == Path(manifest['runtime_module']).resolve()]
    if len(module_records) != 1 or module_records[0]['sha256'] != manifest['runtime_module_sha256']:
        raise ValueError('Grouping module provenance mismatch')
    return {r['sample']: int(r['quantified_study_groups']) for r in records}


def run(manifest, task):
    names = preflight(manifest, task)
    folder = Path(task['output']).parent
    folder.mkdir(parents=True, exist_ok=True)
    try:
        print(f"Starting {task['id']}: {len(names)} acquisitions; unchanged simplified grouping", flush=True)
        spec = importlib.util.spec_from_file_location('reviewed_study', manifest['runtime_module'])
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.group_searches(task['search_root'], task['output'], q_threshold=0.01,
                              config_name=task['config_name'], acquisition_names=names)
        counts = validate_output(manifest, task)
        write_json(folder/'COMPLETE.json', dict(status='PASS', arm=task['id'], acquisitions=len(counts),
                   median_quantified_groups=statistics.median(counts.values()),
                   module_sha256=manifest['runtime_module_sha256'], slurm_job=os.environ.get('SLURM_JOB_ID'),
                   completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        print((folder/'COMPLETE.json').read_text(), flush=True)
    except Exception as error:
        write_json(folder/'FAILED.json', dict(status='FAIL', arm=task['id'], error=str(error), traceback=traceback.format_exc()))
        raise


def collect(manifest):
    out = Path(manifest['output_root'])/'comparison'
    out.mkdir(exist_ok=True)
    arms, failures = {}, {}
    for task in manifest['tasks']:
        try:
            preflight(manifest, task)
            done = json.loads((Path(task['output']).parent/'COMPLETE.json').read_text())
            if done['status'] != 'PASS':
                raise ValueError('Task did not pass')
            arms[task['id']] = validate_output(manifest, task)
        except Exception as error:
            failures[task['id']] = str(error)
    if failures:
        write_json(out/'STATUS.json', dict(status='FAIL', failures=failures))
        raise ValueError(f'Comparison incomplete: {failures}')
    names = Path(manifest['samples_file']).read_text().splitlines()
    medians = {arm: statistics.median(counts.values()) for arm, counts in arms.items()}
    ratios = {}
    for arm, counts in arms.items():
        if arm == 'FL':
            continue
        paired = [arms['FL'][s]/counts[s] for s in names if counts[s] > 0]
        ratios[arm] = dict(ratio_of_medians=medians['FL']/medians[arm] if medians[arm] else None,
                          median_paired_ratio=statistics.median(paired) if paired else None,
                          paired_ratio_positive_denominator_n=len(paired),
                          zero_denominator_n=sum(counts[s] == 0 for s in names),
                          FastaLake_higher_n=sum(arms['FL'][s] > counts[s] for s in names),
                          equal_n=sum(arms['FL'][s] == counts[s] for s in names))
    with (out/'matched_protein_group_counts.tsv').open('w') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['sample', *[f'groups_{arm}' for arm in arms]])
        writer.writerows([s, *[counts[s] for counts in arms.values()]] for s in names)
    result = dict(status='PASS', acquisitions=len(names), median_quantified_groups=medians,
                  comparisons=ratios, module_sha256=manifest['runtime_module_sha256'],
                  population='Same acquisition roster used for study inference separately in every strategy',
                  zeros='Included in count medians and paired wins; paired ratios use positive comparator counts only',
                  endpoint='Positive quantified groups under direct canonical-peptide assignment; groups do not correspond one-to-one across strategies')
    write_json(out/'SUMMARY.json', result)
    write_json(out/'STATUS.json', dict(status='PASS', acquisitions=len(names)))
    lines = ['# Matched simplified protein-group comparison', '', f"All strategies use the same {len(names)} acquisitions.", '',
             '| Strategy | Median quantified groups | FL / comparator medians |', '|---|---:|---:|']
    for arm, value in medians.items():
        ratio = ratios.get(arm, {}).get('ratio_of_medians')
        lines.append(f"| {arm} | {value:,.1f} | {ratio:.4f} |" if ratio is not None else f"| {arm} | {value:,.1f} | — |")
    lines += ['', result['endpoint'] + '.', result['zeros'] + '.', '', 'See SUMMARY.json for paired ratios, denominator counts and paired wins.',
              'This reports post-search quantitative units; it does not establish unique protein/organism identity or calibrated group FDR.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['run', 'collect', 'preflight'])
    parser.add_argument('manifest', type=Path)
    parser.add_argument('--arm')
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    if args.action == 'collect':
        collect(manifest)
    elif args.action == 'preflight':
        for task in manifest['tasks']:
            print(task['id'], len(preflight(manifest, task)), 'PASS')
    else:
        run(manifest, next(t for t in manifest['tasks'] if t['id'] == args.arm))
