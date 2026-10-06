"""Pilot-gated CAMPI control on the exact broad-UniProt production inputs."""
import collections
import csv
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time

import run_experiment as m
import audit_fusion as audit

ROOT = Path(__file__).resolve().parent
PROV = json.loads((ROOT / 'inputs/PROVENANCE.json').read_text())
SEEDS = (20260907, 20260908)

def verified_inputs():
    m.verify_code()
    for r in PROV['files']:
        p = Path(r['path']); s = p.stat()
        assert (s.st_size, s.st_mtime_ns) == (r['bytes'], r['mtime_ns']), ('input changed', str(p))
        if 'sha256' in r: assert m.sha(p) == r['sha256'], ('input hash changed', str(p))

def preflight():
    verified_inputs()
    out = ROOT / 'preflight'; out.mkdir(exist_ok=False)
    with (out / 'fixture_tests.log').open('w') as f:
        subprocess.run([sys.executable, '-m', 'unittest', '-v', 'test_experiment'],
                       cwd=ROOT, stdout=f, stderr=subprocess.STDOUT, check=True)
    print('Existing fusion-construction and entrapment-only-retention fixtures passed.', flush=True)
    # Reproduce unmodified production selection before releasing the new intervention.
    r = m.samples('CAMPI_MOCK')[0]
    dest = out / 'unmodified_selection'
    m.run([m.binary('parsimony'), '--strategy', 'razor', '--tiebreak', 'hash-acc',
           '--fasta', r['candidates'], '--peptides', r['denovo_csv'], '--output-dir', dest,
           '--sample', r['sample_id'], '--min-length', '9', '--max-length', '50',
           '--normalize-il', '-t', os.environ.get('SLURM_CPUS_PER_TASK', '8')], dest)
    generated = dest / (r['sample_id'] + '_razor.fasta')
    expected = collections.Counter(seq for _, seq in m.fasta(r['db_fasta']))
    observed = collections.Counter(seq for _, seq in m.fasta(generated))
    assert expected == observed, 'Unmodified production protein-sequence multiset did not reproduce'
    print('Production S01 selection reproduced:', sum(observed.values()), 'entries.', flush=True)
    hashes = {}
    # Full reservoir hash runs on the compute node; its 136-GB contents are never staged on login.
    for p in [PROV['reservoir']['path']] + [r[k] for r in m.samples('CAMPI_MOCK') for k in ('denovo_csv', 'mzml')]:
        hashes[p] = m.sha(p)
    verified_inputs()
    m.save(out / 'PASS.json', dict(status='PASS', job_id=os.getenv('SLURM_JOB_ID'),
           baseline_sample=r['sample_id'], baseline_protein_entries=sum(observed.values()),
           baseline_sequence_multiset_reproduced=True, fixture_tests_passed=True,
           input_sha256=hashes, original_reservoir_sequences=375333654,
           config_sha256=m.sha(ROOT/'config.json'), source_manifest_sha256=m.sha(ROOT/'SOURCE_SHA256.json')))
    print('Preflight PASS.', flush=True)

def prepare(index):
    verified_inputs()
    gate = json.loads((ROOT / 'preflight/PASS.json').read_text()); assert gate['status'] == 'PASS'
    seed = SEEDS[index]
    m.prepare('CAMPI_MOCK', seed)
    path = ROOT / 'runs_v2/CAMPI_MOCK' / str(seed) / 'PREPARED.json'
    prepared = json.loads(path.read_text())
    prepared['original_reservoir_sha256'] = gate['input_sha256'][PROV['reservoir']['path']]
    prepared['production_arm'] = PROV['benchmark_root']
    m.save(path, prepared)

def sample(flat_index):
    verified_inputs()
    sample_index, seed_index = divmod(flat_index, 2)
    seed = SEEDS[seed_index]
    r = m.samples('CAMPI_MOCK')[sample_index]
    m.execute_sample(sample_index, seed)
    base = ROOT / 'runs_v2/CAMPI_MOCK' / str(seed) / 'samples' / r['sample_id']
    original_cfg = json.loads(Path(r['sage_config']).read_text())
    actual_cfg = json.loads((base / 'search_input/sage.json').read_text())
    for cfg in [original_cfg, actual_cfg]:
        cfg['database'].pop('fasta'); cfg.pop('output_directory')
    assert actual_cfg == original_cfg, 'Search settings differ beyond FASTA/output path'
    complete = json.loads((base/'COMPLETE.json').read_text())
    assert complete['razor']['both_members_retained_in_every_group']
    task = dict(index=flat_index, cohort='CAMPI_MOCK', sample=r['sample_id'], seed=seed,
                directory=str(base), expected_results_sha256=complete['results_sha256'])
    audit.source_commit = lambda: 'frozen-sha256:' + m.sha(ROOT/'audit_fusion.py')
    audit.run_task(task, base/'independent_audit')
    independent = json.loads((base/'independent_audit/COMPLETE.json').read_text())
    assert independent['independent_accounting_matches'] and independent['independent_final_pairing_matches']
    assert any(x['n_canonical'] > 0 for x in complete['summary'])
    verified_inputs()
    m.save(base/'VALIDATED.json', dict(status='PASS', sample=r['sample_id'], seed=seed,
           job_id=os.getenv('SLURM_JOB_ID'), chemistry_and_search_settings_match_production=True,
           source_inputs_unchanged=True, complete_sha256=m.sha(base/'COMPLETE.json'),
           independent_audit_sha256=m.sha(base/'independent_audit/COMPLETE.json')))
    print('Search and independent audit PASS:', r['sample_id'], seed, flush=True)

def collect():
    verified_inputs()
    m.collect()
    rows, completions = [], []
    for seed in SEEDS:
        for r in m.samples('CAMPI_MOCK'):
            base = ROOT/'runs_v2/CAMPI_MOCK'/str(seed)/'samples'/r['sample_id']
            v = json.loads((base/'VALIDATED.json').read_text()); assert v['status'] == 'PASS'
            assert v['complete_sha256'] == m.sha(base/'COMPLETE.json')
            assert v['independent_audit_sha256'] == m.sha(base/'independent_audit/COMPLETE.json')
            c = json.loads((base/'independent_audit/COMPLETE.json').read_text())
            assert c['outputs']['decomposition.tsv'] == m.sha(base/'independent_audit/decomposition.tsv')
            with (base/'independent_audit/decomposition.tsv').open() as f:
                rows.extend(csv.DictReader(f, delimiter='\t'))
            completions.append({'path':str(base/'VALIDATED.json'),'sha256':m.sha(base/'VALIDATED.json')})
    assert len(completions) == 30 and len(rows) == 450
    m.table(ROOT/'results/independent_per_acquisition.tsv', rows)
    groups = collections.defaultdict(list)
    for row in rows:
        groups[(int(row['seed']),float(row['q_threshold']),int(row['min_residues']),int(row['max_residues']))].append(row)
    pooled = []
    for (seed,q,low,high), records in sorted(groups.items()):
        assert len(records) == 15
        total = sum(int(x['n_canonical']) for x in records)
        numerator = sum(int(x['paired_numerator']) for x in records)
        unknown = sum(int(x['n_unpaired_accepted']) for x in records)
        individual = [float(x['paired_diagnostic']) for x in records if x['paired_diagnostic']]
        pooled.append(dict(seed=seed,q_threshold=q,min_residues=low,max_residues=high,
                      acquisitions=15,n_canonical=total,paired_numerator=numerator,n_unpaired_accepted=unknown,
                      pooled_diagnostic=(numerator+unknown)/total if total else None,
                      paired_component=numerator/total if total else None,
                      unpaired_component=unknown/total if total else None,
                      median_acquisition_diagnostic=statistics.median(individual) if individual else None))
    m.table(ROOT/'results/pooled_by_seed.tsv', pooled)
    m.save(ROOT/'results/VALIDATION.json', dict(status='PASS', expected_searches=30,validated_searches=30,
           reservoir_sequences=375333654, completion_receipts=completions,
           outputs={p.name:m.sha(p) for p in (ROOT/'results').glob('*.tsv')},
           scientific_review='pending; matched control diagnostics do not establish target-only whole-procedure FDR'))
    lines=['# Broad-reservoir CAMPI fusion-control results','',
           'Thirty searches and independent pairing/accounting audits passed. Scientific review remains pending.',
           '', '| Seed | Length | Pooled diagnostic | Paired component | Unpaired component |',
           '|---|---|---:|---:|---:|']
    for r in pooled:
        if r['q_threshold'] == .01:
            lines.append(f"| {r['seed']} | {r['min_residues']}–{r['max_residues']} | {100*r['pooled_diagnostic']:.3f}% | {100*r['paired_component']:.3f}% | {100*r['unpaired_component']:.3f}% |")
    (ROOT/'results/REPORT.md').write_text('\n'.join(lines)+'\n')
    print('All 30 matched-configuration searches and audits complete.', flush=True)

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'preflight': preflight()
    elif mode == 'prepare': prepare(int(os.environ['SLURM_ARRAY_TASK_ID']))
    elif mode == 'sample': sample(int(os.environ['SLURM_ARRAY_TASK_ID']))
    elif mode == 'collect': collect()
    else: raise ValueError(mode)
