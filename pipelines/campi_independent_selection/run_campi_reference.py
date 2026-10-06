"""CAMPI acquisition A/B comparison from a common reference reservoir.

Developmental biological diagnostic; neither a true-FDR measurement nor a
replacement for the 233-million-sequence production-reservoir experiment.
"""
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[1]
FROZEN = PROJECT / 'analyses/fusion_entrapment_2026-09-07'
sys.path.insert(0, str(FROZEN))
from run_experiment import fasta, split_and_pairs, read_scores, estimates


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def run(argv, directory):
    directory.mkdir(parents=True, exist_ok=True)
    command = list(map(str, argv))
    start = time.time()
    with (directory / 'stdout.log').open('w') as stdout, (directory / 'stderr.log').open('w') as stderr:
        result = subprocess.run(command, stdout=stdout, stderr=stderr)
    record = dict(argv=command, returncode=result.returncode, seconds=time.time() - start)
    (directory / 'COMMAND.json').write_text(json.dumps(record, indent=2))
    if result.returncode:
        raise RuntimeError(f'Failed: {directory}')
    print(directory.relative_to(ROOT), round(record['seconds'], 2), flush=True)


def main():
    mode = sys.argv[1]
    assert mode in ('ordinary', '20260908', '20260909')
    source_manifest = json.loads((ROOT / 'CAMPI_SOURCE_SHA256.json').read_text())
    for name, value in source_manifest.items():
        assert sha(ROOT / name) == value
    cfg = json.loads((FROZEN / 'config.json').read_text())
    bins = {}
    for key in ('extractor', 'parsimony', 'sage', 'builder'):
        spec = cfg['binaries'][key]
        assert sha(spec['path']) == spec['sha256'], key
        bins[key] = spec['path']
    with (FROZEN / 'inputs/CAMPI_MOCK_samples.tsv').open() as stream:
        samples = {r['sample_id']: r for r in csv.DictReader(stream, delimiter='\t')
                   if r['sample_id'] in ('S01', 'S02')}
    assert set(samples) == {'S01', 'S02'}
    out = ROOT / 'campi_ref_v1' / mode
    out.mkdir(parents=True, exist_ok=False)
    source = PROJECT / 'projects/campi_reference_2026/0_refs/SIHUMI_DB1UNIPROT.faa'
    targets = [(a, s) for a, s in fasta(source)
               if not a.lower().startswith(('rev_', 'reverse_', 'decoy_'))]
    reference = out / 'reference.fasta'
    reference.write_text(''.join(f'>{a}\n{s}\n' for a, s in targets))
    reservoir = reference
    if mode != 'ordinary':
        reservoir = out / 'linked_reservoir.fasta'
        run([bins['builder'], reference, reservoir, mode, len(targets)], out / 'build_fusion')
    provenance = dict(source=str(source), source_sha256=sha(source), reference_targets=len(targets),
                      source_code=source_manifest, frozen_binaries=cfg['binaries'], mode=mode,
                      scope='fixed common SIHUMIx reference; developmental biological diagnostic',
                      known_truth='community compatibility only; not spectrum truth',
                      samples={s: {k: r[k] for k in ('denovo_csv', 'mzml', 'sage_config')} for s, r in samples.items()})
    for sample, r in samples.items():
        provenance['samples'][sample]['mzml_sha256'] = sha(r['mzml'])
        provenance['samples'][sample]['denovo_sha256'] = sha(r['denovo_csv'])
    assert provenance['samples']['S01']['mzml_sha256'] != provenance['samples']['S02']['mzml_sha256']
    (out / 'PROVENANCE.json').write_text(json.dumps(provenance, indent=2))
    threads = str(int(os.environ.get('SLURM_CPUS_PER_TASK', '4')))
    selections = {}
    for sample, r in samples.items():
        base = out / 'selection' / sample
        pool = base / 'predictions'
        pool.mkdir(parents=True)
        (pool / (sample + '_predictions.csv')).symlink_to(r['denovo_csv'])
        current = reservoir
        for stage in ('1_cohort', '3_acquisition'):
            dest = base / stage
            target = dest / 'selected.fasta'
            run([bins['extractor'], '--peptides-dir', pool, '--database', current, '--output', target,
                 '--output-hits', dest / 'hits.tsv', '--top-percent', '0.30', '--min-length', '9',
                 '--max-length', '50', '--normalize-il', 'true', '-t', threads], dest)
            current = target
        dest = base / '4_razor'
        run([bins['parsimony'], '--strategy', 'razor', '--tiebreak', 'hash-acc', '--fasta', current,
             '--peptides', r['denovo_csv'], '--output-dir', dest, '--sample', sample,
             '--min-length', '9', '--max-length', '50', '--normalize-il', '-t', threads], dest)
        selections[sample] = dest / (sample + '_razor.fasta')
    matrix, results = [], []
    # Create each search database once; fixed and each acquisition-selected set.
    databases = {'fixed': reservoir, **selections}
    pairs_by_db = {}
    if mode != 'ordinary':
        for key, linked in list(databases.items()):
            dest = out / 'search_inputs' / key
            dest.mkdir(parents=True)
            pairs, counts = split_and_pairs(linked, dest)
            pairs_by_db[key] = pairs
            databases[key] = dest / 'search.fasta'
    for evaluated, r in samples.items():
        independent = 'S02' if evaluated == 'S01' else 'S01'
        for arm, key in [('fixed', 'fixed'), ('self', evaluated), ('independent', independent)]:
            dest = out / 'searches' / (evaluated + '_' + arm)
            dest.mkdir(parents=True)
            search_config = json.loads(Path(r['sage_config']).read_text())
            search_config['database']['fasta'] = str(databases[key])
            search_config['output_directory'] = str(dest)
            search_config['mzml_paths'] = [r['mzml']]
            config_path = dest / 'sage.json'
            config_path.write_text(json.dumps(search_config, indent=2))
            run([bins['sage'], config_path, '--write-pin'], dest)
            result_path = dest / 'results.sage.tsv'
            record = dict(evaluated=evaluated, arm=arm, selection_acquisition=None if key == 'fixed' else key,
                          evaluated_spectra_enter_selection=(arm == 'self'), database_sha256=sha(databases[key]),
                          result=str(result_path), result_sha256=sha(result_path))
            matrix.append(record)
            if mode != 'ordinary':
                scores, classes = read_scores(result_path)
                for q in (.001, .005, .01, .02, .05):
                    for low in (5, 9):
                        results.append(dict(evaluated=evaluated, arm=arm, seed=int(mode),
                                            **estimates(scores, classes, pairs_by_db[key], q, low)))
    (out / 'SEARCH_MATRIX.json').write_text(json.dumps(matrix, indent=2))
    if results:
        with (out / 'diagnostics.tsv').open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(results[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(results)
    (out / 'COMPLETE.json').write_text(json.dumps(dict(searches=len(matrix), calibration_established=False), indent=2))


if __name__ == '__main__':
    main()
