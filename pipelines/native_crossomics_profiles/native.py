"""Expanded cohort native geometry with frozen sequence hashing and statistics."""
from pathlib import Path
import gc
import importlib.util
import json
import sys
import numpy as np
import pandas as pd
from scipy import sparse

C = Path(__file__).resolve().parent
D = C.parent
R = Path(__import__('os').environ['FASTALAKE_DATA_DIR'])
N = C / 'native'
N.mkdir(exist_ok=True)
source = Path(json.loads((R / 'review/v15_integration/analyses/public_only_consistency_2026-09-29/CROSSOMICS_SOURCES.json').read_text())['native'])
spec = importlib.util.spec_from_file_location('frozen_native', source)
native = importlib.util.module_from_spec(spec)
spec.loader.exec_module(native)
old = source.parent / 'results/CAPSCAN'
original = json.loads((old / 'COMPLETE.json').read_text())
assert original['status'] == 'PASS' and original['code_sha256'] == native.sha(source)
read = native.read
mode = sys.argv[1]
assert mode in ['prepare', 'pilot', 'full']


def scenario_design(name):
    design = read(C / 'expanded_design.tsv', index_col=0)
    design['specimen'] = design.aliquot_main_id.astype(int).astype(str)
    membership = read(D / 'species/scenarios.tsv')
    ids = set(membership.loc[membership.scenario.eq(name), 'sample_id'])
    design = design.loc[design.index.isin(ids)]
    og = read(C / 'variants_full/CAPSCAN_OG_microbial_single_og.tsv.gz', index_col=0)
    og = og.drop(index='FunctionallyUnassigned', errors='ignore').T.reindex(design.index)
    valid = set(design.loc[og.sum(1) > 0, 'specimen'])
    met = read(C / 'metabolomics/CAPSCAN_metabolite_intensity.tsv', index_col=0)
    met.index = met.index.astype(str)
    met = met.loc[met.index.isin(valid)]
    dspec = design.groupby('specimen').first().loc[met.index]
    return design, dspec, met


def prepare():
    assert not (N / 'INPUTS.json').exists()
    assert json.loads((C / 'geometry_expanded/COMPLETE.json').read_text())['status'] == 'PASS'
    runs = read(native.MG / 'manifests/runs.tsv')
    runs.specimen = runs.specimen.astype(str)
    dna = read(native.MG / 'aggregate/metaG_OG.tsv', index_col=0)
    good = set(dna.index[dna.sum(1) > 0])
    populations = json.loads((C / 'geometry_expanded/COMPLETE.json').read_text())['populations']
    scenarios = {}
    allrows = {}
    for pop in populations:
        name, people = pop['scenario'], pop['participants']
        _, dm, _ = scenario_design(name)
        dm = dm[dm.subject_id.isin(people) & dm.group.isin(native.REG)]
        rows = []
        for r in runs[runs.specimen.isin(dm.index) & runs.run_accession.isin(good)].itertuples():
            folder = native.MG / 'results' / r.run_accession
            done = folder / ('genes.complete.json' if (folder / 'genes.complete.json').exists() else 'all.complete.json')
            row = dict(id=r.run_accession, specimen=r.specimen, participant=int(dm.loc[r.specimen, 'subject_id']),
                region=dm.loc[r.specimen, 'group'], tpm=native.rec(folder / 'salmon_quant/quant.sf'),
                fasta=native.rec(folder / 'genes.faa'), completion=native.rec(done))
            assert row['participant'] == int(r.subject_id) and row['region'] == r.group
            if row['id'] in allrows:
                assert row == allrows[row['id']]
            allrows[row['id']] = row
            rows.append(row['id'])
        assert {(allrows[x]['participant'], allrows[x]['region']) for x in rows} == {(p, r) for p in people for r in native.REG}
        scenarios[name] = dict(participants=people, sources=rows)
    oldmanifest = json.loads((source.parent / 'INPUT_MANIFEST.json').read_text())
    assert original['manifest_sha256'] == native.sha(source.parent / 'INPUT_MANIFEST.json')
    baseline = oldmanifest['items']['CAPSCAN']
    assert scenarios['current331']['sources'] == [r['id'] for r in baseline]
    for r in baseline:
        assert allrows[r['id']] == r
    for name in ['distance_matrices.npz', 'correlations.tsv']:
        assert native.sha(old / name) == original['outputs'][name]
    native.dump(N / 'INPUTS.json', dict(source=str(source), source_sha256=native.sha(source),
        worker_sha256=native.sha(__file__), old=str(old), scenarios=scenarios, rows=allrows,
        original_DNA_roster_reproduced=True))
    print('PREPARED', len(allrows), 'DNA sources;', len(scenarios), 'scenarios', flush=True)


def cached_source(co, row, out):
    # Reuse immutable source-level vectors only after exact source-record and cache-hash checks.
    shared = N / 'source_cache'
    folder = shared / 'sources'
    folder.mkdir(parents=True, exist_ok=True)
    rec = folder / (row['id'] + '.json')
    if not rec.exists():
        previous = old / 'sources' / rec.name
        if previous.exists():
            record = json.loads(previous.read_text())
            assert record['source'] == row and record['status'] == 'PASS'
            for suffix in ['.json', '.npz']:
                p = old / 'sources' / (row['id'] + suffix)
                (folder / p.name).symlink_to(p)
    for field in ['tpm', 'fasta', 'completion']:
        native.check(row[field])
    return source_vector(co, row, shared)


source_vector = native.source_vector
native.source_vector = cached_source


def expanded_other_layers(name, ci, out, matrix):
    design, dspec, met = scenario_design(name)
    frames = {}
    for key, frame in [('Non-host protein', matrix.loc[:, ~matrix.columns.str.contains('_HUMAN', regex=False)]),
                       ('Host protein', matrix.loc[:, matrix.columns.str.contains('_HUMAN', regex=False)])]:
        byspec = native.relative(frame.reindex(design.index)).groupby(design.specimen).median()
        ix = byspec.index.intersection(dspec.index)
        f = byspec.loc[ix].groupby([dspec.loc[ix, 'subject_id'], dspec.loc[ix, 'group']]).mean().loc[ci]
        assert np.isfinite(f.values).all() and (f.sum(1) > 0).all()
        frames[key] = sparse.csr_matrix(f.values)
        pd.Series(f.columns, name='representative').to_csv(out / (key.replace(' ', '_') + '_features.tsv.gz'), sep='\t', index=False)
    frames['Metabolites'] = np.log10(met.where(met > 0)).groupby([dspec.subject_id, dspec.group]).mean().loc[ci].values
    keep = dspec.apply(lambda r: (r.subject_id, r.group) in set(ci), axis=1)
    dspec.loc[keep, ['subject_id', 'group']].to_csv(out / 'protein_metabolite_specimen_roster.tsv', sep='\t')
    return frames


def run():
    m = json.loads((N / 'INPUTS.json').read_text())
    assert m['source_sha256'] == native.sha(source) and m['worker_sha256'] == native.sha(__file__)
    if mode == 'pilot':
        control = N / 'pilot'
        control.mkdir(exist_ok=True)
        ci = pd.MultiIndex.from_product([m['scenarios']['current331']['participants'], native.REG], names=['participant', 'region'])
        layers = native.other_layers('CAPSCAN', ci, control)
        saved = np.load(old / 'distance_matrices.npz')
        for label, threshold in [('primary_prevalence10', .1), ('sensitivity_all_features', 0.)]:
            for key, x in layers.items():
                dist, _ = native.profile_distance(x, key, threshold, capscan=True)
                assert np.allclose(dist, saved[label + '__' + key], atol=1e-8, rtol=1e-10)
        new = [r for r in m['rows'].values() if r['participant'] == 15]
        chosen = [next(iter(m['rows'].values())), max(new, key=lambda r: r['fasta']['bytes'])]
        for row in chosen:
            native.source_vector('CAPSCAN', row, control)
        native.dense_control()
        native.dump(control / 'COMPLETE.json', dict(status='PASS', source_ids=[r['id'] for r in chosen],
            existing_protein_metabolite_distances_reproduced=True, sparse_dense_control=True,
            worker_sha256=native.sha(__file__)))
        return
    pilot = json.loads((N / 'pilot/COMPLETE.json').read_text())
    assert pilot['status'] == 'PASS' and pilot['worker_sha256'] == native.sha(__file__)
    base = D / 'expanded_grouping/full/cohorts/CAPSCAN_noclust/reaggregated'
    matrix = read(base / 'simplified_study_group_matrix.tsv.gz', index_col=0).fillna(0).T
    assert np.isfinite(matrix.values).all() and (matrix.values >= 0).all()
    summaries = []
    for name, scenario in m['scenarios'].items():
        out = N / name
        out.mkdir(exist_ok=True)
        assert not (out / 'COMPLETE.json').exists()
        ci = pd.MultiIndex.from_product([scenario['participants'], native.REG], names=['participant', 'region'])
        rows = [m['rows'][s] for s in scenario['sources']]
        layers, geneqc = native.gene_layers('CAPSCAN', rows, out, ci)
        layers.update(expanded_other_layers(name, ci, out, matrix))
        archive, quality, results = {}, {}, []
        for label, threshold in [('primary_prevalence10', .1), ('sensitivity_all_features', 0.)]:
            ds, quality[label] = {}, {}
            for key, x in layers.items():
                ds[key], quality[label][key] = native.profile_distance(x, key, threshold, capscan=True)
                archive[label + '__' + key] = ds[key]
            if name == 'current331':
                previous = np.load(old / 'distance_matrices.npz')
                assert np.allclose(ds['DNA'], previous[label + '__DNA'], atol=1e-8, rtol=1e-10)
                assert np.allclose(ds['Metabolites'], previous[label + '__Metabolites'], atol=1e-8, rtol=1e-10)
            t = native.correlations(ds, 'CAPSCAN')
            t.insert(1, 'analysis', label)
            t['scenario'] = name
            results.append(t)
        result = pd.concat(results, ignore_index=True)
        result.to_csv(out / 'correlations.tsv', sep='\t', index=False)
        summaries.append(result)
        np.savez_compressed(out / 'distance_matrices.npz', **archive)
        native.dump(out / 'feature_quality.json', quality)
        native.dump(out / 'COMPLETE.json', dict(status='PASS', scenario=name, participants=scenario['participants'],
            n_profiles=len(ci), gene_profile=geneqc, all_observed_primary=True,
            scope='Expanded assembly-inclusive native cohort sensitivity; no cached cohort distances reused.'))
        print('SCENARIO PASS', name, flush=True)
        del layers
        gc.collect()
    pd.concat(summaries).to_csv(N / 'correlations_all_scenarios.tsv', sep='\t', index=False)
    native.dump(N / 'COMPLETE.json', dict(status='PASS', scenarios=len(summaries),
        worker_sha256=native.sha(__file__), input_manifest_sha256=native.sha(N / 'INPUTS.json'),
        scope='Native cohort sensitivities; requires scientific review, final cohort decision and manuscript integration.'))


if mode == 'prepare':
    prepare()
else:
    run()
