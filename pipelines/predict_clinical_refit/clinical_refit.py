"""Frozen-roster PREDICT compartment sensitivity; original inference reproduced first."""
from pathlib import Path
from itertools import combinations
import json, hashlib, os, sys
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu

H = Path(__file__).resolve().parent
ROOT = Path(__import__('os').environ['FASTALAKE_DATA_DIR'])
OLD = ROOT/'MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/PREDICT_repaired__simplified/source_data/Figure4'
NEW = H/'profiles/PREDICT_noclust'
OUTPUTS = ['visit1_patient_profiles.tsv', 'paired_antibiotic_profiles.tsv',
           'clinical_paired_statistics.tsv', 'clinical_group_tests.tsv',
           'clinical_group_mean_intervals.tsv']

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8<<20), b''): h.update(block)
    return h.hexdigest()

def load(path): return json.loads(Path(path).read_text())
def dump(path, data): Path(path).write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')

def bh(p):
    p = np.asarray(p, dtype=float)
    assert np.isfinite(p).all()
    order = np.argsort(p); ranked = p[order]*len(p)/np.arange(1, len(p)+1)
    out = np.empty_like(p); out[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return out

def derive(ac, classes, out):
    # Transcribed from frozen derive_clinical_pc.py, retaining ordering/seeds/tails.
    metrics = classes + ['pg_parsimony']
    base = ac[ac.visit.eq(1)].groupby(['patient_id', 'clinical_group'])[metrics].mean().reset_index()
    assert len(base) == base.patient_id.nunique() == 91
    states = ac[ac.antibiotic_state.ne('Unknown')].groupby(['patient_id', 'antibiotic_state'])
    paired = states[classes].mean().join(states.pg_parsimony.median()).join(states.size().rename('n_acquisitions')).reset_index()
    eligible = paired.groupby('patient_id').antibiotic_state.nunique()
    paired = paired[paired.patient_id.isin(eligible[eligible.eq(2)].index)]
    assert len(paired) == 84 and paired.patient_id.nunique() == 42
    stats = []
    for idx, (metric, scale, unit) in enumerate([('Human', 100, 'percentage points'), ('pg_parsimony', 1, 'quantified simplified groups')]):
        wide = paired.pivot(index='patient_id', columns='antibiotic_state', values=metric)[['Off', 'On']]*scale
        delta = (wide.On-wide.Off).to_numpy()
        rng = np.random.default_rng(20260907)
        boots = delta[rng.integers(0, len(delta), (50000, len(delta)))].mean(axis=1)
        lo, hi = np.quantile(boots, [.025, .975])
        rng = np.random.default_rng(20260907+idx); count = 0
        for _ in range(100):
            null = ((rng.integers(0, 2, (2000, len(delta)))*2-1)*delta).mean(axis=1)
            count += int((abs(null) >= abs(delta.mean())-1e-12).sum())
        stats.append(dict(metric=metric, n_patients=42, mean_off=wide.Off.mean(), mean_on=wide.On.mean(),
                          mean_on_minus_off=delta.mean(), CI95_low=lo, CI95_high=hi, units=unit,
                          p_signflip_MC=(count+1)/200001))
    stats = pd.DataFrame(stats); stats['q_two_clinical_endpoints'] = bh(stats.p_signflip_MC)
    order = ['SDC', 'UDC', 'pre-ACLF', 'ACLF']
    vectors = {g: base.loc[base.clinical_group.eq(g), 'Human'].to_numpy()*100 for g in order}
    kg = kruskal(*vectors.values())
    dump(out/'clinical_group_global.json', dict(test='Kruskal-Wallis', statistic=float(kg.statistic), p=float(kg.pvalue)))
    tests = []
    for a, b in combinations(order, 2):
        test = mannwhitneyu(vectors[a], vectors[b], alternative='two-sided', method='asymptotic', use_continuity=True)
        tests.append(dict(group_a=a, group_b=b, n_a=len(vectors[a]), n_b=len(vectors[b]), p=float(test.pvalue)))
    tests = pd.DataFrame(tests); tests['q_BH_six'] = bh(tests.p)
    intervals = []
    for i, g in enumerate(order):
        vals = vectors[g]; rng = np.random.default_rng(20260908+i)
        lo, hi = np.quantile(vals[rng.integers(0, len(vals), (50000, len(vals)))].mean(axis=1), [.025, .975])
        intervals.append(dict(group=g, n=len(vals), mean=vals.mean(), median=np.median(vals), CI95_low=lo, CI95_high=hi))
    for name, frame in zip(OUTPUTS, [base, paired, stats, tests, pd.DataFrame(intervals)]):
        frame.to_csv(out/name, sep='\t', index=False)
    return stats

def compare_file(old, new, columns=None):
    a = pd.read_csv(old, sep='\t'); b = pd.read_csv(new, sep='\t')
    if columns is not None: a, b = a[columns], b[columns]
    assert list(a.columns) == list(b.columns) and a.shape == b.shape, str(old)
    errors = {}
    for c in a.columns:
        if pd.api.types.is_numeric_dtype(a[c]):
            assert np.allclose(a[c], b[c], rtol=1e-10, atol=1e-12, equal_nan=True), (str(old), c)
            errors[c] = float(np.nanmax(abs(a[c]-b[c])))
        else: assert a[c].equals(b[c]), (str(old), c)
    return errors

def main(mode):
    assert mode in ['baseline', 'refit']
    manifest = load(H/'CLINICAL_INPUTS.json')
    for path, expected in manifest['files'].items(): assert sha(path) == expected, path
    code = load(H/'CLINICAL_CODE.json')
    for name, expected in code.items(): assert sha(H/name) == expected, name
    ac = pd.read_csv(Path(os.environ['FASTALAKE_CONTROLLED_DIR'])/'acquisition_accounting.tsv', sep='\t')
    assert len(ac) == ac.sample_id.nunique() == 306
    out = H/('clinical_baseline' if mode == 'baseline' else 'clinical_refit')
    out.mkdir(exist_ok=False)
    checks = {}
    classes = ['Human', 'Microbial', 'Mixed', 'Unknown']
    if mode == 'refit':
        gate = load(H/'clinical_baseline/COMPLETE.json')
        assert gate['status'] == 'PASS' and gate['code_manifest_sha256'] == sha(H/'CLINICAL_CODE.json')
        assert gate['input_manifest_sha256'] == sha(H/'CLINICAL_INPUTS.json')
        tax = pd.read_csv(NEW/'taxonomy_signal.tsv', sep='\t').set_index('sample')
        assert tax.index.is_unique and set(tax.index) == set(ac.sample_id)
        tax = tax.loc[ac.sample_id]
        assert np.allclose(ac.total_LFQ_intensity, tax.total_LFQ, rtol=1e-9, atol=1e-5)
        classes = ['Human', 'Microbial_species', 'Microbial_above_species', 'Other_annotated_origins', 'Unresolved_origin']
        labels = ['Human-assigned', 'Microbial species', 'Microbial above species', 'Other annotated origins', 'Unresolved origin']
        ac = ac.drop(columns=['Human', 'Microbial', 'Mixed', 'Unknown'])
        for cls, label in zip(classes, labels):
            ac[cls] = tax[label].to_numpy()/tax.total_LFQ.to_numpy() if label in tax else 0.
        assert np.allclose(ac[classes].sum(axis=1), 1, rtol=1e-9, atol=1e-12)
        checks['signal_and_acquisition_roster_match_A33'] = True
        checks['composition_labels'] = dict(zip(classes, labels))
    ac.to_csv(out/'acquisition_accounting.tsv', sep='\t', index=False)
    stats = derive(ac, classes, out)
    if mode == 'baseline':
        for name in OUTPUTS: checks[name] = compare_file(OLD/name, out/name)
        old = load(OLD/'clinical_group_global.json'); new = load(out/'clinical_group_global.json')
        assert old['test'] == new['test']
        assert all(np.isclose(old[k], new[k], rtol=1e-10, atol=1e-12) for k in ['statistic', 'p'])
        checks['global_test_reproduced'] = True
    else:
        checks['paired_population_and_depth'] = compare_file(OLD/'paired_antibiotic_profiles.tsv', out/'paired_antibiotic_profiles.tsv', ['patient_id', 'antibiotic_state', 'pg_parsimony', 'n_acquisitions'])
        checks['baseline_population_and_depth'] = compare_file(OLD/'visit1_patient_profiles.tsv', out/'visit1_patient_profiles.tsv', ['patient_id', 'clinical_group', 'pg_parsimony'])
        before = pd.read_csv(OLD/'clinical_paired_statistics.tsv', sep='\t')
        before.merge(stats, on='metric', suffixes=('_original', '_replacement'), validate='one_to_one').to_csv(out/'paired_method_comparison.tsv', sep='\t', index=False)
        before = pd.read_csv(OLD/'clinical_group_tests.tsv', sep='\t')
        before.merge(pd.read_csv(out/'clinical_group_tests.tsv', sep='\t'), on=['group_a', 'group_b'], suffixes=('_original', '_replacement'), validate='one_to_one').to_csv(out/'group_method_comparison.tsv', sep='\t', index=False)
    import scipy
    dump(out/'COMPLETE.json', dict(status='PASS', mode=mode, utc=datetime.now(timezone.utc).isoformat(),
        job_id=os.environ.get('SLURM_JOB_ID'), code_manifest_sha256=sha(H/'CLINICAL_CODE.json'), input_manifest_sha256=sha(H/'CLINICAL_INPUTS.json'),
        checks=checks, acquisitions=306, paired_patients=42, baseline_patients=91,
        scope='Compartment-classification sensitivity on identical accepted signal. Original simplified-group counts retained; no regrouping, feature selection, EDS or cohort endpoint closure.',
        inference='Original two-sided 200000 sign flips, 50000 bootstrap, BH two paired endpoints and BH six clinical pairwise comparisons.',
        environment=dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
        outputs={p.name: sha(p) for p in out.iterdir() if p.is_file()}))
    print('PASS', mode, flush=True)

if __name__ == '__main__': main(sys.argv[1])
