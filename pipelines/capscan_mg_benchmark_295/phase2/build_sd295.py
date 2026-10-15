#!/usr/bin/env python3
"""Phase 2, step 2: write the 295-acquisition versions of every CapScan-matched Source Data table into
phase2/repo/source_data (a copy of release_prep/fastalake-paper_v28). The 282 existing rows are kept
byte-for-byte (lines copied); the 13 new rows are appended, computed with the same definitions:

  Fig4/same_spectrum_benchmark.tsv               +26  (5-50, 9-50; FL_public vs MG, frozen Confidence reader)
  SFig06/same_spectrum_benchmark_1_capscan.tsv   +52  (FL_public/MG and FL_cohortMG/MG, 5 and 9)
  supplementary_methods_6_7/confidence_ordinary_all_pairs.tsv  +468 (13 x 36 gate rows)
  Fig4/capscan_dna_yield_by_acquisition.tsv      +13  (5-50; assembly_proteins = sanitizer n_kept; DNA not recorded
                                                      for participant-15 capsule specimens -> empty)
  SFig09/capscan_acquisitions_9_50.tsv, capscan_assembly_dna_yield_gain_by_acquisition.tsv, a08_benchmark282_{old,new}.tsv  +13
  SData1/v28d_capscan_mg_pairs.tsv               +13
  Fig4/crossfit_per_acquisition.tsv, supplementary_methods_6_7/crossfit_per_acquisition.tsv  +13  (--crossfit)
     per-row definitions copied from public347 run/analyze_final.py and run/recovery.py; odd-scan gains from
     spectrum_split347 audit (FULL/HALF_SELF/HALF_DISJOINT) over mg295 MG_ODD, as iteration3/disjoint/analyze.py.
"""
import csv, json, sys, math
from pathlib import Path
import numpy as np
M = Path(__file__).resolve().parents[1]; FINAL = M.parent
SD = M / 'phase2/repo/source_data'
P = M / 'phase2/per_acq13'
LAB = Path('LAB_ROOT')
SAN = {r['sample']: r for r in csv.DictReader(open(LAB / 'canonical/publication_data/sources/assemblies_prefixed/CAPSCAN_sanitize_report.tsv'), delimiter='\t')}
acq = list(csv.DictReader(open(M / 'inputs/ACQ13.tsv'), delimiter='\t'))
NEW = [a['sample'] for a in acq]; A = {a['sample']: a for a in acq}
pairs = {(r['sample'], r['test_arm'], int(r['min_length'])): r for r in csv.DictReader(open(P / 'pairs.tsv'), delimiter='\t')}
conf = [r for r in csv.DictReader(open(P / 'confidence_rows.tsv'), delimiter='\t') if r['sample'] in NEW]
assert len(conf) == 13 * 36


def append(rel, rows, expect_old):
    p = SD / rel
    lines = p.read_text().splitlines(keepends=True)
    hdr = lines[0].rstrip('\n').split('\t')
    assert sum(1 for l in lines[1:] if 'T760' in l) == expect_old or expect_old is None, (rel, sum(1 for l in lines[1:] if 'T760' in l))
    assert not any(s in ''.join(lines) for s in NEW), (rel, 'already contains a new sample')
    assert lines[-1].endswith('\n')
    with open(p, 'a', newline='') as f:
        w = csv.DictWriter(f, fieldnames=hdr, delimiter='\t', lineterminator='\n', extrasaction='raise')
        for r in rows:
            assert set(r) == set(hdr), (rel, set(hdr) ^ set(r))
            w.writerow(r)
    print('appended', len(rows), 'to', rel)


def pr(s, arm, L): return pairs[(s, arm, L)]
def ratio(a, b): return a / b


def same_spectrum():
    rows = []
    for s in NEW:
        for L in (5, 9):
            r = pr(s, 'FL_public', L); t, f, o = int(r['n_test']), int(r['n_reference']), int(r['overlap'])
            rows.append(dict(cohort='CAPSCAN', acquisition_id=s, participant_id=15, sampling_site=A[s]['region'], min_length=L,
                             n_fl=t, n_mg=f, overlap=o, count_ratio=repr(t / f), mg_recovery_fraction=repr(o / f)))
    append('Fig4/same_spectrum_benchmark.tsv', rows, 564)
    rows = []
    for s in NEW:
        for arm in ('FL_public', 'FL_cohortMG'):
            for L in (5, 9):
                r = pr(s, arm, L); t, f, o = int(r['n_test']), int(r['n_reference']), int(r['overlap'])
                rows.append(dict(cohort='CAPSCAN', sample=s, subject_id=15, group=A[s]['region'], aliquot='pellet', test_arm=arm, reference_arm='MG',
                                 min_length=L, n_test=t, n_reference=f, overlap=o, gained=t - o, lost=f - o, count_ratio=repr(t / f),
                                 reference_recovery=repr(o / f), test_recovery=repr(o / t)))
    append('SFig06/same_spectrum_benchmark_1_capscan.tsv', rows, 1790)
    rows = []
    for r in conf:
        n_fl, n_mg = int(r['n_fl']), int(r['n_mg'])
        base = next(x for x in conf if x['sample'] == r['sample'] and x['gate'] == 'q_only' and x['q'] == r['q'] and x['min_length'] == r['min_length'] and x['max_length'] == r['max_length'])
        rows.append(dict(r, gain_pct=repr(100 * (n_fl / n_mg - 1)) if n_mg else 'nan',
                         fl_retained=repr(n_fl / int(base['n_fl'])) if int(base['n_fl']) else 'nan',
                         mg_retained=repr(n_mg / int(base['n_mg'])) if int(base['n_mg']) else 'nan'))
    append('supplementary_methods_6_7/confidence_ordinary_all_pairs.tsv', rows, 10152)


def dna_tables():
    def row(s, L, empty):
        a = pr(s, 'FL_public', L); a5 = pr(s, 'FL_public', 5)
        fl, mg = int(a['n_test']), int(a['n_reference'])
        return dict(sample=s, region=A[s]['region'], subject=15, biosample=A[s]['biosample'], assembly_proteins=int(SAN[A[s]['biosample']]['n_kept']),
                    dna_ng=empty, volume_ul=empty, FL_peptides=fl, MG_peptides=mg, gain_percent=round(100 * (fl / mg - 1), 3),
                    MG_identity_recovery=round(int(a5['overlap']) / int(a5['n_reference']), 4))
    append('Fig4/capscan_dna_yield_by_acquisition.tsv', [row(s, 5, '') for s in NEW], 282)
    append('SFig09/capscan_assembly_dna_yield_gain_by_acquisition.tsv', [row(s, 5, 'nan') for s in NEW], 282)
    append('SFig09/capscan_acquisitions_9_50.tsv', [row(s, 9, '') for s in NEW], 282)
    b = []
    for s in NEW:
        a5, a9 = pr(s, 'FL_public', 5), pr(s, 'FL_public', 9)
        b.append(dict(sample=s, group=A[s]['region'], n_fl_5=a5['n_test'], n_mg_5=a5['n_reference'], rec_5=repr(int(a5['overlap']) / int(a5['n_reference'])),
                      n_fl_9=a9['n_test'], n_mg_9=a9['n_reference'], rec_9=repr(int(a9['overlap']) / int(a9['n_reference']))))
    append('SFig09/a08_benchmark282_old.tsv', b, 282)
    append('SFig09/a08_benchmark282_new.tsv', b, 282)
    append('SData1/v28d_capscan_mg_pairs.tsv', [dict(acquisition_id=s, specimen_id_A16='NA', specimen_id_biology=f"main:{A[s]['workbook_Main']}",
                                                     specimen_id_crossomics=A[s]['workbook_Main'], participant_id=15, region=A[s]['region'],
                                                     n_test=pr(s, 'FL_public', 5)['n_test'], n_reference=pr(s, 'FL_public', 5)['n_reference']) for s in NEW], 282)


def crossfit():
    ro = {r['sample']: r for r in csv.DictReader(open(M / 'disjoint/roster_mg295.tsv'), delimiter='\t')}
    aud = {}
    for r in csv.DictReader(open(FINAL / 'public347/spectrum_split347/audit/FULL_SOURCE_DATA.tsv'), delimiter='\t'):
        aud[(r['sample'], r['arm'])] = r
    sup_rows, fig_rows = [], []
    sup_hdr = (SD / 'supplementary_methods_6_7/crossfit_per_acquisition.tsv').read_text().split('\n', 1)[0].split('\t')
    for s in NEW:
        j = json.loads((M / f'crossfit/outputs/CAPSCAN/{s}/COMPLETE.json').read_text())
        oj = json.loads((M / f'disjoint/outputs/CAPSCAN/{s}/COMPLETE_odd.json').read_text()); assert oj['odd_mzml_identical']
        r = ro[s]; i = int(r['i'])
        d = dict(i=i, cohort='CAPSCAN', sample=s, group=r['group'], subject_id='15',
                 fold_ms2_total=sum(f['fold_ms2'] for f in j['folds'].values()),
                 fl_fasta_sequences_median=float(np.median([f['fl_fasta_sequences'] for f in j['folds'].values()])),
                 lda_fallback=';'.join(a for a in ['FL', 'FULL', 'MG'] if j[f'{a}_pooled']['rescoring'] != 'lda'))
        for L in ['5', '9']:
            k = f'peptides_{L}_50'
            for arm in ['FL', 'FULL', 'MG']:
                for how in ['pooled', 'union']:
                    d[f'{arm}_{how}_{L}'] = j[f'{arm}_{how}'][k]
                for x in (0, 1):
                    d[f'{arm}_pert{x + 1}_{L}'] = j[f'{arm}_pooled']['solver_perturbation'][x][k]
                d[f'{arm}_pooled_naive_{L}'] = j[f'{arm}_pooled_naive_diagnostic'][k]
            p = pr(s, 'FL_public', int(L))
            d[f'frozen_fl_{L}'] = int(p['n_test']); d[f'frozen_mg_{L}'] = int(p['n_reference'])
            g = lambda a, b: 100 * (a / b - 1) if b > 0 else float('nan')
            d[f'gain_cf_pooled_{L}'] = g(d[f'FL_pooled_{L}'], d[f'MG_pooled_{L}'])
            d[f'gain_cf_union_{L}'] = g(d[f'FL_union_{L}'], d[f'MG_union_{L}'])
            d[f'gain_cf_pert1_{L}'] = g(d[f'FL_pert1_{L}'], d[f'MG_pert1_{L}'])
            d[f'gain_cf_pert2_{L}'] = g(d[f'FL_pert2_{L}'], d[f'MG_pert2_{L}'])
            d[f'gain_full_pooled_{L}'] = g(d[f'FULL_pooled_{L}'], d[f'MG_pooled_{L}'])
            d[f'gain_full_union_{L}'] = g(d[f'FULL_union_{L}'], d[f'MG_union_{L}'])
            d[f'gain_frozen_{L}'] = g(d[f'frozen_fl_{L}'], d[f'frozen_mg_{L}'])
            d[f'cf_vs_full_pooled_{L}'] = g(d[f'FL_pooled_{L}'], d[f'FULL_pooled_{L}'])
            d[f'cf_vs_full_union_{L}'] = g(d[f'FL_union_{L}'], d[f'FULL_union_{L}'])
            d[f'calib_full_minus_frozen_{L}'] = d[f'gain_full_pooled_{L}'] - d[f'gain_frozen_{L}']
            d[f'mg_pooled_vs_frozen_{L}'] = g(d[f'MG_pooled_{L}'], d[f'frozen_mg_{L}'])
            mg_odd = oj['MG_ODD'][k]
            for arm, c in [('FULL', 'gain_odd_full'), ('HALF_SELF', 'gain_odd_self'), ('HALF_DISJOINT', 'gain_odd_disjoint')]:
                d[f'{c}_{L}'] = 100 * (int(aud[(s, arm)][k]) / mg_odd - 1)
        d['cluster_order'] = 6   # participant 15's existing cluster_order in this table
        assert set(d) == set(sup_hdr), set(d) ^ set(sup_hdr)
        sup_rows.append({k: (repr(v) if isinstance(v, float) else v) for k, v in d.items()})
        S = {a: set((M / f'crossfit/outputs/CAPSCAN/{s}/{a}_POOLED_ACCEPTED.txt').read_text().split()) for a in ['FL', 'MG']}
        assert len(S['FL']) == d['FL_pooled_5'] and len(S['MG']) == d['MG_pooled_5']
        rec = {}
        for L in (5, 9):
            T = {a: {x for x in v if L <= len(x) <= 50} for a, v in S.items()}
            rec[L] = 100 * len(T['FL'] & T['MG']) / len(T['MG'])
        fig_rows.append(dict(crossfit_index=i, cohort='CAPSCAN', acquisition_id=s, participant_id=15, participant_sort_rank=6, sampling_site=r['group'],
                             lda_fallback=d['lda_fallback'], FL_pooled_5=d['FL_pooled_5'], MG_pooled_5=d['MG_pooled_5'], frozen_fl_5=d['frozen_fl_5'],
                             frozen_mg_5=d['frozen_mg_5'], FL_pooled_9=d['FL_pooled_9'], MG_pooled_9=d['MG_pooled_9'], frozen_fl_9=d['frozen_fl_9'],
                             frozen_mg_9=d['frozen_mg_9'], recovery_cf_5_pct=repr(rec[5]), recovery_cf_9_pct=repr(rec[9])))
    append('supplementary_methods_6_7/crossfit_per_acquisition.tsv', sup_rows, 282)
    append('Fig4/crossfit_per_acquisition.tsv', fig_rows, 282)
    fb = [x for x in sup_rows if x['lda_fallback']]
    print('LDA fallbacks among the 13:', [(x['sample'], x['lda_fallback']) for x in fb])
    json.dump(dict(lda_fallbacks_13=[(x['sample'], x['lda_fallback']) for x in fb]), open(M / 'phase2/CROSSFIT13_NOTES.json', 'w'))


if __name__ == '__main__':
    steps = sys.argv[1:] or ['same', 'dna']
    if 'same' in steps: same_spectrum()
    if 'dna' in steps: dna_tables()
    if 'crossfit' in steps: crossfit()
