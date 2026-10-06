"""Collect worker parts -> per_acquisition.tsv, summary tables, checks. usage: python collect.py"""
import glob, json
from pathlib import Path
import numpy as np, pandas as pd

O = Path(__file__).resolve().parent
SF7 = O.parent.parent / 'len950/SFig07/work'
rows = [r for f in sorted(glob.glob(str(O / 'parts/part_*.json'))) for r in json.load(open(f))]
tasks = json.loads((O / 'tasks.json').read_text())
assert len(rows) == len(tasks) == 701, len(rows)
d = pd.DataFrame(rows).sort_values(['cohort', 'acquisition']).reset_index(drop=True)
assert d.ms_level_total_equals_spectrumList.all() and d.mzml_same_file_all_arms.all()
ARMS = ['FLpub', 'MG', 'FLcmg', 'FLown', 'FLunion']
PART = {'FLpub': 'MG', 'MG': 'FLpub', 'FLcmg': 'FLunion', 'FLown': 'FLunion'}
for a in ARMS:
    assert (d.filter(like=f'{a}_peptide_protein_pairs_unplaced').fillna(0) == 0).all().all()
    d[f'{a}_ids_per_1000_ms2_9_50'] = 1000 * d[f'{a}_peptides_9_50'] / d['ms2_spectra']
    d[f'{a}_ids_per_1000_ms2_5_50'] = 1000 * d[f'{a}_peptides_5_50'] / d['ms2_spectra']
# FLunion's paired-candidates column depends on its partner (FLcmg in PREDICT/CapScan, FLown in HSTOOL)
d['FLunion_median_scored_candidates_paired'] = d['FLunion_median_scored_candidates_paired_vs_FLcmg'].fillna(d['FLunion_median_scored_candidates_paired_vs_FLown'])
for a, b in PART.items():
    d[f'{a}_median_scored_candidates_paired'] = d[f'{a}_median_scored_candidates_paired_vs_{b}']
# checks against existing recounts
chk = {}
c = pd.concat([pd.read_csv(f, sep='\t') for f in sorted(glob.glob(str(SF7 / 'clin/part_*.tsv')))]).set_index('sample')
m = d[d.FLcmg_peptides_5_50.notna()].set_index('acquisition')
chk['clinical_own_common_rows_compared'] = len(m)
chk['clinical_own_common_counts_match_len950_SFig07'] = bool(all((c.loc[m.index, f'original_{L}'] == m[f'FLcmg_peptides_{L}_50']).all() and (c.loc[m.index, f'common_{L}'] == m[f'FLunion_peptides_{L}_50']).all() for L in (5, 9)))
h = pd.read_csv(SF7 / 'hstool/identity_per_acquisition_len.tsv', sep='\t'); h = h[h.arm == 'union_individual']
m = d[d.cohort == 'HSTOOL'].set_index('acquisition')
chk['hstool_counts_match_len950_SFig07'] = bool(all((h[h.min_length == L].set_index('sample').loc[m.index, 'original'] == m[f'FLown_peptides_{L}_50']).all() and (h[h.min_length == L].set_index('sample').loc[m.index, 'test'] == m[f'FLunion_peptides_{L}_50']).all() for L in (5, 9)))
e = pd.read_csv(O.parent / 'E1_per_acquisition.tsv', sep='\t').set_index('sample'); m = d[d.MG_peptides_5_50.notna()].set_index('acquisition')
chk['scored_candidates_match_task_E1'] = bool((e.loc[m.index, 'median_a'] == m.FLpub_median_scored_candidates_paired).all() and (e.loc[m.index, 'median_b'] == m.MG_median_scored_candidates_paired).all() and (e.loc[m.index, 'paired_spectra'] == m.paired_spectra_FLpub_MG).all())
chk['flpub_mg_sets_reproduce_A10_A08'] = bool(d.FLpub_reproduces_existing_sets.all() and d.loc[d.MG_peptides_5_50.notna(), 'MG_reproduces_existing_sets'].all())
chk['psm_rows_read_total'] = int(d.filter(like='_psm_rows_read').sum().sum())
chk['ms2_spectra_total'] = int(d.ms2_spectra.sum())
print(json.dumps(chk, indent=1)); assert all(v for k, v in chk.items() if isinstance(v, bool))

MET = ['db_proteins', 'db_residues', 'median_scored_candidates_paired', 'median_scored_candidates_all_rank1', 'peptides_5_50', 'peptides_9_50',
       'proteins_with_9_50_peptide', 'peptides_per_protein_mean', 'peptides_per_protein_median', 'coverage_median', 'ids_per_1000_ms2_9_50', 'ids_per_1000_ms2_5_50']
cols = ['cohort', 'acquisition', 'flpub_source', 'ms2_spectra'] + [f'{a}_{k}' for a in ARMS for k in MET] + \
       [f'paired_spectra_{a}_{b}' for a, b in [('FLpub', 'MG'), ('FLcmg', 'FLunion'), ('FLown', 'FLunion')]] + \
       [f'{a}_psm_rows_read' for a in ARMS] + ['mzml'] + [f'{a}_{k}' for a in ARMS for k in ('results', 'fasta')]
d[cols].to_csv(O / 'per_acquisition.tsv', sep='\t', index=False, float_format='%.6g')

# cohort medians per arm (strata: PREDICT 296; CapScan 282 matched; CapScan FL_public all 347; CapScan own/common 331; HSTOOL 58)
S = []
strata = [('PREDICT', d[d.cohort == 'PREDICT'], ['FLpub', 'MG', 'FLcmg', 'FLunion']),
          ('CAPSCAN_matched282', d[(d.cohort == 'CAPSCAN') & d.MG_peptides_5_50.notna()], ['FLpub', 'MG', 'FLcmg', 'FLunion']),
          ('CAPSCAN_all347', d[d.cohort == 'CAPSCAN'], ['FLpub']),
          ('CAPSCAN_own_common331', d[(d.cohort == 'CAPSCAN') & d.FLcmg_peptides_5_50.notna()], ['FLcmg', 'FLunion']),
          ('HSTOOL', d[d.cohort == 'HSTOOL'], ['FLown', 'FLunion'])]
for name, x, arms in strata:
    for a in arms:
        r = {'stratum': name, 'arm': a, 'acquisitions': len(x), 'ms2_spectra': x.ms2_spectra.median()}
        r.update({k: x[f'{a}_{k}'].median() for k in MET}); S.append(r)
S = pd.DataFrame(S); S.to_csv(O / 'summary_medians.tsv', sep='\t', index=False, float_format='%.6g')
# paired ratios with sign counts
P = []
for name, x, a, b in [('PREDICT', strata[0][1], 'FLpub', 'MG'), ('CAPSCAN_matched282', strata[1][1], 'FLpub', 'MG'),
                      ('PREDICT', strata[0][1], 'FLunion', 'FLcmg'), ('CAPSCAN_own_common331', strata[3][1], 'FLunion', 'FLcmg'),
                      ('HSTOOL', strata[4][1], 'FLunion', 'FLown')]:
    for k in MET:
        u, v = x[f'{a}_{k}'], x[f'{b}_{k}']; ok = u.notna() & v.notna() & (v != 0)
        P.append({'stratum': name, 'ratio': f'{a}/{b}', 'metric': k, 'n': int(ok.sum()), 'median_ratio': float(np.median(u[ok] / v[ok])),
                  'n_higher': int((u[ok] > v[ok]).sum()), 'n_lower': int((u[ok] < v[ok]).sum()), 'n_equal': int((u[ok] == v[ok]).sum())})
P = pd.DataFrame(P); P.to_csv(O / 'summary_paired_ratios.tsv', sep='\t', index=False, float_format='%.6g')
(O / 'CHECKS.json').write_text(json.dumps(chk, indent=1) + '\n')
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(S.to_string()); print(P.to_string())
