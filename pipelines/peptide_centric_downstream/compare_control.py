"""Baseline control: the copied chain on the frozen 331 inputs must equal the A33 / v25 outputs exactly.
usage: compare_control.py <control_run_dir> [<control_gated_dir>]"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
A33 = Path('LAB_ROOT/MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19')
SD = Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/figure_sd/source_data/package/FastaLake_Source_Data/Figure5')
C = Path(sys.argv[1]); report = {}
def tsv_equal(a, b, sort=None):
    x, y = pd.read_csv(a, sep='\t'), pd.read_csv(b, sep='\t')
    if sort: x, y = x.sort_values(sort).reset_index(drop=True), y.sort_values(sort).reset_index(drop=True)
    assert list(x.columns) == list(y.columns), (a, list(x.columns), list(y.columns))
    pd.testing.assert_frame_equal(x, y, check_exact=True)
    return 'IDENTICAL'
for rule in ['simplified', 'original']:
    a, b = C/f'CAPSCAN__{rule}/source_data/Figure5', A33/f'results/CAPSCAN__{rule}/source_data/Figure5'
    za, zb = np.load(a/'function_species.npz'), np.load(b/'function_species.npz')
    for k in ['ko', 'species', 'samples']: assert np.array_equal(za[k], zb[k]), (rule, k)
    assert np.array_equal(za['values'], zb['values']), (rule, 'values')
    report[f'{rule}/function_species.npz'] = 'IDENTICAL'
    for f, s in [('contributor_cells.tsv', ['KO', 'participant']), ('contributor_region_cells.tsv', ['KO', 'participant', 'region']),
                 ('contributor_agreement_cells.tsv', ['KO', 'participant', 'region']), ('contributor_agreement_donors.tsv', ['participant', 'region']),
                 ('contributor_agreement_well_resolved.tsv', ['participant', 'region']), ('contributor_mosaic.tsv', ['KO', 'participant']),
                 ('mosaic_donor_regions.tsv', None), ('mosaic_donors.tsv', None), ('mosaic_abundance.tsv', None)]:
        report[f'{rule}/{f}'] = tsv_equal(a/f, b/f, s)
    sa = json.loads((a/'SUMMARY.json').read_text())['contributors']; sb = json.loads((b/'SUMMARY.json').read_text())['contributors']
    assert sa == sb, rule; report[f'{rule}/SUMMARY contributors'] = 'IDENTICAL'
report['matched_capsule_stool_donors.tsv'] = tsv_equal(C/'results/matched_capsule_stool_donors.tsv', A33/'results/matched_capsule_stool_donors.tsv')
for f in ['participant_leaders_resolved.tsv', 'participant_leaders_strict_both_methods.tsv']:
    report[f] = tsv_equal(C/'render/tables'/f, A33/'render/tables'/f)
    report[f + ' vs v25 Source Data'] = tsv_equal(C/'render/tables'/f, SD/f)
report['mosaic_donor_regions vs v25 Source Data'] = tsv_equal(C/'CAPSCAN__simplified/source_data/Figure5/mosaic_donor_regions.tsv', SD/'mosaic_donor_regions.tsv')
n = json.loads((C/'NUMBERS.json').read_text())
exp = dict(res_o=0.410498, res_s=0.44, str_o=0.567, str_s=0.58208)
got = dict(res_o=n['fig5f_leaders']['resolved']['original'], res_s=n['fig5f_leaders']['resolved']['simplified'],
           str_o=n['fig5f_leaders']['strict_both_methods']['original'], str_s=n['fig5f_leaders']['strict_both_methods']['simplified'])
for k in exp: assert round(got[k]*100, 1) == round(exp[k]*100, 1), (k, got[k])
assert n['fig5e_mosaic_direct']['same_leader'] == 657 and n['fig5e_mosaic_direct']['comparisons'] == 1344
assert n['fig5e_mosaic_direct']['persistence'] == 98 and n['fig5e_mosaic_direct']['persistence_profiles'] == 336
assert n['sf12_direct']['eligible_cells'] == 20352 and n['sf12_direct']['eligible_functions'] == 1416 and round(100*n['sf12_direct']['median_resolved_fraction'], 1) == 29.4
report['v25 printed values'] = 'REPRODUCED (41.0/44.0, 56.7/58.2, 657/1344, 98/336, 20,352 cells, 1,416 functions, 29.4%)'
if len(sys.argv) > 2:
    g = Path(sys.argv[2])
    for stem in ['KO_species', 'KO_species_original']:
        za, zb = np.load(g/f'CAPSCAN_{stem}.npz'), np.load(A33/f'inputs/gated_joint/CAPSCAN_{stem}.npz')
        for k in ['ko', 'species', 'samples', 'values']: assert np.array_equal(za[k], zb[k]), (stem, k)
        report[f'gated_joint {stem}.npz'] = 'IDENTICAL'
(C/'CONTROL.json').write_text(json.dumps(dict(status='PASS', checks=report), indent=2) + '\n')
print(json.dumps(report, indent=1)); print('CONTROL PASS')
