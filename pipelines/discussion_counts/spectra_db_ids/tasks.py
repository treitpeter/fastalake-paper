"""Build the task list: one row per acquisition, with every arm's existing Sage output, searched FASTA and mzML.
No search is run. Sources (all existing):
  PREDICT 296 / CapScan 282 paired FL_public + MG: sequence_accumulation_2026-09-29/manifest.json (same searches as A16 / Fig 4 / Supp Fig 6d,e)
  CapScan FL_public unmatched: production A08 results (49 of 331) and public347/A08_347 results (16 participant-15 acquisitions)
  FL own (cohort-MG FastaLake, noclust) and FL study-union FASTA, one search per acquisition: common-FASTA manifest
      (MM_review_analysis/end_of_day_2026-09-15/common_cache/comparison/inputs/manifest.json; recounted in len950/SFig07, array 6072755)
  HSTOOL 58: A21 manifest_revision2 arms 'original' (FL own) and 'union_individual' (FL study union, separate searches)
Searched FASTA and mzML are read from each search's own results.json (database.fasta, mzml_paths).
usage: python tasks.py   -> tasks.json"""
import csv, json, os
from pathlib import Path

OUT = Path(__file__).resolve().parent
LAB = Path(os.environ['FASTALAKE_DATA_DIR'])  # LAB_ROOT
FINAL = OUT.parent.parent
ACC = LAB / 'review/v15_integration/analyses/sequence_accumulation_2026-09-29/manifest.json'
A08 = LAB / 'MM_review_analysis/A08_4_5_capscan_public_only/results'
A08N = FINAL / 'public347/A08_347/results'
ROSTER = FINAL.parent / 'coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/roster/fixed_roster_347.tsv'
COMMON = LAB / 'MM_review_analysis/end_of_day_2026-09-15/common_cache/comparison/inputs/manifest.json'
A21 = LAB / 'MM_review_analysis/A21_7_hstool_quant_comparison/inputs/manifest_revision2.json'


def arm(results, fasta_expected=None, expected=None):
    results = Path(results); cfg = json.loads((results.parent / 'results.json').read_text())
    fa = cfg['database']['fasta']; mz = cfg['mzml_paths']; assert len(mz) == 1, results
    assert results.exists() and Path(fa).exists() and Path(mz[0]).exists(), results
    if fasta_expected is not None:
        assert os.path.realpath(fa) == os.path.realpath(fasta_expected), (fa, fasta_expected)
    return {'results': str(results), 'fasta': fa, 'mzml': mz[0], 'mzml_real': os.path.realpath(mz[0]), 'expected': expected}


rows = {}
acc = json.loads(ACC.read_text())['rows']
for r in acc:
    k = (r['cohort'], r['sample']); a = r['arms']
    rows[k] = {'cohort': r['cohort'], 'acquisition': r['sample'], 'arms': {
        'FLpub': arm(a['FL_public']['results'], a['FL_public']['fasta'], {L: a['FL_public']['expected_sets'][L]['path'] for L in ('5', '9')}),
        'MG': arm(a['MG']['results'], a['MG']['fasta'], {L: a['MG']['expected_sets'][L]['path'] for L in ('5', '9')})}}
roster = [x['sample_id'] for x in csv.DictReader(open(ROSTER), delimiter='\t')]
assert len(roster) == 347 == len(set(roster))
for s in roster:
    if ('CAPSCAN', s) in rows:
        continue
    d = A08 / s if (A08 / s / 'COMPLETE.json').exists() else A08N / s
    src = 'A08_production' if d.parent == A08 else 'A08_347'
    rows[('CAPSCAN', s)] = {'cohort': 'CAPSCAN', 'acquisition': s, 'flpub_source': src, 'arms': {
        'FLpub': arm(d / 'search/results.sage.tsv', None, {'5': str(d / 'peptide_q01.txt'), '9': str(d / 'length9_50.txt')})}}
for k, r in rows.items():
    if k[0] == 'CAPSCAN':
        r.setdefault('flpub_source', 'A08_production')
# own (cohort-MG FL) and study-union FASTA searches, one per acquisition
for t in json.loads(COMMON.read_text())['tasks']:
    k = (t['cohort'], t['sample'])
    if k not in rows:
        continue  # PREDICT: the 10 acquisitions without a matched metagenome are outside the 296 roster
    rc = json.loads(Path(t['common_receipt']).read_text()); assert rc['status'] == 'PASS'
    by = {Path(x['path']).name: x['path'] for x in rc['outputs']}
    rows[k]['arms']['FLcmg'] = arm(t['baseline_psm']['path'])
    rows[k]['arms']['FLunion'] = arm(by['results.sage.tsv'], rc['union']['path'])
# HSTOOL
m = json.loads(A21.read_text()); recs = {a['name']: {x['id']: x for x in a['records']} for a in m['arms']}
assert len(m['samples']) == 58
for s in m['samples']:
    o, u = recs['original'][s], recs['union_individual'][s]
    rows[('HSTOOL', s)] = {'cohort': 'HSTOOL', 'acquisition': s, 'arms': {
        'FLown': arm(o['psm']['path'], o['fasta']['path']), 'FLunion': arm(u['psm']['path'], u['fasta']['path'])}}
rows = sorted(rows.values(), key=lambda r: (r['cohort'], r['acquisition']))
# one mzML per acquisition: every arm must point to the same file
for r in rows:
    reals = {a['mzml_real'] for a in r['arms'].values()}
    r['mzml_same_file_all_arms'] = len(reals) == 1
    r['mzml'] = sorted(reals)[0]
from collections import Counter
print('rows', len(rows), Counter(r['cohort'] for r in rows))
print('arm sets', Counter((r['cohort'], tuple(sorted(r['arms']))) for r in rows))
print('acquisitions whose arms point to different mzML files', sum(not r['mzml_same_file_all_arms'] for r in rows))
(OUT / 'tasks.json').write_text(json.dumps(rows, indent=1) + '\n')
