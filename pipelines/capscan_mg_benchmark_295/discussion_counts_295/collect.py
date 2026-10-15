"""Collect: control (282 rerun) against the archived outputs, then the 295 values. Read-only on inputs."""
import csv, filecmp, json
from pathlib import Path

D = Path(__file__).resolve().parent
OLD = D.parent / 'discussion_counts'
ACC = Path('LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29')
SUB = 'review/v15_integration/analyses/sequence_accumulation_2026-09-29'
T = {t: D / 'trees' / t / SUB for t in ['c282', 'r295']}
tsv = lambda p: list(csv.DictReader(open(p), delimiter='\t'))


def vals(out):
    A = json.loads((out / 'A_RESULT.json').read_text()); E = {r['cohort']: r for r in tsv(out / 'E1_summary.tsv')}
    c = {(r['cohort'], r['arm']): r for r in A['curve']}
    v = {}
    for co in ['PREDICT', 'CAPSCAN']:
        v[co] = {'specimens': c[(co, 'MG')]['specimens'], 'MG_union': c[(co, 'MG')]['fixed_kfinal'], 'FL_union': c[(co, 'FL_public')]['fixed_kfinal'],
                 'fold': c[(co, 'MG')]['fixed_kfinal'] / c[(co, 'FL_public')]['fixed_kfinal'], 'acquisitions_E1': int(E[co]['acquisitions']),
                 'E1_FL_9_50': float(E[co]['cohort_median_of_median_a_9_50']), 'E1_MG_9_50': float(E[co]['cohort_median_of_median_b_9_50'])}
    return v


res = {'checks': {}}
ck = res['checks']
# control: rerun 282 trees reproduce archived worker outputs byte for byte (except COMPLETE.json job id)
for f in ['accumulation.tsv', 'per_specimen.tsv', 'per_acquisition.tsv', 'header_audit.tsv']:
    ck[f'c282 full/CAPSCAN/{f} byte-identical to archive'] = filecmp.cmp(T['c282'] / 'full/CAPSCAN' / f, ACC / 'full/CAPSCAN' / f, shallow=False)
for f in ['per_acquisition.tsv', 'cohort_median_acquisition.tsv', 'cohort_median_participant.tsv']:
    ck[f'c282 candidates/full/{f} byte-identical to archive'] = filecmp.cmp(T['c282'] / 'candidates/full' / f, ACC / 'candidates/full' / f, shallow=False)
for f in ['A_accumulation_readout.tsv', 'A_per_specimen_medians.tsv', 'A_per_specimen_ids.tsv', 'E1_per_acquisition.tsv', 'E1_summary.tsv']:
    ck[f'out_c282/{f} byte-identical to discussion_counts/{f}'] = filecmp.cmp(D / 'out_c282' / f, OLD / f, shallow=False)
c282 = vals(D / 'out_c282'); old = vals(OLD); r295 = vals(D / 'out_r295')
ck['c282 values == archived values'] = c282 == old
ck['c282 CAPSCAN MG union == 28,547,993'] = c282['CAPSCAN']['MG_union'] == 28547993
ck['c282 CAPSCAN FL union == 416,703'] = c282['CAPSCAN']['FL_union'] == 416703
ck['c282 CAPSCAN E1 MG 9-50 == 3322.75'] = c282['CAPSCAN']['E1_MG_9_50'] == 3322.75
ck['c282 CAPSCAN E1 FL 9-50 == 277.5'] = c282['CAPSCAN']['E1_FL_9_50'] == 277.5
ck['r295 PREDICT values == archived'] = r295['PREDICT'] == old['PREDICT']
ck['r295 CAPSCAN specimens == 263, acquisitions == 295'] = (r295['CAPSCAN']['specimens'], r295['CAPSCAN']['acquisitions_E1']) == (263, 295)
# 295: every 282 / PREDICT per-acquisition row unchanged; exactly 13 new
e282 = {r['sample']: r for r in tsv(D / 'out_c282/E1_per_acquisition.tsv')}; e295 = {r['sample']: r for r in tsv(D / 'out_r295/E1_per_acquisition.tsv')}
ck['r295 E1 rows: 578 old rows identical, 13 new'] = all(e295[s] == e282[s] for s in e282) and len(e295) == 591
p282 = {(r['sample'], r['arm']): r for r in tsv(T['c282'] / 'full/CAPSCAN/per_acquisition.tsv')}
p295 = {(r['sample'], r['arm']): r for r in tsv(T['r295'] / 'full/CAPSCAN/per_acquisition.tsv')}
ck['r295 accumulation per_acquisition: 564 old rows identical, 26 new'] = all(p295[k] == p282[k] for k in p282) and len(p295) == 590
cp = json.loads((T['r295'] / 'full/CAPSCAN/COMPLETE.json').read_text())
ck['r295 worker COMPLETE PASS 295/263/15'] = (cp['status'], cp['acquisitions'], cp['specimens'], cp['participants']) == ('PASS', 295, 263, 15)
ck['r295 MG union > 282 MG union and FL union > 282 FL union'] = r295['CAPSCAN']['MG_union'] > c282['CAPSCAN']['MG_union'] and r295['CAPSCAN']['FL_union'] > c282['CAPSCAN']['FL_union']
res.update({'archived': old, 'control_282': c282, 'roster_295': r295, 'all_checks_pass': all(ck.values())})
(D / 'COLLECT.json').write_text(json.dumps(res, indent=1) + '\n'); print(json.dumps(res, indent=1))
