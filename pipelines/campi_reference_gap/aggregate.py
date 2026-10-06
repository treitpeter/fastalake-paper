#!/usr/bin/env python3
"""Validate and aggregate frozen per-acquisition CAMPI gap audits."""
import csv, json, os, sys
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
RUN = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(RUN / 'frozen'))
from lib import dump, meta, verify, table, sha

def rows(p):
    with open(p) as f: return list(csv.DictReader(f, delimiter='\t'))

m = json.loads((RUN / 'inputs/aggregate_manifest.json').read_text())
for x in m['inputs'] + m['code']: verify(x)
original = json.loads((RUN / 'inputs/manifest.json').read_text())
stage_rows, member_rows, pairs = [], [], []
pooled = defaultdict(set)
labels = {}
for task in original['tasks']:
    sample = task['sample']; out = RUN / 'results' / sample
    c = json.loads((out / 'COMPLETE.json').read_text())
    assert c['status'] == 'PASS' and c['sample'] == sample
    for x in c['outputs']: verify(x)
    gaps = rows(out / 'reference_gap.tsv')
    assert len(gaps) == c['reference_peptides']
    assert dict(Counter(r['stage'] for r in gaps)) == c['stages']
    stage_rows += rows(out / 'stage_counts.tsv')
    member_rows += rows(out / 'member_counts.tsv')
    sets = {}
    for r in rows(out / 'peptide_membership.tsv'):
        members = frozenset(filter(None, r['members'].split(';')))
        assert len(members) == int(r['n_members'])
        assert r['peptide'] not in labels or labels[r['peptide']] == members
        labels[r['peptide']] = members
        for arm in r['arms'].split(';'):
            sets.setdefault(arm, set()).add(r['peptide'])
            pooled[arm].add(r['peptide'])
    assert {a:len(v) for a,v in sets.items()} == c['counts']
    for lo in (5,9):
        ref = {p for p in sets['REF'] if len(p) >= lo}
        for arm, peps in sets.items():
            p = {s for s in peps if len(s) >= lo}
            pairs.append(dict(sample=sample, arm=arm, min_length=lo, peptides=len(p),
                reference=len(ref), overlap=len(p&ref), gained=len(p-ref), lost=len(ref-p),
                count_ratio=len(p)/len(ref), reference_recovery=len(p&ref)/len(ref)))
out = RUN / 'results'
table(out / 'stage_counts_all.tsv', stage_rows)
table(out / 'member_counts_all.tsv', member_rows)
table(out / 'paired_reference_recovery.tsv', pairs)
summary, members_summary = [], []
members = [r['member'] for r in rows(original['members']['path'])]
for lo in (5,9):
    ref = {p for p in pooled['REF'] if len(p) >= lo}
    for arm, allp in pooled.items():
        p = {s for s in allp if len(s) >= lo}
        pr = [r for r in pairs if r['arm']==arm and r['min_length']==lo]
        summary.append(dict(arm=arm, min_length=lo, acquisitions=len(pr),
            pooled_peptides=len(p), pooled_reference=len(ref), pooled_overlap=len(p&ref),
            pooled_reference_recovery=len(p&ref)/len(ref),
            median_acquisition_count_ratio=median(r['count_ratio'] for r in pr),
            median_acquisition_reference_recovery=median(r['reference_recovery'] for r in pr),
            known_member_specific=sum(len(labels[s])==1 for s in p),
            known_member_shared=sum(len(labels[s])>1 for s in p),
            absent_eight_member_reference=sum(not labels[s] for s in p)))
        for member in members:
            values=[int(r['specific_peptides']) for r in member_rows
                if r['arm']==arm and int(r['min_length'])==lo and r['member']==member]
            assert len(values)==15
            members_summary.append(dict(arm=arm,min_length=lo,member=member,
                pooled_specific_peptides=sum(labels[s]=={member} for s in p),
                acquisitions_ge1=sum(v>=1 for v in values), acquisitions_ge10=sum(v>=10 for v in values),
                median_specific_peptides=median(values)))
table(out / 'recovery_summary.tsv', summary)
table(out / 'member_summary.tsv', members_summary)
stage_summary=[]
for lo in (5,9):
    for stage in sorted({r['stage'] for r in stage_rows}):
        rr=[r for r in stage_rows if int(r['min_length'])==lo and r['stage']==stage]
        stage_summary.append(dict(min_length=lo,stage=stage,
            peptide_acquisition_observations=sum(int(r['peptides']) for r in rr),
            no_direct_prediction_observations=sum(int(r['no_direct_prediction']) for r in rr)))
table(out / 'stage_summary.tsv', stage_summary)
report=['# CAMPI five-arm recovery and residual gap','',
    'All 15 acquisition audits and their hashed outputs passed. Peptides are unmodified I/L-canonical targets at peptide q ≤ 0.01; 9–50 is a nested sensitivity population, not a new FDR calibration.', '',
    'The reference contains the eight named CAMPI proteomes plus 116 other entries. Member-specific means unique among these eight proteomes, not globally species-specific. Presence outside this reference does not by itself establish an error.', '',
    '| Arm | Pooled peptides | Pooled REF recovery | Median paired REF recovery |',
    '|---|---:|---:|---:|']
for r in summary:
    if r['min_length']==5: report.append(f"| {r['arm']} | {r['pooled_peptides']:,} | {r['pooled_reference_recovery']:.2%} | {r['median_acquisition_reference_recovery']:.2%} |")
report += ['', 'Gap stages describe the observed pipeline state, not a causal intervention. An absent direct prediction can coexist with recovery through another peptide in the same protein. Only rank-0 exported predictions were available; do not describe these as all beams.', '',
    'Stage totals count peptide–acquisition observations. Pooled identity recovery instead counts distinct peptide sequences across the cohort. The final-digest check uses the actual Sage KR/P digest with at most two missed cleavages and no initiator-methionine clipping.', '',
    'Scientific interpretation, figure review and DOCX integration remain open.']
(out / 'REPORT.md').write_text('\n'.join(report)+'\n')
files=['stage_counts_all.tsv','member_counts_all.tsv','paired_reference_recovery.tsv',
       'recovery_summary.tsv','member_summary.tsv','stage_summary.tsv','REPORT.md']
dump(out / 'AGGREGATE_VALIDATION.json',dict(status='PASS',acquisitions=15,arms=5,
    job_id=os.environ.get('SLURM_JOB_ID'), manifest_sha256=sha(RUN/'inputs/aggregate_manifest.json'),
    outputs=[meta(out/n) for n in files], scientific_review='PENDING'))
print('CAMPI aggregate PASS', flush=True)
