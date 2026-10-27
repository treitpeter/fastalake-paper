#!/usr/bin/env python3
"""282 control: every output of the relocated chain must equal the archived production output.
Files that embed absolute paths/timings are compared after removing only those fields (listed in the report).
Usage: compare_control.py <phase3 item2 dir>"""
import gzip, hashlib, json, sys
from pathlib import Path
I = Path(sys.argv[1]); C = I / 'runs/ctl282'
FLR = Path('LAB_ROOT')
OG = FLR / 'simplified_grouping_results/CAPSCAN_two_arm_282_2026-09-11_immediate'
OE = FLR / 'reviews/supplementary_completion_2026-09-13/assignment_exports'
OH = FLR / 'Hartout_edits/H09_evidence_ablation/full_CAPSCAN'
REPO = I.parents[1] / 'phase2/repo_base/source_data/SFig08/clinical_fixed_dictionary_ablation.tsv'
def h(p):
    x = hashlib.sha256(); op = gzip.open if str(p).endswith('.gz') else open
    with op(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''): x.update(b)
    return x.hexdigest()
res = []
def cmp(a, b, label, mode='content'):
    if mode == 'content': ok = h(a) == h(b)
    else: ok = mode(a, b)
    res.append(dict(item=label, new=str(a), archived=str(b), identical=bool(ok), mode='decompressed sha256' if mode == 'content' else 'json minus path/timing fields'))
def jdrop(keys):
    def f(a, b):
        A = json.loads(Path(a).read_text()); B = json.loads(Path(b).read_text())
        for k in keys: A.pop(k, None); B.pop(k, None)
        return A == B
    return f
for arm in ['FL', 'metaG']:
    for f in sorted((OG / 'cohorts' / arm / 'direct').iterdir()):
        if f.name == 'input_manifest.json':
            A = json.loads((C / 'cohorts' / arm / 'direct' / f.name).read_text()); B = json.loads(f.read_text())
            ok = [x['sha256'] for x in A[:-1]] == [x['sha256'] for x in B[:-1]] and A[-1]['sha256'] == B[-1]['sha256']
            res.append(dict(item=f'grouping {arm} input_manifest sha256 list', identical=ok, mode='sha256 sequence (paths differ by relocation)'))
        else:
            cmp(C / 'cohorts' / arm / 'direct' / f.name, f, f'grouping {arm} {f.name}')
cmp(C / 'comparison/SUMMARY.json', OG / 'comparison/SUMMARY.json', 'comparison SUMMARY.json')
cmp(C / 'comparison/matched_protein_group_counts.tsv', OG / 'comparison/matched_protein_group_counts.tsv', 'comparison matched_protein_group_counts.tsv')
PATHY = {'TASK.json', 'INPUT_RECEIPTS.json', 'SEARCH_INPUTS.json', 'original_input_manifest.json', 'PROGRESS.json'}
for t in ['CAPSCAN_282_FL', 'CAPSCAN_282_metaG']:
    for f in sorted((OE / t).iterdir()):
        if f.name in PATHY: continue
        if f.name == 'COMPLETE.json': cmp(C / 'export' / t / f.name, f, f'export {t} COMPLETE.json', jdrop(['elapsed_seconds']))
        else: cmp(C / 'export' / t / f.name, f, f'export {t} {f.name}')
    A = json.loads((C / 'export' / t / 'INPUT_RECEIPTS.json').read_text()); B = json.loads((OE / t / 'INPUT_RECEIPTS.json').read_text())
    res.append(dict(item=f'export {t} INPUT_RECEIPTS sha256 list', identical=[x['sha256'] for x in A] == [x['sha256'] for x in B], mode='sha256 sequence'))
for f in ['acquisition_summary.tsv', 'group_ablation.tsv.gz', 'KO_ablation_equal_split.tsv.gz']:
    cmp(C / 'h09/full_CAPSCAN' / f, OH / f, f'H09 {f}')
cmp(C / 'values/clinical_fixed_dictionary_ablation_282.tsv', REPO, 'SFig08 source table clinical_fixed_dictionary_ablation.tsv (repo_base)')
ok = all(r['identical'] for r in res)
(I / 'controls/CONTROL_282.json').write_text(json.dumps(dict(status='PASS' if ok else 'FAIL', n_checks=len(res), checks=res), indent=1) + '\n')
print('CONTROL', 'PASS' if ok else 'FAIL', len(res), 'checks'); [print(r['identical'], r['item']) for r in res]
