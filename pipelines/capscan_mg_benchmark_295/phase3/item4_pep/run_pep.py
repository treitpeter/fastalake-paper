#!/usr/bin/env python3
"""Item 4 (cross-fitted PEP pooling) for the 13 participant-15 pairs, with the C1 producer imported UNCHANGED
(c1_copy/cf_pep.py, byte-identical to coauthor_comment_check_v25_2026-10-05/compute/C1_pep_primary/crossfit/cf_pep.py).
The only thing set from outside is cf_pep.IT4 (the cross-fit output root read by run_one); design code (sage_lda, sage_q)
is the iteration4 copy that cf_pep imports, byte-identical to mg295/crossfit/design.

Modes
  control <i> [<i> ...]  282-roster rows, read from iteration4 (C1 default); result must equal C1 per_acq/<i>.json exactly
                         (keys job/code_sha256/pep_model excluded). Writes control/<i>.json + control/<i>.cmp.json.
  pool1311               row 335 (1311) read from the mg295 rerun (public347 347-pooled caches); compared with C1 335.json
                         (diagnostic for the fold-cache pool, item 6). Writes control/pool1311.json.
  new <k>                k-th of the 13 (roster_mg295 rows 578+k), read from mg295/crossfit/outputs. Writes per_acq13/<i>.json.
"""
import sys, json, os, hashlib, csv
from pathlib import Path
H = Path(__file__).resolve().parent
sys.path.insert(0, str(H / 'c1_copy'))
import cf_pep as m
C1 = Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C1_pep_primary/crossfit')
MG295 = H.parents[1]
IT4_DEFAULT = m.IT4
SKIP = {'job', 'code_sha256', 'pep_model'}
assert m.VARIANT == json.loads((C1 / 'PEP_MODEL.json').read_text())

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def cmp(rec, i):
    ref = json.loads((C1 / 'per_acq' / f'{i:03d}.json').read_text())
    a = {k: v for k, v in rec.items() if k not in SKIP}; b = {k: v for k, v in ref.items() if k not in SKIP}
    a = json.loads(json.dumps(a)); diff = [k for k in sorted(set(a) | set(b)) if a.get(k) != b.get(k)]
    return dict(i=i, equal=not diff, differing_keys=diff, ref_job=ref.get('job'))

def stamp(rec):
    rec['job'] = os.getenv('SLURM_JOB_ID'); rec['code_sha256'] = sha(m.__file__); rec['wrapper_sha256'] = sha(__file__)
    rec['pep_model'] = m.VARIANT; rec['it4_root'] = str(m.IT4); return rec

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'control':
        rows = m.roster(); (H / 'control').mkdir(exist_ok=True)
        for i in map(int, sys.argv[2:]):
            m.IT4 = IT4_DEFAULT
            rec = m.run_one(i, rows[i]); c = cmp(rec, i)
            (H / 'control' / f'{i:03d}.json').write_text(json.dumps(stamp(rec), indent=1))
            (H / 'control' / f'{i:03d}.cmp.json').write_text(json.dumps(c, indent=1)); print('CONTROL', c, flush=True)
    elif mode == 'pool1311':
        rows = m.roster(); m.IT4 = MG295 / 'crossfit'
        rec = m.run_one(335, rows[335]); c = cmp(rec, 335)
        (H / 'control' / 'pool1311.json').write_text(json.dumps(dict(record=stamp(rec), cmp_vs_iteration4=c), indent=1)); print('POOL1311', c, flush=True)
    elif mode == 'new':
        k = int(sys.argv[2]); i = 578 + k
        r = list(csv.DictReader(open(MG295 / 'disjoint/roster_mg295.tsv'), delimiter='\t'))[i]
        assert int(r['i']) == i and r['cohort'] == 'CAPSCAN' and r['subject_id'] == '15', r
        m.IT4 = MG295 / 'crossfit'
        f = H / 'per_acq13' / f'{i:03d}.json'
        rec = m.run_one(i, dict(cohort=r['cohort'], sample=r['sample'], subject_id=r['subject_id'], group=r['group']))
        f.write_text(json.dumps(stamp(rec), indent=1)); print('done', i, r['sample'], {a: rec[a]['baseline_equals_frozen'] for a in m.ARMS}, flush=True)
