"""Build run directories for the Fig. 6 chain with the C2 fixed roster (331 frozen + 16 participant-15 capsules).

Per-acquisition stages (A26 revision2, A20, A27) get only the 16 new participant-15 tasks; their frozen code
is copied byte-identical. The collectors (A24 CAPSCAN, A28 revision4 full_CAPSCAN) get two run directories:
  baseline331 : the original 331 CapScan tasks unchanged (must reproduce the v25 outputs exactly), and
  fixed347    : the same 331 tasks plus the 16 participant-15 tasks (the only change).
The tryptic two-recounter check (A10 FL_cohortMG_5_50.txt, absent for participant 15 because A10 ran on the
331 roster) is fed by an independent recount of the same tryptic results.sage.tsv using the A05 frozen reader.
"""
import os
import json, shutil, sys, pathlib, hashlib, copy
import pandas as pd
F6 = pathlib.Path(sys.argv[1]).resolve()
MM = pathlib.Path(os.environ['FASTALAKE_DATA_DIR'], 'MM_review_analysis')
CAP = pathlib.Path(os.environ['FASTALAKE_DATA_DIR'], 'canonical/publication_data/arms/noclust_2026-08-07/cohorts/CAPSCAN')
RAW = pathlib.Path(os.environ['FASTALAKE_CAPSCAN_MZML_DIR'])
PRED = pathlib.Path(os.environ['FASTALAKE_DATA_DIR'], 'projects/capscan_unified/mzml_combined/alphanovo_inference')
sys.path.insert(0, str(F6/'A05/frozen'))
from lib import meta
from search import recount  # A05 frozen tryptic/unspecific recounter (independent of A20/A27 readers)
J = lambda p: json.loads(pathlib.Path(p).read_text())
def W(p, x): p = pathlib.Path(p); p.parent.mkdir(parents=True, exist_ok=True); p.write_text(json.dumps(x, indent=2) + '\n')
def sha(p): return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

p15 = pd.read_csv(F6.parent/'roster/participant15_added.tsv', sep='\t')
fixed = pd.read_csv(F6.parent/'roster/fixed_roster_347.tsv', sep='\t')
assert len(p15) == 16 and len(fixed) == 347

def copy_frozen(src_run, dst_run, entries):
    (dst_run/'frozen').mkdir(parents=True, exist_ok=False)
    new = []
    for x in entries:
        p = pathlib.Path(x['path']); assert sha(p) == x['sha256'], p
        q = dst_run/'frozen'/p.name; shutil.copy2(p, q); assert sha(q) == x['sha256']; new.append(meta(q))
    for extra in (src_run/'frozen').glob('*.py'):  # helper modules imported but not listed
        if not (dst_run/'frozen'/extra.name).exists(): shutil.copy2(extra, dst_run/'frozen'/extra.name)
    return new

# ---------- independent tryptic recount for participant 15 (replaces the A10 file in the two-recounter check)
base = {}
for s in p15.sample_id:
    _, sets, _ = recount(str(CAP/'5_sage/razor_hash-acc'/s/'results.sage.tsv'))
    f = F6/'baseline_tryptic_recount'/s/'FL_tryptic_5_50.txt'; f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text(''.join(p + '\n' for p in sorted(sets['peptide_q01']))); base[s] = meta(f)
# control: the same recounter reproduces the A10 file for three frozen acquisitions
for s in ['20210513_EXPL0_PeTr_SA_T760_1249', '20211017_EXPL0_PeTr_SA_T760_2126', '20210513_EXPL0_PeTr_SA_T760_1250']:
    _, sets, _ = recount(str(CAP/'5_sage/razor_hash-acc'/s/'results.sage.tsv'))
    a10 = set((MM/'A10_4_3_4_5_identity_recovery/results/CAPSCAN'/s/'FL_cohortMG_5_50.txt').read_text().splitlines())
    assert sets['peptide_q01'] == a10, ('recounter does not reproduce A10', s)

def design(r): return dict(participant=str(int(r.subject_id)), region=r.region, aliquot='pellet', specimen_id=f'main:{int(r.workbook_Main)}')
def route(sub, s):
    d = CAP/sub/'razor_hash-acc'/s
    return dict(search=meta(d/'results.sage.tsv'), config=meta(d/'results.json'))
A05R = lambda s: str(F6/'A05/results'/(s + '__unspecific')/'COMPLETE.json')

# ---------- A26 revision2 (de novo thresholds), 16 tasks
src = MM/'A26_4_7_denovo_thresholds'; old = J(src/'inputs/manifest_revision2.json'); run = F6/'A26'
code = copy_frozen(src, run, old['code'])
tasks = [dict(cohort='CAPSCAN', sample=r.sample_id, prediction=meta(PRED/(r.sample_id + '_predictions.csv')),
              fasta=meta(CAP/'4_parsimony/razor_hash-acc'/r.sample_id/'parsimony.fasta'), design=design(r)) for r in p15.itertuples()]
m = {k: v for k, v in old.items() if k not in ('tasks', 'code')}; m.update(tasks=tasks, code=code, note='C2: 16 participant-15 tasks; frozen code byte-identical')
W(run/'inputs/manifest_revision2.json', m)

# ---------- A20 terminal routes, 16 tasks
src = MM/'A20_4_7_terminal_routes'; old = J(src/'inputs/manifest.json'); run = F6/'A20'
code = copy_frozen(src, run, old['frozen_code'])
tasks = []
for r in p15.itertuples():
    s = r.sample_id; fa = meta(CAP/'4_parsimony/razor_hash-acc'/s/'parsimony.fasta'); pr = meta(PRED/(s + '_predictions.csv'))
    routes = dict(tryptic=route('5_sage', s), semi=route('7_semitryptic', s))
    tasks.append(dict(cohort='CAPSCAN', sample=s, raw=meta(RAW/(s + '.mzML'), False), fasta=fa, prediction=pr, design=design(r), routes=routes,
                      baseline_peptides=base[s], unspecific_receipt=A05R(s),
                      inputs=[fa, pr, base[s], routes['tryptic']['search'], routes['tryptic']['config'], routes['semi']['search'], routes['semi']['config']]))
m = {k: v for k, v in old.items() if k not in ('tasks', 'frozen_code')}; m.update(tasks=tasks, frozen_code=code, note='C2: 16 participant-15 tasks')
W(run/'inputs/manifest.json', m)

# ---------- A27 histidine per-acquisition (available = idx 0-15, unspecific = idx 16-31)
src = MM/'A27_4_7_histidine_defence'; old = J(src/'inputs/manifest.json'); run = F6/'A27'
code = copy_frozen(src, run, old['code'])
tasks = [dict(t, pilot=False, inputs=t['inputs']) for t in m['tasks']]  # same task records as A20 (fields identical)
tasks = [{k: t[k] for k in ['cohort', 'sample', 'raw', 'fasta', 'prediction', 'design', 'routes', 'baseline_peptides', 'unspecific_receipt', 'inputs', 'pilot']} for t in tasks]
m27 = {k: v for k, v in old.items() if k not in ('tasks', 'code')}; m27.update(tasks=tasks, code=code, note='C2: 16 participant-15 tasks')
W(run/'inputs/manifest.json', m27)

# ---------- A24 terminal inference: baseline331 and fixed347
src = MM/'A24_4_7_terminal_inference'; old = J(src/'inputs/manifest.json')
orig_cap = [t for t in old['tasks'] if t['cohort'] == 'CAPSCAN']; assert len(orig_cap) == 331
assert {t['sample'] for t in orig_cap} | set(p15.sample_id) == set(fixed.sample_id)
new24 = [dict(cohort='CAPSCAN', sample=r.sample_id, design=design(r),
              denovo_receipt=str(F6/'A26/results_revision2/CAPSCAN'/r.sample_id/'COMPLETE.json'),
              routes_receipt=str(F6/'A20/results/CAPSCAN'/r.sample_id/'COMPLETE.json')) for r in p15.itertuples()]
for name, tt in [('A24_baseline331', orig_cap), ('A24_fixed347', orig_cap + new24)]:
    run = F6/name; code = copy_frozen(src, run, old['code'])
    mm = {k: v for k, v in old.items() if k not in ('tasks', 'code')}; mm.update(tasks=tt, code=code, note=f'C2 {name}')
    W(run/'inputs/manifest.json', mm)

# ---------- A28 revision4 full_CAPSCAN: baseline331 and fixed347
src = MM/'A28_4_7_histidine_inference'; old = J(src/'inputs/manifest_revision4.json')
orig_cap = [t for t in old['tasks'] if t['cohort'] == 'CAPSCAN']; assert len(orig_cap) == 331
new28 = [dict(cohort='CAPSCAN', sample=r.sample_id, fasta=meta(CAP/'4_parsimony/razor_hash-acc'/r.sample_id/'parsimony.fasta'), design=design(r),
              receipts=dict(available=str(F6/'A27/results/available/CAPSCAN'/r.sample_id/'COMPLETE.json'),
                            unspecific=str(F6/'A27/results/unspecific/CAPSCAN'/r.sample_id/'COMPLETE.json'))) for r in p15.itertuples()]
for name, tt in [('A28_baseline331', orig_cap), ('A28_fixed347', orig_cap + new28)]:
    run = F6/name; code = copy_frozen(src, run, [x if isinstance(x, dict) else meta(x) for x in old['code']])
    mm = {k: v for k, v in old.items() if k not in ('tasks', 'code')}; mm.update(tasks=tt, code=code, note=f'C2 {name}')
    W(run/'inputs/manifest_revision4.json', mm)
for d in ['A26', 'A20', 'A27', 'A24_baseline331', 'A24_fixed347', 'A28_baseline331', 'A28_fixed347']: (F6/d/'logs').mkdir(exist_ok=True)
print('PREPARED chain; P15 tasks', len(p15))
