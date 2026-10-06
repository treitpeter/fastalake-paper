#!/usr/bin/env python3
"""Per-task receipt after frozen/search.py has written COMPLETE.json (PASS).
Records input sha256, tool versions, peptide counts at q<=0.01 for 5-50 and 9-50
residues, and /usr/bin/time resources; then removes the private mzML copy
(its sha256 stays recorded in COMPLETE.json and here)."""
import json, os, pathlib, re, sys
RUN = pathlib.Path(sys.argv[1]).resolve(); idx = int(sys.argv[2])
sys.path.insert(0, str(RUN / 'frozen'))
from lib import dump, sha  # noqa: E402
m = json.loads((RUN / 'inputs/manifest.json').read_text()); t = m['tasks'][idx]
case = RUN / 'results' / t['id']
c = json.loads((case / 'COMPLETE.json').read_text())
assert c['status'] == 'PASS' and c['sample'] == t['sample'] and c['arm'] == t['arm']
n = {lo: sum(1 for l in open(case / f'accepted_{lo}_50.txt') if l.strip()) for lo in (5, 9)}
assert n[9] == c['accepted_length9_50']
def tv(p):
    """Parse peak RSS, wall time and exit status from a /usr/bin/time -v log."""
    s = p.read_text(errors='replace') if p.exists() else ''
    g = lambda k: (re.search(re.escape(k) + r':\s*(.+)', s) or [None, None])[1]
    return dict(max_rss_kb=int(g('Maximum resident set size (kbytes)') or 0), elapsed=g('Elapsed (wall clock) time (h:mm:ss or m:ss)'), exit_status=g('Exit status'))
def first(p, pat):
    """First line of a file matching a regular expression, or None."""
    if not p.exists(): return None
    for l in p.read_text(errors='replace').splitlines():
        if re.search(pat, l): return l.strip()
mz = case / (t['sample'] + '.mzML'); removed = None
if mz.exists():
    assert sha(mz) == t['raw']['sha256']
    mz.unlink(); removed = dict(path=str(mz), sha256=t['raw']['sha256'])
dump(case / 'RECEIPT.json', dict(
    status='PASS', id=t['id'], cohort=t['cohort'], arm=t['arm'], sample=t['sample'],
    inputs=dict(fasta=t['fasta'], raw=t['raw'], source_config=t['source_config']),
    frozen_code={pathlib.Path(x['path']).name: x['sha256'] for x in m['frozen_code']},
    versions=dict(msfragger=dict(jar_sha256=m['jar']['sha256'], banner=first(case / 'msfragger.out', r'MSFragger version')),
                  percolator=dict(sha256=m['percolator']['sha256'], banner=first(case / 'percolator.err', r'Percolator version')),
                  java=m['java']['sha256'], seed=m['seed']),
    target_entries=c['target_entries'], labels=c['labels'],
    peptides_q01=dict(len5_50=n[5], len9_50=n[9]),
    resources=dict(msfragger=tv(case / 'msfragger.err'), percolator=tv(case / 'percolator.err')),
    slurm=dict(job=os.environ.get('SLURM_JOB_ID'), array_job=os.environ.get('SLURM_ARRAY_JOB_ID'), task=os.environ.get('SLURM_ARRAY_TASK_ID'), node=os.environ.get('SLURMD_NODENAME')),
    mzml_copy_removed=removed))
print('RECEIPT', t['id'], n[5], n[9], flush=True)
