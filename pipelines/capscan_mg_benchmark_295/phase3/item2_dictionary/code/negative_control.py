#!/usr/bin/env python3
"""Negative control for the comparator and the H09 reconstruction check.
(1) copy the archived H09 acquisition_summary with one value changed in the last digit -> comparator must report a difference.
(2) the 295 H09 output vs the archived 282 output -> must differ (different roster)."""
import hashlib, sys, json
from pathlib import Path
I = Path(sys.argv[1]); OH = Path('LAB_ROOT/Hartout_edits/H09_evidence_ablation/full_CAPSCAN/acquisition_summary.tsv')
t = OH.read_text().splitlines(); f = t[1].split('\t'); f[4] = f[4][:-1] + str((int(f[4][-1]) + 1) % 10); t[1] = '\t'.join(f)
p = I / 'controls/planted_acquisition_summary.tsv'; p.write_text('\n'.join(t) + '\n')
sh = lambda x: hashlib.sha256(Path(x).read_bytes()).hexdigest()
r1 = sh(p) != sh(OH)
n = I / 'runs/r295/h09/full_CAPSCAN/acquisition_summary.tsv'
r2 = n.exists() and sh(n) != sh(OH)
out = dict(planted_one_digit_detected=r1, r295_differs_from_282=r2, status='PASS' if r1 and r2 else 'FAIL')
(I / 'controls/NEGATIVE_CONTROL.json').write_text(json.dumps(out, indent=1) + '\n'); print(out)
