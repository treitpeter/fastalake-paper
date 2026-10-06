#!/usr/bin/env python3
"""Reconvert a corrupt mzML from its hash-identified original instrument RAW.

Never alter historical inputs. Gate reuse on full scan counts and prespecified
intact-spectrum comparison against the historical vendor-centroided data.
"""
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from pyteomics import mzml

MM = Path(__file__).resolve().parents[1]
ROOT = MM / 'recovery_2026-09-15/raw_recovery'
m = json.loads((ROOT / 'manifest.json').read_text())
t = m['tasks'][int(sys.argv[1])]
out = ROOT / t['sample']
out.mkdir(exist_ok=True)
assert not (out / 'COMPLETE.json').exists()

def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def dump(name, x):
    (out / name).write_text(json.dumps(x, indent=2) + '\n')

assert digest(__file__) == m['worker_sha256']
assert digest(m['converter']['path']) == m['converter']['sha256']
old = Path(t['mzml'])
raw = Path(t['raw'])
with old.open('rb') as f:
    header = f.read(16000).decode()
expected_sha1 = re.search(r'name="SHA-1" value="([0-9a-f]{40})"', header)[1]
expected_scans = int(re.search(r'<spectrumList count="(\d+)"', header)[1])
raw_sha1 = digest(raw, 'sha1')
assert raw_sha1 == expected_sha1, 'Original instrument RAW does not match the mzML source identity'
identity = dict(source_mzml=str(old), source_mzml_sha256=digest(old), instrument_raw=str(raw), instrument_raw_sha1=raw_sha1,
                expected_spectra=expected_scans, converter=m['converter'], worker_sha256=m['worker_sha256'])
dump('SOURCE_IDENTITY.json', identity)
env = dict(os.environ, DOTNET_BUNDLE_EXTRACT_BASE_DIR=str(ROOT / 'dotnet_extract'))
(ROOT / 'dotnet_extract').mkdir(exist_ok=True)
new = out / (t['sample'] + '.mzML')
assert not new.exists(), 'Preserve any previous conversion attempt'
argv = [m['converter']['path'], '-i', str(raw), '-b', str(new), '-f', '2']
with (out / 'convert.out').open('w') as so, (out / 'convert.err').open('w') as se:
    subprocess.run(argv, stdout=so, stderr=se, env=env, check=True)
dump('CONVERSION.json', {'argv': argv, 'scope': 'New vendor-centroided conversion; not yet approved as equivalent', 'job_id': os.environ.get('SLURM_JOB_ID')})

def snapshot(sp):
    ions = sp.get('precursorList', {}).get('precursor', [{}])[0].get('selectedIonList', {}).get('selectedIon', [{}])[0]
    return {'mz': np.asarray(sp['m/z array'], np.float32), 'intensity': np.asarray(sp['intensity array'], np.float32),
            'precursor_mz': float(ions.get('selected ion m/z', float('nan'))), 'charge': ions.get('charge state')}

# Select every 257th native scan from the valid prefix, independently of peptide/H outcomes.
reference = {}
old_count = 0
old_error = None
try:
    with mzml.MzML(str(old), use_index=False, iterative=True, decode_binary=True) as reader:
        for sp in reader:
            old_count += 1
            scan = int(re.search(r'scan=(\d+)', sp['id'])[1])
            if sp.get('ms level') == 2 and scan % 257 == 0:
                reference[scan] = snapshot(sp)
except Exception as exc:
    old_error = {'class': type(exc).__name__, 'message': str(exc)}
assert old_error is not None, 'Expected historical parser failure was not reproduced'
checks = []
count = 0
ms2 = 0
with mzml.MzML(str(new), use_index=False, iterative=True, decode_binary=True) as reader:
    for sp in reader:
        count += 1
        mz = np.asarray(sp['m/z array'])
        intensity = np.asarray(sp['intensity array'])
        assert len(mz) == len(intensity) and np.isfinite(mz).all() and np.isfinite(intensity).all() and (intensity >= 0).all()
        ms2 += sp.get('ms level') == 2
        scan = int(re.search(r'scan=(\d+)', sp['id'])[1])
        if scan in reference:
            a, b = reference[scan], snapshot(sp)
            checks.append({'scan': scan, 'mz_float32_exact': np.array_equal(a['mz'], b['mz']),
                           'intensity_float32_exact': np.array_equal(a['intensity'], b['intensity']),
                           'precursor_within_1ppm': abs(a['precursor_mz'] - b['precursor_mz']) / a['precursor_mz'] <= 1e-6,
                           'charge_exact': a['charge'] == b['charge']})
validation = {'old_spectra_read_before_error': old_count, 'old_error': old_error, 'expected_spectra': expected_scans,
              'reconverted_spectra': count, 'reconverted_MS2': ms2, 'sampled_reference_scans': len(reference), 'checks': checks}
dump('VALIDATION.json', validation)
assert count == expected_scans, 'Reconverted scan roster differs'
assert len(checks) == len(reference) and len(checks) >= 30
assert all(all(v for k, v in row.items() if k != 'scan') for row in checks), 'Vendor conversion equivalence not established'
dump('COMPLETE.json', dict(status='PASS', scope='Source reconversion only; downstream searches must be rerun and validated separately',
                           sample=t['sample'], identity=identity, new_mzml=str(new), new_mzml_sha256=digest(new),
                           validation=str(out / 'VALIDATION.json'), spectra=count, MS2=ms2, matched_intact_MS2=len(checks),
                           job_id=os.environ.get('SLURM_JOB_ID')))
print('SOURCE RECONVERSION PASS', t['sample'], count, len(checks), flush=True)
