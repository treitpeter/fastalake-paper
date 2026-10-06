"""Revalidate wave-2 tasks whose receipt FAILED only because csv's default field limit (131072) was exceeded while
reading results.sage.tsv (very long protein lists). Sage exited 0 and wrote all outputs; the search is NOT rerun.
Repeats every post-search check of run_search.py with csv.field_size_limit raised. Usage: revalidate_fieldlimit.py <idx>"""
import sys, json, csv, re, time, os, shutil
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_search import sha, streamsha, clean
csv.field_size_limit(1 << 30)
O = Path(__file__).resolve().parent.parent
idx = int(sys.argv[1]); p = O/'receipts'/('%03d.json' % idx)
old = json.loads(p.read_text())
assert old['status'] == 'FAIL' and old['error'] == "Error('field larger than field limit (131072)')", old.get('error')
archive = p.with_name('%03d.failed_field_limit.json' % idx); assert not archive.exists(); shutil.copy2(p, archive)
t = json.loads((O/'tasks.json').read_text())[idx]; dest = Path(old['output'])
assert sha(t['source']) == t['source_sha'] and sha(t['config']) == t['config_sha'] and sha(t['binary']) == t['binary_sha']
cfg = json.loads((dest/'config.json').read_text())
expected = clean(json.loads(Path(t['source']).read_text())); expected['database']['enzyme']['min_len'] = 7
assert clean(cfg) == expected and clean(json.loads((dest/'results.json').read_text())) == expected
assert streamsha(cfg['database']['fasta']) == t['fasta_sha']
assert all(streamsha(k) == v for k, v in old['raw_hashes'].items())
accepted = 0; files = set()
with (dest/'results.sage.tsv').open() as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['rank'] == '1' and r['label'] == '1' and float(r['peptide_q']) <= .01:
            seq = re.sub(r'\[[^\]]*\]|[^A-Z]', '', r['peptide']); assert 7 <= len(seq) <= 50
            accepted += 1; files.add(r.get('filename', ''))
assert accepted > 0 and (dest/'results.sage.pin').stat().st_size > 0
assert (dest/'lfq.tsv').exists() == bool(cfg['quant']['lfq'])
rec = dict(old)
if t['kind'] == 'joint':
    from joint_lfq_audit import audit_joint
    rec['lfq_audit'] = audit_joint(dest, cfg)
    assert len(files) == len(cfg['mzml_paths']), (len(files), len(cfg['mzml_paths']))
names = ['results.json', 'results.sage.tsv', 'lfq.tsv'] + (['lfq_target_decoy_audit.tsv'] if t['kind'] == 'joint' else [])
rec.update(status='PASS', revalidated_utc_epoch=time.time(), revalidation_job=os.getenv('SLURM_JOB_ID'), prior_failed_receipt=str(archive),
           prior_failed_receipt_sha256=sha(archive), revalidation='Existing Sage outputs (exit 0); csv field_size_limit raised in the validator only; search not rerun; all post-search checks repeated',
           accepted_psms=accepted, n_acquisitions=len(cfg['mzml_paths']), resolved_settings_match=True, outputs={n: streamsha(dest/n) for n in names})
rec.pop('error', None); rec.pop('traceback', None)
tmp = p.with_suffix('.tmp'); tmp.write_text(json.dumps(rec, indent=2) + '\n'); tmp.replace(p)
print('PASS revalidated', idx, accepted)
