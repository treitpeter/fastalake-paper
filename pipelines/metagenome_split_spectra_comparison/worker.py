"""MG comparator on the disjoint-half spectra (D7, round 3).

For roster row i:
  1. Rebuild the odd-scan mzML with the FROZEN split-half code (pipeline.odd_mzml, imported
     unchanged and hash-checked) from the exact mzML the FastaLake split searches used, and
     require its SHA-256 to equal the archived PARITY.json output_sha256.
  2. MG_ODD: the FastaLake HALF_DISJOINT config.json, changing only database.fasta to the
     frozen MG database (output_directory changed to this folder). Same Sage binary.
  3. MG_FULL (sanity): the frozen MG results.json config, unchanged apart from output_directory
     and the dropped version field (as in the split pipeline), on the frozen full mzML.
     The canonical 5-50 and 9-50 accepted sets must equal the frozen MG sets.
  4. Count with the frozen pipeline.accepted() (rank 1, target, peptide_q<=0.01, mods stripped,
     I/L collapsed). MG_FULL uses a copy of the counter without the odd-scan assertion.
"""
import sys, json, csv, re, hashlib, importlib.util, subprocess, gzip, shutil, os, time
from pathlib import Path

W = Path(__file__).resolve().parent
SH = Path('LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/split_half_repair')
PIPE_SHA = 'd983b1be57abbf8edee1a237fbaff633ee8353ff4e56098c1df3854078d511c4'
CONF_MAN = Path('LAB_ROOT/review/v15_integration/ready_integration_2026-09-28/source_additions/Review_v16/Confidence/lineage/manifest.json')


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


assert sha(SH / 'pipeline.py') == PIPE_SHA, 'frozen pipeline changed'
spec = importlib.util.spec_from_file_location('splitpipe', SH / 'pipeline.py')
P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)
SAGE = P.SAGE


def count_any(path, out):
    """pipeline.accepted() without the odd-scan assertion (full-spectra sanity search)."""
    acc = set(); nrows = 0; ndecoy = 0
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            nrows += 1; ndecoy += int(r['label']) == -1
            if int(r['rank']) == 1 and int(r['label']) == 1 and float(r['peptide_q']) <= .01:
                pep = re.sub(r'\[[^\]]*\]', '', r['peptide']).replace('I', 'L')
                assert pep.isalpha() and pep.isupper(), pep
                if 5 <= len(pep) <= 50:
                    acc.add(pep)
    assert nrows > 0 and ndecoy > 0
    out.write_text('\n'.join(sorted(acc)) + '\n')
    return dict(psms=nrows, decoy_psms=ndecoy, peptides_5_50=len(acc), peptides_9_50=sum(len(p) >= 9 for p in acc)), acc


def search(cfg, d):
    d.mkdir(parents=True, exist_ok=True)
    cfg['output_directory'] = str(d)
    P.dump(d / 'config.json', cfg)
    t = time.time()
    P.command([SAGE, d / 'config.json', '--write-pin'], d / 'sage')
    return time.time() - t


def tidy(d):
    for n in ['results.sage.pin']:
        (d / n).unlink(missing_ok=True)
    src = d / 'results.sage.tsv'
    with open(src, 'rb') as a, gzip.open(str(src) + '.gz', 'wb', compresslevel=4) as b:
        shutil.copyfileobj(a, b)
    src.unlink()


def main(i, arms):
    rows = list(csv.DictReader(open(W / os.environ.get('FL_ROSTER', 'roster.tsv')), delimiter='\t'))
    r = rows[i]
    out = W / os.environ.get('FL_OUTDIR', 'outputs') / r['cohort'] / r['sample']; out.mkdir(parents=True, exist_ok=True)
    res = dict(i=i, cohort=r['cohort'], sample=r['sample'], sage=str(SAGE), sage_sha256=sha(SAGE), pipeline_sha256=PIPE_SHA)
    if 'odd' in arms:
        scratch = out / 'odd_only.mzML'
        t = time.time()
        P.odd_mzml(r['mzml_split'], r['csv'], scratch, out / 'PARITY.json')
        res['odd_split_seconds'] = time.time() - t
        got = sha(scratch)
        res['odd_mzml_sha256'] = got; res['odd_mzml_identical'] = got == r['odd_mzml_sha256']
        assert got == r['odd_mzml_sha256'], 'odd-scan mzML differs from the archived FastaLake split'
        cfg = json.loads(Path(r['fl_disjoint_config']).read_text())
        assert cfg['mzml_paths'][0].endswith('odd_only.mzML')
        res['fl_disjoint_fasta'] = cfg['database']['fasta']
        cfg['database']['fasta'] = r['mg_fasta']; cfg['mzml_paths'] = [str(scratch)]
        res['mg_odd_seconds'] = search(cfg, out / 'MG_ODD')
        res['MG_ODD'], _ = P.accepted(out / 'MG_ODD/results.sage.tsv', out / 'MG_ODD/ACCEPTED.txt')
        tidy(out / 'MG_ODD'); scratch.unlink()
    if 'full' in arms:
        cfg = json.loads(Path(r['mg_frozen_config']).read_text()); cfg.pop('version', None)
        res['mg_full_seconds'] = search(cfg, out / 'MG_FULL')
        res['MG_FULL'], acc = count_any(out / 'MG_FULL/results.sage.tsv', out / 'MG_FULL/ACCEPTED.txt')
        man = json.loads(CONF_MAN.read_text())
        o = next(x for x in man['ordinary'] if x['cohort'] == r['cohort'] and x['sample'] == r['sample'])
        for L in ['5', '9']:
            exp = set(Path(o['expected_sets']['MG'][L]['path']).read_text().split())
            mine = {p for p in acc if len(p) >= int(L)}
            res[f'mg_full_equals_frozen_{L}_50'] = mine == exp
            res[f'mg_full_frozen_{L}_50'] = len(exp)
            res[f'mg_full_only_{L}'] = len(mine - exp); res[f'frozen_only_{L}'] = len(exp - mine)
        tidy(out / 'MG_FULL')
    res['utc'] = P.now()
    P.dump(out / ('COMPLETE.json' if set(arms) == {'odd', 'full'} else f"COMPLETE_{'_'.join(arms)}.json"), res)
    print(json.dumps(res), flush=True)


if __name__ == '__main__':
    # usage: worker.py <list file>; SLURM_ARRAY_TASK_ID indexes the list of roster rows
    ids = open(W / sys.argv[1]).read().split()
    main(int(ids[int(os.environ['SLURM_ARRAY_TASK_ID'])]), ['odd', 'full'])
