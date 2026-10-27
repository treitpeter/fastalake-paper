"""D11 10-fold spectrum cross-fit: shared paths, frozen-code loading, fold rule, I/O helpers.

Frozen producer: split_half_repair/pipeline.py (sha256 d983b1be...), imported unchanged.
Task inputs: split_half_repair/mzml_recovery_2026-09-27/MANIFEST.json (710 tasks; five PREDICT
mzML paths replaced, as used by every split-half and round-3 search).
Roster: iteration3/disjoint/roster.tsv (578 matched FL/MG acquisitions; sha256 33eb88f8...).
"""
import csv, hashlib, importlib.util, json, re, contextlib, os
from pathlib import Path
from lxml import etree as ET

DES = Path(__file__).resolve().parent
CF = DES.parent
REVIEW = Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03')  # public347: production review root (roster, iteration3), copy lives elsewhere
SH = Path('LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/split_half_repair')
REC = SH / 'mzml_recovery_2026-09-27'
PIPE_SHA = 'd983b1be57abbf8edee1a237fbaff633ee8353ff4e56098c1df3854078d511c4'
ROSTER = Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/disjoint/roster_mg295.tsv')  # mg295: production 578-row roster + 13 participant-15 MG-matched rows (578-590)
ROSTER_SHA = '580795333e612c6a071b24fa20e10c96c4de7e78a230f498879014c5455c354c'
CONF_MAN = Path('LAB_ROOT/review/v15_integration/ready_integration_2026-09-28/source_additions/Review_v16/Confidence/lineage/manifest.json')

K = 10
FOLD_SEED = 'FastaLake-D11-crossfit-20261004'
csv.field_size_limit(10_000_000)


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


assert sha(SH / 'pipeline.py') == PIPE_SHA, 'frozen pipeline changed'
_spec = importlib.util.spec_from_file_location('splitpipe', SH / 'pipeline.py')
P = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(P)
MAN = json.loads(Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/public347/crossfit_capscan347/inputs/MANIFEST_347.json').read_text())  # mg295: the public347 manifest, read in place  # public347: production manifest, CAPSCAN tasks repointed to A08_347 + 16 appended
for _c in MAN['code']:
    if Path(_c['path']).name in ('pipeline.py', 'fasta_extractor', 'parsimony_engine', 'sage'):
        assert sha(_c['path']) == _c['sha256'], _c['path']
assert str(P.SAGE) == [c['path'] for c in MAN['code'] if c['path'].endswith('/sage')][0]
M = P.M


def fold_of(scan, k=K):
    """Deterministic fold of a native scan number: SHA-256 of '<seed>:<scan>', first 8 bytes, mod k."""
    if k == 1:
        return 0
    return int.from_bytes(hashlib.sha256(f'{FOLD_SEED}:{int(scan)}'.encode()).digest()[:8], 'big') % k


def roster():
    assert sha(ROSTER) == ROSTER_SHA
    return list(csv.DictReader(open(ROSTER), delimiter='\t'))


def task_for(row):
    t = MAN['tasks'][int(row['split_index'])]
    assert t['cohort'] == row['cohort'] and t['sample'] == row['sample'], (t['sample'], row['sample'])
    assert t['mzml'] == row['mzml_split'] and t['csv'] == row['csv']
    return t


def dump(p, j):
    P.dump(p, j)


def write_partition_csv(src, dests, route):
    """Copy of the frozen P.split() row handling, generalised to any routing.

    dests: {key: path}; route(scan) -> iterable of keys the row is written to.
    Same DictReader/DictWriter formatting and the same assertions (integer spec_idx, unique scan,
    beam_rank 0). Returns per-key row counts and per-key sorted scan lists.
    """
    for p in dests.values():
        Path(p).parent.mkdir(parents=True, exist_ok=True)
    counts = {k: 0 for k in dests}; ids = set()
    with contextlib.ExitStack() as st:
        fs = {k: st.enter_context(open(p, 'w', newline='')) for k, p in dests.items()}
        with open(src) as f:
            rd = csv.DictReader(f)
            ws = {k: csv.DictWriter(v, fieldnames=rd.fieldnames) for k, v in fs.items()}
            for w in ws.values():
                w.writeheader()
            for row in rd:
                n = float(row['spec_idx']); assert n.is_integer() and n >= 0
                scan = int(n); assert scan not in ids; ids.add(scan)
                assert float(row.get('beam_rank', '0')) == 0, 'not top-ranked export'
                for k in route(scan):
                    if k in ws:
                        ws[k].writerow(row); counts[k] += 1
    return counts, ids


def csv_scans(path):
    with open(path) as f:
        return {int(float(r['spec_idx'])) for r in csv.DictReader(f)}


def fold_mzml(src, csvpath, dests, assign, receipt_path):
    """Frozen P.odd_mzml() generalised from 'keep odd MS2' to 'MS2 of fold k' for several outputs
    in one pass. dests: {k: path}. assign(scan) -> fold. Every MS1 spectrum goes to every output.
    With dests={1: p} and assign=scan%2 the output is the frozen odd-scan file (tested byte-for-byte).
    """
    pred = {}
    with open(csvpath) as f:
        for row in csv.DictReader(f):
            pred[int(float(row['spec_idx']))] = (float(row['prec_mz']), int(float(row['prec_charge'])))
    count = {k: 0 for k in dests}; ms2 = {k: 0 for k in dests}
    seen = set(); matched = 0; ms1 = 0; run_attrs = None; root_attrs = None; root_ns = None; headers = []; slist_attrs = None
    for event, e in ET.iterparse(src, events=['start', 'end'], huge_tree=True):
        if event == 'start':
            if e.tag == M + 'mzML': root_attrs = dict(e.attrib); root_ns = e.nsmap
            if e.tag == M + 'run': run_attrs = dict(e.attrib)
            if e.tag == M + 'spectrumList': slist_attrs = dict(e.attrib)
            continue
        if e.tag == M + 'spectrum':
            cv = {x.get('accession'): x.get('value') for x in e.iter(M + 'cvParam')}
            level = int(cv['MS:1000511']); scan = int(re.search(r'scan=(\d+)', e.get('id')).group(1))
            if level == 2:
                assert scan not in seen; seen.add(scan)
                if scan in pred:
                    mz, ch = pred[scan]; assert abs(float(cv['MS:1000744']) - mz) < .001, (scan, 'm/z')
                    assert int(cv.get('MS:1000041', '0')) == ch, (scan, 'charge'); matched += 1
                f = assign(scan)
                if f in count:
                    count[f] += 1; ms2[f] += 1
            else:
                ms1 += level == 1
                for k in count:
                    count[k] += 1
            e.clear()
            while e.getprevious() is not None: del e.getparent()[0]
        elif e.getparent() is not None and e.getparent().tag == M + 'mzML' and e.tag != M + 'run':
            headers.append(ET.tostring(e))
    assert matched == len(pred) and set(pred) <= seen
    assert root_attrs is not None and run_attrs is not None and slist_attrs is not None
    written = {k: 0 for k in dests}; written_scans = {k: [] for k in dests}
    with contextlib.ExitStack() as st:
        lists = {}
        for k, path in dests.items():
            xf = st.enter_context(ET.xmlfile(str(path), encoding='utf-8'))
            xf.write_declaration()
            st.enter_context(xf.element(M + 'mzML', root_attrs, nsmap=root_ns))
            for h in headers: xf.write(ET.fromstring(h))
            st.enter_context(xf.element(M + 'run', run_attrs))
            sa = dict(slist_attrs); sa['count'] = str(count[k])
            st.enter_context(xf.element(M + 'spectrumList', sa))
            lists[k] = xf
        for _, e in ET.iterparse(src, events=['end'], tag=M + 'spectrum', huge_tree=True):
            level = int(next(x.get('value') for x in e.findall(M + 'cvParam') if x.get('accession') == 'MS:1000511'))
            scan = int(re.search(r'scan=(\d+)', e.get('id')).group(1))
            targets = list(dests) if level != 2 else [assign(scan)]
            for k in targets:
                if k in lists:
                    e.set('index', str(written[k])); lists[k].write(e); written[k] += 1
                    if level == 2: written_scans[k].append(scan)
            e.clear()
            while e.getprevious() is not None: del e.getparent()[0]
    assert written == count
    for k in dests:
        assert all(assign(s) == k for s in written_scans[k])
    rec = dict(source=str(src), source_bytes=Path(src).stat().st_size, predictions=len(pred), predictions_matched=matched,
               original_ms2=len(seen), retained_ms1=ms1,
               outputs={str(k): dict(path=str(p), ms2=ms2[k], spectra=written[k], sha256=sha(p), bytes=Path(p).stat().st_size) for k, p in dests.items()},
               native_ids_preserved=True, fold_seed=FOLD_SEED, utc=P.now())
    dump(receipt_path, rec)
    return rec, written_scans


def select(sample, csvpath, db, out):
    """Frozen P.select() with only the output folder changed (P.select hard-codes the frozen tree)."""
    pool = out / 'queries'; pool.mkdir(parents=True, exist_ok=True)
    dest = pool / (sample + '_predictions.csv')
    if dest.is_symlink() or dest.exists(): dest.unlink()
    dest.symlink_to(csvpath)
    P.extract(pool, db, out / 'retrieval')
    par = out / 'parsimony'; par.mkdir(exist_ok=True)
    P.command([P.PAR, '--strategy', 'razor', '--tiebreak', 'hash-acc', '--fasta', out / 'retrieval/selected.fasta', '--peptides', csvpath,
               '--output-dir', par, '--sample', sample, '--min-length', '9', '--max-length', '50', '--normalize-il', '-t', '8'], par / 'select')
    choices = list(par.glob('*.fasta')); assert len(choices) == 1, [str(p) for p in choices]
    return choices[0]


def sage_search(config_path, fasta, mzml, d):
    """Frozen per-acquisition Sage configuration; only database.fasta, mzml_paths, output_directory changed."""
    d.mkdir(parents=True, exist_ok=True)
    cfg = json.loads(Path(config_path).read_text()); cfg.pop('version', None)
    cfg['database']['fasta'] = str(fasta); cfg['mzml_paths'] = [str(mzml)]; cfg['output_directory'] = str(d)
    dump(d / 'config.json', cfg)
    P.command([P.SAGE, d / 'config.json', '--write-pin'], d / 'sage')
    return d / 'results.sage.tsv'


def strip(pep):
    p = re.sub(r'\[[^\]]*\]', '', pep).replace('I', 'L')
    assert p.isalpha() and p.isupper(), pep
    return p


def accepted_set(psms, qcol='peptide_q'):
    """Frozen pipeline.accepted() rule on a DataFrame: rank 1, target, q <= 0.01, mods stripped, I/L, 5-50."""
    m = (psms['rank'] == 1) & (psms['label'] == 1) & (psms[qcol] <= 0.01)
    acc = {strip(p) for p in psms.loc[m, 'peptide']}
    return {p for p in acc if 5 <= len(p) <= 50}
