"""Independent Sage input-coverage audit (read-only).

For every Sage search in scope, record the mzML that Sage was configured to read,
whether that mzML ends with a closing </mzML>/</indexedmzML> tag, the declared
spectrumList count, and the largest scan number that actually appears in
results.sage.tsv. Sage 0.14.6 exits 0 when it hits an XML error part-way through
a file, so the highest scored scan relative to the declared count is the
independent signal for partial ingestion.

Scope (argv[1]):
  all       primary public-only 637 (current), archived pre-repair public PREDICT 306,
            assembly-inclusive canonical cohorts, public-only route searches
  preflight 2 positive (archived truncated) + 2 negative (repaired / ordinary) searches
Output: TSVs + SUMMARY.json in argv[2].
"""
import csv, glob, json, os, re, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

FL = Path('LAB_ROOT')
POC = FL / 'review/v15_integration/analyses/public_only_consistency_2026-09-29'
CAN = FL / 'canonical/publication_data'
SCAN_RE = re.compile(r'scan=(\d+)')
COUNT_RE = re.compile(rb'<spectrumList[^>]*count="(\d+)"')
SOFT_RE = re.compile(rb'<software id="([^"]+)" version="([^"]+)"')
csv.field_size_limit(10**9)


def searches(scope):
    out = []
    prim = json.loads((POC / 'PRIMARY_INPUTS.json').read_text())
    for t in prim['tasks']:
        d = Path(t['inputs']['config']['path']).parent
        out.append(dict(lineage='primary_public_current', cohort=t['cohort'], sample=t['sample'],
                        repaired=t.get('raw_repair_applied', False), dir=str(d)))
        if t['cohort'] == 'PREDICT' and t.get('archived_public_config'):
            out.append(dict(lineage='archived_public_prerepair', cohort='PREDICT', sample=t['sample'],
                            repaired=t.get('raw_repair_applied', False),
                            dir=str(Path(t['archived_public_config']).parent)))
    for coh in ['PREDICT', 'CAPSCAN', 'HSTOOL', 'CAMPI_STOOL', 'CAMPI_MOCK', 'SIHUMI', 'SALIVA']:
        for rj in sorted(glob.glob(str(CAN / 'cohorts' / coh / '5_sage/razor_hash-acc/*/results.json'))):
            d = Path(rj).parent
            out.append(dict(lineage='canonical_assembly_inclusive', cohort=coh, sample=d.name, repaired=False, dir=str(d)))
    for rj in sorted(glob.glob(str(POC / 'search_routes/results/*/search/results.json'))):
        d = Path(rj).parent
        coh, sample, route = d.parent.name.split('__')
        out.append(dict(lineage='public_route_' + route, cohort=coh, sample=sample, repaired=None, dir=str(d)))
    if scope == 'preflight':
        pos = [s for s in out if s['lineage'] == 'archived_public_prerepair' and s['sample'].endswith(('_115', '_181'))]
        neg = [s for s in out if s['lineage'] == 'primary_public_current' and s['sample'].endswith(('_115', '_114'))]
        out = pos + neg
    return out


def mzml_info(path):
    p = Path(path)
    r = dict(mzml=path, exists=p.exists())
    if not r['exists']:
        return r
    st = p.stat()
    r['bytes'] = st.st_size
    with open(p, 'rb') as f:
        head = f.read(400_000)
        f.seek(max(0, st.st_size - 1024))
        tail = f.read()
    m = COUNT_RE.search(head)
    r['declared_spectra'] = int(m.group(1)) if m else None
    r['software'] = ';'.join(f'{a.decode()}={b.decode()}' for a, b in SOFT_RE.findall(head)[:3])
    t = tail.rstrip()
    r['closing_tag'] = t.endswith(b'</indexedmzML>') or t.endswith(b'</mzML>')
    r['nul_in_tail'] = b'\x00' in tail
    return r


def scan_results(s):
    d = Path(s['dir'])
    r = dict(s)
    try:
        cfg = json.loads((d / 'results.json').read_text())
        r['mzml_cfg'] = cfg['mzml_paths'][0] if len(cfg['mzml_paths']) == 1 else '|'.join(cfg['mzml_paths'])
        r['sage_version'] = cfg.get('version')
        r['n_mzml'] = len(cfg['mzml_paths'])
    except Exception as e:  # noqa
        r['error'] = f'results.json: {e}'
        return r
    tsv = d / 'results.sage.tsv'
    if not tsv.exists():
        gz = d / 'results.sage.tsv.gz'
        r['error'] = 'results.sage.tsv missing' + (' (gz present)' if gz.exists() else '')
        return r
    maxscan, rows, peps = 0, 0, set()
    with open(tsv, newline='') as f:
        rd = csv.reader(f, delimiter='\t')
        h = next(rd)
        si, li, qi, pi = h.index('scannr'), h.index('label'), h.index('peptide_q'), h.index('peptide')
        fi = h.index('filename') if r.get('n_mzml', 1) > 1 else None
        for row in rd:
            rows += 1
            m = SCAN_RE.search(row[si])
            if m:
                v = int(m.group(1))
                if v > maxscan:
                    maxscan = v
            if row[li] == '1' and float(row[qi]) <= 0.01:
                peps.add(row[pi])
    r.update(rows=rows, max_scan=maxscan, target_peptides_q01_modified=len(peps))
    return r


def main():
    scope, outdir = sys.argv[1], Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)
    ss = searches(scope)
    workers = int(os.environ.get('SLURM_CPUS_PER_TASK', '2'))
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(scan_results, ss, chunksize=4))
    mz = sorted({p for r in res for p in str(r.get('mzml_cfg', '')).split('|') if p})
    with ProcessPoolExecutor(workers) as ex:
        info = {i['mzml']: i for i in ex.map(mzml_info, mz, chunksize=8)}
    for r in res:
        i = info.get(r.get('mzml_cfg'), {}) if r.get('n_mzml') == 1 else {}
        r.update({k: i.get(k) for k in ['bytes', 'declared_spectra', 'closing_tag', 'nul_in_tail', 'software']})
        if r.get('declared_spectra') and r.get('max_scan'):
            r['max_scan_fraction'] = round(r['max_scan'] / r['declared_spectra'], 4)
        frac = r.get('max_scan_fraction')
        r['flag'] = bool((r.get('closing_tag') is False) or (frac is not None and frac < 0.9) or r.get('error'))
    cols = ['lineage', 'cohort', 'sample', 'repaired', 'flag', 'closing_tag', 'nul_in_tail', 'declared_spectra',
            'max_scan', 'max_scan_fraction', 'rows', 'target_peptides_q01_modified', 'bytes', 'software',
            'sage_version', 'n_mzml', 'mzml_cfg', 'dir', 'error']
    with open(outdir / 'search_coverage.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for r in res:
            w.writerow(r)
    with open(outdir / 'mzml_integrity.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, ['mzml', 'exists', 'bytes', 'declared_spectra', 'closing_tag', 'nul_in_tail', 'software'],
                           delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for k in mz:
            w.writerow(info[k])
    summ = {}
    for r in res:
        k = r['lineage'] + '|' + r['cohort']
        s = summ.setdefault(k, dict(n=0, flagged=[], errors=0, versions=set()))
        s['n'] += 1
        s['versions'].add(str(r.get('sage_version')))
        if r.get('error'):
            s['errors'] += 1
        if r['flag']:
            s['flagged'].append(dict(sample=r['sample'], frac=r.get('max_scan_fraction'), closing=r.get('closing_tag'),
                                     peps=r.get('target_peptides_q01_modified'), err=r.get('error')))
    for v in summ.values():
        v['versions'] = sorted(v['versions'])
    bad_mz = [i for i in info.values() if i.get('closing_tag') is False or not i.get('exists')]
    json.dump(dict(scope=scope, n_searches=len(res), n_mzml=len(mz), mzml_without_closing_tag=bad_mz,
                   by_lineage=summ), open(outdir / 'SUMMARY.json', 'w'), indent=1, default=str)
    print(json.dumps(dict(n=len(res), n_mzml=len(mz), bad_mzml=len(bad_mz),
                          flagged={k: len(v['flagged']) for k, v in summ.items()}), indent=1))


if __name__ == '__main__':
    main()
