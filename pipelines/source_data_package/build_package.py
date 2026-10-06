"""Build the round-2 Source Data package from SELECTION.tsv.

Usage:
  build_package.py SELECTION.tsv OUT_ROOT WORK_DIR KEY_DIR

For every selected member: read the original bytes, apply the special transform
for that package path (specials.py), then
  (1) path cleaning: absolute paths and internal working-directory paths are
      reduced to file names (URLs and package-relative paths are kept),
  (2) PREDICT participant and tube codes -> random pseudonyms (key in KEY_DIR,
      outside the package); rows that were sorted by study code are re-sorted
      by pseudonym so that row order does not reveal the key,
  (3) CapScan capsule device IDs -> removed (device_id column dropped; IDs inside
      composite strings replaced by random capsule labels); swallow/recovery
      dates and times dropped,
  (4) xlsx cell comments removed and document properties neutralised.
Writes the files under OUT_ROOT/FastaLake_Source_Data/<package_path> and
MANIFEST_build.tsv, REDACTION_LOG.tsv, ORDER_LOG.tsv in WORK_DIR.
The residual gate is a separate script (gate.py) run on the finished package.
"""
import csv, gzip, hashlib, io, json, os, re, sys, zipfile, secrets
sys.path.insert(0, os.path.dirname(__file__))
import selection as S
import specials as SP

PKG = 'FastaLake_Source_Data'
csv.field_size_limit(10**9)

# ----------------------------------------------------------------- path rules
R0 = '/f[s]/gpfs41/lv12/fileset02/pool/pool-mann-projects/Peter/FastaLake/'
PATH_RULES = [
    (re.compile(r'/f[s]/gpfs41/lv12/fileset02/pool/pool-mann-projects/Peter/FastaLake/'), 'FastaLake/'),
    (re.compile(r'/f[s]/pool/pool-mann-projects/Peter/FastaLake/'), 'FastaLake/'),
    (re.compile(r'/f[s]/gpfs41/lv12/fileset02/pool/pool-mann-projects/'), '<storage>/'),
    (re.compile(r'/f[s]/pool/pool-mann-projects/'), '<storage>/'),
    (re.compile(r'/f[s]/hom[e]/[A-Za-z0-9_.-]+/'), '<home>/'),
    (re.compile(r'/f[s]/[A-Za-z0-9_.-]+(/[A-Za-z0-9_.-]+)*'), '<cluster-path>'),
    (re.compile(r'(?<![A-Za-z0-9_.<>-])/tmp/[^\s"\'\t,;]*'), '<tmp-path>'),
    (re.compile(r'(?<![A-Za-z0-9_.<>-])/scratch/[^\s"\'\t,;]*'), '<scratch-path>'),
    (re.compile(r'(?<![A-Za-z0-9_.<>-])/(Users|home)/[A-Za-z0-9_.-]+/'), '<home>/'),
]
RAW_PATH = re.compile(r'/f[s]/(?:gpfs41/lv12/fileset02/)?pool/pool-mann-projects/(?!Peter/FastaLake/)[^\s"\'\t,;|<>]*')
WIN_URI = re.compile(r'file:///[A-Za-z]:\\[^\s"\'\t|]*')
WIN_PATH = re.compile(r'(^|[^A-Za-z])([A-Za-z]:\\[^\s"\'\t|]+)')
TEXT_SUBS = [
    (re.compile(r'\bPeter confirmation\b'), 'author confirmation'),
    (re.compile(r'\(Peter, (\d{4}-\d{2}-\d{2})\)'), r'(author confirmation, \1)'),
    (re.compile(r'\bPeter confirms\b'), 'The authors confirm'),
    (re.compile(r'permutation_p_from_codex_full'), 'permutation_p_full'),
    (re.compile(r'matches_codex_all_observed'), 'matches_independent_all_observed'),
    (re.compile(r'\b([A-Za-z0-9_.-]+)\.sbatch\b'), r'\1 job script'),
]
EMAIL = re.compile(r'[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}')
# path-like tokens: no whitespace, quotes or separators; must contain '/'
PATH_TOKEN = re.compile(r'[^\s"\'\t,;|<>()\[\]{}=]*/[^\s"\'\t,;|()\[\]{}=]*')
INTERNAL_ROOTS = {'FastaLake', '<storage>', '<home>', '<cluster-path>', '<tmp-path>', '<scratch-path>', 'canonical',
                  'MM_review_analysis', 'analyses', 'projects', 'review', 'manuscript', 'legacy', 'worktrees', 'publication_data',
                  'publication', 'source_data', 'code', 'scratchpad', 'atelier', 'external_runs', 'v26', 'paper_repo', '..', '.'}
INTERNAL_MARK = re.compile(r'(00_PACKAGE|cleanup_2026|comment_audit|review_decisions|MM_review_analysis|fable|jess_integration|claude_review|'
                           r'passage_review|review_package|biology_refit|manuscript_completion|analysis_followup|nine_point|reappraisal|'
                           r'manuscript_final_pass|main_methods_pass|mm_v11|methods_coordination|Masters_Thesis|publication_data|'
                           r'final_revision|thread_selection_math|recovery_2026|/review/|^review/)', re.I)

def clean_path_token(tok):
    if re.match(r'(?i)(https?|ftp|doi):', tok) or tok.startswith('//') or tok.startswith('www.'):
        return tok
    segs = [s for s in tok.split('/') if s != '']
    if not segs:
        return tok
    if segs[0] in INTERNAL_ROOTS or segs[0].startswith('00_PACKAGE') or INTERNAL_MARK.search(tok):
        # keep a trailing ':member' or '!sheet' spec of the last segment
        return segs[-1]
    return tok

def clean_paths(s):
    return PATH_TOKEN.sub(lambda m: clean_path_token(m.group(0)), s)

ID_HEADER = re.compile(r'(participant|patient|donor|subject|individual|person|^pid$|_pid$)', re.I)
DROP_COLS = {'device_id', 'swallow_date', 'swallow_time', 'recover_date', 'recover_time'}
JOBKEY = re.compile(r'"(job_?ids?|slurm[a-z_]*|array_job[a-z_]*)"\s*:', re.I)
JOBKEY_NAME = re.compile(r'^(job_?ids?|slurm[a-z_]*|array_job[a-z_]*)$', re.I)
POSITIONAL = re.compile(r'(axes|order|columns|MATRIX|_values\.npz|heatmap)', re.I)
FORBIDDEN_NAME = re.compile(r'(\.log$|SUBMI|squeue|sacct|HANDOFF|RESUME|CLOSEOUT|QUEUE|failed_|before_|\.docx$|_DRAFT|PASSAGE_CHANGES|jobids|\.py$|\.sh$|\.sbatch$|\.ipynb$|\.R$)', re.I)


def load_or_make_key(path, items, prefix, width):
    """Random pseudonyms; the key file is reused if it already covers exactly these items."""
    if os.path.exists(path):
        k = {}
        for r in csv.DictReader(open(path), delimiter='\t'):
            k[r['study_code']] = r['package_code']
        if set(k) == set(items):
            return k
    rng = secrets.SystemRandom()
    nums = list(range(1, len(items) + 1)); rng.shuffle(nums)
    k = {it: f'{prefix}{n:0{width}d}' for it, n in zip(sorted(items), nums)}
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write('study_code\tpackage_code\n')
        for a in sorted(k): f.write(f'{a}\t{k[a]}\n')
    return k


class Redactor:
    def __init__(self, part, tube, dev):
        self.part, self.tube, self.dev = part, tube, dev
        alt = '|'.join(sorted(map(re.escape, self.part), key=len, reverse=True))
        self.part_token = re.compile(r'(?<![A-Za-z0-9.])(' + alt + r')(?:\.0)?(?![A-Za-z0-9]|\.\d)')
        talt = '|'.join(sorted(map(re.escape, self.tube), key=len, reverse=True))
        self.tube_token = re.compile(r'(?<![A-Za-z0-9])(' + talt + r')(?![A-Za-z0-9])')
        dalt = '|'.join(sorted(map(re.escape, self.dev), key=len, reverse=True)) or 'NO_DEVICE_IDS'
        self.dev_token = re.compile(r'(?<![A-Za-z0-9])(' + dalt + r')(?![A-Za-z0-9])')
        self.log = []
        self.keep_email = False

    def paths(self, s):
        n = 0
        s, k = RAW_PATH.subn(lambda m: m.group(0).rstrip('/').split('/')[-1], s); n += k
        s, k = WIN_URI.subn('<raw-data>', s); n += k
        s, k = WIN_PATH.subn(lambda m: m.group(1) + m.group(2).rstrip('\\').split('\\')[-1], s); n += k
        for rx, rep in TEXT_SUBS:
            s, k = rx.subn(rep, s); n += k
        for rx, rep in PATH_RULES:
            s, k = rx.subn(rep, s); n += k
        s2 = clean_paths(s)
        if s2 != s: n += 1
        s = s2
        if not self.keep_email:
            s, k = EMAIL.subn('[e-mail address removed]', s); n += k
        return s, n

    def ids_in_string(self, s, whole_cell_ok):
        s, nd = self.dev_token.subn(lambda m: self.dev[m.group(1)], s)
        s2, nt = self.tube_token.subn(lambda m: self.tube[m.group(1)], s)
        nt += nd
        stripped = s2.strip()
        m = re.fullmatch(r'(\d{6})(\.0)?', stripped)
        if m and m.group(1) in self.part:
            if whole_cell_ok:
                return self.part[m.group(1)], nt, 1, 0
            return s2, nt, 0, 1
        if re.search(r'[A-Za-z_|:;/,\- ]', s2):
            s3, np_ = self.part_token.subn(lambda mm: self.part[mm.group(1)], s2)
            return s3, nt, np_, 0
        return s2, nt, 0, 0

    def code_of(self, raw):
        """Return (kind, original_code) if the raw cell is a study code, else None."""
        v = (raw or '').strip()
        if v.endswith('.0'): v = v[:-2]
        if v in self.part: return v
        if v in self.tube: return v
        return None


def _ascending_fraction(seq):
    k = [(0, int(x)) if x.isdigit() else (1, x) for x in seq]
    return sum(1 for i in range(len(k) - 1) if k[i] < k[i + 1]) / max(1, len(k) - 1)


def resort_rows(name, header, rows, R, order_log):
    """rows: list of lists of ORIGINAL cell strings.  If the first-appearance order of the study codes in
    a column is mostly ascending (>=0.75 of successive pairs; >=10 distinct codes), the coded rows are
    re-sorted by pseudonym (stable); non-coded rows keep their positions.  Participant columns first."""
    if POSITIONAL.search(name):
        return rows, False
    cols = sorted(range(len(header)), key=lambda c: (0 if ID_HEADER.search(header[c]) else 1, c))
    for ci in cols:
        coded = [(i, R.code_of(r[ci] if ci < len(r) else '')) for i, r in enumerate(rows)]
        coded = [(i, c) for i, c in coded if c]
        first, seen = [], set()
        for _, c in coded:
            if c not in seen: seen.add(c); first.append(c)
        if len(first) < 10 or _ascending_fraction(first) < 0.75:
            continue
        pseud = {c: (R.part.get(c) or R.tube.get(c)) for c in seen}
        idx = [i for i, _ in coded]
        new_order = sorted(coded, key=lambda t: (pseud[t[1]], t[0]))
        out = list(rows)
        for slot, (i, _) in zip(idx, new_order):
            out[slot] = rows[i]
        order_log.append([name, header[ci], 'resorted_by_pseudonym', len(idx)])
        return out, True
    return rows, False


def drop_cols(header):
    return [i for i, h in enumerate(header) if h.strip() in DROP_COLS]


def sha(b):
    return hashlib.sha256(b).hexdigest()


def read_member(sid, member):
    if sid in ('AUTH', 'ABS'):
        return open(member, 'rb').read()
    kind, p = S.SOURCES[sid]
    if kind == 'zip':
        return zipfile.ZipFile(p).read(member)
    with open(os.path.join(p, member), 'rb') as f:
        return f.read()


def origin_label(sid, member):
    if sid == 'AUTH':
        return 'authored for this package (round 2)'
    if sid == 'ABS':
        return S.abs_label(member) + ': ' + member.split('/')[-1]
    lab = S.PACKET_LABEL[sid]
    parts = member.split('/')
    b = parts[-1]
    m = re.match(r'([0-9a-f]{12})_(?:[0-9a-f]{12}_)?(.*)$', b)
    if m:
        b = m.group(1) + '_' + m.group(2).split('__')[-1]
    mem = '/'.join(parts[:-1] + [b])
    mem = clean_paths(mem) if INTERNAL_MARK.search(mem) else mem
    return f'{lab}: {mem}'


def text_table(name, s, R, order_log, stats):
    lines = s.split('\n')
    hdr_i = next((i for i, l in enumerate(lines) if l and not l.startswith('#')), None)
    if hdr_i is None:
        return s
    header = lines[hdr_i].split('\t')
    body_idx = [i for i in range(hdr_i + 1, len(lines)) if lines[i] != '' and not lines[i].startswith('#')]
    rows = [lines[i].split('\t') for i in body_idx]
    rows, moved = resort_rows(name, header, rows, R, order_log)
    dc = drop_cols(header)
    if dc:
        stats['dropped_cols'] += len(dc)
    keep = [i for i in range(len(header)) if i not in dc]
    out_lines = list(lines)
    out_lines[hdr_i] = '\t'.join(R.ids_in_string(header[i], False)[0] for i in keep)
    for slot, cells in zip(body_idx, rows):
        new = []
        for ci in keep:
            if ci >= len(cells):
                continue
            c = cells[ci]
            hdr = header[ci] if ci < len(header) else ''
            c2, nt, np_, amb = R.ids_in_string(c, bool(ID_HEADER.search(hdr)))
            stats['tubes'] += nt; stats['participants'] += np_; stats['ambiguous'] += amb
            if amb: R.log.append([name, 'ambiguous_numeric_cell_not_replaced', hdr, c])
            new.append(c2)
        # cells beyond the header width are kept
        new += [R.ids_in_string(c, False)[0] for c in cells[len(header):]]
        out_lines[slot] = '\t'.join(new)
    return '\n'.join(out_lines)


def transform(name, data, R, order_log):
    low = name.lower()
    tot = {'paths': 0, 'tubes': 0, 'participants': 0, 'ambiguous': 0, 'dropped_cols': 0}
    if low.endswith(('.tsv', '.csv', '.txt', '.md', '.json', '.html')):
        try:
            s = data.decode('utf-8')
        except UnicodeDecodeError:
            s = data.decode('latin-1')
        s1, tot['paths'] = R.paths(s)
        if low.endswith('.json') and JOBKEY.search(s1):
            obj = json.loads(s1)
            def scrub(o):
                if isinstance(o, dict):
                    return {k: scrub(v) for k, v in o.items() if not JOBKEY_NAME.match(k)}
                if isinstance(o, list):
                    return [scrub(v) for v in o]
                return o
            s1 = json.dumps(scrub(obj), indent=1, ensure_ascii=False) + '\n'
        if low.endswith('.tsv'):
            s2 = text_table(name, s1, R, order_log, tot)
        else:
            s2, nt = R.dev_token.subn(lambda m: R.dev[m.group(1)], s1)
            s2, nt2 = R.tube_token.subn(lambda m: R.tube[m.group(1)], s2)
            s2, np_ = R.part_token.subn(lambda m: R.part[m.group(1)], s2)
            tot['tubes'] += nt + nt2; tot['participants'] += np_
        if s2 == s:
            return data, 'byte-identical copy', tot
        return s2.encode('utf-8'), 'text rewritten (paths/codes/columns)', tot
    if low.endswith('.gz'):
        inner = gzip.decompress(data)
        new, note, st = transform(name[:-3], inner, R, order_log)
        for k in tot: tot[k] += st.get(k, 0)
        if new == inner:
            return data, 'byte-identical copy', tot
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode='wb', mtime=0, compresslevel=6, filename='') as g:
            g.write(new)
        return buf.getvalue(), 'gzip content rewritten (re-compressed, mtime 0)', tot
    if low.endswith('.xlsx'):
        import openpyxl
        zin = zipfile.ZipFile(io.BytesIO(data))
        has_comments = any(re.search(r'comments\d*\.xml$|threadedComment', n) for n in zin.namelist())
        core = zin.read('docProps/core.xml').decode('utf-8', 'replace') if 'docProps/core.xml' in zin.namelist() else ''
        wb = openpyxl.load_workbook(io.BytesIO(data))
        changed = 0
        for ws in wb.worksheets:
            rows = [list(r) for r in ws.iter_rows()]
            if not rows:
                continue
            header = [str(c.value) if c.value is not None else '' for c in rows[0]]
            # comments
            for r in rows:
                for c in r:
                    if c.comment is not None:
                        c.comment = None; changed += 1
            # re-sort rows sorted by study code
            vals = [[('' if c.value is None else str(c.value)) for c in r] for r in rows[1:]]
            new_vals, moved = resort_rows(name + '!' + ws.title, header, vals, R, order_log)
            if moved:
                raw = [[c.value for c in r] for r in rows[1:]]
                pos = {id(v): i for i, v in enumerate(vals)}
                perm = [pos[id(v)] for v in new_vals]
                for ri, src in enumerate(perm):
                    for ci, c in enumerate(rows[ri + 1]):
                        c.value = raw[src][ci] if ci < len(raw[src]) else None
                changed += 1
            for ri, r in enumerate(rows):
                for ci, c in enumerate(r):
                    v = c.value
                    hdr = header[ci] if ci < len(header) else ''
                    if isinstance(v, str):
                        v2, n = R.paths(v)
                        v3, nt, np_, amb = R.ids_in_string(v2, bool(ID_HEADER.search(hdr)) and ri > 0)
                        tot['paths'] += n; tot['tubes'] += nt; tot['participants'] += np_; tot['ambiguous'] += amb
                        if amb: R.log.append([name + '!' + ws.title, 'ambiguous_numeric_cell_not_replaced', hdr, v])
                        if v3 != v:
                            c.value = v3; changed += 1
                    elif ri > 0 and isinstance(v, (int, float)) and not isinstance(v, bool) and float(v).is_integer() and str(int(v)) in R.part:
                        if ID_HEADER.search(hdr):
                            c.value = R.part[str(int(v))]; changed += 1; tot['participants'] += 1
                        else:
                            tot['ambiguous'] += 1; R.log.append([name + '!' + ws.title, 'ambiguous_numeric_cell_not_replaced', hdr, v])
            dc = drop_cols(header)
            for ci in sorted(dc, reverse=True):
                ws.delete_cols(ci + 1); changed += 1; tot['dropped_cols'] += 1
        personal = re.search(r'<dc:creator>(?!openpyxl)[^<]+</dc:creator>|<cp:lastModifiedBy>[^<]+</cp:lastModifiedBy>', core)
        if not changed and not has_comments and not personal:
            return data, 'byte-identical copy', tot
        wb.properties.creator = 'FastaLake'; wb.properties.lastModifiedBy = 'FastaLake'
        out = io.BytesIO(); wb.save(out)
        what = []
        if changed: what.append(f'{changed} cell/column edits')
        if has_comments: what.append('cell comments removed')
        if personal: what.append('document properties neutralised')
        return out.getvalue(), 'xlsx rewritten (' + '; '.join(what) + '; re-saved with openpyxl)', tot
    if low.endswith('.npz'):
        import numpy as np
        z = np.load(io.BytesIO(data), allow_pickle=False)
        arrs = {}; changed = 0
        for k in z.files:
            a = z[k]
            k2, n1 = R.paths(k)
            k3 = R.part_token.sub(lambda m: R.part[m.group(1)], R.tube_token.sub(lambda m: R.tube[m.group(1)], k2))
            if a.dtype.kind in 'US':
                flat = a.astype(str).ravel().tolist()
                new = []
                for x in flat:
                    x2, n = R.paths(x); tot['paths'] += n
                    x3, nt, np_, amb = R.ids_in_string(x2, True)
                    tot['tubes'] += nt; tot['participants'] += np_
                    new.append(x3)
                if new != flat:
                    a = np.array(new).reshape(a.shape); changed += 1
            if k3 != k: changed += 1
            arrs[k3] = a
        if not changed:
            return data, 'byte-identical copy', tot
        compressed = any(i.compress_type == zipfile.ZIP_DEFLATED for i in zipfile.ZipFile(io.BytesIO(data)).infolist())
        out = io.BytesIO()
        (np.savez_compressed if compressed else np.savez)(out, **arrs)
        return out.getvalue(), 'npz string arrays/keys rewritten (re-saved with numpy)', tot
    if low.endswith('.pdf'):
        return data, 'byte-identical copy', tot
    raise SystemExit(f'unhandled type {name}')


def load_codes():
    data = read_member('Z16', 'current/Supplementary_Data_1_Sample_rosters.xlsx')
    import openpyxl
    wb = openpyxl.load_workbook(io.BytesIO(data), read_only=True)
    parts, tubes = set(), set()
    for ws in wb.worksheets:
        if not ws.title.startswith('PREDICT'): continue
        rows = ws.iter_rows(values_only=True); hdr = [str(h) for h in next(rows)]
        for r in rows:
            for h, v in zip(hdr, r):
                if v is None: continue
                v = str(v).strip()
                if v.endswith('.0'): v = v[:-2]
                if 'participant_id' in h and re.fullmatch(r'17[12]\d{3}', v): parts.add(v)
                if ('tube' in h or h == 'specimen_id') and re.fullmatch(r'17_[A-Z]_\d{3}_[A-Z0-9]{3}', v): tubes.add(v)
    return parts, tubes


def collect_capscan(rows):
    """Pre-pass: device IDs and collection dates from every selected table with those columns."""
    import openpyxl
    devs, dates = set(), set()
    for r in rows:
        low = r['member'].lower()
        if not low.endswith(('.tsv', '.tsv.gz', '.xlsx')):
            continue
        data = read_member(r['source'], r['member'])
        if low.endswith('.gz'):
            data = gzip.decompress(data)
        if low.endswith('.xlsx'):
            zf = zipfile.ZipFile(io.BytesIO(data))
            blob = b''.join(zf.read(n) for n in zf.namelist() if n.startswith('xl/') and n.endswith('.xml'))
            if b'device_id' not in blob and b'swallow_date' not in blob:
                continue
            wb = openpyxl.load_workbook(io.BytesIO(data), read_only=True)
            for ws in wb.worksheets:
                it = ws.iter_rows(values_only=True)
                try: hdr = [str(h) if h is not None else '' for h in next(it)]
                except StopIteration: continue
                for row in it:
                    for h, v in zip(hdr, row):
                        if v is None: continue
                        if h == 'device_id': devs.add(str(v).strip())
                        if h in ('swallow_date', 'recover_date'): dates.add(str(v).strip()[:10])
            continue
        if b'device_id' not in data[:20000] and b'swallow_date' not in data[:20000]:
            continue
        s = data.decode('utf-8', 'replace').split('\n')
        hdr = s[0].split('\t')
        for line in s[1:]:
            c = line.split('\t')
            for h, v in zip(hdr, c):
                if h == 'device_id' and v.strip() not in ('', 'nan', 'NA'): devs.add(v.strip())
                if h in ('swallow_date', 'recover_date') and v.strip(): dates.add(v.strip()[:10])
    devs = {d for d in devs if re.fullmatch(r'[0-9A-Fa-f]{6,12}', d)}
    return devs, dates


def main():
    sel_path, out_root, work, keydir = sys.argv[1:5]
    pk = os.path.join(out_root, PKG)
    if os.path.exists(pk) and os.listdir(pk):
        raise SystemExit(f'{pk} not empty; refusing to overwrite')
    os.makedirs(pk, exist_ok=True); os.makedirs(keydir, exist_ok=True)
    rows = list(csv.DictReader(open(sel_path), delimiter='\t'))
    for r in rows:
        if FORBIDDEN_NAME.search(r['package_path']):
            raise SystemExit(f'forbidden name selected: {r["package_path"]}')
    parts, tubes = load_codes()
    devs, dates = collect_capscan(rows)
    part = load_or_make_key(os.path.join(keydir, 'PREDICT_PARTICIPANT_KEY.tsv'), parts, 'PREDICT_P', 3)
    tube = load_or_make_key(os.path.join(keydir, 'PREDICT_TUBE_KEY.tsv'), tubes, 'PREDICT_T', 4)
    dev = load_or_make_key(os.path.join(keydir, 'CAPSCAN_DEVICE_KEY.tsv'), devs, 'capsule_', 3)
    with open(os.path.join(work, 'CAPSCAN_DATES_internal.txt'), 'w') as f:
        f.write('\n'.join(sorted(dates)) + '\n')
    R = Redactor(part, tube, dev)
    manifest = []; redlog = []; order_log = []
    for r in rows:
        rel = r['package_path']
        data = read_member(r['source'], r['member']); olabel = origin_label(r['source'], r['member'])
        pre, snote = SP.apply(rel, data)
        R.keep_email = rel.startswith('Reference_snapshots/')
        if rel == 'Annotation_verification/native_row_verification.tsv.gz':
            # ROUND3: 9.4 million rows of public UHGG gene identifiers, seed orthologues and alignment scores; no paths
            # or study codes can occur. Shipped unmodified (the gate scans its decompressed content).
            new, note, st = pre, 'byte-identical copy (large public-identifier table, scanned by the gate)', {'paths': 0, 'tubes': 0, 'participants': 0, 'ambiguous': 0, 'dropped_cols': 0}
        elif rel.startswith('Reference_snapshots/'):
            # third-party page captures are shipped unmodified
            new, note, st = pre, 'byte-identical copy (third-party page capture, unmodified)', {'paths': 0, 'tubes': 0, 'participants': 0, 'ambiguous': 0, 'dropped_cols': 0}
        else:
            new, note, st = transform(rel, pre, R, order_log)
        if snote:
            note = snote + ('' if note == 'byte-identical copy' else '; ' + note)
        dst = os.path.join(pk, rel); os.makedirs(os.path.dirname(dst), exist_ok=True)
        if os.path.exists(dst):
            raise SystemExit(f'duplicate destination {rel}')
        with open(dst, 'wb') as f: f.write(new)
        manifest.append([rel, len(new), sha(new), olabel, sha(data), note, r['use']])
        redlog.append([rel, st['paths'], st['tubes'], st['participants'], st['ambiguous'], st['dropped_cols'], note])
    with open(os.path.join(work, 'REDACTION_LOG.tsv'), 'w') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['package_path', 'path_edits', 'tube_or_device_codes', 'participant_codes', 'ambiguous_cells', 'dropped_columns', 'transformation'])
        w.writerows(redlog)
        for x in R.log: w.writerow(['#'] + [str(y)[:200] for y in x])
    with open(os.path.join(work, 'ORDER_LOG.tsv'), 'w') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(['file', 'column', 'action', 'rows']); w.writerows(order_log)
    with open(os.path.join(work, 'MANIFEST_build.tsv'), 'w') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['package_path', 'bytes', 'sha256', 'origin', 'origin_sha256', 'transformation', 'use'])
        w.writerows(manifest)
    print(json.dumps({'files': len(manifest), 'bytes': sum(m[1] for m in manifest),
                      'rewritten_files': sum(1 for m in manifest if m[5] != 'byte-identical copy'),
                      'participants_in_key': len(part), 'tubes_in_key': len(tube), 'devices_in_key': len(dev),
                      'capscan_dates': len(dates), 'resorted_tables': len(order_log)}))


if __name__ == '__main__':
    main()
