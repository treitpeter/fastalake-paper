"""Rebuild Figure 2 with panels b, d, e, f, g at 9-50 residues (len950, 2026-10-06).

Usage (Slurm only):
    python build_fig2_len950.py control   # PLOT_MIN_LENGTH = 5; must reproduce expected/Fig2_v28f.png
    python build_fig2_len950.py len950    # PLOT_MIN_LENGTH = 9; writes Figure2_len950.pdf/.png

Step 0 is ``panels_len950.build_canvas`` (copy of figures/Fig2/panels.py; only the length
filter of b, d, e, f, g changed). Steps 1-6 repeat the PDF edits of
figures/Fig2/ASSEMBLY.md, ported from the original scripts named there; no step replots data:
  1 build_figures_r12.py, figure 2 branch + labels_and_small_fixes + trim_blank_bands
  2 rebuild_fig2f.py (orange <10 frames, 'Members with / 10 or more peptides', n/8 counts,
    'Fewer than 10' deleted from the content stream)
  3 fix_current_figures.py, figure 2 branch (g axis title, e legend keys, 180 -> 183 mm)
  4 build_FAR15_figures.py, figure 2 branch (REF -> Reference overlays)
  5 build_FAR16_figures.py, figure 2 branch (16 species labels italic, 0.20 shear)
  6 v26b/fig_clean: PyMuPDF text redaction of four panel a strings
Generalisations: the stream regexes of steps 2 and 3 also accept 'Tf 0 0 Td[<' (this
environment's matplotlib writes a Td after Tf); step 2 removes whatever number of <10 frames and n/8 counts the canvas has
(the original asserted the 5-50 counts: 6 frames, 8/8 x4 + 6/8).
"""
import collections
import hashlib
import json
import math
import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
LAB = Path(os.environ['FASTALAKE_DATA_DIR'])
REPO = LAB / 'release_prep/fastalake-paper_v28'
FONTDIR = Path(os.environ.get('FASTALAKE_FONT_DIR', ''))
os.environ.setdefault('FASTALAKE_FONT_DIR', str(FONTDIR))
sys.path[:0] = [str(HERE), str(REPO)]

import fitz  # noqa: E402
import numpy as np  # noqa: E402
import pymupdf as f  # noqa: E402

from style import fl_style as fl  # noqa: E402
import panels_len950 as panels  # noqa: E402

MODE = sys.argv[1]
assert MODE in ('control', 'len950')
panels.PLOT_MIN_LENGTH = 5 if MODE == 'control' else 9
OUT = HERE / 'work' / MODE
OUT.mkdir(parents=True, exist_ok=True)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()  # noqa: E731
log = {'mode': MODE, 'plot_min_length': panels.PLOT_MIN_LENGTH}


def spans_of(page):
    return [s for b in page.get_text('dict')['blocks'] for l in b.get('lines', []) for s in l['spans']]


# ---------------------------------------------------------------- step 0: canvas
fl.init()
canvas = panels.build_canvas(REPO / 'source_data', REPO / 'figures/Fig2/artwork/bottle_61.pdf', OUT)
d = canvas['d']
rc = panels.recount_membership(REPO / 'source_data', d)          # asserts tables == membership recount
log['step0'] = dict(pdf=str(canvas['pdf']), sha256=sha(canvas['pdf']), recount_checks=rc['checks'])
if MODE == 'control':
    ref = REPO / 'figures/Fig2/expected/CAMPI_Figure2_review.png'
    from PIL import Image, ImageChops
    diff = ImageChops.difference(Image.open(canvas['png']).convert('RGB'), Image.open(ref).convert('RGB')).getbbox()
    log['step0']['png_equal_expected_canvas'] = diff is None


# ---------------------------------------------------------------- step 1
def remove_text(p, strings):
    found = []
    for b in p.get_text('dict')['blocks']:
        for l in b.get('lines', []):
            for s in l['spans']:
                if s['text'] in strings:
                    p.add_redact_annot(f.Rect(s['bbox']), fill=False, cross_out=False); found.append(s)
    if found:
        p.apply_redactions(images=0, graphics=0)
    return found


def labels_and_small_fixes(doc, k):
    p = doc[0]; ss = spans_of(p)
    changes = [s for s in ss if re.fullmatch('[a-i]', s['text']) and s['size'] >= 5]
    for s in changes:
        p.add_redact_annot(f.Rect(s['bbox']), fill=False, cross_out=False)
    if changes:
        p.apply_redactions(images=0, graphics=0)
    for key, file in [('flbold', 'soehne-halbfett.ttf'), ('flreg', 'soehne-buch.ttf')]:
        p.insert_font(fontname=key, fontfile=str(FONTDIR / file))
    for s in changes:
        p.insert_text(list(s['origin']), s['text'], fontsize=8, fontname='flbold', color=(.15, .15, .15))
    return doc


def trim_blank_bands(doc, threshold=15, pad=5):
    p = doc[0]; pix = p.get_pixmap(dpi=72, alpha=False); a = np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, 3)
    white = np.all(a > 249, axis=(1, 2)); cuts = []; start = None
    for y, v in enumerate(list(white) + [False]):
        if v and start is None:
            start = y
        elif not v and start is not None:
            if y - start > threshold:
                cuts.append((start + pad, y - pad))
            start = None
    if not cuts:
        return doc, []
    out = f.open(); dest = out.new_page(width=p.rect.width, height=p.rect.height - sum(b - a for a, b in cuts)); srcy = 0; dsty = 0
    for a, b in cuts + [(p.rect.height, p.rect.height)]:
        if a > srcy:
            rect = f.Rect(0, srcy, p.rect.width, a); dest.show_pdf_page(f.Rect(0, dsty, p.rect.width, dsty + rect.height), doc, 0, clip=rect); dsty += rect.height
        srcy = b
    return out, cuts


src = canvas['pdf']
dd = f.open(src); p = dd[0]
found = remove_text(p, {'d–f', 'g', 'h', 'i'})
for s in found:
    p.insert_text(s['origin'], {'d–f': 'd', 'g': 'f', 'h': 'g', 'i': 'h'}[s['text']], fontname='hebo', fontsize=8)
p.insert_text((260, 295), 'e', fontname='hebo', fontsize=8)
tags = remove_text(p, {'No exported prediction', 'Prediction exists'})
assert len(tags) == 2
yy = tags[0]['origin'][1]
p.draw_rect(f.Rect(44, yy - 7, 192, yy + 2), color=None, fill=(1, 1, 1))
p.insert_font(fontname='legendreg', fontfile=str(FONTDIR / 'soehne-buch.ttf'))
for xx, lab, col in [(47.6, 'No exported prediction', (0, .447, .698)), (135, 'Prediction exists', (.902, .624, 0))]:
    p.draw_rect(f.Rect(xx, yy - 4, xx + 5.7, yy), color=None, fill=col)
    p.insert_text((xx + 10, yy), lab, fontname='legendreg', fontsize=5.7, color=(.306, .306, .306))
dd = labels_and_small_fixes(dd, 2)
dd, cuts = trim_blank_bands(dd)
s1 = OUT / 'step1_Figure2.pdf'; dd.save(s1, garbage=4, deflate=True)
log['step1'] = dict(pdf=str(s1), sha256=sha(s1), blank_bands_removed=cuts, page=list(f.open(s1)[0].rect))

# ---------------------------------------------------------------- step 2
d2 = fitz.open(s1); p = d2[0]
ORANGE = (0.902, 0.624, 0)
HEAT = fitz.Rect(119.1, 444.8, 414.0, 518.7)


def spans2(page):
    return [(s['text'], fitz.Rect(s['bbox'])) for s in spans_of(page)]


mt = d['mt']
n_frames_expected = 2 * (int((mt.to_numpy() < 10).sum()) + 1)      # cells < 10 plus the key box; PyMuPDF reports 2 drawings per frame (5-50: 6, as asserted by rebuild_fig2f.py)
counts_text = [f'{int((mt[a] >= 10).sum())}/8' for a in panels.ARMS]
before = spans2(p); before_numbers = sorted(t for t, r in before if r.intersects(HEAT) and re.fullmatch(r'[\d,]+', t))
assert len(before_numbers) == 40, len(before_numbers)
n_orange = lambda pg: len([x for x in pg.get_drawings() if x.get('color') and all(abs(a - b) < 0.02 for a, b in zip(x['color'], ORANGE)) and x.get('fill') is None])  # noqa: E731
assert n_orange(p) == n_frames_expected, (n_orange(p), n_frames_expected)
pat = re.compile(r'q 1 G [\d. ]+ re W n \.9019608 \.62352946 0 RG 1\.2 w [\d.]+ [\d.]+ m(?: [\d.]+ [\d.]+ l)+ h S Q ?')
edited = {}
for x in range(1, d2.xref_length()):
    if not d2.xref_is_stream(x):
        continue
    s = d2.xref_stream(x)
    if not s or b'.9019608 .62352946 0 RG 1.2 w' not in s:
        continue
    t = s.decode('latin1'); t2, n = pat.subn('', t); assert n > 0; edited[x] = n
    d2.update_stream(x, t2.encode('latin1'))
assert edited
tmp1 = OUT / '_s2a.pdf'; d2.save(tmp1, garbage=0); d2 = fitz.open(tmp1); p = d2[0]
assert n_orange(p) == 0


def hexpat(txt):
    return r'(?:>-?[\d.]+<)?'.join('%04X' % ord(c) for c in txt)


TXT = ['Members with', '10 or more peptides', 'Fewer than 10'] + counts_text
removed = []
for x in list(edited):
    t = d2.xref_stream(x).decode('latin1')
    for txt in sorted(set(TXT)):
        pat2 = re.compile(r'q 1 0 0 1 [\d.-]+ [\d.-]+ cm [^Q]*?BT/F1 [\d.]+ Tf(?: 0 0 Td)?\[<' + hexpat(txt) + r'>\]TJ ET Q ?')
        t, n = pat2.subn('', t); removed += [txt] * n
    d2.update_stream(x, t.encode('latin1'))
assert sorted(removed) == sorted(TXT), removed
tmp2 = OUT / '_s2b.pdf'; d2.save(tmp2, garbage=0); d2 = fitz.open(tmp2); p = d2[0]
after = spans2(p); after_numbers = sorted(t for t, r in after if r.intersects(HEAT) and re.fullmatch(r'[\d,]+', t))
assert after_numbers == before_numbers
assert collections.Counter(t for t, _ in before) - collections.Counter(t for t, _ in after) == collections.Counter(TXT)
assert not (collections.Counter(t for t, _ in after) - collections.Counter(t for t, _ in before))
assert n_orange(p) == 0 and len(p.get_image_info()) == 2
nb, na = len(fitz.open(s1)[0].get_drawings()), len(p.get_drawings())
assert nb - na == n_frames_expected, (nb, na)
s2 = OUT / 'step2_Figure2.pdf'; d2.save(s2, garbage=4, deflate=True); tmp1.unlink(); tmp2.unlink()
log['step2'] = dict(pdf=str(s2), sha256=sha(s2), streams_edited=edited, frames_removed=n_frames_expected, text_removed=removed)

# ---------------------------------------------------------------- step 3
BUCH = str(FONTDIR / 'soehne-buch.ttf')
doc = fitz.open(s2); p = doc[0]
target = 'FL peptides (thousands)'
hp = r'(?:>-?[\d.]+<)?'.join('%04X' % ord(c) for c in target)
pat = re.compile(r'q [\d. \-]+ cm [^Q]*?BT/F\d+ [\d.]+ Tf(?: 0 0 Td)?\[<' + hp + r'>\]TJ ET Q ?')
nrem = 0
for x in range(1, doc.xref_length()):
    if not doc.xref_is_stream(x):
        continue
    data = doc.xref_stream(x)
    if not data:
        continue
    text = data.decode('latin1'); new, n = pat.subn('', text)
    if n:
        doc.update_stream(x, new.encode('latin1')); nrem += n
assert nrem == 1
legend = fitz.Rect(273, 284.5, 447, 294)
p.add_redact_annot(legend, fill=(1, 1, 1))
p.apply_redactions(images=0, graphics=0)
p.draw_rect(legend, color=None, fill=(1, 1, 1))
p.insert_font(fontname='SohneFix', fontfile=BUCH)
for x, txt, c in [(275, 'One member', (0, .447, .698)), (335, 'Several members', (.5, .5, .5)), (410, 'No member', (.85, .85, .85))]:
    p.draw_rect(fitz.Rect(x, 287, x + 5, 291), color=None, fill=c)
    p.insert_text((x + 9, 291.4), txt, fontname='SohneFix', fontsize=5.6, color=(.27, .27, .27))
p.insert_text((75, 618), target, fontname='SohneFix', fontsize=6.3, rotate=90, color=(.27, .27, .27))
nums = lambda pg: collections.Counter(re.findall(r'(?<![A-Za-z])\d[\d,.]*', pg.get_text()))  # noqa: E731
assert nums(p) == nums(fitz.open(s2)[0])
scaled = fitz.open(); width = 183 / 25.4 * 72
page = scaled.new_page(width=width, height=p.rect.height * width / p.rect.width); page.show_pdf_page(page.rect, doc, 0)
s3 = OUT / 'step3_Figure2.pdf'; scaled.save(s3, garbage=4, deflate=True)
assert min(s['size'] for s in spans_of(fitz.open(s3)[0]) if s['text'].strip()) >= 5
log['step3'] = dict(pdf=str(s3), sha256=sha(s3))

# ---------------------------------------------------------------- step 4 (FAR15, figure 2 branch)
fitz.TOOLS.set_small_glyph_heights(True)
doc = fitz.open(s3); page = doc[0]
spans = []
for b in page.get_text('dict')['blocks']:
    for line in b.get('lines', []):
        for sp in line['spans']:
            spans.append(dict(sp, direction=line['dir']))
xref_fix = [e[0] for e in page.get_fonts(full=True) if e[4] == 'SohneFix']
assert len(set(xref_fix)) == 1           # the original read xref 39, which is this font
fontbuf = doc.extract_font(xref_fix[0])[3]
fonts = {'regular': fitz.Font(fontbuffer=fontbuf)}
page.insert_font(fontname='RenamedRegular', fontbuffer=fontbuf)
edits = []
for sp in spans:
    old = sp['text']; box = fitz.Rect(sp['bbox']); cx = (box.x0 + box.x1) / 2
    if not re.search(r'\bREF\b', old):
        continue
    new = re.sub(r'\bREF\b', 'Reference', old)
    new = new.replace('FL+Reference+MG', 'FL + Reference + MG').replace('FL+Reference', 'FL + Reference')
    size = sp['size']; align = 'center'; x = cx; y = sp['origin'][1]; rotate = 0
    kind = 'italic' if 'Oblique' in sp['font'] else 'regular'
    if box.y0 < 145 and old not in ['REF proteomes']:
        align = 'left'; x = box.x0
    if 252 < box.y0 < 261:
        if new == 'FL + Reference':
            new = 'FL +\nReference'
        if new == 'FL + Reference + MG':
            new = 'FL +\nReference\n+ MG'
    if 340 < box.y0 < 400:
        align = 'right'; x = box.x1
    if old == 'REF recovered':
        x = 490
    if old.startswith('Unrecovered REF'):
        new = old.replace('REF', 'reference')
    if 'no community composition' in old:
        new = old.replace('CAMPI REF', 'CAMPI reference proteomes')
    if old.startswith('REF and MG were'):
        new = old.replace('REF', 'Reference proteomes')
    edits.append(dict(old=old, new=new, old_bbox=list(box), size=size, align=align, x=x, y=y, rotate=rotate, kind=kind, color=sp['color']))
    page.draw_rect(box + (-.3, -.3, .3, .3), color=None, fill=(1, 1, 1), overlay=True)
for kind, font in fonts.items():
    page.insert_font(fontname='Renamed' + kind.title(), fontbuffer=font.buffer)
for e in edits:
    font = fonts[e['kind']]; color = tuple(((e['color'] >> shift) & 255) / 255 for shift in [16, 8, 0])
    for n, line in enumerate(e['new'].split('\n')):
        if not line:
            continue
        w = font.text_length(line, fontsize=e['size']); x = e['x']; y = e['y'] + n * e['size'] * 1.15
        if e['align'] == 'center':
            x -= w / 2
        elif e['align'] == 'right':
            x -= w
        page.insert_text((x, y), line, fontname='Renamed' + e['kind'].title(), fontsize=e['size'], color=color, rotate=e['rotate'])
s4 = OUT / 'step4_Figure2.pdf'; doc.save(s4, garbage=4, deflate=True)
log['step4'] = dict(pdf=str(s4), sha256=sha(s4), labels=[(e['old'], e['new']) for e in edits])

# ---------------------------------------------------------------- step 5 (FAR16, figure 2 branch)
doc = fitz.open(s4); page = doc[0]
spans = spans_of(page)
taxa = re.compile(r'^[A-Z]\.\s+[a-z]+$')
font = fitz.Font(fontfile=str(FONTDIR / 'soehne-buch.ttf'))
page.insert_font(fontname='SpeciesItalic', fontbuffer=font.buffer)
ed5 = []
for sp in spans:
    label = sp['text']
    if not taxa.fullmatch(label):
        continue
    box = fitz.Rect(sp['bbox'])
    if box.y0 < 150:
        continue
    size = sp['size']; w = font.text_length(label, fontsize=size)
    x = box.x1 - w; y = sp['origin'][1]
    page.draw_rect(box + (-.35, -.25, .35, .25), color=None, fill=(1, 1, 1), overlay=True)
    point = fitz.Point(x, y)
    color = tuple(((sp['color'] >> shift) & 255) / 255 for shift in [16, 8, 0])
    page.insert_text(point, label, fontname='SpeciesItalic', fontsize=size, color=color, morph=(point, fitz.Matrix(1, 0, .20, 1, 0, 0)))
    ed5.append(label)
assert len(ed5) == 16, len(ed5)
s5 = OUT / 'step5_Figure2.pdf'; doc.save(s5, garbage=4, deflate=True)
log['step5'] = dict(pdf=str(s5), sha256=sha(s5), italic_labels=len(ed5))

# ---------------------------------------------------------------- step 6 (v26b fig_clean: four panel a strings)
REMOVE = ['8-member gut model · 15 acquisitions',
          'FL: selection from 375 M gut-catalogue and UniProtKB proteins;',
          'no community composition, CAMPI reference proteomes or specimen MG supplied.',
          'Reference proteomes and MG were added to the selected FL database after parsimony.']
doc = fitz.open(s5); page = doc[0]
hits = [s for s in spans_of(page) if s['text'].strip() in REMOVE]
assert {s['text'].strip() for s in hits} == set(REMOVE), [s['text'] for s in hits]
for s in hits:
    page.add_redact_annot(fitz.Rect(s['bbox']), fill=False, cross_out=False)
page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE)
left = page.get_text()
assert not any(r in left for r in REMOVE)
stem = 'Figure2_control' if MODE == 'control' else 'Figure2_len950'
final_pdf = (OUT if MODE == 'control' else HERE) / f'{stem}.pdf'
doc.save(final_pdf, garbage=4, deflate=True)
final_png = final_pdf.with_suffix('.png')
with fitz.open(final_pdf) as q:
    pix = q[0].get_pixmap(dpi=300, alpha=False); pix.set_dpi(300, 300); pix.save(final_png)
    log['final'] = dict(pdf=str(final_pdf), pdf_sha256=sha(final_pdf), png=str(final_png), png_sha256=sha(final_png),
                        png_px=[pix.width, pix.height], page_pt=list(q[0].rect)[2:])
if MODE == 'control':
    from PIL import Image
    a = np.asarray(Image.open(final_png).convert('RGB')).astype(int)
    b = np.asarray(Image.open(REPO / 'figures/Fig2/expected/Fig2_v28f.png').convert('RGB')).astype(int)
    log['control_vs_v28f'] = dict(shape_equal=a.shape == b.shape)
    if a.shape == b.shape:
        dif = np.any(a != b, axis=2)
        ys, xs = np.nonzero(dif)
        log['control_vs_v28f'].update(n_diff_px=int(dif.sum()), max_abs=int(np.abs(a - b).max()),
                                      bbox=[int(xs.min()), int(ys.min()), int(xs.max()), int(ys.max())] if dif.any() else None)
(OUT / 'BUILD_LOG.json').write_text(json.dumps(log, indent=2, default=str) + '\n')
print(json.dumps({k: v for k, v in log.items() if k in ('mode', 'final', 'control_vs_v28f')}, indent=1, default=str))
print('step0 equal expected canvas:', log['step0'].get('png_equal_expected_canvas'))
