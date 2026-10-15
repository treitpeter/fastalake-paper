"""Figure 2 composite: the PDF edits that turn the producer canvas into the v30 figure.

Moved from the v30 producer ``build_fig2_len950.py`` (see PROVENANCE.md, ASSEMBLY.md), steps 1-6,
logic unchanged; only paths are parameters. No step replots data.

  1 panel letters to v28d order, legend tags of d, Söhne panel letters, trim blank bands
  2 remove the orange <10 frames, 'Members with / 10 or more peptides', n/8 counts and
    'Fewer than 10' from panel f
  3 panel g axis title, panel e legend keys, scale to 183 mm width
  4 'REF' -> 'Reference' overlays
  5 16 species labels italic (0.20 shear)
  6 remove four panel a strings

Fonts: the Söhne files in ``FASTALAKE_FONT_DIR`` (``soehne-buch.ttf``, ``soehne-halbfett.ttf``).
"""
import collections
import hashlib
import os
import re
from pathlib import Path

import fitz  # noqa: E402
import numpy as np  # noqa: E402
import pymupdf as f  # noqa: E402

import panels

#: Composite width (mm); NBT double column (ASSEMBLY.md step 3).
COMPOSITE_WIDTH_MM = 183
#: Raster resolution of the composite PNG.
COMPOSITE_DPI = 300

# v33 (Fig2_v31): the d and e colour keys sit in one row directly above their axes,
# left-aligned with the axes' left spine; entries spaced by measured text width.
PT = 72 / 25.4
AX_TOP = 112 * PT                   # top of the d and e axes (canvas ax(..., 112, ..., 27))
D_AX_X0, E_AX_X0 = 43 * PT, 114 * PT
KEY_GAP_PT = 4.5                    # key baseline above the axes top
KEY_ITEM_GAP_PT = 8.0               # space between one key entry's text and the next swatch
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()  # noqa: E731


def assemble(canvas, out_dir, font_dir=None, stem='Figure2'):
    """Apply steps 1-6 to the canvas and write ``<stem>.pdf`` and a 300 dpi ``<stem>.png``.

    Parameters
    ----------
    canvas : dict
        Output of :func:`panels.build_canvas` (``pdf`` path and inputs ``d``).
    out_dir : pathlib.Path
        Directory for the step PDFs and the composite (``build/Fig2``).
    font_dir : pathlib.Path, optional
        Söhne font directory; defaults to ``FASTALAKE_FONT_DIR``.
    stem : str
        Output file stem.

    Returns
    -------
    dict
        Step PDFs with sha256, the composite ``pdf`` and ``png`` paths and their sha256.
    """
    FONTDIR = Path(font_dir or os.environ['FASTALAKE_FONT_DIR'])
    OUT = Path(out_dir)
    d = canvas['d']
    log = {}
    def spans_of(page):
        return [s for b in page.get_text('dict')['blocks'] for l in b.get('lines', []) for s in l['spans']]


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
    _lf = f.Font(fontfile=str(FONTDIR / 'soehne-buch.ttf')); xx = D_AX_X0; yyk = AX_TOP - KEY_GAP_PT
    for lab, col in [('No exported prediction', (0, .447, .698)), ('Prediction exists', (.902, .624, 0))]:
        p.draw_rect(f.Rect(xx, yyk - 4, xx + 5.7, yyk), color=None, fill=col)
        p.insert_text((xx + 10, yyk), lab, fontname='legendreg', fontsize=5.7, color=(.306, .306, .306))
        xx += 10 + _lf.text_length(lab, fontsize=5.7) + KEY_ITEM_GAP_PT
    dd = labels_and_small_fixes(dd, 2)
    dd, cuts = trim_blank_bands(dd)
    s1 = OUT / 'step1_Figure2.pdf'; dd.save(s1, garbage=4, deflate=True)
    log['step1'] = dict(pdf=s1.name, sha256=sha(s1), blank_bands_removed=cuts, page=list(f.open(s1)[0].rect))

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
    log['step2'] = dict(pdf=s2.name, sha256=sha(s2), streams_edited=edited, frames_removed=n_frames_expected, text_removed=removed)

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
    _ef = fitz.Font(fontfile=BUCH); x = E_AX_X0; yk = AX_TOP - KEY_GAP_PT
    for txt, c in [('One member', (0, .447, .698)), ('Several members', (.5, .5, .5)), ('No member', (.85, .85, .85))]:
        p.draw_rect(fitz.Rect(x, yk - 4.4, x + 5, yk - .4), color=None, fill=c)
        p.insert_text((x + 9, yk), txt, fontname='SohneFix', fontsize=5.6, color=(.27, .27, .27))
        x += 9 + _ef.text_length(txt, fontsize=5.6) + KEY_ITEM_GAP_PT
    p.insert_text((75, 618), target, fontname='SohneFix', fontsize=6.3, rotate=90, color=(.27, .27, .27))
    nums = lambda pg: collections.Counter(re.findall(r'(?<![A-Za-z])\d[\d,.]*', pg.get_text()))  # noqa: E731
    assert nums(p) == nums(fitz.open(s2)[0])
    scaled = fitz.open(); width = COMPOSITE_WIDTH_MM / 25.4 * 72
    page = scaled.new_page(width=width, height=p.rect.height * width / p.rect.width); page.show_pdf_page(page.rect, doc, 0)
    s3 = OUT / 'step3_Figure2.pdf'; scaled.save(s3, garbage=4, deflate=True)
    assert min(s['size'] for s in spans_of(fitz.open(s3)[0]) if s['text'].strip()) >= 5
    log['step3'] = dict(pdf=s3.name, sha256=sha(s3))

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
    log['step4'] = dict(pdf=s4.name, sha256=sha(s4), labels=[(e['old'], e['new']) for e in edits])

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
    log['step5'] = dict(pdf=s5.name, sha256=sha(s5), italic_labels=len(ed5))

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

    final_pdf = OUT / f'{stem}.pdf'
    doc.save(final_pdf, garbage=4, deflate=True)
    final_png = final_pdf.with_suffix('.png')
    with fitz.open(final_pdf) as q:
        pix = q[0].get_pixmap(dpi=COMPOSITE_DPI, alpha=False); pix.set_dpi(COMPOSITE_DPI, COMPOSITE_DPI); pix.save(final_png)
    log['final'] = dict(pdf=final_pdf.name, pdf_sha256=sha(final_pdf), png=final_png.name, png_sha256=sha(final_png))
    return log
