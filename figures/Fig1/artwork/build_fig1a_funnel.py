"""Figure 1 panel a (v30): rebuild the FastaLake column of the composite as one funnel.

Repository port of the v30 producer ``build_fig1a_funnel.py`` (LAB_ROOT path in
PROVENANCE.md). The drawing logic below the path block is unchanged; only the input and
output paths were made repository-relative and the font directory is read from
``FASTALAKE_FONT_DIR``.

Inputs (all in ``figures/Fig1/artwork/``)
    Figure1_v28f.pdf               round-8 composite, the base page (read-only)
    Figure1_v27_art_source.pdf     text-free source of the clipped icon and peptide bars
    Figure1_v29_intermediate.pdf   v29 build, used only for the outside-column assertion
Outputs (``FASTALAKE_OUT_DIR``, default ``build/Fig1/funnel``)
    Figure1.pdf, Figure1.png (300 dpi), _build_receipt.json

Usage: ``FASTALAKE_FONT_DIR=<fonts> python figures/Fig1/artwork/build_fig1a_funnel.py``
(about one minute, 8 GB). The committed ``artwork/Figure1_v30.pdf`` is the producer's output;
the notebook rebuilds it with this script and compares the two at 300 dpi.

Original producer docstring follows.

Figure 1 panel a, FastaLake column (right third) rebuilt as one funnel. Base: fig1a_perfecter/Figure1.pdf (round 8, read-only).

Approved conception (Peter, 6 Oct): MS/MS -> De novo; De novo -> AlphaNovo peptide predictions (down) and -> scrambled
peptides (right); one straight arrow from each box into the funnel rim. Funnel: top layer = gut sequence reservoir, then
Exact match, Candidates this acquisition, Parsimony (glyph). Below: two straight arrows through one Database search box to
the two result boxes (kept exactly as in round 8). Real/Scrambled table re-laid so its rows sit on the funnel bands.
No step badges, no 'DDA', no bent arrows, no lookup box.

What it does
  1. Shows the round-8 page unchanged, then paints the column background (row-wise gradient sampled from the page at
     x = 236 pt) over the FastaLake column from y 23 pt (below the title) to y 307.6 pt (top of the result boxes). The title
     (already centred on the funnel axis, 350.5 pt) and the two result boxes are left as they are.
  2. Clips three pieces of artwork from a text-free copy of the round-8 page: the MS/MS spectrum icon (on the step-1 box
     fill) and the peptide bars of the two query boxes (on their box fills). Boxes, funnel, glyph, table grid and arrows are
     drawn as vectors; every label is typed in Soehne (real text).
  3. Asserts 0 changed pixels at 300 dpi outside the FastaLake column (x < 230 pt or y >= 372 pt) against the round-8 PDF.
Run on Slurm (producer: srun -c2 --mem=8G -t 30, interpreter python build_fig1a_funnel.py).

v29 (6 Oct 2026, Peter approved): band 3 label "Candidates, this acquisition" -> "Per-acquisition match" (bands now read
Reservoir [label unchanged] / Exact match / Per-acquisition match / Parsimony; numbers only in the table); one orange #d94d00
for all orange text, arrows and table values (round 8 had #e85d0f for "0 peptides" and #ff4d00 for the table).

Fig1_v2 (6 Oct 2026, Peter's requested changes, FastaLake column only):
  1. AlphaNovo and scrambled boxes side by side on one level above the rim, both shrunk to QW pt; one straight orange
     arrow AlphaNovo -> scrambled (replaces De novo -> scrambled).
  2. MS/MS and De novo centred on the column axis X (title / funnel axis); straight vertical MS/MS -> De novo -> AlphaNovo.
     The rim (116 pt) is narrower than the two boxes, so the pair cannot be centred on X with the De novo arrow landing on
     AlphaNovo. Solution: both rim arrows leave their box at the same offset E from the box's left edge, which makes the
     rim entries symmetric (X -/+ (QW + G)/2); E is set so the De novo arrow lands in the middle half of the AlphaNovo top.
  3. One vertical arrow from each box into the rim (navy, orange).
  4. Band 2 label "Exact match ·" / "study-supported proteins" (two lines, same style).
  5. Table Scrambled / Study-supported cell "–" -> "5,841" (same face, size and x as the 5,841 below it); the scrambled
     result-box spectrum line art recoloured #e85d0f -> #d94d00 (covered with the box fill and redrawn at the same geometry).
Funnel, table x, Database search and result boxes keep the v29 geometry (the rim centre CY is pinned to its v29 value)."""
import hashlib, json, math, os
from pathlib import Path
import fitz
import numpy as np

ART = Path(__file__).resolve().parent                # figures/Fig1/artwork
REPO = ART.parents[2]
OUT = Path(os.environ.get('FASTALAKE_OUT_DIR', REPO / 'build' / 'Fig1' / 'funnel'))
OUT.mkdir(parents=True, exist_ok=True)
BASE = ART / 'Figure1_v28f.pdf'                      # = producer input final/fig1a_perfecter/Figure1.pdf (same sha256)
FD = Path(os.environ['FASTALAKE_FONT_DIR'])          # Soehne fonts (producer: the lab style-guide directory)
FONTS = {'H': FD / 'soehne-halbfett.ttf', 'B': FD / 'soehne-buch.ttf'}
FNAME = {'H': 'FLfHalb', 'B': 'FLfBuch'}
CAP = 0.72
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
hexc = lambda h: tuple(int(h[i:i + 2], 16) / 255 for i in (1, 3, 5))
rgb = lambda t: tuple(v / 255 for v in t)
NAVY, ORANGE = hexc('#0a2d5e'), hexc('#d94d00')      # authentic arm / scrambled arm (as briefed)
LW, HL, HW = 1.0131829, 4.37, 2.525                  # the column-1/2 arrow style (shaft, head length, head half-width)
TXT = NAVY

src = fitz.open(BASE); sp = src[0]
W, H = sp.rect.width, sp.rect.height
font = {k: fitz.Font(fontfile=str(v)) for k, v in FONTS.items()}
tw = lambda t, f, s: font[f].text_length(t, s)

def spans(page):
    for b in page.get_text('dict')['blocks']:
        for l in b.get('lines', []):
            for s in l['spans']:
                if s['text'].strip(): yield s

# ---------------- base renders, background, text-free art copy ----------------
K = 300 / 72
bpm = sp.get_pixmap(dpi=300, alpha=False)
def px(x, y): return bpm.pixel(int(x * K), int(y * K))
ROWS = {j: bpm.pixel(int(236 * K), j) for j in range(int(20 * K), int(312 * K))}
def bg(y): return rgb(ROWS[min(max(int(y * K), int(20 * K)), int(312 * K) - 1)])
BAND1_FILL, BAND2_FILL = rgb(px(256, 76)), rgb(px(256, 105))      # round-8 step-1 / step-2 box fills
# Artwork pieces are clipped from fig_v27/Figure1.pdf (the shallow page round 8 clipped them from). Clipping from the round-8
# page itself is impractical: its resources are nested five overlay passes deep and PyMuPDF re-scans them on every placement.
ART_SRC = ART / 'Figure1_v27_art_source.pdf'   # = producer input final/fig_v27/Figure1.pdf (same sha256)
art = fitz.open(); art.insert_pdf(fitz.open(ART_SRC)); ap = art[0]
_apm = ap.get_pixmap(dpi=300, alpha=False)
NAVYBOX_FILL = rgb(_apm.pixel(int(256.5 * K), int(110.0 * K)))         # v27 AlphaNovo box fill (= round-8 navy box)
SCRBOX_FILL = rgb(_apm.pixel(int(416.4 * K), int(110.0 * K)))          # v27 scrambled box fill
ap.add_redact_annot(ap.rect, fill=False)
ap.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE, text=fitz.PDF_REDACT_TEXT_REMOVE)
assert not list(spans(ap))
def flat_edges(r):          # a clip must sit on flat fill: its 1-px border is one colour
    pm = ap.get_pixmap(dpi=300, clip=r, alpha=False)
    a = np.frombuffer(pm.samples, np.uint8).reshape(pm.height, pm.width, 3).astype(int)
    e = np.concatenate([a[0], a[-1], a[:, 0], a[:, -1]])
    return int(np.abs(e - e[0]).max())
ICON = fitz.Rect(361.9, 167.3, 398.6, 186.3)        # v27 MS/MS spectrum on funnel band 1 (the clip round 8 used)
BARS_N = fitz.Rect(256.5, 110.2, 323.3, 116.8)      # v27 AlphaNovo-box peptide bars (ink 257.28-322.56 x 110.76-116.22)
BARS_S = fitz.Rect(416.4, 110.2, 483.2, 116.8)      # v27 scrambled-box peptide bars (ink 417.16-482.44)
for nm, r in (('icon', ICON), ('navy bars', BARS_N), ('scrambled bars', BARS_S)):
    assert flat_edges(r) <= 2, (nm, flat_edges(r))
# the clip background must equal the fill the new box is drawn with
for r, f in ((BARS_N, NAVYBOX_FILL), (BARS_S, SCRBOX_FILL)):
    assert max(abs(a - b) for a, b in zip(_apm.pixel(int(r.x0 * K) + 1, int(r.y0 * K) + 1), [round(v * 255) for v in f])) <= 2, (r, f)
def ink(r):                 # ink bbox of a clip (pt), at 1200 dpi
    pm = ap.get_pixmap(dpi=1200, clip=r, alpha=False)
    a = np.frombuffer(pm.samples, np.uint8).reshape(pm.height, pm.width, 3).astype(int)
    m = (np.abs(a - a[0, 0]).sum(2) > 60); ys, xs = np.nonzero(m); s = 72 / 1200
    return fitz.Rect(r.x0 + xs.min() * s, r.y0 + ys.min() * s, r.x0 + (xs.max() + 1) * s, r.y0 + (ys.max() + 1) * s)
ICON_INK, BN_INK, BS_INK = ink(ICON), ink(BARS_N), ink(BARS_S)

# ---------------- layout (pt) ----------------
X = 350.5                     # funnel axis = title centre = Database search / result-box symmetry axis (round 8)
G = 12.365                    # every vertical arrow gap (= MS/MS -> De novo gap set by the scrambled-box clearance)
R = 3.04                      # box corner radius (round-8 query and result boxes)
# MS/MS box (step-1 fill), De novo box (step-2 fill)
MS_LBL = ('MS/MS', 'H', 5.741); DN_LBL = ('De novo', 'H', 5.741)
PAD = 6.0
MS_W = PAD + tw(*MS_LBL) + 5.0 + ICON_INK.width + PAD
MSBOX = fitz.Rect(X - MS_W / 2, 26.0, X + MS_W / 2, 47.0)               # v2: centred on the axis
DN_W = tw(*DN_LBL) + 2 * 9.0
DNBOX = fitz.Rect(X - DN_W / 2, 59.36, X + DN_W / 2, 76.36)             # v2: centred on the axis (v29 y)
DN_CY = (DNBOX.y0 + DNBOX.y1) / 2
assert abs(DNBOX.y0 - MSBOX.y1 - G) < 0.01
QW, QH = 70.0, 35.73                                 # v2: query boxes shrunk (were 87.13 wide); v29 height
E = 10.0                                             # rim arrows leave each box E pt right of its left edge
D = (QW + G) / 2                                     # -> rim entries at X -/+ D (symmetric)
QL = X - D - E                                       # AlphaNovo box left edge
ALPHB = fitz.Rect(QL, DNBOX.y1 + G, QL + QW, DNBOX.y1 + G + QH)
SCRB = fitz.Rect(ALPHB.x1 + G, ALPHB.y0, ALPHB.x1 + G + QW, ALPHB.y1)
XN, XO = ALPHB.x0 + E, SCRB.x0 + E                   # navy / orange rim arrow x
assert abs(XN - (X - D)) < 1e-9 and abs(XO - (X + D)) < 1e-9
assert ALPHB.x0 + QW / 4 <= X <= ALPHB.x1 - QW / 4, 'De novo arrow outside the middle half of the AlphaNovo top'
# funnel
RX, RY, RIM_SW = 58.0, 8.3, 0.8
rim_top = lambda x, cy: cy - RY * math.sqrt(1 - ((x - X) / RX) ** 2) - RIM_SW / 2
CY = ALPHB.y1 + G + (RY * math.sqrt(1 - (42.0 / RX) ** 2) + RIM_SW / 2)   # rim centre pinned to v29 (arms then at X -/+ 42)
assert abs(CY - 142.95) < 0.01, CY
WB = 30.0                                            # half-width at the funnel bottom (wall end)
RESULT_TOP = 307.72                                  # round-8 result boxes (kept)
DBH = 16.21
DBBOX = fitz.Rect(X - 66.0, RESULT_TOP - G - DBH, X + 66.0, RESULT_TOP - G)
XB = 6.0                                             # the arms leave the funnel bottom at X -/+ XB
RB = RY * WB / RX
YB = DBBOX.y0 - G - RB * math.sqrt(1 - (XB / WB) ** 2)
half = lambda y: RX + (WB - RX) * (y - CY) / (YB - CY)
ry_of = lambda w: RY * w / RX
BAND_H0 = 35.0
YK = [CY, CY + BAND_H0]
for i in range(2): YK.append(YK[-1] + (YB - CY - BAND_H0) / 3)
YK.append(YB)
COL = [rgb((176, 226, 246)), rgb((177, 213, 232)), rgb((155, 201, 225)), rgb((130, 194, 200))]  # reservoir, old bands 3,4,5
CAP_COL = rgb((129, 179, 174))                                                                     # old bowl bottom
RIM_FILL, RIM_STROKE = rgb((220, 243, 255)), hexc('#1fa7df')                                       # reservoir cylinder top
# visible band spans on the axis (arc bottoms): band k from top[k] to top[k+1]
VTOP = [CY + RY] + [YK[k] + ry_of(half(YK[k])) for k in (1, 2, 3)] + [YB + RB]

# ---------------- compose ----------------
out = fitz.open(); out.insert_pdf(src)              # base page, unchanged
# old column text under the cover is removed (text-only redaction; line art and images untouched), so the text layer
# holds only what is visible. The title (y < 23) and the result-box text (y > 307.6) stay.
OLD_COL_TXT = [s_ for s_ in spans(sp) if s_['bbox'][0] > 230 and s_['bbox'][1] > 23 and s_['bbox'][3] < 307.6]
for s_ in OLD_COL_TXT:
    r_ = fitz.Rect(s_['bbox']); out[0].add_redact_annot(fitz.Rect(r_.x0 + .2, r_.y0 + 1.0, r_.x1 - .2, r_.y1 - 1.0), fill=False)
# v29: the round-8 '0 peptides' result-box title is #e85d0f; it is removed here and retyped in #d94d00 at the same origin
ZERO = [s_ for s_ in spans(sp) if s_['text'].strip() == '0 peptides' and s_['bbox'][1] > 307]
assert len(ZERO) == 1 and '%06x' % ZERO[0]['color'] == 'e85d0f', ZERO
ZERO = ZERO[0]
_zr = fitz.Rect(ZERO['bbox']); out[0].add_redact_annot(fitz.Rect(_zr.x0 + .2, _zr.y0 + 1.0, _zr.x1 - .2, _zr.y1 - 1.0), fill=False)
out[0].apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE, text=fitz.PDF_REDACT_TEXT_REMOVE)
ov = fitz.open(); pg = ov.new_page(width=W, height=H)   # overlay: everything new is drawn here
WRITERS = {}                                          # one TextWriter per colour (page.insert_text re-scans the deep base resources per call)

def cover(r):
    sh = pg.new_shape(); y = r.y0
    while y < r.y1:
        y2 = min(y + 0.5, r.y1)
        sh.draw_rect(fitz.Rect(r.x0, y, r.x1, min(y2 + 0.6, r.y1))); sh.finish(color=None, fill=bg((y + y2) / 2), width=0)
        y = y2
    sh.commit()
COVER = fitz.Rect(236.0, 23.0, 510.5, 307.6)
cover(COVER)

TEXT = []                                             # (text, rect) of every typed label, for overlap checks
def put(t, f, s, x, y, col=TXT):
    WRITERS.setdefault(tuple(col), fitz.TextWriter(pg.rect)).append((x, y), t, font=font[f], fontsize=s)
    TEXT.append((t, fitz.Rect(x, y - CAP * s, x + tw(t, f, s), y)))
def put_c(t, f, s, cx, y, col=TXT): put(t, f, s, cx - tw(t, f, s) / 2, y, col)
def place(clip, dx, dy): pg.show_pdf_page(clip + (dx, dy, dx, dy), art, 0, clip=clip)
def box(r, fill, stroke=None, sw=0): pg.draw_rect(r, color=stroke, fill=fill, width=sw, radius=R / min(r.width, r.height))

put(ZERO['text'], 'H' if 'Halbfett' in ZERO['font'] else 'B', ZERO['size'], ZERO['origin'][0], ZERO['origin'][1], ORANGE)
# v2: scrambled result-box spectrum (#e85d0f line art, nested in the round-8 resources) is covered with the box fill and
# redrawn in #d94d00. Geometry = the round-8 line art (baseline 0.608 pt, sticks 0.540 pt, butt caps), measured at 1200 dpi:
# ink x 365.36-451.52, baseline 347.26-347.92, peak top 329.98.
SPEC = fitz.Rect(364.6, 329.2, 452.3, 348.7)
_spm = sp.get_pixmap(dpi=1200, clip=SPEC, alpha=False)
_sa = np.frombuffer(_spm.samples, np.uint8).reshape(_spm.height, _spm.width, 3).astype(int)
SPEC_FILL = tuple(int(v) for v in _sa[0, 0])
assert SPEC_FILL == (255, 241, 232) and all(tuple(_sa[i, j]) == SPEC_FILL for i in (0, -1) for j in (0, -1)), SPEC_FILL
pg.draw_rect(SPEC, color=None, fill=rgb(SPEC_FILL), width=0)
SPEC_BASE_Y, SPEC_DX = 347.59, 36.96
STICKS = [(331.13, 302.46), (336.89, 302.69), (342.66, 301.98), (348.42, 302.69), (354.19, 285.56), (359.95, 302.69),
          (365.72, 302.22), (371.49, 302.69), (377.25, 302.46), (383.02, 302.69), (388.78, 301.98), (394.55, 302.69),
          (400.31, 302.46), (406.08, 302.69), (411.84, 302.69)]       # round-8 stick x / top in the spectrum's own frame
for x_, t_ in STICKS:
    pg.draw_line((x_ + SPEC_DX, t_ + 44.42), (x_ + SPEC_DX, SPEC_BASE_Y), color=ORANGE, width=0.5403642, lineCap=0)
pg.draw_line((365.39, SPEC_BASE_Y), (451.51, SPEC_BASE_Y), color=ORANGE, width=0.6079097, lineCap=0)
# MS/MS: label left, spectrum right, both centred vertically
box(MSBOX, BAND1_FILL)
mcy = (MSBOX.y0 + MSBOX.y1) / 2
put(MS_LBL[0], MS_LBL[1], MS_LBL[2], MSBOX.x0 + PAD, mcy + CAP * MS_LBL[2] / 2)
place(ICON, MSBOX.x1 - PAD - ICON_INK.x1, mcy - (ICON_INK.y0 + ICON_INK.y1) / 2)
# De novo
box(DNBOX, BAND2_FILL)
put_c(DN_LBL[0], DN_LBL[1], DN_LBL[2], X, DN_CY + CAP * DN_LBL[2] / 2)
assert TEXT[-1][1] in DNBOX and fitz.Rect(TEXT[-2][1]) in MSBOX, 'label outside its box'
# query boxes: one internal layout (round-8 scrambled box: title baselines +6.01 / +12.77, count +19.86, bars ink top +24.92),
# v2: bars scaled to BARW wide, and the whole block centred vertically in the box
QT1, QT2, QCNT, QBAR = 6.01, 12.77, 19.86, 24.92
BARW = QW - 2 * 7.0
def place_scaled(clip, ink_r, cx, ytop, w):
    sc = w / ink_r.width
    x0 = cx - (ink_r.x0 - clip.x0 + ink_r.width / 2) * sc; y0 = ytop - (ink_r.y0 - clip.y0) * sc
    pg.show_pdf_page(fitz.Rect(x0, y0, x0 + clip.width * sc, y0 + clip.height * sc), art, 0, clip=clip)
    return fitz.Rect(cx - w / 2, ytop, cx + w / 2, ytop + ink_r.height * sc)
QBLK_TOP, QBLK_BOT = QT1 - CAP * 5.741, QBAR + BN_INK.height * BARW / BN_INK.width
QOFF = (QH - (QBLK_BOT - QBLK_TOP)) / 2 - QBLK_TOP                  # shift so the content block is centred
BARS_DRAWN = []
for b_, fill, title, col, cnt, bars, bink in ((ALPHB, NAVYBOX_FILL, ('AlphaNovo peptide', 'predictions'), NAVY, '24,560 unique peptides', BARS_N, BN_INK),
                                              (SCRB, SCRBOX_FILL, ('Composition-matched', 'scrambled peptides'), ORANGE, '36,588 unique peptides', BARS_S, BS_INK)):
    box(b_, fill)
    cx = (b_.x0 + b_.x1) / 2
    put_c(title[0], 'H', 5.741, cx, b_.y0 + QOFF + QT1, col); put_c(title[1], 'H', 5.741, cx, b_.y0 + QOFF + QT2, col)
    put_c(cnt, 'B', 5.404, cx, b_.y0 + QOFF + QCNT)
    BARS_DRAWN.append(place_scaled(bars, bink, cx, b_.y0 + QOFF + QBAR, BARW))
    assert max(tw(t, 'H', 5.741) for t in title) <= b_.width - 2 * 4.5 and tw(cnt, 'B', 5.404) <= b_.width - 2 * 4.5, title

# funnel: paint each band from its top boundary down to the bottom, top band first
def arc(cx, cy, rx, ry, t0, t1, n=60):
    return [fitz.Point(cx + rx * math.cos(t0 + (t1 - t0) * i / n), cy + ry * math.sin(t0 + (t1 - t0) * i / n)) for i in range(n + 1)]
def region(yk, col):
    w = half(yk)
    pts = arc(X, yk, w, ry_of(w), math.pi, 0)                  # lower arc of the boundary, left -> right
    pts += arc(X, YB, WB, RB, 0, math.pi)                        # right wall down, bottom arc right -> left, left wall up
    pg.draw_polyline(pts, color=None, fill=col, closePath=True, width=0)
for k in range(4): region(YK[k], COL[k])
# bowl bottom: the old funnel's darker cap
CAPY = YB - 3.0                                          # cap = band between the arc at YB-3 and the bottom arc
w_c = half(CAPY)
pg.draw_polyline(arc(X, CAPY, w_c, ry_of(w_c), math.pi, 0) + arc(X, YB, WB, RB, 0, math.pi), color=None, fill=CAP_COL, closePath=True, width=0)
# rim = top of the reservoir layer (reservoir cylinder colours)
pg.draw_oval(fitz.Rect(X - RX, CY - RY, X + RX, CY + RY), color=RIM_STROKE, fill=RIM_FILL, width=RIM_SW)

# band labels
def lines_c(lines, ytop, ybot, pitch=6.3):
    # block from the first line's cap top to the last baseline, centred between ytop and ybot
    b0 = (ytop + ybot) / 2 - (pitch * (len(lines) - 1) + CAP * lines[0][2]) / 2 + CAP * lines[0][2]
    for i, (t, f, s) in enumerate(lines):
        put_c(t, f, s, X, b0 + i * pitch)
lines_c([('Gut sequence reservoir', 'H', 5.741), ('233 M proteins', 'H', 5.741), ('UHGG · GMGC · GMSC', 'B', 5.083),
         ('+ human Swiss-Prot', 'B', 5.083)], VTOP[0], VTOP[1])
lines_c([('Exact match ·', 'H', 5.404), ('study-supported proteins', 'H', 5.404)], VTOP[1], VTOP[2])   # v2
lines_c([('Per-acquisition match', 'H', 5.404)], VTOP[2], VTOP[3])          # v29: band labels are the actions
# parsimony band: label + glyph as one centred group
GW, GH = 18.0, 11.9
pl = ('Parsimony', 'H', 5.404)
grp = tw(*pl) + 4.0 + GW
gx0 = X - grp / 2
pcy = (VTOP[3] + VTOP[4]) / 2
put(pl[0], pl[1], pl[2], gx0, pcy + CAP * pl[2] / 2)
GX = gx0 + tw(*pl) + 4.0; GY = pcy - GH / 2
# glyph: three proteins (bars) carrying peptide ticks; bars A and B cover all four peptides, bar C (peptides 2-3 only) is redundant
TICK = [2.0, 6.5, 11.0, 15.5]                          # peptide x positions in the glyph
BARS = [(0.0, 13.0, [0, 1, 2], hexc('#2a7ab9')),       # A
        (5.0, 18.0, [2, 3], hexc('#3a9e5b')),          # B
        (1.0, 13.0, [1, 2], hexc('#d3dade'))]          # C (redundant, grey)
BH, BG_, TX = 1.8, 2.7, 0.55               # bar height, bar gap, tick overhang (ticks of two rows never touch)
for i, (a, b, ticks, col) in enumerate(BARS):
    y0 = GY + i * (BH + BG_) + TX
    pg.draw_rect(fitz.Rect(GX + a, y0, GX + b, y0 + BH), color=None, fill=col, width=0, radius=0.5)
    tc = NAVY if i < 2 else hexc('#9aa5ae')
    for t in ticks:
        pg.draw_line((GX + TICK[t], y0 - TX), (GX + TICK[t], y0 + BH + TX), color=tc, width=0.55, lineCap=0)
GLYPH = fitz.Rect(GX, GY, GX + GW, GY + 3 * BH + 2 * BG_ + 2 * TX)
assert BG_ - 2 * TX >= 1.5 and abs(GLYPH.height - GH) < 0.01
assert all(any(t in b[2] for b in BARS[:2]) for t in range(4))       # A + B cover every peptide

# Database search box (column-1/2 style: white fill, navy 0.675 pt stroke, 2.03 pt corners)
pg.draw_rect(DBBOX, color=NAVY, fill=(1, 1, 1), width=0.6754552, radius=2.03 / DBBOX.height)
put_c('Database search', 'B', 6.079, X, (DBBOX.y0 + DBBOX.y1) / 2 + CAP * 6.079 / 2)

# comparison table: round-8 table, rows re-laid on the funnel bands (x unchanged)
TAB_TXT = [s for s in spans(sp) if s['bbox'][0] > 400 and 176 < s['bbox'][1] < 279]
assert len(TAB_TXT) == 18, len(TAB_TXT)
OLD_Y = [189.81, 211.41, 233.02, 254.62, 276.23]
TAB_X0, TAB_X1, TAB_VX, TAB_COL, TAB_W = 445.42, 507.42, 476.41, hexc('#cfd9de'), 0.27450090646743774
NEW_Y = VTOP
for s in TAB_TXT:
    y = s['origin'][1]
    if y < OLD_Y[0]: dy = NEW_Y[0] - OLD_Y[0]
    else:
        r_ = max(i for i in range(4) if y > OLD_Y[i])
        dy = (NEW_Y[r_] + NEW_Y[r_ + 1]) / 2 - (OLD_Y[r_] + OLD_Y[r_ + 1]) / 2
    f = 'H' if 'Halbfett' in s['font'] else 'B'
    c_ = '#%06x' % s['color']
    t_, x_ = s['text'], s['origin'][0]
    if t_ == '–':                                     # v2: scrambled study-supported = its candidates (worked acquisition only)
        q_ = [q for q in TAB_TXT if q['text'] == '5,841'][0]
        assert '%06x' % s['color'] == 'ff4d00' and abs(s['size'] - q_['size']) < 1e-6 and q_['font'] == s['font']
        t_, x_ = '5,841', q_['origin'][0]
    put(t_, f, s['size'], x_, y + dy, col=ORANGE if c_ == '#ff4d00' else hexc(c_))   # v29: one orange
pg.draw_rect(fitz.Rect(TAB_X0, NEW_Y[0], TAB_X1, NEW_Y[4]), color=TAB_COL, width=TAB_W)
pg.draw_line((TAB_VX, NEW_Y[0]), (TAB_VX, NEW_Y[4]), color=TAB_COL, width=TAB_W)
for y_ in NEW_Y[1:4]: pg.draw_line((TAB_X0, y_), (TAB_X1, y_), color=TAB_COL, width=TAB_W)

# ---------------- arrows: one style, all straight ----------------
def arrow(p0, p1, col):
    p0, p1 = fitz.Point(p0), fitz.Point(p1)
    assert abs(p0.x - p1.x) < 1e-9 or abs(p0.y - p1.y) < 1e-9, 'arrow not straight'
    v = p1 - p0; u = v / abs(v); n = fitz.Point(-u.y, u.x); base = p1 - u * HL
    pg.draw_line(p0, base + u * 0.2, color=col, width=LW, lineCap=0)
    pg.draw_polyline([p1, base + n * HW, base - n * HW], color=None, fill=col, closePath=True, width=0)
    return dict(p0=[round(p0.x, 2), round(p0.y, 2)], p1=[round(p1.x, 2), round(p1.y, 2)], length=round(abs(v), 2))
A = {}
A['MS/MS -> De novo'] = arrow((X, MSBOX.y1), (X, DNBOX.y0), NAVY)
A['De novo -> AlphaNovo'] = arrow((X, DNBOX.y1), (X, ALPHB.y0), NAVY)
A['AlphaNovo -> scrambled'] = arrow((ALPHB.x1, (ALPHB.y0 + ALPHB.y1) / 2), (SCRB.x0, (SCRB.y0 + SCRB.y1) / 2), ORANGE)
A['AlphaNovo -> rim'] = arrow((XN, ALPHB.y1), (XN, rim_top(XN, CY)), NAVY)
A['scrambled -> rim'] = arrow((XO, SCRB.y1), (XO, rim_top(XO, CY)), ORANGE)
yb = lambda x: YB + RB * math.sqrt(1 - ((x - X) / WB) ** 2)
A['funnel -> Database search (navy)'] = arrow((X - XB, yb(X - XB)), (X - XB, DBBOX.y0 - 0.34), NAVY)
A['funnel -> Database search (orange)'] = arrow((X + XB, yb(X + XB)), (X + XB, DBBOX.y0 - 0.34), ORANGE)
OUTA_CX, OUTS_CX = (243.16 + 342.08) / 2, (358.97 + 457.89) / 2       # round-8 result boxes
assert DBBOX.x0 + 6 < OUTA_CX and OUTS_CX < DBBOX.x1 - 6
A['Database search -> 29,663'] = arrow((OUTA_CX, DBBOX.y1 + 0.34), (OUTA_CX, RESULT_TOP), NAVY)
A['Database search -> 0'] = arrow((OUTS_CX, DBBOX.y1 + 0.34), (OUTS_CX, RESULT_TOP), ORANGE)

for col, wr in WRITERS.items(): wr.write_text(pg, color=col)
# attach the overlay to the base page as one Form XObject (no resource re-scan of the base page)
out.insert_pdf(ov)
op1 = out[1]
_cont = op1.read_contents()
_rt, _rv = out.xref_get_key(op1.xref, 'Resources')
_res = _rv if _rt == 'dict' else _rv
FX = out.get_new_xref()
out.update_object(FX, f'<< /Type /XObject /Subtype /Form /BBox [0 0 {W} {H}] /Resources {_res} >>')
out.update_stream(FX, _cont)
p0 = out[0]
_owner, _key = p0.xref, 'Resources'                  # walk indirect objects down to the XObject dictionary
_t, _v = out.xref_get_key(_owner, _key)
if _t == 'xref': _owner, _key = int(_v.split()[0]), ''
_t, _v = out.xref_get_key(_owner, (_key + '/' if _key else '') + 'XObject')
if _t == 'xref': out.xref_set_key(int(_v.split()[0]), 'FLfunnelv2', f'{FX} 0 R')
else: out.xref_set_key(_owner, (_key + '/' if _key else '') + 'XObject/FLfunnelv2', f'{FX} 0 R')
CX_ = out.get_new_xref(); out.update_object(CX_, '<< >>'); out.update_stream(CX_, b'Q q /FLfunnelv2 Do Q')
CQ_ = out.get_new_xref(); out.update_object(CQ_, '<< >>'); out.update_stream(CQ_, b'q')   # base content wrapped in q ... Q
_ct, _cv = out.xref_get_key(p0.xref, 'Contents')
_old = _cv.strip('[]') if _ct == 'array' else _cv
out.xref_set_key(p0.xref, 'Contents', f'[{CQ_} 0 R {_old} {CX_} 0 R]')
out.delete_page(1)
out.save(OUT / 'Figure1.pdf', garbage=4, deflate=True); out.close()

# ---------------- checks ----------------
o = fitz.open(OUT / 'Figure1.pdf'); op = o[0]
opm = op.get_pixmap(dpi=300, alpha=False); opm.save(OUT / 'Figure1.png')
A_ = np.frombuffer(opm.samples, np.uint8).reshape(opm.height, opm.width, 3)
B_ = np.frombuffer(bpm.samples, np.uint8).reshape(bpm.height, bpm.width, 3)
row0, colL = int(372 * K), int(230 * K)
mask = np.ones(A_.shape[:2], bool); mask[:row0, colL:] = False
DIFF_OUTSIDE = int(((A_ != B_).any(axis=2) & mask).sum())
assert DIFF_OUTSIDE == 0, DIFF_OUTSIDE
# result boxes unchanged (below the arrow tips)
_d = (A_[int(309 * K):row0, colL:] != B_[int(309 * K):row0, colL:]).any(axis=2)
_zx0, _zy0, _zx1, _zy1 = [int(v * K) for v in (_zr.x0 - 1, _zr.y0 - 1, _zr.x1 + 1, _zr.y1 + 1)]
_m = np.ones_like(_d); _m[_zy0 - int(309 * K):_zy1 - int(309 * K), _zx0 - colL:_zx1 - colL] = False
_sx0, _sy0, _sx1, _sy1 = [int(v * K) for v in (SPEC.x0 - 1, SPEC.y0 - 1, SPEC.x1 + 1, SPEC.y1 + 1)]
_m[_sy0 - int(309 * K):_sy1 - int(309 * K), _sx0 - colL:_sx1 - colL] = False      # v2: spectrum recoloured
rb_ = int((_d & _m).sum()); ZERO_DIFF = int((_d & ~_m).sum())
assert rb_ == 0 and ZERO_DIFF > 0, (rb_, ZERO_DIFF)
# v2: the redrawn spectrum has the round-8 shape: ink masks at 1200 dpi agree except at antialiased edges
def _inkmask(page):
    pm = page.get_pixmap(dpi=1200, clip=SPEC, alpha=False)
    a = np.frombuffer(pm.samples, np.uint8).reshape(pm.height, pm.width, 3).astype(int)
    return np.abs(a - np.array(SPEC_FILL)).sum(2) > 150
_mo = _inkmask(sp); _mn = _inkmask(op)
SPEC_XOR, SPEC_INK = int((_mo ^ _mn).sum()), int(_mo.sum())
assert SPEC_XOR <= 0.05 * SPEC_INK, (SPEC_XOR, SPEC_INK)
OLD_ORANGE_PX = int((A_[:row0, colL:] == np.array([232, 93, 15])).all(2).sum())
assert OLD_ORANGE_PX == 0, OLD_ORANGE_PX
# v2: pixel diff against the v29 Fig1 build outside the FastaLake column
V29 = ART / 'Figure1_v29_intermediate.pdf'      # = producer input v29/figures/Fig1/Figure1.pdf (same sha256)
_vpm = fitz.open(V29)[0].get_pixmap(dpi=300, alpha=False)
V_ = np.frombuffer(_vpm.samples, np.uint8).reshape(_vpm.height, _vpm.width, 3)
DIFF_OUTSIDE_V29 = int(((A_ != V_).any(axis=2) & mask).sum())
assert DIFF_OUTSIDE_V29 == 0, DIFF_OUTSIDE_V29
_dv = (A_ != V_).any(axis=2); _ys, _xs = np.nonzero(_dv)
DIFF_V29_BBOX_PT = [round(float(v) / K, 1) for v in (_xs.min(), _ys.min(), _xs.max() + 1, _ys.max() + 1)]
# text: every label is real text, sizes 5-7 pt, nothing removed is left
col_sp = [s for s in spans(op) if s['bbox'][0] > 230 and s['bbox'][1] < 372]
col_txt = [s['text'].strip() for s in col_sp]
for t in ['MS/MS', 'De novo', 'AlphaNovo peptide', 'predictions', '24,560 unique peptides', 'Composition-matched', 'scrambled peptides',
          '36,588 unique peptides', 'Gut sequence reservoir', '233 M proteins', 'UHGG · GMGC · GMSC', '+ human Swiss-Prot', 'Exact match ·', 'study-supported proteins',
          'Per-acquisition match', 'Parsimony', 'Database search', '29,663 peptides', '0 peptides',
          'Accepted peptides (authentic)', 'Accepted peptides (scrambled)', 'Real', 'Scrambled', '18 M', '6.67 M', '5,841', '12,606', '525', 'FastaLake']:
    assert t in col_txt, t
for bad in ['MS/MS (DDA)', 'DDA', 'Study-level lookup', 'Pool queries across 58 acquisitions', '18 M study-supported proteins',
            'De novo predictions', '1', '2', '3', '4', '5', 'Exact-match', 'Sample-', 'specific', 'FASTA']:
    assert bad not in col_txt, bad
assert col_txt.count('233 M') == 2 and col_txt.count('5,841') == 2 and '–' not in col_txt
assert col_txt.count('Candidates,') == 1 and col_txt.count('this acquisition') == 1   # table row label only
ORANGES = sorted({'%06x' % s['color'] for s in col_sp if s['color'] != 0x0a2d5e and s['color'] != 0x126ca4})
assert ORANGES == ['d94d00'], ORANGES
assert all(5.0 <= s['size'] <= 7.0 for s in col_sp if s['text'].strip() != 'FastaLake'), sorted({round(s['size'], 2) for s in col_sp})
# no two labels overlap; no label touches a box edge, the glyph, an arrow, the funnel rim outline
tr = [(t, r) for t, r in TEXT]
for i, (t1, r1) in enumerate(tr):
    for t2, r2 in tr[i + 1:]:
        assert not (r1 + (-0.5, -0.5, 0.5, 0.5)).intersects(r2), (t1, t2)
    assert not r1.intersects(GLYPH + (-1, -1, 1, 1)), t1
    for br in BARS_DRAWN: assert not r1.intersects(br + (-0.8, -0.8, 0.8, 0.8)), t1
for name, a in A.items():
    (x0, y0), (x1, y1) = a['p0'], a['p1']
    seg = fitz.Rect(min(x0, x1) - LW / 2 - HW, min(y0, y1), max(x0, x1) + LW / 2 + HW, max(y0, y1)) if x0 == x1 else \
          fitz.Rect(min(x0, x1), y0 - HW, max(x0, x1), y0 + HW)
    for t, r in tr: assert not seg.intersects(r), (name, t)
    assert a['length'] >= HL + 3, (name, a['length'])
for t, r in tr:    # band labels inside the funnel walls with >= 3 pt margin; table labels clear of the walls by >= 5 pt
    if VTOP[0] - 1 < r.y0 and r.y1 < VTOP[4] + 1 and r.x0 < 400:
        for y in (r.y0, r.y1): assert r.x0 > X - half(y) + 3 and r.x1 < X + half(y) - 3, (t, r)
    if r.x0 > 400 and VTOP[0] < r.y1 < VTOP[4]:
        assert r.x0 > X + half(r.y0) + 5, (t, r.x0, X + half(r.y0))
# boxes do not overlap each other; clearance >= 5 pt from the panel edges (233.3 / 512.45 pt)
BOXES = {'MS/MS': MSBOX, 'De novo': DNBOX, 'AlphaNovo': ALPHB, 'scrambled': SCRB, 'Database search': DBBOX,
         'funnel': fitz.Rect(X - RX, CY - RY, X + RX, YB + RB), 'table': fitz.Rect(405.0, NEW_Y[0] - 12, TAB_X1, NEW_Y[4])}
ks = list(BOXES)
for i, a in enumerate(ks):
    for b in ks[i + 1:]:
        if {a, b} == {'funnel', 'table'}: continue
        assert not BOXES[a].intersects(BOXES[b]), (a, b)
    assert BOXES[a].x0 >= 233.3 + 5 and BOXES[a].x1 <= 512.45 - 5, a
rec = dict(base=BASE.name, base_sha256=sha(BASE), output_pdf_sha256=sha(OUT / 'Figure1.pdf'), png_px=[opm.width, opm.height],
           pixels_changed_outside_column_300dpi=DIFF_OUTSIDE, column_region_pt=[230, 0, round(W, 2), 372],
           arrows=A, arrow_style=dict(line_pt=LW, head_len_pt=HL, head_width_pt=2 * HW),
           layout=dict(axis_x=X, arm_x=[XN, XO], rim=[X, round(CY, 2), RX, RY], band_boundaries_y=[round(v, 2) for v in YK],
                       band_visible_on_axis=[round(v, 2) for v in VTOP], boxes={k: [round(v, 2) for v in r] for k, r in BOXES.items()}),
           fonts=sorted({s['font'] for s in col_sp}), result_box_pixels_changed_outside_0_peptides=rb_, zero_peptides_title_and_spectrum_changed_px=ZERO_DIFF,
           spectrum_inkmask_xor_px_1200dpi=[SPEC_XOR, SPEC_INK], old_orange_e85d0f_px_in_column=OLD_ORANGE_PX,
           v29_build=V29.name, v29_sha256=sha(V29), pixels_changed_outside_column_vs_v29_300dpi=DIFF_OUTSIDE_V29, diff_vs_v29_bbox_pt=DIFF_V29_BBOX_PT, text_colours=ORANGES, sizes=sorted({round(s['size'], 3) for s in col_sp}))
(OUT / '_build_receipt.json').write_text(json.dumps(rec, indent=1, ensure_ascii=False) + '\n')
print(json.dumps(rec['layout'], indent=1)); print(json.dumps(A)); print('outside-column changed px:', DIFF_OUTSIDE, 'vs v29:', DIFF_OUTSIDE_V29, 'v29 diff bbox', DIFF_V29_BBOX_PT, 'spec xor', SPEC_XOR, SPEC_INK)
for nm, r in ((('full', fitz.Rect(0, 0, W, 372)), ('col', fitz.Rect(228, 0, W, 372))) if (OUT / '_work').is_dir() else ()):   # review crops, only if _work/ exists
    op.get_pixmap(dpi=300 if nm == 'full' else 450, clip=r, alpha=False).save(OUT / f'_work/check_{nm}.png')
