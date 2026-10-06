"""Figure 3 at 9-50 residues: rebuild panels b-e with MIN_LENGTH = 9 and replace them in the v28f base.

Base: LAB/Figure3_v28f.pdf (= final/fig3_labels/Figure3.pdf; panel b legend "Same series" /
"Different donor"). Only the b-e band (page y 79.92-364.80 pt = 300-dpi rows 333-1519) is changed:
1. panels_len950.panels_bcde (repo panels.py with MIN_LENGTH 9; nothing else changed) draws b-e.
2. The band is covered with white and the new b-e PDF placed where the assembly chain placed the
   producer canvas (repo panels.composite_rects()['bcde']).
3. Text is re-typeset exactly as build_FAR21.py did (Soehne Buch 7 pt axis labels, 6.5 pt ticks,
   Soehne Halbfett 8 pt letters; old span covered with white), and the b legend is relabelled exactly
   as fig3_labels/relabel.py did ("Same series", "Different donor").
Usage: build_fig3_len950.py MIN_LENGTH OUTSTEM   (MIN_LENGTH 5 = control that must reproduce the base)
"""
import os  # release: paths come from environment variables
import hashlib, json, re, sys
from pathlib import Path
import fitz
import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panels_len950 as P  # noqa: E402

LAB = Path(os.environ['FASTALAKE_DATA_DIR'])
BASE = LAB / 'Figure3_v28f.pdf'
BASE_PNG = LAB / 'review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/fig3_labels/Figure3.png'
SD = LAB / 'release_prep/fastalake-paper_v28/source_data/Fig3/transfer_pool_comparisons.tsv'
STYLE = Path(os.environ.get('FASTALAKE_FONT_DIR', ''))
REG, BOLD = STYLE / 'soehne-buch.ttf', STYLE / 'soehne-halbfett.ttf'
BAND = fitz.Rect(0, 333 * 72 / 300, P.COMPOSITE_PAGE_PT[0], 1520 * 72 / 300)
LEGEND_OLD = {'Same-design donor (T1)', 'Designated donor control (T4)'}
RELABEL_BOX = fitz.Rect(54.0, 186.8, 156.0, 204.7)          # relabel.py
RELABEL = [(193.862, 'Same series'), (202.560, 'Different donor')]  # relabel.py
RELABEL_X = 54.517
INK = (0x4e / 255, 0x4f / 255, 0x51 / 255)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main(min_length, outstem):
    outstem = Path(outstem)
    outstem.parent.mkdir(parents=True, exist_ok=True)
    P.MIN_LENGTH = int(min_length)
    comp = pd.read_csv(SD, sep='\t')
    med = P.panels_bcde(comp, outstem.parent / (outstem.name + '_panels_bcde'))
    panel_pdf = outstem.parent / (outstem.name + '_panels_bcde.pdf')
    src = fitz.open(panel_pdf)
    pw, ph = src[0].rect.width, src[0].rect.height
    rect = fitz.Rect(*P.composite_rects()['bcde'])
    s = min(rect.width / pw, rect.height / ph)
    ox = rect.x0 + (rect.width - pw * s) / 2
    oy = rect.y0 + (rect.height - ph * s) / 2
    tr = lambda b: fitz.Rect(ox + b[0] * s, oy + b[1] * s, ox + b[2] * s, oy + b[3] * s)

    doc = fitz.open(BASE)
    page = doc[0]
    page.draw_rect(BAND, color=None, fill=(1, 1, 1), overlay=True)
    page.show_pdf_page(rect, src, 0, clip=None)
    page.insert_font(fontname='L950Reg', fontfile=str(REG))
    page.insert_font(fontname='L950Bold', fontfile=str(BOLD))
    font, bf = fitz.Font(fontfile=str(REG)), fitz.Font(fontfile=str(BOLD))
    spans = [dict(sp, direction=l['dir']) for b in src[0].get_text('dict')['blocks']
             for l in b.get('lines', []) for sp in l['spans']]
    edits = []
    for sp in spans:
        text = sp['text']
        box = tr(sp['bbox'])
        assert BAND.y0 <= box.y0 and box.y1 <= BAND.y1, (text, box)
        if text in LEGEND_OLD:
            continue
        letter = text in ['b', 'c', 'd', 'e']
        numeric = bool(re.fullmatch(r'\d+(?:\.\d+)?', text))
        size = 8 if letter else 6.5 if numeric else 7
        f, fname = (bf, 'L950Bold') if letter else (font, 'L950Reg')
        length = f.text_length(text, fontsize=size)
        x, y = ox + sp['origin'][0] * s, oy + sp['origin'][1] * s
        rotate = 0
        if sp['direction'] == (0.0, -1.0):
            rotate = 90
            y = (box.y0 + box.y1) / 2 + length / 2
        elif letter:
            x = box.x0
        elif numeric:
            x = box.x1 - length
        else:
            x = (box.x0 + box.x1) / 2 - length / 2
        page.draw_rect(box + (-.2, -.2, .2, .2), color=None, fill=(1, 1, 1), overlay=True)
        color = tuple(((sp['color'] >> v) & 255) / 255 for v in [16, 8, 0])
        page.insert_text((x, y), text, fontname=fname, fontsize=size, color=color, rotate=rotate)
        edits.append(dict(text=text, bbox=[round(v, 3) for v in box], size=size))
    page.draw_rect(RELABEL_BOX, color=None, fill=(1, 1, 1), overlay=True)
    page.insert_font(fontname='L950Buch', fontfile=str(REG))
    for y, t in RELABEL:
        page.insert_text((RELABEL_X, y), t, fontname='L950Buch', fontsize=6.5, color=INK)
    doc.save(str(outstem) + '.pdf', garbage=4, deflate=True)
    pix = fitz.open(str(outstem) + '.pdf')[0].get_pixmap(dpi=300, alpha=False)
    pix.save(str(outstem) + '.png')

    new = np.array(Image.open(str(outstem) + '.png').convert('RGB')).astype(int)
    base = np.array(Image.open(BASE_PNG).convert('RGB')).astype(int)
    assert new.shape == base.shape, (new.shape, base.shape)
    diff = (new != base).any(axis=2)
    r0, r1 = 333, 1520
    outside = diff.copy()
    outside[r0:r1, :] = False
    rec = dict(min_length=P.MIN_LENGTH, base_pdf=str(BASE), base_pdf_sha256=sha(BASE), base_png=str(BASE_PNG),
               base_png_sha256=sha(BASE_PNG), out_pdf_sha256=sha(str(outstem) + '.pdf'),
               out_png_sha256=sha(str(outstem) + '.png'), replaced_band_rows_300dpi=[r0, r1 - 1],
               replaced_band_pt=list(BAND), placement_rect_pt=list(rect), scale=s,
               differing_pixels_outside_band=int(outside.sum()), differing_pixels_inside_band=int(diff[r0:r1].sum()),
               retyped_spans=len(edits),
               plotted_medians={'transfer': {'|'.join(k): v for k, v in med['transfer'].items()},
                                'pool': {'|'.join(k): v for k, v in med['pool'].items()}})
    Image.fromarray((np.where(diff[..., None], [255, 0, 0], new)).astype(np.uint8)).save(str(outstem) + '_diff.png')
    Path(str(outstem) + '_receipt.json').write_text(json.dumps(rec, indent=1) + '\n')
    print(json.dumps({k: v for k, v in rec.items() if k != 'plotted_medians'}, indent=1))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
