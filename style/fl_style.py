"""Canonical FastaLake figure style — 1:1 with the AlphaNovo Figure 5 look.

Imports the same Mann-lab maker modules the AlphaNovo figures use
(a lab styleguide; its three used symbols are reimplemented in
_house.py for this repository): Soehne font, ink #4e4f51, 0.6 line weights, NO bold anywhere,
Okabe-Ito colorblind-safe palette, NBT widths, TrueType-embedded multi-format save.

This module is also a DROP-IN SUPERSET of the old manuscript/figures/extras/_fl_style.py
API (FL, METAG, REF, HOST, ACCENT, INK, MM, init_style, panel_letter, despine, hnum,
read_tsv, title_block, grid_y, fl_cmap, save), so a v16.67 figure adapts to the standard
by swapping its inline style block for `import fl_style as fl` — the semantic colour
names keep working but now resolve to the Okabe-Ito / Soehne standard.

Native API:  fl.init(); plt.figure(figsize=fl.mm(183,120)); fl.style(ax);
             fl.panel(ax,'a'); fl.save(fig, path)

NBT HARD RULE (enforced here): figures carry NO subpanel titles and NO floating
narrative/method text — only axis/tick labels, direct data-value labels, legends
and panel letters; titles/interpretation live in the caption. `panel()` draws the
LETTER only (any `title=` is dropped with a warning), `title_block()` is a disabled
no-op, and `save()` runs `enforce_nbt()` to strip any stray on-canvas title.
"""

import os
import sys
import warnings
from pathlib import Path

import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

# ---- paper-repository adaptation (the only change from the canonical copy) ----
# The canonical module imports two lab styleguide modules from FL_STYLE_DIR. Here the
# three symbols it uses (COLORBLIND_SAFE_4, lighten_color, init_plotting) are
# reimplemented in _house.py with identical behaviour, so the only external input is
# the licensed Söhne typeface, read from FASTALAKE_FONT_DIR (FL_STYLE_DIR is accepted
# as an alias). Without the fonts every colour, size and line is unchanged and text
# renders in the matplotlib default face, so text pixels differ from the published
# panels. See style/README.md.
try:
    from . import _house as _sp  # noqa: E402  (imported as style.fl_style)
except ImportError:
    import _house as _sp  # noqa: E402  (imported as fl_style with style/ on sys.path)
_cp = _sp

_STYLE = Path(os.environ.get("FASTALAKE_FONT_DIR") or os.environ.get("FL_STYLE_DIR") or "")
_HAVE_STYLE = bool(str(_STYLE)) and (_STYLE / "soehne-buch.otf").is_file()
if not _HAVE_STYLE:
    warnings.warn(
        "FASTALAKE_FONT_DIR is unset or lacks soehne-buch.otf: using the default "
        "typeface. Colours, sizes and layout are unchanged; text pixels will differ "
        "from the published panels.",
        stacklevel=2,
    )

# ---- palette (Okabe-Ito colorblind-safe, the AlphaNovo standard) -----------
INK = "#4e4f51"
GREY = "#7f7f7f"
MUTE = "#8a8a8a"
GRID = "#E9E9E9"
BLUE, ORANGE, GREEN, VERM = _cp.COLORBLIND_SAFE_4  # #0072B2 #E69F00 #009E73 #D55E00
# full Okabe-Ito 8-colour categorical set (colourblind-safe). Use these ONLY — no
# off-palette hues anywhere. cat(n) hands back the first n in a good categorical order.
SKY = "#56B4E9"  # sky blue
YELLOW = "#F0E442"  # yellow (weak on white; use late)
PURPLE = "#CC79A7"  # reddish purple
BLACK = "#000000"
OKABE_ITO = [BLUE, ORANGE, GREEN, VERM, SKY, PURPLE, YELLOW, BLACK]


def cat(n):
    """First n Okabe-Ito categorical colours (n<=8)."""
    return OKABE_ITO[:n]


# semantic names (drop-in for _fl_style.py); all resolve to Okabe-Ito hues
FL = BLUE  # FastaLake (hero)
METAG = ORANGE  # metagenome / external comparator
REF = GREY  # reference / neutral (grey is a neutral, not a category)
HOST = VERM  # host / accent
ACCENT = VERM
MM = 1 / 25.4  # inches per millimetre (legacy figs do width*MM)

# 🔴 NATURE COMMUNICATIONS COLUMN WIDTHS, verified 2026-08-10 against their how-to-submit
# page: "A single column width measures 88 mm and a double column width measures 180 mm."
# The project previously built at the NBT 89/183. Use these constants, never a literal --
# a retyped width is how two figures end up at different scales and nobody notices.
SINGLE_MM = 88.0
DOUBLE_MM = 180.0
# NComms: "The thinnest lines in the final figure should be no smaller than one point wide."
MIN_LINE_PT = 1.0

# Nature Communications: labelling must be legible at final size; the working floor this
# project has used is 5 pt. Measured 2026-08-11 across the 57 composed panels, EIGHT sat
# below it -- 4.6 pt in Fig2 b/c and Fig6 d/e, 4.9 pt in Fig1 c/e/g and Fig2 e/g. Every one
# came from matplotlib scaling a label down to fit, which no rcParam prevents because the
# rcParam sets the REQUESTED size and the autoscaler runs afterwards.
MIN_TYPE_PT = 5.0

# The house title size. Panels built on different days topped out at 7, 8, 9 and 10 pt, so
# the same kind of label was 7 pt in one panel and 9 pt in the one beside it. rcParams set
# this for anything drawn through them; `normalise_type()` catches the explicit overrides.
TITLE_PT = 7.0

# 🔴 THE ONE PRINCIPLED EXCEPTION TO THE 7 pt CEILING, and it exists to satisfy the FLOOR.
#
# A log axis is labelled by matplotlib's LogFormatterSciNotation, which writes `$10^{3}$`
# through mathtext -- and mathtext draws a superscript at 0.7x the surrounding size. At the
# house 7 pt that puts the EXPONENT at 4.9 pt, below the 5 pt floor, and no font setting
# reaches it: the ratio is internal to mathtext, not an rcParam. Measured 2026-08-11, that
# was eight of the eleven remaining floor violations across the composed panels, in six
# different panels, and it is the same defect CLAUDE.md already records for a check that
# passed a panel putting 4.55 pt glyphs on the page.
#
# 5.0 / 0.7 = 7.15, so a log axis is labelled at 7.2 pt and its exponents land at 5.04.
# The alternative -- rewriting log ticks as plain numbers -- turns 10^9 into 1000000000 and
# is worse. The base label is 0.2 pt above the house size, on log axes only.
LOG_TICK_PT = 7.2
_MATHTEXT_RATIO = 0.7

PALETTE = {
    "host": VERM,
    "microbial": BLUE,
    "metag": PURPLE,
    "metat": SKY,
    "metap": VERM,
    "accent": ORANGE,
    "ink": INK,
    "grey": GREY,
}

fl_cmap = LinearSegmentedColormap.from_list("fl", [_sp.lighten_color(FL, 0.6), FL, "#0a3550"])
# Okabe-anchored scales so sequential/diverging maps stay inside the palette:
seq_cmap = LinearSegmentedColormap.from_list("fl_seq", ["#ffffff", BLUE, "#0a3550"])  # low->high
div_cmap = LinearSegmentedColormap.from_list("fl_div", [BLUE, "#f2f2f2", VERM])  # -/0/+
lighten = _sp.lighten_color


def mm(w_mm, h_mm):
    return (w_mm * MM, h_mm * MM)


def _register_soehne_weights():
    """Make the bold and semibold Söhne faces AVAILABLE to matplotlib. Idempotent.

    🔴 THIS CHANGES NO POLICY. It fixes a defect where asking for bold produced regular
    and said nothing.

    `shared_functions_plotting.init_plotting()` registers exactly ONE file --
    `soehne-buch.otf` -- as the family `soehne`, and then sets `font.family = "soehne"`.
    `soehne-kraftig.otf` (semibold) and `soehne-halbfett.otf` (bold) sit unused in the
    same directory. matplotlib's font resolver, asked for family `soehne` at
    `weight="bold"`, finds only the one registered entry, scores it as the best available
    match, and returns it. So `fontweight="bold"` rendered REGULAR, with no warning and
    no visible failure -- the caller sees text, just not the weight it asked for.

    That matters now because NComms requires bold lower-case panel letters
    ("Figures divided into parts should be labelled with a lower-case bold a, b, ..."),
    and `fl_style.panel()` already passes `fontweight="bold"` on the strength of it --
    against a font manager that could not deliver it.

    ⚠️ THE BOLD POLICY IS UNRESOLVED AND IS NOT DECIDED HERE. `publication/figures/
    compose.py` stamps the composed figure's panel letters at `fontweight="normal"` and
    its docstring says "No bold anywhere, per the style rule"; FIGURE_STANDARD.md §2 says
    the panel letter is the one place bold IS allowed. Those two contradict and Peter has
    not ruled. This function only makes the weight OBTAINABLE, so that whichever way he
    rules, the code does what it says. compose.py is untouched.

    Registered at the END of the ttflist, not the front: `init_plotting` inserts `buch`
    at position 0 and the regular face must stay the default for an unweighted request.
    """
    from matplotlib import font_manager as _fm

    faces = {"soehne-kraftig.otf": "semibold", "soehne-halbfett.otf": "bold"}
    registered, missing = [], []
    have = {(e.name, e.weight) for e in _fm.fontManager.ttflist}
    for fname, weight in faces.items():
        path = _STYLE / fname
        if not path.is_file():
            missing.append(fname)
            continue
        if ("soehne", weight) in have:
            continue
        _fm.fontManager.ttflist.append(
            _fm.FontEntry(fname=str(path), name="soehne", weight=weight)
        )
        registered.append(f"{fname} as soehne/{weight}")
    if missing:
        warnings.warn(
            "Söhne bold faces not found in FASTALAKE_FONT_DIR: "
            + ", ".join(missing)
            + ". fontweight='bold' will silently render REGULAR.",
            stacklevel=2,
        )
    return registered


def init():
    """Register Soehne + lab rcParams (Okabe-Ito ink, >=1 pt lines, no bold). Idempotent."""
    _sp.init_plotting(font_path=str(_STYLE / "soehne-buch.otf"))
    if _HAVE_STYLE:
        _register_soehne_weights()
    mpl.rcParams.update(
        {
            "font.size": 7,
            "axes.titleweight": "normal",
            "axes.titlesize": 7,
            "legend.fontsize": 6.5,
            "xtick.major.size": 2.6,
            "ytick.major.size": 2.6,
            # 🔴 LINE-WEIGHT FLOOR, NComms: no line thinner than 1 pt in the final figure.
            # The project built at 0.6 pt throughout, which is below the stated floor.
            # Set here so every figure inherits it; builders passing an explicit lw=0.6
            # still override and are swept separately.
            "lines.linewidth": MIN_LINE_PT,
            "axes.linewidth": MIN_LINE_PT,
            "grid.linewidth": MIN_LINE_PT,
            "xtick.major.width": MIN_LINE_PT,
            "ytick.major.width": MIN_LINE_PT,
            "xtick.minor.width": MIN_LINE_PT,
            "ytick.minor.width": MIN_LINE_PT,
            "patch.linewidth": MIN_LINE_PT,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            # 🔴 NOT "tight". Changed 2026-08-11, and this is the root cause of a defect
            # that had reached 17 panels.
            #
            # `bbox_inches="tight"` CROPS the saved file to the drawn content, so a panel
            # authored at 180 mm lands at 175.98 mm whenever its outermost few millimetres
            # happen to be empty. Nothing errors. The panel is simply the wrong size, by an
            # amount that depends on the DATA -- a legend one character shorter moves it --
            # and compose.py then refuses it, or worse, a repair pads it back and the
            # padding has to be re-run after every rebuild.
            #
            # The authored figsize IS the specification here: this project draws every
            # panel at its true placed width so that 6.5 pt type is 6.5 pt on paper. A save
            # mode that silently renegotiates that width defeats the entire arrangement.
            #
            # fig_layout.apply_type() has undone this line since 2026-08-03 for exactly
            # this reason, with a comment explaining it. That fix only ever covered panels
            # drawn through fig_layout.Canvas; everything drawn straight onto a matplotlib
            # figure and saved through fl.save() stayed exposed. Fixing it at the source
            # covers both, and fig_layout's override is now a no-op that agrees.
            "savefig.bbox": None,
            "savefig.pad_inches": 0.0,
            "svg.fonttype": "none",
        }
    )
    # Mathtext must use Soehne too. Added 2026-07-28 (Peter) after Figure 5 panel d
    # was found embedding DejaVuSans: its axis label uses mathtext (`log$_{10}$`),
    # and matplotlib renders mathtext through `mathtext.fontset`, which is
    # INDEPENDENT of font.family and defaults to 'dejavusans'. Setting the family
    # alone never touched it, so any mathtext in any figure silently fell back.
    # 7 of the 25 figures the v23.17 build places were affected.
    #
    # Blast radius is exactly the defective glyphs: this can only change spans that
    # are currently rendering in the WRONG font, so it cannot alter correct type.
    # Each figure picks it up when its own freeze window re-renders.
    mpl.rcParams.update(
        {
            "mathtext.fontset": "custom",
            "mathtext.rm": "soehne",
            "mathtext.it": "soehne",
            "mathtext.bf": "soehne",
            "mathtext.sf": "soehne",
            "mathtext.tt": "soehne",
            "mathtext.cal": "soehne",
            "mathtext.default": "regular",
        }
    )


init_style = init  # legacy alias


def style(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    # width=MIN_LINE_PT, not 0.6. tick_params writes the artist directly, so this call
    # overrode the xtick/ytick.major.width rcParams init() sets -- the same "the floor is
    # set, then quietly put back" pattern fig_layout.apply_type() had. Fixed 2026-08-11.
    ax.tick_params(length=2.6, width=MIN_LINE_PT, labelsize=7)
    ax.xaxis.label.set_size(7)
    ax.yaxis.label.set_size(7)
    return ax


def panel(ax, letter, title="", fontsize=8):
    """Panel LETTER only, left-aligned, lab ink (non-bold). NBT rule: NO subpanel
    titles on the canvas — captions live under the figure. A non-empty `title` is
    dropped with a warning (signature kept for backward compatibility)."""
    if title:
        warnings.warn(
            f"fl_style.panel: subpanel title {title!r} dropped "
            "(NBT: titles belong in the caption, not on the figure)",
            stacklevel=2,
        )
    # 🔴 BOLD, lower case. NComms: "Figures divided into parts should be labelled with a
    # lower-case bold a, b, and so on, in the same type size as used elsewhere in the
    # figure." This is the ONLY place bold is permitted -- the no-bold house rule holds
    # everywhere else in the figure.
    ax.set_title(
        str(letter).strip().lower(), loc="left", fontsize=fontsize, color=INK,
        fontweight="bold", pad=4
    )


def _is_title_text(s):
    """A bare panel letter ('a', 'b)', 'c.') is allowed; longer text is a subpanel
    title that must move to the caption under the NBT rule."""
    s = (s or "").strip().rstrip(").:").strip()
    return len(s) > 2


def enforce_nbt(fig):
    """Strip any on-canvas subpanel title / figure suptitle before saving (NBT rule).
    Bare panel letters, data labels and legends are preserved. Returns (and warns) the
    list of stripped strings so a builder can move them into the caption."""
    stripped = []
    for ax in fig.axes:
        for t in (ax.title, getattr(ax, "_left_title", None), getattr(ax, "_right_title", None)):
            if t is not None and _is_title_text(t.get_text()):
                stripped.append(t.get_text().strip())
                t.set_text("")
    sup = getattr(fig, "_suptitle", None)
    if sup is not None and sup.get_text().strip():
        stripped.append(sup.get_text().strip())
        sup.set_text("")
    if stripped:
        warnings.warn(
            "fl_style.enforce_nbt stripped on-canvas title(s) "
            "(move to caption): " + " | ".join(stripped[:5]) + (" …" if len(stripped) > 5 else ""),
            stacklevel=2,
        )
    return stripped


def panel_letter(ax, s, x=-0.16, y=1.06, fontsize=9):
    """Legacy positioned panel letter — now NON-BOLD, lab ink (AlphaNovo standard)."""
    ax.text(
        x,
        y,
        s,
        transform=ax.transAxes,
        fontsize=fontsize,
        fontweight="normal",
        color=INK,
        va="bottom",
        ha="right",
    )


def despine(ax, keep=("left", "bottom")):
    for sp in ("top", "right", "left", "bottom"):
        ax.spines[sp].set_visible(sp in keep)


def grid_y(ax):
    # lw=MIN_LINE_PT, not 0.7. A grid rule is a line in the final figure and the NComms
    # floor does not carve out decoration. It stays visually quiet because the WEIGHT is
    # not what makes it recede -- GRID is #E9E9E9, a 9% grey. Raised 2026-08-11 with the
    # rest of the sweep, so no single figure differs from the other five.
    ax.grid(axis="y", color=GRID, lw=MIN_LINE_PT, zorder=0)
    ax.set_axisbelow(True)


def hnum(n):
    n = float(n)
    if n >= 1e6:
        return f"{n / 1e6:.1f}M" if n < 1e8 else f"{n / 1e6:.0f}M"
    if n >= 1e3:
        return f"{n / 1e3:.1f}k".replace(".0k", "k") if n < 1e5 else f"{int(n):,}"
    return f"{int(n):,}"


def read_tsv(path):
    with open(path) as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip() and not ln.startswith("#")]
    hdr, data = rows[0], rows[1:]
    return hdr, [dict(zip(hdr, r)) for r in data]


def title_block(fig, title, subtitle=None, x=0.02, y=0.99):
    """DISABLED under the NBT rule — on-canvas figure titles and subtitles are
    narrative text that belongs in the caption, not on the figure. Kept as a no-op
    (with a warning) for backward compatibility; draws nothing."""
    warnings.warn(
        f"fl_style.title_block is disabled (NBT: figure title/subtitle "
        f"{title!r} belong in the caption). Nothing drawn.",
        stacklevel=2,
    )
    return


def table(
    ax,
    col_labels,
    rows,
    col_widths=None,
    header_color=None,
    zebra=True,
    fontsize=7,
    highlight_cols=None,
    align=None,
):
    """Render a clean NBT-style table on `ax` (no vertical rules, thin horizontal
    header + bottom rules, optional zebra shading, Soehne ink). rows = list of lists
    of strings. highlight_cols = set of column indices to tint (e.g. the differing
    field). Returns the drawn height in axes fraction. Overlap-safe: wraps nothing —
    caller must size columns; long cells are the caller's responsibility."""
    ax.set_axis_off()
    ncol = len(col_labels)
    if col_widths is None:
        col_widths = [1.0 / ncol] * ncol
    x0 = [sum(col_widths[:i]) for i in range(ncol)]
    align = align or ["left"] * ncol
    n = len(rows)
    row_h = 1.0 / (n + 1.3)
    y_hdr = 1.0 - row_h
    hc = header_color or INK
    # header rules
    ax.plot([0, 1], [1.0, 1.0], color=INK, lw=MIN_LINE_PT, clip_on=False)
    ax.plot([0, 1], [y_hdr, y_hdr], color=INK, lw=MIN_LINE_PT, clip_on=False)
    for j, lab in enumerate(col_labels):
        ha = align[j]
        xx = x0[j] + (
            0.008 if ha == "left" else col_widths[j] - 0.008 if ha == "right" else col_widths[j] / 2
        )
        ax.text(
            xx,
            y_hdr + row_h * 0.32,
            lab,
            fontsize=fontsize,
            color=hc,
            ha=ha,
            va="center",
            fontweight="normal",
        )
    for i, row in enumerate(rows):
        y = y_hdr - (i + 1) * row_h
        if zebra and i % 2 == 1:
            ax.add_patch(
                mpl.patches.Rectangle(
                    (0, y), 1, row_h, facecolor="#f4f5f6", edgecolor="none", zorder=0
                )
            )
        if highlight_cols:
            # pale warm wash (Okabe-Ito orange blended ~86% toward white) — an elegant
            # column emphasis, not a heavy saturated bar.
            _o = mpl.colors.to_rgb(ORANGE)
            _tint = tuple(c + 0.86 * (1.0 - c) for c in _o)
            for j in highlight_cols:
                ax.add_patch(
                    mpl.patches.Rectangle(
                        (x0[j], y),
                        col_widths[j],
                        row_h,
                        facecolor=_tint,
                        edgecolor="none",
                        zorder=0,
                    )
                )
        for j, cell in enumerate(row):
            ha = align[j]
            xx = x0[j] + (
                0.008
                if ha == "left"
                else col_widths[j] - 0.008
                if ha == "right"
                else col_widths[j] / 2
            )
            ax.text(
                xx, y + row_h * 0.42, str(cell), fontsize=fontsize, color=INK, ha=ha, va="center"
            )
    _foot = y_hdr - n * row_h
    ax.plot([0, 1], [_foot, _foot], color=INK, lw=MIN_LINE_PT, clip_on=False)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return 1.0


def column_check(pdf_path, tol_mm=0.5):
    """Warn if a written panel is not at a Nature Communications column width.

    🔴 THE POINT IS *WHEN* THIS RUNS. compose.py already refuses a panel that is neither
    88 nor 180 mm, and that is the hard gate -- but it fires at composition, hours and a
    SLURM chain away from the builder that produced the wrong size, and it names the file
    rather than the code. Checking at WRITE time puts the complaint in the log of the job
    that caused it, next to the panel's own name.

    A warning, deliberately, not an exception. A builder mid-development legitimately
    writes panels at odd sizes, and turning that into a crash would push people to stop
    using fl.save() -- which is where the check lives. compose.py stays the gate.

    Returns (width_mm, height_mm, ok) so a caller can assert if it wants to.
    """
    try:
        from pypdf import PdfReader

        box = PdfReader(str(pdf_path)).pages[0].mediabox
    except Exception:  # noqa: BLE001 - a missing/odd PDF is the caller's problem, not this
        return None, None, True
    w_mm = float(box.width) / (72.0 / 25.4)
    h_mm = float(box.height) / (72.0 / 25.4)
    near = min((SINGLE_MM, DOUBLE_MM), key=lambda t: abs(t - w_mm))
    if abs(near - w_mm) <= tol_mm:
        return w_mm, h_mm, True
    warnings.warn(
        f"{os.path.basename(str(pdf_path))} is {w_mm:.2f} mm wide, not {SINGLE_MM:.0f} or "
        f"{DOUBLE_MM:.0f} mm (nearest {near:.0f}). compose.py will refuse it. If it is a "
        f"near miss, the cause is almost always a bbox that cropped the canvas -- REDRAW "
        f"at the placed width; do not scale and do not pad.",
        stacklevel=3,
    )
    return w_mm, h_mm, False


def clip_check(fig, name="figure", tol_mm=0.2):
    """Warn if anything drawn falls OUTSIDE the canvas and will be cut off.

    🔴 THIS IS THE OTHER HALF OF DROPPING `bbox_inches="tight"`, AND IT IS THE DANGEROUS
    HALF. Found 2026-08-11 during the layout pass, twice, independently.

    `tight` did two things, and only one of them was the defect. It CROPPED a canvas whose
    edges were empty -- that is the bug, it made panel widths depend on the data. But it
    also GREW the canvas to take in artists drawn outside the figure bounds, and plenty of
    builders had come to rely on that without knowing it. `fl_style.panel_letter()` draws
    at x = -0.16 in axes fraction, right-aligned; under matplotlib's default left margin
    that anchor sits 0.2 mm from the left edge of a 180 mm canvas, so the letter hangs off
    it and `tight` quietly pulled the edge out to include it.

    Remove the trim and those artists are simply CLIPPED AWAY. The first Fig3 rebuild came
    out at a correct 180.00 mm with NO PANEL LETTERS AT ALL, and nothing said so -- a
    MediaBox check passes either way, and so does every size check in this repository. It
    was caught by a human looking at a PNG.

    That is precisely the failure mode `canonical/CLAUDE.md` names as the one to ask about:
    not "did it error?" but "would anything have told me?". Nothing would have. So the
    check goes next to the save, where the loss happens.

    Exposed artists are not only panel letters: rotated tick labels, a legend placed with
    `bbox_to_anchor` outside the axes, an `annotate` with `annotation_clip=False`.

    THE FIX WHEN THIS FIRES IS NOT TO PUT THE TRIM BACK. Give the content room inside the
    authored canvas -- `layout="constrained"`, or explicit `fig.subplots_adjust(...)`. The
    canvas size is the specification; the drawing fits it, not the other way round.

    A warning, not an exception: a builder may legitimately clip a deliberate bleed. But it
    now has to be a decision rather than a surprise.
    """
    try:
        fig.canvas.draw()
        tb = fig.get_tightbbox(fig.canvas.get_renderer())
    except Exception:  # noqa: BLE001 - a backend without a renderer must not break saving
        return True
    w_in, h_in = fig.get_size_inches()
    tol = tol_mm / 25.4
    over = {
        "left": -tb.x0,
        "bottom": -tb.y0,
        "right": tb.x1 - w_in,
        "top": tb.y1 - h_in,
    }
    bad = {k: v * 25.4 for k, v in over.items() if v > tol}
    if not bad:
        return True
    warnings.warn(
        f"{name}: content extends past the canvas and WILL BE CLIPPED -- "
        + ", ".join(f"{k} by {v:.2f} mm" for k, v in sorted(bad.items()))
        + ". Panel letters drawn at negative axes coordinates are the usual cause. Give it "
        "room inside the authored size (layout='constrained', or subplots_adjust); do NOT "
        "restore bbox_inches='tight', which fixes this by changing the panel's width.",
        stacklevel=3,
    )
    return False


def normalise_type(fig, name="figure", verbose=True):
    """Bring one figure onto the house style, immediately before it is written.

    WHY AT SAVE TIME AND NOT IN THE BUILDERS. Thirty-odd builders across four trees draw
    the panels in this paper, written on different days, some of them hand-tuned. Setting
    rcParams fixes only what is drawn THROUGH the rcParams: an explicit `fontsize=9`, a
    `lw=0.6`, or matplotlib's own autoscaling of a tick label to make it fit, all sail past
    them. Measured 2026-08-11 across the 57 composed panels, that left:

      * type at 4.6 and 4.9 pt, below the 5 pt floor, in eight panels
      * largest type at 7 pt in the original panels, 8 pt in Fig1 and Fig4, 9 pt in Fig6
        and 10 pt in Fig2 -- so an axis title was 7 pt in Figure 6 panel c and 9 pt in
        panel g, in the same figure
      * lines at 0.6 pt against the >= 1 pt floor, where a builder passed lw= explicitly

    Fixing that builder by builder means finding all of them and keeping them in step
    forever. Doing it here means every panel that goes through `save()` is correct by
    construction, and a new builder cannot reintroduce the problem.

    WHAT IT DOES NOT DO. It does not shrink anything, and it does not touch the panel
    letter -- that is stamped by `compose.py` on the composed figure, at its own size, and
    is the one glyph allowed to be larger. Type between the floor and the title size is
    left exactly as authored: this is a floor and a ceiling, not a re-typesetting.
    """
    raised_type, capped_type, raised_line, fixed_ticks = 0, 0, 0, 0

    # 🔴 TICK LABELS FIRST, AND THROUGH tick_params -- NOT through the Text objects.
    # matplotlib REGENERATES tick label Text objects during draw, from the locator and the
    # axis's stored label size, so setting `.set_fontsize()` on the ones `findobj` returns
    # is thrown away and the panel is written at the old size. That is exactly what
    # happened on the first attempt at this function: it reported "0 labels raised" and the
    # re-rendered panels still carried 4.9 pt ticks. The 4.9 is 7 pt x 0.694, matplotlib's
    # `x-small` -- a relative size a builder asked for, which no absolute rcParam bounds.
    for ax in fig.get_axes():
        for axis, which in ((ax.xaxis, "x"), (ax.yaxis, "y")):
            labels = axis.get_ticklabels()
            if not labels:
                continue
            sz = labels[0].get_fontsize()
            if sz is None:
                continue
            # A tick label containing `$` is mathtext, and its superscript renders at
            # _MATHTEXT_RATIO of the label size. Bound the SUPERSCRIPT by the floor, not
            # the label -- bounding the label is what let 4.9 pt exponents through.
            is_math = any("$" in t.get_text() for t in labels)
            lo = MIN_TYPE_PT / _MATHTEXT_RATIO if is_math else MIN_TYPE_PT
            hi = LOG_TICK_PT if is_math else TITLE_PT
            if sz < lo:
                ax.tick_params(axis=which, labelsize=lo)
                fixed_ticks += 1
            elif sz > hi:
                ax.tick_params(axis=which, labelsize=hi)
                fixed_ticks += 1

    # 🔴 TICK LABELS ARE EXCLUDED HERE -- the loop above already owns them, and letting this
    # one see them UNDID it: a log axis clamped to 7.15 pt so its exponent would clear the
    # floor was immediately capped back to 7.0 as "above the house size", and the panel was
    # written with 4.9 pt exponents again. Two passes with overlapping scope, the second
    # reverting the first, and nothing in the report showed it -- it printed both "1 tick
    # axis clamped" and "16 capped to 7 pt" and looked like it had worked.
    _ticklabels = set()
    for ax in fig.get_axes():
        for axis in (ax.xaxis, ax.yaxis):
            _ticklabels.update(id(t) for t in axis.get_ticklabels(which="both"))

    for t in fig.findobj(match=lambda o: hasattr(o, "get_fontsize") and hasattr(o, "get_text")):
        if id(t) in _ticklabels:
            continue
        try:
            sz = t.get_fontsize()
        except Exception:                                          # noqa: BLE001
            continue
        if sz is None or not t.get_text():
            continue
        if sz < MIN_TYPE_PT:
            t.set_fontsize(MIN_TYPE_PT)
            raised_type += 1
        elif sz > TITLE_PT:
            t.set_fontsize(TITLE_PT)
            capped_type += 1

    for o in fig.findobj(match=lambda o: hasattr(o, "get_linewidth")):
        try:
            lw = o.get_linewidth()
        except Exception:                                          # noqa: BLE001
            continue
        # A collection returns a sequence; a zero width is deliberate (invisible) and is
        # left alone -- raising it to 1 pt would draw a line nobody asked for.
        if lw is None or isinstance(lw, (list, tuple)) or lw == 0:
            continue
        if 0 < lw < MIN_LINE_PT:
            o.set_linewidth(MIN_LINE_PT)
            raised_line += 1

    if verbose and (raised_type or capped_type or raised_line or fixed_ticks):
        print(
            f"  [style] {name}: {raised_type} label(s) raised to {MIN_TYPE_PT:.0f} pt, "
            f"{capped_type} capped to {TITLE_PT:.0f} pt, {fixed_ticks} tick axis/axes "
            f"clamped, {raised_line} line(s) raised to {MIN_LINE_PT:.0f} pt"
        )
    return raised_type, capped_type, raised_line, fixed_ticks


def save(fig, basepath, dpi=300, formats=("png", "pdf", "svg")):
    """Write a panel at its AUTHORED size, in every requested format.

    🔴 DOES NOT DELEGATE TO `_sp.save_figure` ANY MORE, and that is the whole change.
    The shared styleguide's `save_figure` hard-codes `bbox_inches="tight"` INSIDE the
    call, so it cannot be overridden by an rcParam or a keyword -- setting
    `savefig.bbox = None` in init() fixes every other save path in the project and does
    nothing here. Panels written through this function were trimmed to their content and
    came out at 175.81-179.22 mm instead of 180, and 83.53-86.90 instead of 88.

    `fig_layout.save_panel()` was added on 2026-08-05 as a way around this for panels
    drawn outside Canvas. It should not have needed to exist; both now do the same thing.

    Kept from the shared helper: white facecolor, and closing the figure.
    """
    enforce_nbt(fig)  # NBT: strip any stray on-canvas title before writing
    base = os.path.splitext(str(basepath))[0]
    # Type and line floors BEFORE the clip check, because raising a 4.6 pt label to 5 pt
    # makes it wider and can push it off the canvas -- so the clip check must see the text
    # at the size that will actually be written, not the size it was authored at.
    normalise_type(fig, os.path.basename(base))
    # BEFORE writing, not after: once the file exists the clipped artist is simply absent
    # and no amount of reading the file back can tell it was ever there.
    clip_check(fig, os.path.basename(base))
    for fmt in formats:
        fig.savefig(
            f"{base}.{fmt}",
            dpi=dpi,
            bbox_inches=None,  # explicit: the authored figsize IS the specification
            pad_inches=0.0,
            facecolor="white",
            edgecolor="none",
        )
    import matplotlib.pyplot as _plt

    _plt.close(fig)
    if "pdf" in formats:
        column_check(base + ".pdf")
    return base + ".pdf", base + ".png"
