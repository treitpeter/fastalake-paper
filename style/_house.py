"""In-repo equivalents of the three lab-styleguide symbols that fl_style.py uses.

The canonical fl_style.py imported ``color_palettes`` and ``shared_functions_plotting``
from a lab styleguide directory. Only three symbols were used. They are reproduced here
with the same behaviour so that the figure style needs no file outside this repository
except the licensed Söhne typeface (read from ``FASTALAKE_FONT_DIR``).
"""

import colorsys
import os

import matplotlib as mpl
import matplotlib.colors as mc
import matplotlib.font_manager as fm

#: Okabe-Ito four-colour palette: blue, orange, green, vermillion.
COLORBLIND_SAFE_4 = ["#0072B2", "#E69F00", "#009E73", "#D55E00"]

#: Ink colour for all text, axes and ticks.
INK = "#4e4f51"


def lighten_color(color, amount=0.5):
    """Lighten a colour in HLS space.

    Parameters
    ----------
    color : str or tuple
        Matplotlib colour name, hex string or RGB tuple.
    amount : float
        Lightness becomes ``1 - amount * (1 - L)``; smaller values give lighter colours.

    Returns
    -------
    tuple of float
        RGB triple in [0, 1].
    """
    c = mc.cnames.get(color, color) if isinstance(color, str) else color
    h, l, s = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(h, 1 - amount * (1 - l), s)


def init_plotting(font_path=None):
    """Set the house rcParams and register Söhne Buch when the font file exists.

    Parameters
    ----------
    font_path : str or None
        Path to ``soehne-buch.otf``. If the file is missing, the family is left at the
        matplotlib default and every other setting is still applied.

    Returns
    -------
    bool
        True if the Söhne face was registered.
    """
    mpl.rcParams["figure.dpi"] = 300
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["svg.fonttype"] = "none"
    have = bool(font_path) and os.path.isfile(font_path)
    if have:
        fe = fm.FontEntry(fname=font_path, name="soehne")
        fm.fontManager.ttflist.insert(0, fe)
        mpl.rcParams["font.family"] = fe.name
    else:
        mpl.rcParams["font.family"] = "sans-serif"
    for key in ("text.color", "axes.labelcolor", "axes.edgecolor", "xtick.color", "ytick.color"):
        mpl.rcParams[key] = INK
    for key in (
        "lines.linewidth",
        "axes.linewidth",
        "xtick.major.width",
        "ytick.major.width",
        "xtick.minor.width",
        "ytick.minor.width",
    ):
        mpl.rcParams[key] = 0.6
    return have
