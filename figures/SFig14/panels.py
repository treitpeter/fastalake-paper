"""Supplementary Figure 14: measured spectra for selected PREDICT histidine-boundary peptides.

``figure_v30`` draws the figure as placed from manuscript v30 on: the same values in the house
style (port of the v30 restyle script ``build_sf14_style.py`` and ``sf_common.py``).

Panels
------
a-d  Four selected PREDICT stool spectra from the non-specific Sage search (HIS_001 to
     HIS_004). Gray: processed spectrum; blue and orange: matched b and y ions. Intensity
     is normalised to the spectrum's maximum processed peak.

The original producer drew the four panels in one matplotlib figure with a figure title
(``HISTIDINE_MEASURED_SPECTRA.{pdf,svg,png}``, 180 dpi). A second step removed the in-image
title from the PDF and rasterised it at 300 dpi; that PNG is the image in the manuscript.
``figure_abcd`` reproduces the first step, ``remove_title`` the second. The drawing code is
moved unchanged from the producer; only paths are repository-relative.

The notebook also re-derives, from the shipped observed peaks, the processed peak lists and
the matched-fragment annotations of each case, and checks the case-selection rule against
the selection tables. These re-derivations port the upstream spectrum-extraction code
(analysis A25) without numba; they are checks, not new identifications.

Source data (``source_data/SFig14/``)
-------------------------------------
PREDICT_spectrum_case_review.tsv
    One row per case: acquisition, compartment, origin, boundary side, peptide, matched
    ions, backbone coverage, q values, bounded-alternative counts, parent protein.
spectra/HIS_00<k>/processed_peaks.tsv
    The plotted spectrum (gray peaks): 200-2000 m/z, at most 150 most intense peaks.
spectra/HIS_00<k>/observed_peaks.tsv.gz
    Every peak of the measured MS2 scan as read from the mzML.
spectra/HIS_00<k>/matched_fragments.tsv
    Matched b and y ions (blue and orange peaks, labels).
spectra/HIS_00<k>/bounded_sequence_alternatives.tsv
    Sequence and modification alternatives present in the local FASTA, with matched ions.
spectra/HIS_00<k>/PROVENANCE.json
    Raw file, FASTA and search configuration identity (file names and sha256) and the
    selected PSM row.
spectra/selected_cases.tsv, all_strata.tsv, complete_candidate_counts.tsv,
spectra/eligible_case_selection_psms.tsv.gz
    The case-selection tables (one stratum per origin x boundary side).
spectra/spectrum_audit.tsv
    Upstream per-case summary (matched ions, backbone coverage, alternatives).
spectra/same_scan_all_search_calls.tsv
    Every Sage PSM on the four selected scans, per search route (A25). Used to re-derive
    ``accepted_in_semi`` of the case table.
"""

import collections
import gzip
import itertools
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------- choices (named constants)
#: matplotlib settings of the producer (build_revision.py line 13), applied on top of the
#: matplotlib defaults. DejaVu Sans is the producer's font; Söhne is not used here.
PRODUCER_RC = {
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "pdf.fonttype": 42, "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 180,
}
#: Producer palette (build_revision.py line 14); b ions COL[0] (blue), y ions COL[1] (orange).
COL = ["#2878a3", "#df8f44", "#619b68", "#b35872", "#8277ac"]
#: Processed-peak colour and line widths (producer lines 283-285).
PEAK_COLOR = "#b5b5b5"
PEAK_LW = 0.6
ION_LW = 0.9
#: Matched ions are labelled when their intensity exceeds 8% of the maximum processed peak
#: (producer line 287).
LABEL_MIN_RELATIVE = 0.08
#: y-axis limit in percent of maximum processed peak (producer line 288).
YLIM = (0, 137)
#: Figure size in inches (producer line 279).
FIGSIZE = (7.1, 8.8)
#: Figure title written by the producer and removed by the title-redaction step.
PRODUCER_TITLE = "Measured PREDICT histidine-boundary spectra"
#: Label displacement in the producer's save(): step in data units, attempts, ceiling.
LABEL_STEP = 10
LABEL_ATTEMPTS = 12
LABEL_CEILING = 120
#: Rasterisation of the title-redacted PDF (prepare_supplement_figures.py line 17).
COMPOSITE_DPI = 300
#: Height (points) under which a text span counts as figure title
#: (inspect_package_inputs.py line 16, prepare_supplement_figures.py line 8).
TITLE_SPAN_MAX_Y = 35
TITLE_SPAN_ASSERT_Y = 20

#: Upstream spectrum-extraction settings. A25 search_revision2.py line 71-72 asserts them
#: against the Sage configuration of the non-specific search and applies them at line 80.
#: Peak window of the plotted spectrum, m/z. TEXT-GAP: the Methods state 200-2,000 only as
#: neutral fragment mass of the Sage index (v28d line 625), not as the peak window of the
#: figure (OPEN_ITEMS.md).
FRAGMENT_MZ = (200, 2000)
#: Most intense peaks kept. TEXT-GAP: Methods, v28d line 625, state 150 peaks only for the
#: healthy-stool production configuration; line 855 says the non-specific searches retained
#: the standard configuration (OPEN_ITEMS.md).
MAX_PEAKS = 150
#: Fragment tolerance, ppm. Methods, v28d lines 625 and 855 (+-20 ppm).
FRAGMENT_PPM = 20
#: Highest fragment charge considered; a precursor of charge z uses fragment charges 1 to
#: min(z - 1, FRAGMENT_CHARGE_MAX), at least 1. Sage max_fragment_charge (A25 line 71).
#: TEXT-GAP: Methods, v28d line 625, state "fragment charge at most 2" only for the
#: healthy-stool configuration (OPEN_ITEMS.md).
FRAGMENT_CHARGE_MAX = 2
#: Producer ppm check tolerance (build_revision.py line 289).
PPM_ATOL = 1e-6
#: Proton mass. The producer's ppm check uses this float64 literal (build_revision.py line
#: 289); Sage and the upstream code use it as float32 (A25 fragment_helpers.py).
PROTON_F64 = 1.0072764
PROTON = np.float32(PROTON_F64)
#: Same-scan acceptance gate of the case audit (audit_predict_spectrum_cases.py line 158):
#: target, peptide q and spectrum q at most 0.01. Methods, v28d line 855.
SAME_SCAN_Q = 0.01
MASS = np.array([71.03711, 0, 103.00919, 115.02694, 129.04259, 147.0684, 57.02146, 137.05891,
                 113.08406, 0, 128.09496, 113.08406, 131.0405, 114.04293, 237.14774, 97.05276,
                 128.05858, 156.1011, 87.03203, 101.04768, 150.95363, 99.06841, 186.07932, 0,
                 163.06332, 0], dtype=np.float32)
TOK = re.compile(r"([A-Z])(?:\[([+-]?[0-9.]+)\])?")
#: Case order and panel letters.
CASES = ["HIS_001", "HIS_002", "HIS_003", "HIS_004"]
LETTERS = "abcd"

WORDS = {1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five", 6: "Six"}


# ---------------------------------------------------------------- loading
def load_case_review(sd):
    """Read the case table the producer iterates over (one row per panel).

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig14``.

    Returns
    -------
    pandas.DataFrame
        Four rows, HIS_001 to HIS_004, in panel order.
    """
    return pd.read_csv(Path(sd) / "PREDICT_spectrum_case_review.tsv", sep="\t")


# ---------------------------------------------------------------- producer drawing code
def _title(ax, letter, text):
    """Producer panel title (build_revision.py line 59).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axis.
    letter : str
        Panel letter.
    text : str
        Title text after the letter.

    Returns
    -------
    None
    """
    ax.set_title(letter + "  " + text, loc="left", pad=9)


def panel_spectrum(ax, i, r, sd):
    """Draw one spectrum panel exactly as the producer did (build_revision.py lines 280-288).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axis.
    i : int
        Panel index (0 = a).
    r : pandas Series
        Case-review row.
    sd : pathlib.Path
        ``source_data/SFig14``.

    Returns
    -------
    None
    """
    cid = r.case_id.replace("PREDICT_", "")
    folder = Path(sd) / "spectra" / cid
    peaks = pd.read_csv(folder / "processed_peaks.tsv", sep="\t")
    frag = pd.read_csv(folder / "matched_fragments.tsv", sep="\t")
    norm = peaks.intensity.max()
    ax.vlines(peaks.mz, 0, peaks.intensity / norm * 100, color=PEAK_COLOR, lw=PEAK_LW)
    for series, g in frag.groupby("series"):
        ax.vlines(g.observed_mz, 0, g.intensity / norm * 100, color=COL[0] if series == "b" else COL[1],
                  lw=ION_LW, label=series + " ions")
        for j, x in enumerate(g.itertuples()):
            if x.intensity / norm > LABEL_MIN_RELATIVE:
                ax.text(x.observed_mz, min(105, x.intensity / norm * 100 + 4 + j % 2 * 6),
                        f"{series}{x.ordinal}" + ("²⁺" if x.charge == 2 else ""),
                        rotation=90, fontsize=7, ha="center", va="bottom")
    ax.set(xlabel="m/z", ylabel="Relative intensity (%)", ylim=YLIM)
    _title(ax, LETTERS[i], cid + "  " + r.peptide + "  (" + r.origin + ", " + r.side + " boundary)")
    ax.legend(frameon=False, loc="upper right", fontsize=7)


def _producer_save(fig, stem):
    """Producer ``save()`` for HISTIDINE_MEASURED_SPECTRA (build_revision.py lines 34-52).

    Adds the figure title, runs tight_layout, moves overlapping ion labels upwards in steps
    of ``LABEL_STEP`` with a gray leader line, and writes pdf, svg and png.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The four-panel figure.
    stem : pathlib.Path
        Output path without suffix.

    Returns
    -------
    None
    """
    fig.suptitle(PRODUCER_TITLE, fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, .96), pad=1.6)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for axis in fig.axes:
        placed = []
        labels = sorted([t for t in axis.texts if re.fullmatch(r"[by]\d+(?:²⁺)?", t.get_text())],
                        key=lambda t: t.get_position())
        for label in labels:
            x, y = label.get_position()
            initial = y
            box = label.get_window_extent(renderer).expanded(1.2, 1.1)
            for attempt in range(LABEL_ATTEMPTS):
                if not any(box.overlaps(b) for b in placed):
                    break
                y += LABEL_STEP
                label.set_position((x, y))
                box = label.get_window_extent(renderer).expanded(1.2, 1.1)
            assert not any(box.overlaps(b) for b in placed), label.get_text()
            assert y < LABEL_CEILING, (label.get_text(), y)
            if y != initial:
                axis.plot([x, x], [initial, y - 2], color=".5", lw=.35)
            placed.append(box)
    for ext in ["pdf", "svg", "png"]:
        fig.savefig(str(stem) + "." + ext)
    plt.close(fig)


def figure_abcd(sd, out_stem):
    """Render panels a-d as one figure, as the producer did, and save pdf/svg/png.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig14``.
    out_stem : pathlib.Path
        Output path without suffix (``build/SFig14/panel_abcd``).

    Returns
    -------
    pathlib.Path
        The png written.
    """
    case = load_case_review(sd)
    with matplotlib.rc_context():
        matplotlib.rcdefaults()
        plt.rcParams.update(PRODUCER_RC)
        fig, axes = plt.subplots(4, 1, figsize=FIGSIZE)
        for i, (_, r) in enumerate(case.iterrows()):
            panel_spectrum(axes[i], i, r, sd)
        _producer_save(fig, out_stem)
    return Path(str(out_stem) + ".png")


def remove_title(pdf_in, pdf_out, png_out, dpi=COMPOSITE_DPI):
    """Remove the in-image figure title and rasterise, as the release step did.

    Ports ``inspect_package_inputs.py`` line 16 (title span = text span with top below
    ``TITLE_SPAN_MAX_Y`` pt) and ``prepare_supplement_figures.py`` lines 7-17 (white
    redaction of that span, check that no numeric label changed, save with garbage=4 and
    deflate, rasterise at 300 dpi).

    Parameters
    ----------
    pdf_in : pathlib.Path
        Producer PDF with title.
    pdf_out, png_out : pathlib.Path
        Composite PDF and PNG.
    dpi : int
        Rasterisation resolution.

    Returns
    -------
    dict
        Removed title text and its bounding box.
    """
    import fitz

    doc = fitz.open(pdf_in)
    p = doc[0]
    before = p.get_text()
    top = [{"text": s["text"], "bbox": s["bbox"], "size": s["size"]}
           for b in p.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]
           if s["bbox"][1] < TITLE_SPAN_MAX_Y]
    title = top[0]
    assert title["bbox"][1] < TITLE_SPAN_ASSERT_Y
    assert title["text"] == PRODUCER_TITLE, title["text"]
    p.add_redact_annot(fitz.Rect(title["bbox"]), fill=(1, 1, 1))
    p.apply_redactions(images=0, graphics=0)
    assert title["text"] not in p.get_text()
    before_nums = collections.Counter(re.findall(r"\d[\d,.]*", before.replace(title["text"], "")))
    after_nums = collections.Counter(re.findall(r"\d[\d,.]*", p.get_text()))
    assert before_nums == after_nums, "numeric label changed"
    doc.save(pdf_out, garbage=4, deflate=True)
    q = fitz.open(pdf_out)[0]
    q.get_pixmap(dpi=dpi, alpha=False).save(png_out)
    return {"removed_title": title["text"], "bbox": [round(v, 3) for v in title["bbox"]]}


# ---------------------------------------------------------------- statistics and checks
def stat_n_cases(sd):
    """Number of spectra shown (legend, v28d line 465: "Four selected PREDICT spectra").

    Statistical unit: PSM (one selected spectrum per case); n = 4 cases from 4 PREDICT
    stool acquisitions.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig14``.

    Returns
    -------
    dict
        ``n`` (int), ``word`` (the count written as a word), ``n_spectrum_folders``,
        ``n_selected`` (rows in selected_cases.tsv).
    """
    case = load_case_review(sd)
    folders = sorted(p.name for p in (Path(sd) / "spectra").glob("HIS_*") if p.is_dir())
    sel = pd.read_csv(Path(sd) / "spectra" / "selected_cases.tsv", sep="\t")
    n = int(case.case_id.nunique())
    return {"n": n, "word": WORDS.get(n, str(n)), "n_spectrum_folders": len(folders),
            "n_selected": int(len(sel)), "case_ids": list(case.case_id)}


def check_fragment_ppm(sd):
    """Producer check (build_revision.py line 289): fragment ppm column equals
    (observed - calculated) / (calculated - proton) x 1e6 for every matched ion.

    Statistical unit: matched fragment ion; n = all rows of the four matched_fragments.tsv.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig14``.

    Returns
    -------
    dict
        Case id to (n ions, passed).
    """
    res = {}
    for cid in CASES:
        frag = pd.read_csv(Path(sd) / "spectra" / cid / "matched_fragments.tsv", sep="\t")
        ok = np.allclose((frag.observed_mz - frag.calculated_mz) / (frag.calculated_mz - PROTON_F64) * 1e6,
                         frag.ppm, atol=PPM_ATOL)
        res[cid] = (int(len(frag)), bool(ok))
    return res


def _parse(s):
    """Residue masses of a Sage peptide string (A25 fragment_helpers.parse, unchanged).

    Parameters
    ----------
    s : str
        Peptide as written by Sage, modifications in brackets (Da).

    Returns
    -------
    tuple of (str, numpy.ndarray)
        Unmodified sequence and float32 residue masses including modifications.
    """
    t = TOK.findall(s)
    assert "".join(a + ("[" + m + "]" if m else "") for a, m in t) == s, s
    return "".join(a for a, m in t), np.array(
        [np.float32(MASS[ord(a) - 65] + np.float32(m or 0)) for a, m in t], dtype=np.float32)


def _match(masses, ints, residues, total, charge):
    """A25 fragment_helpers.match without the numba decorator (same float32 arithmetic).

    Parameters
    ----------
    masses : numpy.ndarray
        Processed peak m/z minus proton (float32, sorted).
    ints : numpy.ndarray
        Processed peak intensities.
    residues : numpy.ndarray
        Residue masses from ``_parse``.
    total : numpy.float32
        Neutral peptide mass (Sage calcmass).
    charge : int
        Precursor charge.

    Returns
    -------
    tuple of int
        b-ion bond mask, y-ion bond mask, matched ion count, eligible ion count.
    """
    n = len(residues)
    b = np.float32(0)
    y = np.float32(total)
    bmask = 0
    ymask = 0
    matched = 0
    eligible = 0
    for i in range(n - 1):
        b = np.float32(b + residues[i])
        y = np.float32(y - residues[i])
        for series in range(2):
            neutral = b if series == 0 else y
            for z in range(1, max(2, min(charge, FRAGMENT_CHARGE_MAX + 1))):
                center = np.float32(neutral / np.float32(z))
                delta = np.float32(np.float32(center * np.float32(FRAGMENT_PPM)) / np.float32(1000000))
                lo = np.float32(center - delta)
                hi = np.float32(center + delta)
                if center + PROTON >= FRAGMENT_MZ[0] and center + PROTON <= FRAGMENT_MZ[1]:
                    eligible += 1
                index = int(np.searchsorted(masses, lo))
                found = False
                while index < len(masses) and masses[index] <= hi:
                    if ints[index] >= 0:
                        found = True
                    index += 1
                if found:
                    matched += 1
                    if series == 0:
                        bmask |= 1 << i
                    else:
                        ymask |= 1 << i
    return bmask, ymask, matched, eligible


def stat_rederive_spectrum(sd, cid):
    """Re-derive one case's processed peaks and matched fragments from its observed peaks.

    Ports A25 search_revision2.py (processed peaks: 200-2000 m/z, 150 most intense, sorted
    by m/z; annotation loop; backbone coverage = bonds covered by any b or y ion / bonds).
    The upstream run used NumPy 1.x, where a float32 ratio times a Python float is
    promoted to float64; under NumPy 2 it stays float32. The ppm line therefore converts
    the float32 ratio to float before multiplying, which reproduces the upstream value.
    Shipped tables are read with ``float_precision="round_trip"`` so that the comparison
    is exact.
    Statistical unit: peak (processed) and fragment ion (matched); one spectrum.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig14``.
    cid : str
        Case id.

    Returns
    -------
    dict
        Counts and booleans: processed peaks equal, annotations equal, matched ions,
        Sage matched_peaks, backbone coverage.
    """
    folder = Path(sd) / "spectra" / cid
    sel = pd.read_csv(Path(sd) / "spectra" / "selected_cases.tsv", sep="\t").set_index("case_id")
    row = sel.loc[cid]
    with gzip.open(folder / "observed_peaks.tsv.gz", "rt") as fh:
        obs = pd.read_csv(fh, sep="\t")
    mz = obs.mz.to_numpy().astype(np.float32)
    intensity = obs.intensity.to_numpy().astype(np.float32)
    k = (mz >= FRAGMENT_MZ[0]) & (mz <= FRAGMENT_MZ[1])
    ma = mz[k]
    it = intensity[k]
    ix = np.argsort(it, kind="stable")[-MAX_PEAKS:]
    ma = ma[ix]
    it = it[ix]
    ix = np.argsort(ma)
    ma = ma[ix]
    it = it[ix]
    shipped = pd.read_csv(folder / "processed_peaks.tsv", sep="\t")
    proc_equal = (len(shipped) == len(ma)
                  and np.array_equal(shipped.mz.to_numpy().astype(np.float32), ma)
                  and np.array_equal(shipped.intensity.to_numpy().astype(np.float32), it))
    masses = np.asarray(ma - PROTON, np.float32)
    seq, res = _parse(row["peptide"])
    total = np.float32(row["calcmass"])
    charge = int(row["charge"])
    bm, ym, n, eligible = _match(masses, it, res, total, charge)
    b = np.float32(0)
    y = total
    ann = []
    for i in range(len(res) - 1):
        b = np.float32(b + res[i])
        y = np.float32(y - res[i])
        for series, neutral, ordinal in [("b", b, i + 1), ("y", y, len(res) - i - 1)]:
            for z in range(1, max(2, min(charge, FRAGMENT_CHARGE_MAX + 1))):
                center = np.float32(neutral / np.float32(z))
                delta = np.float32(center * np.float32(FRAGMENT_PPM) / np.float32(1000000))
                hit = np.flatnonzero((masses >= np.float32(center - delta)) & (masses <= np.float32(center + delta)))
                if len(hit):
                    j = hit[np.argmax(it[hit])]
                    ann.append(dict(series=series, ordinal=ordinal, charge=z, observed_mz=float(ma[j]),
                                    calculated_mz=float(center + PROTON), intensity=float(it[j]),
                                    ppm=float((masses[j] - center) / center) * 1e6))
    ann = pd.DataFrame(ann)
    frag = pd.read_csv(folder / "matched_fragments.tsv", sep="\t", float_precision="round_trip")
    cols = ["series", "ordinal", "charge", "observed_mz", "calculated_mz", "intensity", "ppm"]
    ann_equal = (len(ann) == len(frag)
                 and ann[cols[:3]].reset_index(drop=True).equals(frag[cols[:3]].reset_index(drop=True))
                 and all(np.array_equal(ann[c].to_numpy(), frag[c].to_numpy()) for c in cols[3:]))
    return {"case_id": cid, "observed_peaks": int(len(obs)), "processed_peaks": int(len(ma)),
            "processed_equal": bool(proc_equal), "matched_ions": int(n), "annotations": int(len(ann)),
            "annotations_equal": bool(ann_equal), "sage_matched_peaks": int(row["matched_peaks"]),
            "backbone_coverage": (bm | ym).bit_count() / (len(seq) - 1), "bonds": len(seq) - 1,
            "eligible_ions": int(eligible)}


def stat_selection(sd):
    """Check the predeclared case-selection rule against the selection tables.

    Rule (A25 search_revision2.py): per stratum (origin x boundary side, PREDICT stool),
    sort acquisitions with at least one candidate peptide by (candidate count, name); take
    the one at index len // 2, falling back to the nearest index only if it has no PSM
    passing both q gates; within it keep the best PSM per canonical peptide, sort by
    (matched_intensity_pct, canonical_peptide, scannr) and take index len // 2.
    Statistical unit: acquisition (choice of acquisition), PSM (choice of spectrum).

    Returns
    -------
    pandas.DataFrame
        One row per case: number of acquisitions with candidates, median index, chosen
        index, acquisition agrees, eligible PSMs, chosen PSM index, PSM agrees.
    """
    sp = Path(sd) / "spectra"
    counts = pd.read_csv(sp / "complete_candidate_counts.tsv", sep="\t")
    sel = pd.read_csv(sp / "selected_cases.tsv", sep="\t")
    elig = pd.read_csv(sp / "eligible_case_selection_psms.tsv.gz", sep="\t")
    strata = pd.read_csv(sp / "all_strata.tsv", sep="\t")
    rows = []
    for _, r in sel.iterrows():
        c = counts[counts.cohort.eq(r.cohort) & counts.compartment.eq(r.compartment)
                   & counts.origin.eq(r.origin) & counts.side.eq(r.side) & counts.candidates.gt(0)]
        choices = sorted(zip(c.candidates, c["sample"]))
        midpoint = len(choices) // 2
        order = sorted(range(len(choices)), key=lambda i: (abs(i - midpoint), i))
        acq_ok = (r.acquisition_median_rank == midpoint and choices[r.acquisition_rank][1] == r["sample"])
        skipped = order.index(r.acquisition_rank)
        g = elig[elig.case_id.eq(r.case_id)].copy()
        g2 = g.sort_values(["matched_intensity_pct", "canonical_peptide", "scannr"]).reset_index(drop=True)
        psm_ok = (list(g2.selection_rank) == sorted(g.selection_rank)
                  and g.canonical_peptide.is_unique and len(g) == r.eligible_peptides
                  and r.selected_rank == len(g) // 2
                  and int(g2.loc[len(g) // 2, "psm_id"]) == int(r.psm_id)
                  and bool(g2.loc[len(g) // 2, "selected"]))
        st = strata[strata.case_id.eq(r.case_id)]
        rows.append(dict(case_id=r.case_id, origin=r.origin, side=r.side,
                         acquisitions_with_candidates=len(choices), median_index=midpoint,
                         chosen_index=int(r.acquisition_rank), fallback_steps=skipped,
                         candidates_in_chosen=int(choices[r.acquisition_rank][0]),
                         acquisition_rule_holds=bool(acq_ok), eligible_psms=int(len(g)),
                         chosen_psm_index=int(r.selected_rank), psm_rule_holds=bool(psm_ok),
                         stratum_status=st.status.iloc[0] if len(st) else "absent"))
    return pd.DataFrame(rows)


def _canonical(s):
    """Peptide without modifications, I replaced by L (audit_predict_spectrum_cases.py line 49).

    Parameters
    ----------
    s : str
        Peptide as written by Sage.

    Returns
    -------
    str
        Canonical sequence.
    """
    return re.sub(r"\[[^\]]*\]", "", s).replace("I", "L")


def stat_case_review(sd, rederived):
    """Compare the case-review table with the upstream audit and the re-derivation.

    Also re-derives the two flag columns of the case table with the rules of the producer,
    audit_predict_spectrum_cases.py lines 158 and 160: ``accepted_in_semi`` is True when a
    same-scan PSM from the route ``semi`` is a target (label 1) with peptide q and spectrum q
    at most ``SAME_SCAN_Q`` and the same canonical peptide; ``N_flanking_H_inferred`` is
    ``side == "N"``. The PREDICT tasks of A25 had only the ``tryptic`` route besides the
    non-specific one (manifest_revision2.json), so the same-scan table has no ``semi`` rows
    and ``accepted_in_semi`` is False for all four cases by construction.

    Statistical unit: case (spectrum); n = 4.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig14``.
    rederived : list of dict
        Output of ``stat_rederive_spectrum`` per case.

    Returns
    -------
    pandas.DataFrame
        One row per case with the review, audit and recomputed matched ions, backbone
        coverage, bounded-alternative counts and the two flag columns.
    """
    case = load_case_review(sd).set_index("case_id")
    calls = pd.read_csv(Path(sd) / "spectra" / "same_scan_all_search_calls.tsv", sep="\t")
    audit = pd.read_csv(Path(sd) / "spectra" / "spectrum_audit.tsv", sep="\t").set_index("case_id")
    rd = {d["case_id"]: d for d in rederived}
    rows = []
    for cid in CASES:
        alt = pd.read_csv(Path(sd) / "spectra" / cid / "bounded_sequence_alternatives.tsv", sep="\t")
        n_orig = int(alt.loc[alt.original, "matched_ions"].iloc[0])
        nonorig = alt[~alt.original]
        c = calls[calls.case_id.eq(cid)]
        same = (c.route.eq("semi") & c.label.eq(1) & c.peptide_q.le(SAME_SCAN_Q)
                & c.spectrum_q.le(SAME_SCAN_Q)
                & c.peptide.map(_canonical).eq(_canonical(case.loc[cid, "peptide"])))
        rows.append(dict(case_id=cid,
                         matched_ions_review=int(case.loc[cid, "matched_ions"]),
                         matched_ions_audit=int(audit.loc[cid, "matched_ions"]),
                         matched_ions_recomputed=rd[cid]["matched_ions"],
                         matched_ions_alternatives_table=n_orig,
                         backbone_review=float(case.loc[cid, "backbone_coverage"]),
                         backbone_recomputed=rd[cid]["backbone_coverage"],
                         alternatives_rows=int(len(alt)),
                         nonoriginal_review=int(case.loc[cid, "nonoriginal_alternatives"]),
                         nonoriginal_recomputed=int(len(nonorig)),
                         nonoriginal_at_least_original_review=int(case.loc[cid, "nonoriginal_matching_at_least_original"]),
                         nonoriginal_at_least_original_recomputed=int(nonorig.at_least_original_matches.sum()),
                         same_scan_routes=";".join(sorted(calls.loc[calls.case_id.eq(cid), "route"].unique())),
                         accepted_in_semi_review=bool(case.loc[cid, "accepted_in_semi"]),
                         accepted_in_semi_recomputed=bool(same.any()),
                         N_flanking_H_inferred_review=bool(case.loc[cid, "N_flanking_H_inferred"]),
                         N_flanking_H_inferred_recomputed=bool(case.loc[cid, "side"] == "N")))
    return pd.DataFrame(rows)


# ============================================================================ v30 house style
# Ported from the v30 restyle helpers (LAB_ROOT/.../final/v30/figures_style/SFig1x/sf_common.py):
# fingerprint/capture/compare prove the restyled figure plots the same values as the producer;
# _v30_finish is the fl_style save sequence. Logic unchanged; fl_style comes from REPO/style.
V30_LETTER_PT = 7.0        # sf_common.LETTER_PT: house panel letters, 7 pt
V30_DPI = 300              # PNG = the PDF rasterised at 300 dpi (sf_common.finish)


def _v30_fl():
    """Return the repository fl_style module with the house rcParams applied (sf_common.house).

    Returns
    -------
    module
        ``style.fl_style`` after ``matplotlib.rcdefaults()`` and ``fl.init()``.
    """
    from style import fl_style as fl
    import matplotlib as _mpl
    _mpl.rcdefaults()
    fl.init()
    return fl


def _v30_is_cbar(ax):
    return hasattr(ax, '_colorbar') or ax.get_label() == '<colorbar>'


def _v30_arr(a):
    a = np.ma.asarray(a)
    return np.ma.filled(a.astype(float), np.nan)


def v30_fingerprint(fig, skip_line=lambda l: False, text_norm=lambda s: s, data_text=lambda t: True):
    """Every plotted array of a figure (sf_common.fingerprint).

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Drawn figure.
    skip_line, text_norm, data_text : callable
        Lines to ignore, text normaliser, and which axes texts are data labels.

    Returns
    -------
    list of dict
        One dict per non-colour-bar axes: limits, in-view ticks and tick labels, images,
        collection offsets, line-collection segments, lines and data texts.
    """
    from matplotlib.collections import LineCollection
    out = []
    for ax in fig.axes:
        if _v30_is_cbar(ax):
            continue
        xl, yl = sorted(ax.get_xlim()), sorted(ax.get_ylim())
        xin = [k for k, v in enumerate(ax.get_xticks()) if xl[0] - 1e-9 <= v <= xl[1] + 1e-9]
        yin = [k for k, v in enumerate(ax.get_yticks()) if yl[0] - 1e-9 <= v <= yl[1] + 1e-9]
        xt, yt = ax.get_xticklabels(), ax.get_yticklabels()
        d = {'xlim': list(ax.get_xlim()), 'ylim': list(ax.get_ylim()),
             'xticks': [float(ax.get_xticks()[k]) for k in xin], 'yticks': [float(ax.get_yticks()[k]) for k in yin],
             'xticklabels': [xt[k].get_text() for k in xin],
             'yticklabels': [yt[k].get_text() for k in yin],
             'images': [], 'offsets': [], 'segments': [], 'lines': [], 'texts': []}
        for im in ax.images:
            d['images'].append({'array': _v30_arr(im.get_array()), 'clim': list(im.get_clim()),
                                'extent': list(im.get_extent())})
        for c in ax.collections:
            if isinstance(c, LineCollection):
                d['segments'].append([np.asarray(s, float) for s in c.get_segments()])
            else:
                d['offsets'].append(np.asarray(c.get_offsets(), float))
        for l in ax.lines:
            if skip_line(l):
                continue
            d['lines'].append(np.column_stack([np.asarray(l.get_xdata(), float), np.asarray(l.get_ydata(), float)]))
        for t in ax.texts:
            if data_text(t):
                d['texts'].append((text_norm(t.get_text()), float(t.get_position()[0])))
        d['texts'].sort()
        out.append(d)
    return out


def v30_capture(run, **fpkw):
    """Run ``run()`` with ``Figure.savefig`` intercepted; return each figure's fingerprint.

    Parameters
    ----------
    run : callable
        Draws and saves one or more figures (the producer).
    **fpkw
        Passed to ``v30_fingerprint``.

    Returns
    -------
    list
        One fingerprint per captured figure; nothing is written by the intercepted saves.
    """
    from matplotlib.figure import Figure
    caught = []
    orig = Figure.savefig

    def fake(self, *a, **k):
        self.canvas.draw()
        if not any(c[0] is self for c in caught):
            caught.append((self, v30_fingerprint(self, **fpkw)))
    Figure.savefig = fake
    try:
        run()
    finally:
        Figure.savefig = orig
    return [c[1] for c in caught]


def _v30_eq(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    return a.shape == b.shape and np.array_equal(a, b, equal_nan=True)


def v30_compare(ref, new, ignore=()):
    """Exact comparison of two fingerprints (sf_common.compare); raises on any difference.

    Parameters
    ----------
    ref, new : list of dict
        Producer and restyled fingerprints.
    ignore : tuple of str
        Keys not compared.

    Returns
    -------
    list of str
        One report line per axes, ending ``IDENTICAL``.
    """
    rep = []
    assert len(ref) == len(new), ('data axes', len(ref), len(new))
    for i, (r, n) in enumerate(zip(ref, new)):
        for k in ['xlim', 'ylim']:
            if k in ignore:
                continue
            assert _v30_eq(r[k], n[k]), (i, k, r[k], n[k])
        for k in ['xticks', 'yticks', 'xticklabels', 'yticklabels']:
            if k in ignore:
                continue
            assert r[k] == n[k], (i, k, r[k], n[k])
        assert len(r['images']) == len(n['images']), (i, 'images')
        for a, b in zip(r['images'], n['images']):
            assert _v30_eq(a['array'], b['array']) and _v30_eq(a['clim'], b['clim']) and _v30_eq(a['extent'], b['extent']), (i, 'image')
        assert len(r['offsets']) == len(n['offsets']), (i, 'n offsets', len(r['offsets']), len(n['offsets']))
        for a, b in zip(r['offsets'], n['offsets']):
            assert _v30_eq(a, b), (i, 'offsets')
        assert len(r['segments']) == len(n['segments']), (i, 'n segment sets')
        for a, b in zip(r['segments'], n['segments']):
            assert len(a) == len(b) and all(_v30_eq(x, y) for x, y in zip(a, b)), (i, 'segments')
        assert len(r['lines']) == len(n['lines']), (i, 'n lines', len(r['lines']), len(n['lines']))
        for a, b in zip(r['lines'], n['lines']):
            assert _v30_eq(a, b), (i, 'line')
        if 'texts' not in ignore:
            assert r['texts'] == n['texts'], (i, 'texts', r['texts'][:5], n['texts'][:5])
        npts = sum(x['array'].size for x in r['images']) + sum(len(x) for x in r['offsets']) \
            + sum(len(s) for x in r['segments'] for s in x) + sum(len(x) for x in r['lines'])
        rep.append(f"axes {i}: images {len(r['images'])}, point sets {len(r['offsets'])}, segment sets "
                   f"{len(r['segments'])}, lines {len(r['lines'])}, data texts {len(r['texts'])}, "
                   f"values compared {npts}: IDENTICAL")
    return rep


def _v30_finish(fl, fig, base, letters, keep_sizes=()):
    """fl_style save sequence of the v30 rebuild (sf_common.finish).

    enforce_nbt, normalise_type, panel letters re-set to ``V30_LETTER_PT``, optional text
    sizes restored, clip_check, PDF save, then the PDF rasterised at ``V30_DPI`` with PyMuPDF.

    Parameters
    ----------
    fl : module
        ``style.fl_style``.
    fig : matplotlib.figure.Figure
        Restyled figure.
    base : pathlib.Path
        Output path without suffix; writes ``.pdf`` and ``.png``.
    letters : list
        Axes carrying panel letters.
    keep_sizes : sequence
        Text artists whose font size must survive ``normalise_type``.

    Returns
    -------
    bool
        ``clip_check`` result.
    """
    import fitz
    base = str(base)
    fl.enforce_nbt(fig)
    sizes = [(t, t.get_fontsize()) for t in keep_sizes]
    fl.normalise_type(fig, Path(base).name, verbose=False)
    for ax in letters:
        ax._left_title.set_fontsize(V30_LETTER_PT)
    for t, s in sizes:
        t.set_fontsize(s)
    ok = fl.clip_check(fig, Path(base).name)
    fig.savefig(f'{base}.pdf', dpi=V30_DPI, bbox_inches=None, pad_inches=0, facecolor='white', edgecolor='none')
    plt.close(fig)
    pix = fitz.open(f'{base}.pdf')[0].get_pixmap(dpi=V30_DPI, alpha=False); pix.set_dpi(V30_DPI, V30_DPI); pix.save(f'{base}.png')
    return ok


# v30 house-style figure (LAB_ROOT/.../final/v30/figures_style/SFig14/build_sf14_style.py)
V30_SIZE_MM = (160, 198.3)       # placed size in the v29b docx
V30_PEAK_COLOR = '#b5b5b5'       # processed peaks, neutral grey
V30_ION_PT, V30_ION2_PT = 7.0, 7.2   # 7.2 pt so the mathtext '2+' superscript is >= 5 pt
V30_H_PAD, V30_W_PAD = 5 / 72, 2 / 72    # constrained-layout pads (inches)
V30_LEADER_COLOR = '.5'


def _v30_norm(s):
    return s.replace('$^{2+}$', '²⁺')


def _v30_is_ion(t):
    return re.fullmatch(r'[by]\d+(?:²⁺)?', _v30_norm(t.get_text())) is not None


def _v30_skip_leader(l):
    return l.get_gid() == 'leader' or (abs(l.get_linewidth() - .35) < 1e-9)


V30_FP = dict(skip_line=_v30_skip_leader, text_norm=_v30_norm, data_text=_v30_is_ion)


def figure_v30(sd, out_stem, ref_stem):
    """Draw Supplementary Figure 14 in the v30 house style and prove its values unchanged.

    Port of ``build_sf14_style.py``: the producer ``figure_abcd`` is run with savefig
    intercepted (reference fingerprint), then the same spectra are drawn with fl_style
    (Söhne, 160 x 198.3 mm, b ions blue, y ions orange, peaks grey, panel titles moved into
    the x-axis label), the producer's label-displacement rule is applied on the new layout,
    and every plotted array is compared exactly with the reference.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig14``.
    out_stem : pathlib.Path
        Output path without suffix; ``.pdf`` and 300 dpi ``.png`` are written.
    ref_stem : pathlib.Path
        Stem passed to the intercepted producer (nothing is written there).

    Returns
    -------
    dict
        ``panels`` (comparison report lines), ``labels_displaced``, ``clip_ok``.
        Statistical unit: spectrum (one PSM per panel), n = 4; no test.
    """
    ref = v30_capture(lambda: figure_abcd(sd, ref_stem), **V30_FP)
    assert len(ref) == 1
    ref = ref[0]
    with matplotlib.rc_context():
        fl = _v30_fl(); INK = fl.INK
        ION_COL = {'b': fl.BLUE, 'y': fl.ORANGE}
        case = load_case_review(sd)
        fig, axes = plt.subplots(4, 1, figsize=fl.mm(*V30_SIZE_MM), layout='constrained')
        fig.get_layout_engine().set(h_pad=V30_H_PAD, w_pad=V30_W_PAD)
        ion_texts = []
        for i, (_, r) in enumerate(case.iterrows()):
            ax = axes[i]; fl.style(ax)
            cid = r.case_id.replace('PREDICT_', ''); folder = sd / 'spectra' / cid
            peaks = pd.read_csv(folder / 'processed_peaks.tsv', sep='\t'); frag = pd.read_csv(folder / 'matched_fragments.tsv', sep='\t')
            normv = peaks.intensity.max()
            ax.vlines(peaks.mz, 0, peaks.intensity / normv * 100, color=V30_PEAK_COLOR, lw=fl.MIN_LINE_PT)
            for series, g in frag.groupby('series'):
                ax.vlines(g.observed_mz, 0, g.intensity / normv * 100, color=ION_COL[series], lw=fl.MIN_LINE_PT, label=series + ' ions')
                for j, x in enumerate(g.itertuples()):
                    if x.intensity / normv > LABEL_MIN_RELATIVE:
                        two = x.charge == 2
                        t = ax.text(x.observed_mz, min(105, x.intensity / normv * 100 + 4 + j % 2 * 6),
                                    f'{series}{x.ordinal}' + ('$^{2+}$' if two else ''), rotation=90,
                                    fontsize=V30_ION2_PT if two else V30_ION_PT, ha='center', va='bottom', color=INK)
                        ion_texts.append(t)
            # v29 panel title content (case ID, peptide, origin, boundary) kept as x-axis label text
            ax.set(xlabel='m/z; ' + cid + ' ' + r.peptide + ' (' + r.origin + ', ' + r.side + ' boundary)',
                   ylabel='Relative intensity (%)', ylim=YLIM)
            ax.legend(frameon=False, loc='upper right', fontsize=7)
            fl.panel(ax, 'abcd'[i])
        # producer label displacement (_producer_save), unchanged rule, on the house layout
        fig.canvas.draw(); fig.set_layout_engine('none')
        renderer = fig.canvas.get_renderer()
        moved = 0
        for axis in axes:
            placed = []
            labels = sorted([t for t in axis.texts if _v30_is_ion(t)], key=lambda t: t.get_position())
            for label in labels:
                x, y = label.get_position(); initial = y
                box = label.get_window_extent(renderer).expanded(1.2, 1.1)
                for attempt in range(LABEL_ATTEMPTS):
                    if not any(box.overlaps(b) for b in placed):
                        break
                    y += LABEL_STEP; label.set_position((x, y)); box = label.get_window_extent(renderer).expanded(1.2, 1.1)
                assert not any(box.overlaps(b) for b in placed), label.get_text()
                assert y < LABEL_CEILING, (label.get_text(), y)
                if y != initial:
                    axis.plot([x, x], [initial, y - 2], color=V30_LEADER_COLOR, lw=fl.MIN_LINE_PT, gid='leader'); moved += 1
                placed.append(box)
        fig.canvas.draw()
        rep = v30_compare(ref, v30_fingerprint(fig, **V30_FP))
        ok = _v30_finish(fl, fig, out_stem, letters=list(axes), keep_sizes=ion_texts)
    assert ok
    return dict(panels=rep, labels_displaced=moved, clip_ok=ok)
