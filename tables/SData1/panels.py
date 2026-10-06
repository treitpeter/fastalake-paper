"""Supplementary Data 1 (Sample_rosters.xlsx): build the roster workbook from source_data.

Two builds are made from ``source_data/SData1/``:

* ``v28d``: the workbook as published with manuscript v28d (34 sheets, CapScan on the
  331-acquisition roster, 710 acquisitions in all). It is checked cell by cell against
  the digests of the published workbook in ``expected/SHEET_DIGESTS_v28d.tsv``, and the
  regenerated ``ANALYSIS_INDEX.tsv`` byte by byte against ``expected/ANALYSIS_INDEX_v28d.tsv``.
* ``v28f``: the same workbook with CapScan on the 347-acquisition roster (participant 15
  added, 726 acquisitions in all). The changes are listed in ``V28F_CHANGES`` and README.md.

The sheet writer (``write_workbook``) is the original producer's ``sheet()`` function and
navigation code (build_workbook.py, see PROVENANCE.md) with rows passed as
lists instead of dicts. Formatting, column widths, tab colours and hyperlinks are unchanged.

Source data read: every ``source_data/SData1/*.tsv`` (see its README.md) and, for v28f,
``source_data/Fig5/design_347.tsv`` (owned by unit Fig5). PREDICT visit columns, the
recorded group code and the antibiotic-pair analysis rows are withheld (CONVENTIONS.md
section 8); with ``FASTALAKE_CONTROLLED_DIR`` set they are merged back in.

No statistical test is computed in this unit. Every function below counts rows.
"""

import csv
import hashlib
import io
import os
import re
from pathlib import Path

#: Roster sizes stated in the v28d text (v28d line 489, 501).
ROSTER_V28D = {"CAMPI": 15, "Healthy stool": 58, "PREDICT": 306, "CapScan": 331}
#: CapScan roster for v28f (participant-15 join fix; v28d line 501 already prints 347).
CAPSCAN_ROSTER_V28F = 347
#: Capsule regions in region order (Methods, v28d line 201: regions 1-4).
CAPSULE_REGIONS = ("duodenum", "jejunum", "ileum", "colon")
REGION_ORDINAL = {"duodenum": 1, "jejunum": 2, "ileum": 3, "colon": 4, "stool": 5}
#: aliquot of the participant-15 intestinal-fluid acquisitions (make_design347.py sets 'pellet').
P15_ALIQUOT = "pellet"

#: Controlled-access file with the withheld PREDICT cells (one row per cell).
CONTROLLED_FILE = "SData1_PREDICT_withheld_cells.tsv"
#: PREDICT columns withheld from the public tables (visit codes and recorded group code).
WITHHELD_COLUMNS = {
    "PREDICT metadata": ("august_visit_id", "august_group", "later_visit_number"),
    "PREDICT mapping changes": ("august_visit_id", "august_group", "later_visit_number"),
    "PREDICT repeated tubes": ("august_visit_id", "august_group", "later_visit_number"),
    "PREDICT stack eligibility": ("recorded_visit_number",),
}
#: PREDICT rows withheld: membership in the antibiotic On/Off pair analysis reveals exposure.
WITHHELD_ROWS = {
    "PREDICT membership": ("analysis_id", "PREDICT_ANTIBIOTIC42"),
    "PREDICT patients": ("analysis_id", "PREDICT_ANTIBIOTIC42"),
}

#: Sheets computed by code rather than read from source_data.
COMPUTED_SHEETS = {"Cohorts", "Sheet directory"}

#: (sheet, description, source_table) in workbook order, as in build_workbook.py.
SHEETS = [
    ("Read me", "How to interpret this workbook", "Derived by build_workbook.py"),
    ("Cohorts", "Four study acquisition rosters; biological denominators are not summed across designs", "Derived by build_workbook.py"),
    ("Analyses", "One row per analysis or experimental-design subset; follow roster_sheet and roster_filter", "Derived by build_workbook.py"),
    ("HS acquisitions", "58 healthy-stool acquisitions with documented relationships and explicit unknown identities", "healthy_campi/HSTOOL_acquisitions.tsv"),
    ("HS relationships", "Author confirmation and documentary links between designs", "healthy_campi/HSTOOL_design_relationships.tsv"),
    ("HS excluded", "39 historical differential-centrifugation acquisitions excluded from final study scope", "healthy_campi/HSTOOL_excluded39.tsv"),
    ("HS CAMPI membership", "Included and excluded acquisition membership for healthy-stool/CAMPI analyses", "healthy_campi/analysis_membership.tsv"),
    ("HS transfer", "Exact target-to-database-source links; biological independence is separately classified", "healthy_campi/HSTOOL_transfer_links.tsv"),
    ("HS pools", "Each target and every database-source acquisition in own/all-twelve/other-eleven pools", "healthy_campi/HSTOOL_pool_links.tsv"),
    ("CAMPI acquisitions", "Retained deposited acquisitions linked to published protocol and fraction labels", "healthy_campi/CAMPI_acquisitions.tsv"),
    ("CAMPI search arms", "75 primary and 45 second-engine acquisition/database searches", "healthy_campi/CAMPI_search_arm_membership.tsv"),
    ("PREDICT metadata", "August and later explicit metadata links side by side", "predict/ACQUISITION_METADATA_LINKS.tsv"),
    ("PREDICT membership", "All 306 acquisitions classified for each PREDICT analysis", "predict/ANALYSIS_MEMBERSHIP.tsv"),
    ("PREDICT patients", "Patient eligibility relative to each documented source population", "predict/PATIENT_MEMBERSHIP.tsv"),
    ("PREDICT mapping changes", "Three later-resolved mappings and the remaining unmapped analytical sample", "predict/METADATA_MAPPING_CHANGES_AND_UNRESOLVED.tsv"),
    ("PREDICT repeated tubes", "Repeated measurements linked to the same recorded tube", "predict/REPEATED_ACQUISITIONS_OF_RECORDED_TUBES.tsv"),
    ("PREDICT MG exclusions", "Ten acquisitions lacking the matched comparator in the producing manifest", "predict/PUBLIC_MG10_EXCLUDED.tsv"),
    ("PREDICT H patients", "Patient membership for each of twelve primary histidine endpoints", "predict/HISTIDINE_ENDPOINT_PATIENT_ROSTER.tsv"),
    ("CAPSCAN H acquisitions", "331 histidine search acquisitions; quantitative zero is retained", "capscan/HISTIDINE_ACQUISITION_ROSTER.tsv"),
    ("CAPSCAN H finite n", "Saved family-specific finite n ranges; no statistical replay", "capscan/HISTIDINE_SAVED_FINITE_N_SUMMARY.tsv"),
    ("CAPSCAN H participants", "Observed finite participant sets for H-source and reciprocal endpoints only", "capscan/HISTIDINE_FINITE_PARTICIPANT_SETS.tsv"),
    ("CAPSCAN acquisitions", "331 acquisition-to-specimen/participant/region links", "capscan/BASE_ACQUISITIONS.tsv"),
    ("CAPSCAN membership", "Included entities; entity_type distinguishes acquisitions, specimens, DNA runs and aggregation units", "capscan/ANALYSIS_MEMBERSHIP.tsv"),
    ("CAPSCAN exclusions", "Documented omissions and source-specific nonmembership", "capscan/EXCLUSIONS.tsv"),
    ("CAPSCAN MS nonmembership", "Every acquisition absent from each MS analysis roster relative to the 331-search design", "capscan/ACQUISITION_NONMEMBERSHIP.tsv"),
    ("CAPSCAN DNA materials", "Upstream DNA specimen categories: 272 gut/stool and 29 saliva", "capscan/DNA301_SOURCE_CATEGORIES.tsv"),
    ("PREDICT 85 to 77", "Eight patients without any complete-stack specimen; matrix row absence is not a known laboratory cause", "predict/EDS85_TO_NATIVE77_PATIENT_EXCLUSIONS.tsv"),
    ("CAPSCAN ID crosswalk", "Exact linked biology/cross-omics specimen namespaces", "capscan/IDENTIFIER_CROSSWALK.tsv"),
    ("CAPSCAN MG pairs", "Matched benchmark pairs preserving A16 and biology specimen labels", "capscan/BENCHMARK282_ACQUISITIONS.tsv"),
    ("CAPSCAN DNA runs", "301 upstream completed DNA runs, including saliva; native subset is in membership", "capscan/DNA301_RUN_SPECIMEN_ROSTER.tsv"),
    ("CAPSCAN DNA eligibility", "Participant-region DNA availability across functional representations", "capscan/DNA_CELL_ELIGIBILITY.tsv"),
    ("CAPSCAN regional pairs", "Four capsule-versus-stool pairs for each of fourteen participants", "capscan/REGIONAL_PAIRED_LINKS.tsv"),
    ("PREDICT stack eligibility", "Recorded assay presence and complete-stack eligibility", "predict/LATER268_ASSAY_ELIGIBILITY.tsv"),
    ("Sheet directory", "Location and original table for each workbook sheet", "Derived by build_workbook.py"),
]

#: Directory rows of the two sheets the Source Data package removed from the published
#: workbook (Sources: absolute lab paths; Author confirmations). The published Sheet
#: directory still lists them; v28d reproduces that, v28f drops them.
REMOVED_SHEET_RECORDS = [
    ("Sources", "94", "Exact inspected inputs; all recorded hashes freshly verified during workbook build", "Derived by build_workbook.py"),
    ("Author confirmations", "1", "Author-confirmed relationships with explicit scope; verbatim record retained in AUTHOR_CONFIRMATIONS.jsonl", "Derived by build_workbook.py"),
]

#: Cohorts sheet text cells, as hard-coded in build_workbook.py; counts are computed.
COHORT_TEXT = {
    "CAMPI": dict(specimens="Defined SIHUMIx mixture; no individual specimen IDs in retained metadata", participants="Not applicable", design="Twelve published workflow labels (S01–S12); S11 contributes four fractions", roster_sheet="CAMPI acquisitions"),
    "Healthy stool": dict(specimens="Original specimen IDs unavailable; within-design relationships documented", participants="Ten between-donor labels plus a shared experimental donor of unknown overlap", design="12 reinjections; 15 preparations; 10 spatial sections; 11 longitudinal samples; 10 between-donor samples", roster_sheet="HS acquisitions"),
    "PREDICT": dict(specimens="300 explicit tube IDs in August metadata; four acquisitions lack tube links", participants="93 recorded August patients; four acquisitions lack participant links", design="Analysis-specific subsets and later metadata are indexed separately", roster_sheet="PREDICT metadata"),
    "CapScan": dict(design="Fourteen participants contribute capsule and stool; one contributes stool only", roster_sheet="CAPSCAN acquisitions"),
}
COHORT_ROSTER_SHEET = {"CAMPI": "CAMPI acquisitions", "Healthy stool": "HS acquisitions",
                       "PREDICT": "PREDICT metadata", "CapScan": "CAPSCAN acquisitions"}

# ---------------------------------------------------------------- v28f (347) constants
#: Upstream file names and sha256 cited in the v28f CapScan membership rows.
V28F_SOURCES = {
    "search": ("design_347.tsv", "30f70233143a370755d1bc4d1dd8f236d94a2fe4b9b27cc10c64e94602461000"),
    "nonempty": ("eligibility_audit.tsv", "3c58ea360554a7631bdef00c21fad5b394ebd0b4fbfab72b6ff537e585ceb36b"),
    "paired": ("own_stool_participants.tsv", "67f22f573d418bdee6e7b2a7ed0bc1ccc2fc140ee5793ab10db2455f6002ff0e"),
}
#: link_basis for the 16 acquisitions added by the participant-15 join fix.
P15_LINK_BASIS = ("Participant-15 join fix: MS pellet ID to source workbook Pellet, Main and "
                  "overview Main; accepted only if unique and Subject and Type agree.")
ZERO_REASON = "Zero total quantitative intensity."
P15_STOOL_ONLY_REASON = "Participant 15 has stool-only proteomics and cannot form the required five-region block."
P15_NO_MG_REASON = ("Participant-15 intestinal-fluid acquisition added in the 347 roster; no "
                    "specimen-matched metagenome database exists for it.")
#: v28f text changes in the computed sheets (sheet, cell, new text).
V28F_COHORT_DESIGN = "All fifteen participants contribute capsule and stool"
V28F_DESCRIPTIONS = {
    "CAPSCAN acquisitions": "347 acquisition-to-specimen/participant/region links",
    "CAPSCAN MS nonmembership": "Every acquisition absent from each MS analysis roster relative to the 347-search design",
    "CAPSCAN regional pairs": "Four capsule-versus-stool pairs for each of fifteen participants",
}
#: CapScan analyses re-indexed on the 347 roster: v28d id -> (v28f id, source key, rule, region order).
V28F_REBUILT = {
    "CAPSCAN_SEARCH_ALL": ("CAPSCAN_SEARCH_ALL", "search", "All 347 rows in the producing CapScan design", "string"),
    "CAPSCAN_QUANT_NONEMPTY": ("CAPSCAN_QUANT_NONEMPTY", "nonempty", "Positive total quantified intensity in the repaired simplified source", "string"),
    "CAPSCAN_REGIONAL_DESCRIPTIVE15": ("CAPSCAN_REGIONAL_DESCRIPTIVE15", "nonempty", "Nonempty acquisitions averaged within specimen, then participant-region", "string"),
    "CAPSCAN_REGIONAL_PAIRED14": ("CAPSCAN_REGIONAL_PAIRED15", "paired", "Participants 1–15 with all five regions; nonempty acquisitions, specimen means then region means", "ordinal"),
}
#: Note appended to the limitations of CapScan analyses not re-indexed on 347.
V28F_NOT_REINDEXED = (" v28f: this row and its membership are on the 331-acquisition roster; "
                      "the 347-acquisition membership is not yet indexed (see OPEN_ITEMS.md).")
V28F_NOT_REINDEXED_IDS = ("CAPSCAN_FUNCTIONAL_NONDNA14", "CAPSCAN_NATIVE12", "CAPSCAN_FUNCTIONAL_DNA12",
                          "CAPSCAN_FUNCTIONAL_NONDNA12_SENS", "CAPSCAN_HISTIDINE_ROUTES")
V28F_DESCRIPTIVE_LIMITATION = "75 observed cells; on the 347 roster this equals the paired five-region design."
V28F_PAIRED_UNIT = "Participant for four paired regions; 75 participant-region profiles"


# ---------------------------------------------------------------- table I/O
def slug(name):
    """File-name form of a sheet name ('PREDICT 85 to 77' -> 'predict_85_to_77')."""
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def sheet_file(name):
    """source_data file name of a published sheet."""
    return f"v28d_{slug(name)}.tsv"


def write_tsv(path, rows, lineterminator="\n"):
    """Write a list of rows (lists of str) as tab-separated text without quoting."""
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator=lineterminator,
                       quoting=csv.QUOTE_NONE, escapechar=None)
        w.writerows(rows)


def read_tsv(path):
    """Read a tab-separated table as a list of rows (lists of str), header first."""
    with open(path, newline="", encoding="utf-8") as fh:
        return [r for r in csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)]


def read_xlsx_values(path):
    """Read every sheet of a workbook as rows of str (empty cells as '').

    Parameters
    ----------
    path : path-like
        Workbook.

    Returns
    -------
    dict
        ``{sheet name: [header, row, ...]}`` in workbook order.
    """
    from openpyxl import load_workbook
    wb = load_workbook(path, read_only=True)
    out = {}
    for ws in wb.worksheets:
        out[ws.title] = [["" if v is None else str(v) for v in r] for r in ws.iter_rows(values_only=True)]
    wb.close()
    return out


def grid_digest(rows):
    """sha256 of a table serialised as tab-joined lines (the comparison key for sheets)."""
    h = hashlib.sha256()
    for r in rows:
        h.update(("\t".join(r) + "\n").encode("utf-8"))
    return h.hexdigest()


# ---------------------------------------------------------------- controlled split / merge
def split_public(name, header, body):
    """Remove withheld PREDICT columns and rows from one sheet.

    Returns
    -------
    tuple
        ``(public_header, public_body, withheld_cells)``; each withheld cell is
        ``(sheet, row_index, column_index, column, value)`` with 0-based data-row index.
    """
    drop_cols = [i for i, h in enumerate(header) if h in WITHHELD_COLUMNS.get(name, ())]
    rule = WITHHELD_ROWS.get(name)
    drop_rows = set()
    if rule:
        ci = header.index(rule[0])
        drop_rows = {i for i, r in enumerate(body) if r[ci] == rule[1]}
    cells = []
    for i, r in enumerate(body):
        for j, h in enumerate(header):
            if i in drop_rows or j in drop_cols:
                cells.append((name, i, j, h, r[j]))
    keep = [j for j in range(len(header)) if j not in drop_cols]
    pub_header = [header[j] for j in keep]
    pub_body = [[r[j] for j in keep] for i, r in enumerate(body) if i not in drop_rows]
    return pub_header, pub_body, cells


def controlled_cells():
    """Withheld cells from FASTALAKE_CONTROLLED_DIR, or None if the file is absent."""
    d = os.environ.get("FASTALAKE_CONTROLLED_DIR")
    p = Path(d) / CONTROLLED_FILE if d else None
    if p is None or not p.is_file():
        print(f"SKIPPED: needs controlled input {CONTROLLED_FILE} (see controlled/CONTROLLED_INPUTS.md); "
              "building the public workbook without the withheld PREDICT cells")
        return None
    rows = read_tsv(p)[1:]
    return [(r[0], int(r[1]), int(r[2]), r[3], r[4]) for r in rows]


def merge_controlled(name, public_rows, cells):
    """Rebuild the full sheet from its public rows and its withheld cells.

    Parameters
    ----------
    name : str
        Sheet name.
    public_rows : list of list of str
        Public table, header first.
    cells : list of tuple
        Withheld cells ``(sheet, row_index, column_index, column, value)``.

    Returns
    -------
    list of list of str
        The sheet as published, header first.
    """
    mine = [c for c in cells if c[0] == name]
    if not mine:
        return public_rows
    if WITHHELD_COLUMNS.get(name):
        names = {c[2]: c[3] for c in mine}
        vals = {(c[1], c[2]): c[4] for c in mine}
        ncols = len(public_rows[0]) + len(names)
        it = iter(public_rows[0])
        header = [names[j] if j in names else next(it) for j in range(ncols)]
        grid = []
        for i, r in enumerate(public_rows[1:]):
            it = iter(r)
            grid.append([vals[(i, j)] if j in names else next(it) for j in range(ncols)])
        return [header] + grid
    header = public_rows[0]
    withheld = {}
    for _, i, j, _, v in mine:
        withheld.setdefault(i, {})[j] = v
    pub = iter(public_rows[1:])
    grid = [[withheld[i][j] for j in range(len(header))] if i in withheld else next(pub)
            for i in range(len(public_rows) - 1 + len(withheld))]
    return [header] + grid


def load_v28d(sd):
    """Read the published data sheets from source_data, merging controlled cells if present.

    Returns
    -------
    tuple
        ``(sheets, complete)``: ``{sheet: rows}`` for every non-computed sheet, and True if
        the withheld PREDICT cells were available.
    """
    cells = controlled_cells()
    sheets = {}
    for name, _, _ in SHEETS:
        if name in COMPUTED_SHEETS:
            continue
        rows = read_tsv(Path(sd) / sheet_file(name))
        sheets[name] = merge_controlled(name, rows, cells) if cells else rows
    return sheets, cells is not None


def as_dicts(rows):
    """Rows (header first) as a list of dicts."""
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def from_dicts(header, dicts):
    """List of dicts as rows (header first)."""
    return [list(header)] + [[d.get(k, "") for k in header] for d in dicts]


# ---------------------------------------------------------------- computed sheets
def cohorts_sheet(sheets, capscan_design=None):
    """Cohorts sheet: one row per cohort, acquisitions counted from the roster sheet.

    Parameters
    ----------
    sheets : dict
        Data sheets of the build.
    capscan_design : str or None
        v28f replacement for the CapScan design text; None keeps the v28d text.

    Returns
    -------
    list of list of str
        Header and four rows (unit: acquisition; CapScan also specimen and participant).
    """
    header = ["cohort", "acquisitions", "specimens", "participants", "design", "roster_sheet"]
    rows = [header]
    for cohort, sheet in COHORT_ROSTER_SHEET.items():
        d = as_dicts(sheets[sheet])
        t = dict(COHORT_TEXT[cohort])
        if cohort == "CapScan":
            t["specimens"] = str(len({r["specimen_id_biology"] for r in d}))
            t["participants"] = str(len({r["participant_id"] for r in d}))
            if capscan_design:
                t["design"] = capscan_design
        rows.append([cohort, str(len(d)), t["specimens"], t["participants"], t["design"], t["roster_sheet"]])
    return rows


def sheet_directory(order, sheets, descriptions, removed_records):
    """Sheet directory rows (sheet, data_rows, description, source_table)."""
    rows = [["sheet", "data_rows", "description", "source_table"]]
    for name, desc, src in order:
        if name == "Sheet directory":
            continue
        rows.append([name, str(len(sheets[name]) - 1), descriptions.get(name, desc), src])
    rows += [list(r) for r in removed_records]
    return rows


# ---------------------------------------------------------------- workbook writer (ported)
def write_workbook(path, order, sheets):
    """Write the workbook; ``sheet()`` formatting and links ported from build_workbook.py.

    Parameters
    ----------
    path : path-like
        Output .xlsx.
    order : list of str
        Sheet names in workbook order.
    sheets : dict
        ``{sheet: rows}`` with header first; every value a str.
    """
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter
    wb = Workbook()
    wb.remove(wb.active)
    for name in order:
        rows = sheets[name]
        keys, values = rows[0], rows[1:]
        assert values, name
        ws = wb.create_sheet(name)
        ws.append(keys)
        for row in values:
            ws.append(row)
        for row in ws:
            for c in row:
                c.data_type = "s"
                c.alignment = Alignment(vertical="top", wrap_text=True)
        for c in ws[1]:
            c.font = Font(name="Calibri", bold=True, color="FFFFFF")
            c.fill = PatternFill("solid", fgColor="24475B")
        for j, key in enumerate(keys, 1):
            maxlen = max([len(key)] + [len(r[j - 1]) for r in values[:100]])
            wide = any(s in key.lower() for s in ["reason", "rule", "purpose", "limit", "note", "interpret", "description"])
            ws.column_dimensions[get_column_letter(j)].width = min(70 if wide else 48, max(18, maxlen + 2))
        ws.freeze_panes = "C2" if len(keys) > 2 else "A2"
        ws.auto_filter.ref = ws.dimensions
        ws.sheet_view.zoomScale = 85
        ws.row_dimensions[1].height = 32
        ws.sheet_properties.pageSetUpPr.fitToPage = True
        ws.page_setup.orientation = "landscape"
        ws.page_setup.paperSize = ws.PAPERSIZE_A3
        ws.page_setup.fitToWidth = 1
        ws.page_setup.fitToHeight = 0
        ws.print_title_rows = "1:1"
    for name in ["Analyses", "Cohorts", "Sheet directory"]:
        ws = wb[name]
        headers = {c.value: c.column for c in ws[1]}
        for key in ["roster_sheet", "exclusion_sheet", "sheet"]:
            if key in headers:
                for row in range(2, ws.max_row + 1):
                    cell = ws.cell(row, headers[key])
                    if cell.value in wb.sheetnames:
                        cell.hyperlink = "#'%s'!A1" % cell.value
                        cell.font = Font(color="0563C1", underline="single")
    for ws in wb:
        ws.sheet_properties.tabColor = ("497E76" if ws.title.startswith("HS") else "517CA4" if ws.title.startswith("PREDICT")
                                        else "98689A" if ws.title.startswith("CAPSCAN") else "B58A43")
    wb.properties.creator = "FastaLake"
    wb.properties.lastModifiedBy = "FastaLake"
    wb.save(path)


def write_analysis_index(path, analyses_rows):
    """ANALYSIS_INDEX.tsv as csv.DictWriter wrote it in build_workbook.py (CRLF line ends)."""
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerows(analyses_rows)


# ---------------------------------------------------------------- validation (ported)
def validate_counts(sheets, complete):
    """Port of build_workbook.py's final checks: every Analyses row's n_MS_acquisitions
    equals the included acquisitions of its membership sheet, without duplicates.

    Parameters
    ----------
    sheets : dict
        Built sheets.
    complete : bool
        False if the antibiotic-pair rows are withheld (that row is then skipped).

    Returns
    -------
    list of tuple
        ``(analysis_id, n_printed, n_counted)`` for every checked row.
    """
    hs = as_dicts(sheets["HS CAMPI membership"])
    pm = as_dicts(sheets["PREDICT membership"])
    cm = as_dicts(sheets["CAPSCAN membership"])
    hacq = as_dicts(sheets["CAPSCAN H acquisitions"])
    out = []
    for r in as_dicts(sheets["Analyses"]):
        aid, n = r["analysis_id"], int(r["n_MS_acquisitions"])
        if r["cohort"] in ("HSTOOL", "CAMPI"):
            found = [x for x in hs if x["analysis_id"] == aid and x["included"] == "yes"]
        elif r["cohort"] == "PREDICT":
            if aid in {v[1] for v in WITHHELD_ROWS.values()} and not complete:
                continue
            found = [x for x in pm if x["analysis_id"] == aid and x["status"] == "included"]
        elif aid == "CAPSCAN_HISTIDINE_ROUTES":
            found = hacq
        else:
            found = [x for x in cm if x["analysis_id"] == aid and x["entity_type"] == "acquisition"]
        assert len(found) == n, (aid, len(found), n)
        assert len({x["acquisition_id"] for x in found}) == n, aid
        out.append((aid, n, len(found)))
    return out


def validate_capscan_entity_counts(sheets):
    """Every CapScan Analyses row's specimen, participant and cell counts equal its membership."""
    cm = as_dicts(sheets["CAPSCAN membership"])
    out = []
    for r in as_dicts(sheets["Analyses"]):
        if r["cohort"] != "CapScan" or r["analysis_id"] == "CAPSCAN_HISTIDINE_ROUTES":
            continue
        for col, et in [("n_recorded_proteomic_specimens", "specimen"), ("n_participants", "participant"),
                        ("n_participant_region_cells", "participant_region")]:
            n = sum(1 for x in cm if x["analysis_id"] == r["analysis_id"] and x["entity_type"] == et)
            assert str(n) == r[col], (r["analysis_id"], col, n, r[col])
            out.append((r["analysis_id"], col, n))
    return out


# ---------------------------------------------------------------- CapScan 347 roster
def capscan_acquisitions_347(v28d_acq, p15_join, fig5_design):
    """CAPSCAN acquisitions on the 347 roster: the 331 published rows plus 16 participant-15 rows.

    Parameters
    ----------
    v28d_acq : list of list of str
        Published 'CAPSCAN acquisitions' sheet (331 rows).
    p15_join : list of list of str
        ``CAPSCAN_participant15_join.tsv`` (A80 join-fix evidence, 16 rows).
    fig5_design : list of list of str
        ``source_data/Fig5/design_347.tsv`` (347 rows; participant, region, specimen_id).

    Returns
    -------
    list of list of str
        Header and 347 rows (unit: acquisition). Asserts agreement with the Fig5 design
        for every acquisition (participant, region, specimen).
    """
    header = v28d_acq[0]
    rows = as_dicts(v28d_acq)
    for r in as_dicts(p15_join):
        main = r["workbook_Main"]
        rows.append(dict(acquisition_id=r["sample_id"], specimen_id_biology=f"main:{main}",
                         participant_id=r["subject_id"], region=r["region"], aliquot=P15_ALIQUOT,
                         aliquot_main_id=main, specimen_id_crossomics=main, link_basis=P15_LINK_BASIS))
    design = {d["sample_id"]: d for d in as_dicts(fig5_design)}
    assert len(rows) == CAPSCAN_ROSTER_V28F == len(design) == len({r["acquisition_id"] for r in rows})
    for r in rows:
        d = design[r["acquisition_id"]]
        assert (d["participant"], d["region"], d["specimen_id"]) == (r["participant_id"], r["region"], r["specimen_id_biology"]), r
    return from_dicts(header, rows)


def membership_rows(analysis_id, acq, source, rule, region_order):
    """Membership rows (acquisition, specimen, participant, participant_region) for one analysis.

    Ordering reproduces the published sheet: acquisitions in roster order, specimens by
    first appearance, participants numerically, cells by (participant, region) as text ('string') or by
    participant number and region order ('ordinal').

    Parameters
    ----------
    analysis_id : str
    acq : list of dict
        Included acquisitions (CAPSCAN acquisitions rows), in roster order.
    source : tuple
        ``(source_file, source_sha256)``.
    rule : str
        inclusion_rule text.
    region_order : {'string', 'ordinal'}

    Returns
    -------
    list of dict
    """
    base = dict(analysis_id=analysis_id, run_accession="", source_file=source[0], source_sha256=source[1], inclusion_rule=rule)
    out = []
    for a in acq:
        out.append(dict(base, entity_type="acquisition", entity_id=a["acquisition_id"], acquisition_id=a["acquisition_id"],
                        specimen_id_biology=a["specimen_id_biology"], specimen_id_crossomics=a["specimen_id_crossomics"],
                        participant_id=a["participant_id"], region=a["region"]))
    seen = {}
    for a in acq:
        seen.setdefault(a["specimen_id_crossomics"], a)
    for sid, a in seen.items():
        out.append(dict(base, entity_type="specimen", entity_id=sid, acquisition_id="",
                        specimen_id_biology=a["specimen_id_biology"], specimen_id_crossomics=sid,
                        participant_id=a["participant_id"], region=a["region"]))
    for p in sorted({a["participant_id"] for a in acq}, key=int):
        out.append(dict(base, entity_type="participant", entity_id=p, acquisition_id="", specimen_id_biology="",
                        specimen_id_crossomics="", participant_id=p, region=""))
    cells = {(a["participant_id"], a["region"]) for a in acq}
    key = (lambda c: (c[0], c[1])) if region_order == "string" else (lambda c: (int(c[0]), REGION_ORDINAL[c[1]]))
    for p, reg in sorted(cells, key=key):
        out.append(dict(base, entity_type="participant_region", entity_id=f"{p}|{reg}", acquisition_id="",
                        specimen_id_biology="", specimen_id_crossomics="", participant_id=p, region=reg))
    return out


def nonmember_rows(analysis_id, roster, included, reason_of):
    """MS nonmembership rows: every roster acquisition not included, with its reason."""
    inc = {a["acquisition_id"] for a in included}
    return [dict(analysis_id=analysis_id, acquisition_id=a["acquisition_id"], specimen_id_biology=a["specimen_id_biology"],
                 specimen_id_crossomics=a["specimen_id_crossomics"], participant_id=a["participant_id"], region=a["region"],
                 disposition_type="not_in_analysis", reason=reason_of(a))
            for a in roster if a["acquisition_id"] not in inc]


def capscan_populations(roster, nonempty_ids):
    """Acquisition populations of the four re-indexed CapScan analyses.

    Parameters
    ----------
    roster : list of dict
        CAPSCAN acquisitions rows (331 or 347).
    nonempty_ids : set of str
        Acquisitions with positive total quantified intensity.

    Returns
    -------
    dict
        ``{'search', 'nonempty', 'paired'}`` -> included acquisitions (unit: acquisition).
        'paired' keeps participants whose nonempty acquisitions cover all five sampling groups.
    """
    nonempty = [a for a in roster if a["acquisition_id"] in nonempty_ids]
    groups = {}
    for a in nonempty:
        groups.setdefault(a["participant_id"], set()).add(a["region"])
    full = {p for p, g in groups.items() if g == set(REGION_ORDINAL)}
    return dict(search=list(roster), nonempty=nonempty, paired=[a for a in nonempty if a["participant_id"] in full])


def regional_pairs(participants, source_file):
    """Capsule-versus-stool pair rows (participant, region, paired_stool_region, source_file)."""
    rows = [["participant_id", "region", "paired_stool_region", "source_file"]]
    for p in sorted(participants, key=int):
        rows += [[p, reg, "stool", source_file] for reg in CAPSULE_REGIONS]
    return rows


def reproduce_v28d_capscan(sheets):
    """Positive control: regenerate the four re-indexed analyses on the published 331 roster.

    Returns
    -------
    dict
        Per analysis: membership and nonmembership rows equal to the published rows (bool).
    """
    roster = as_dicts(sheets["CAPSCAN acquisitions"])
    mem = as_dicts(sheets["CAPSCAN membership"])
    non = as_dicts(sheets["CAPSCAN MS nonmembership"])
    excl = [x for x in as_dicts(sheets["CAPSCAN exclusions"]) if x["analysis_id"] == "CAPSCAN_QUANT_NONEMPTY"]
    zero = {x["entity_id"] for x in excl if x["entity_type"] == "acquisition"}
    pops = capscan_populations(roster, {a["acquisition_id"] for a in roster} - zero)
    res = {}
    for aid, (_, key, _, order) in V28F_REBUILT.items():
        pub = [x for x in mem if x["analysis_id"] == aid]
        src = (pub[0]["source_file"], pub[0]["source_sha256"])
        gen = membership_rows(aid, pops[key], src, pub[0]["inclusion_rule"], order)
        reason = lambda a: ZERO_REASON if a["acquisition_id"] in zero else P15_STOOL_ONLY_REASON
        gnon = nonmember_rows(aid, roster, pops[key], reason) if key != "search" else []
        pnon = [x for x in non if x["analysis_id"] == aid]
        res[aid] = dict(membership_equal=gen == pub, nonmembership_equal=gnon == pnon, n=len(pops[key]))
    pairs = regional_pairs({a["participant_id"] for a in pops["paired"]}, "regional_host_tests.tsv")
    res["CAPSCAN regional pairs"] = dict(equal=pairs == sheets["CAPSCAN regional pairs"])
    return res


def build_v28f(v28d, p15_join, fig5_design, eligibility):
    """Apply the v28f CapScan changes to the v28d sheets.

    Parameters
    ----------
    v28d : dict
        Data sheets of the v28d build.
    p15_join, fig5_design, eligibility : list of list of str
        ``CAPSCAN_participant15_join.tsv``, ``source_data/Fig5/design_347.tsv``,
        ``CAPSCAN_eligibility_audit_347.tsv``.

    Returns
    -------
    dict
        Data sheets of the v28f build.
    """
    s = {k: [list(r) for r in v] for k, v in v28d.items()}
    acq = capscan_acquisitions_347(v28d["CAPSCAN acquisitions"], p15_join, fig5_design)
    s["CAPSCAN acquisitions"] = acq
    xw = v28d["CAPSCAN ID crosswalk"][0]
    s["CAPSCAN ID crosswalk"] = from_dicts(xw, as_dicts(acq))
    roster = as_dicts(acq)
    elig = {r["acquisition"]: r["current_positive"] for r in as_dicts(eligibility)}
    fig5_pos = {r["sample_id"]: r["og_profile_positive"] for r in as_dicts(fig5_design)}
    assert set(elig) == set(fig5_pos) and elig == fig5_pos, "eligibility audit disagrees with Fig5 design"
    nonempty = {a for a, v in elig.items() if v == "True"}
    zero = set(elig) - nonempty
    pops = capscan_populations(roster, nonempty)

    mem_header = v28d["CAPSCAN membership"][0]
    old = as_dicts(v28d["CAPSCAN membership"])
    new_mem, done = [], set()
    for r in old:
        aid = r["analysis_id"]
        if aid in V28F_REBUILT:
            if aid not in done:
                nid, key, rule, order = V28F_REBUILT[aid]
                new_mem += membership_rows(nid, pops[key], V28F_SOURCES[key], rule, order)
                done.add(aid)
            continue
        new_mem.append(r)
    s["CAPSCAN membership"] = from_dicts(mem_header, new_mem)

    non_header = v28d["CAPSCAN MS nonmembership"][0]
    new_non, done = [], set()
    reason = lambda a: ZERO_REASON if a["acquisition_id"] in zero else P15_STOOL_ONLY_REASON
    p15_rows = [a for a in roster if a["link_basis"] == P15_LINK_BASIS]
    old_non = as_dicts(v28d["CAPSCAN MS nonmembership"])
    for k, r in enumerate(old_non):
        aid = r["analysis_id"]
        if aid in V28F_REBUILT:
            if aid not in done:
                nid, key, _, _ = V28F_REBUILT[aid]
                new_non += nonmember_rows(nid, roster, pops[key], reason)
                done.add(aid)
            continue
        new_non.append(r)
        last_of_block = k + 1 == len(old_non) or old_non[k + 1]["analysis_id"] != aid
        if aid == "CAPSCAN_PUBLIC_MG_282" and last_of_block:
            new_non += [dict(analysis_id=aid, acquisition_id=a["acquisition_id"], specimen_id_biology=a["specimen_id_biology"],
                             specimen_id_crossomics=a["specimen_id_crossomics"], participant_id=a["participant_id"],
                             region=a["region"], disposition_type="not_in_analysis", reason=P15_NO_MG_REASON) for a in p15_rows]
    s["CAPSCAN MS nonmembership"] = from_dicts(non_header, new_non)

    ex_header = v28d["CAPSCAN exclusions"][0]
    s["CAPSCAN exclusions"] = from_dicts(ex_header, [r for r in as_dicts(v28d["CAPSCAN exclusions"])
                                                     if r["analysis_id"] != "CAPSCAN_REGIONAL_PAIRED14"])
    s["CAPSCAN regional pairs"] = regional_pairs({a["participant_id"] for a in pops["paired"]}, "four_region_tests.tsv")

    an_header = v28d["Analyses"][0]
    mem = as_dicts(s["CAPSCAN membership"])
    count = lambda aid, et: str(sum(1 for x in mem if x["analysis_id"] == aid and x["entity_type"] == et))
    rows = []
    for r in as_dicts(v28d["Analyses"]):
        aid = r["analysis_id"]
        if aid in V28F_REBUILT:
            nid = V28F_REBUILT[aid][0]
            r = dict(r, analysis_id=nid, n_MS_acquisitions=count(nid, "acquisition"),
                     n_recorded_proteomic_specimens=count(nid, "specimen"), n_participants=count(nid, "participant"),
                     n_participant_region_cells=count(nid, "participant_region"),
                     roster_filter=r["roster_filter"].replace(aid, nid), exclusion_filter=r["exclusion_filter"].replace(aid, nid),
                     inclusion_rule=V28F_REBUILT[aid][2], provenance=V28F_SOURCES[V28F_REBUILT[aid][1]][0])
            if aid == "CAPSCAN_REGIONAL_DESCRIPTIVE15":
                r["limitations"] = V28F_DESCRIPTIVE_LIMITATION
            if aid == "CAPSCAN_REGIONAL_PAIRED14":
                r["analysis_unit"] = r["aggregation"] = V28F_PAIRED_UNIT
        elif aid in V28F_NOT_REINDEXED_IDS:
            r = dict(r, limitations=r["limitations"] + V28F_NOT_REINDEXED)
        rows.append(r)
    s["Analyses"] = from_dicts(an_header, rows)
    return s


def discrepancies(table, roster_rows):
    """Recount the specimen-number discrepancy census against a roster.

    Returns
    -------
    tuple
        ``(rows, n_files, n_in_roster)``; rows carry in_analysis_roster recomputed
        (unit: acquisition).
    """
    ids = {r["acquisition_id"] for r in as_dicts(roster_rows)}
    d = as_dicts(table)
    for r in d:
        r["in_analysis_roster"] = "yes" if r["acquisition"] in ids else "no"
    rows = from_dicts(table[0], d)
    return rows, len(d), sum(r["in_analysis_roster"] == "yes" for r in d)


# ---------------------------------------------------------------- numbers (v28d lines 489, 501, 503)
def stat_numbers(v28d, v28f, p15_join, disc_v28f, complete):
    """Every number of v28d lines 489, 501 and 503 that the rosters determine.

    Parameters
    ----------
    v28d, v28f : dict
        Built sheets.
    p15_join : list of list of str
    disc_v28f : tuple
        Output of ``discrepancies`` on the 347 roster.
    complete : bool
        True if the withheld PREDICT rows are present (needed for the antibiotic-pair count).

    Returns
    -------
    dict
        ``{key: value}``; a value is None if it needs a controlled input that is absent.
    """
    n = {}
    coh = {r[0]: r for r in cohorts_sheet(v28f)[1:]}
    n["L489_total_acquisitions_v28d"] = sum(int(r[1]) for r in cohorts_sheet(v28d)[1:])
    n["L489_total_acquisitions_v28f"] = sum(int(r[1]) for r in coh.values())
    campi = as_dicts(v28f["CAMPI acquisitions"])
    n["L501_campi_acquisitions"] = len(campi)
    n["L501_campi_S11_fractions"] = sum(r["publication_workflow_id"] == "S11" for r in campi)
    n["L501_campi_laboratories"] = len({r["submitting_team"] for r in campi})
    hs = as_dicts(v28f["HS acquisitions"])
    n["L501_hs_acquisitions"] = len(hs)
    for design in ("TechnicalReplicates", "WorkflowReplicates", "OnePiece", "Intra", "Interindividual"):
        n[f"L501_hs_{design}"] = sum(r["design"] == design for r in hs)
    n["L501_hs_excluded_dc"] = len(as_dicts(v28f["HS excluded"]))
    pm = as_dicts(v28f["PREDICT metadata"])
    n["L501_predict_acquisitions"] = len(pm)
    linked = [r for r in pm if r["august_participant_id"]]
    n["L501_predict_linked_acquisitions"] = len(linked)
    n["L501_predict_patients"] = len({r["august_participant_id"] for r in linked})
    mem = as_dicts(v28f["PREDICT membership"])
    n["L501_predict_mg_comparisons"] = sum(r["analysis_id"] == "PREDICT_PUBLIC_MG296" and r["status"] == "included" for r in mem)
    n["L501_predict_visit1"] = sum(r["analysis_id"] == "PREDICT_VISIT1_91" and r["status"] == "included" for r in mem)
    pat = as_dicts(v28f["PREDICT patients"])
    inc = lambda aid: sum(r["analysis_id"] == aid and r["status"] == "included" for r in pat)
    n["L501_predict_antibiotic_pairs"] = inc("PREDICT_ANTIBIOTIC42") if complete else None
    n["L501_predict_eds85"] = inc("PREDICT_EDS85_PRIMARY")
    n["L501_predict_native77"] = inc("PREDICT_NATIVE77")
    n["L501_capscan_acquisitions"] = int(coh["CapScan"][1])
    n["L501_capscan_specimens"] = int(coh["CapScan"][2])
    n["L501_capscan_participants"] = int(coh["CapScan"][3])
    acq = as_dicts(v28f["CAPSCAN acquisitions"])
    fluid = {r["participant_id"] for r in acq if r["region"] in CAPSULE_REGIONS}
    stool = {r["participant_id"] for r in acq if r["region"] == "stool"}
    n["L503_participants_fluid_and_stool"] = len(fluid & stool)
    n["L503_p15_fluid_acquisitions"] = sum(r["subject_id"] == "15" and r["region"] in CAPSULE_REGIONS for r in as_dicts(p15_join))
    n["L503_discrepant_acquisitions"] = disc_v28f[1]
    n["L503_discrepant_in_roster"] = disc_v28f[2]
    return n


# ---------------------------------------------------------------- builds and checks
def build_v28d(sd):
    """All 34 sheets of the published (v28d) workbook from source_data.

    Returns
    -------
    tuple
        ``(sheets, complete)``; see ``load_v28d``.
    """
    sheets, complete = load_v28d(sd)
    sheets["Cohorts"] = cohorts_sheet(sheets)
    sheets["Sheet directory"] = sheet_directory(SHEETS, sheets, {}, REMOVED_SHEET_RECORDS)
    return sheets, complete


def build_v28f_all(v28d, sd, repo):
    """All 34 sheets of the v28f workbook (CapScan on 347 acquisitions)."""
    data = {k: v for k, v in v28d.items() if k not in COMPUTED_SHEETS}
    sheets = build_v28f(data, read_tsv(Path(sd) / "CAPSCAN_participant15_join.tsv"),
                        read_tsv(Path(repo) / "source_data" / "Fig5" / "design_347.tsv"),
                        read_tsv(Path(sd) / "CAPSCAN_eligibility_audit_347.tsv"))
    sheets["Cohorts"] = cohorts_sheet(sheets, V28F_COHORT_DESIGN)
    sheets["Sheet directory"] = sheet_directory(SHEETS, sheets, V28F_DESCRIPTIONS, [])
    return sheets


def check_v28d(xlsx, analysis_index, expected, complete, receipt):
    """Compare the rebuilt v28d workbook with the published one and write a JSON receipt.

    Every sheet is read back from ``xlsx`` and its digest compared with
    ``expected/SHEET_DIGESTS_v28d.tsv`` (full digest with the controlled cells, public
    digest without). The two computed sheets are also compared cell by cell with their
    published copies, and ``ANALYSIS_INDEX.tsv`` byte by byte.

    Returns
    -------
    dict
        The receipt (also written to ``receipt``).
    """
    import json
    expected = Path(expected)
    back = read_xlsx_values(xlsx)
    exp = {r[0]: r for r in read_tsv(expected / "SHEET_DIGESTS_v28d.tsv")[1:]}
    sheets = {}
    for name, _, _ in SHEETS:
        rows = back[name]
        if complete:
            same = grid_digest(rows) == exp[name][3]
        else:
            h, b, _ = split_public(name, rows[0], rows[1:])
            same = grid_digest([h] + b) == exp[name][4]
        sheets[name] = "IDENTICAL" if same else "DIFFERENT"
    cell_diffs = {}
    for name in sorted(COMPUTED_SHEETS):
        pub = read_tsv(expected / f"{slug(name)}_v28d.tsv")
        got = back[name]
        cell_diffs[name] = [dict(row=i, column=pub[0][j], published=pub[i][j], rebuilt=got[i][j])
                            for i in range(max(len(pub), len(got))) for j in range(len(pub[0]))
                            if i >= len(got) or i >= len(pub) or pub[i][j] != got[i][j]]
    ai_same = Path(analysis_index).read_bytes() == (expected / "ANALYSIS_INDEX_v28d.tsv").read_bytes()
    rec = dict(unit="SData1", build="v28d", mode="complete" if complete else "public (withheld PREDICT cells absent)",
               sheet_order_identical=list(back) == [s[0] for s in SHEETS], sheets=sheets,
               computed_sheet_cell_differences=cell_diffs, analysis_index_byte_identical=ai_same,
               result="IDENTICAL" if ai_same and all(v == "IDENTICAL" for v in sheets.values()) else "DIFFERENT")
    Path(receipt).write_text(json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    return rec


def diff_builds(a, b):
    """Per-sheet change summary from build ``a`` to build ``b`` (rows added, removed, kept)."""
    out = {}
    for name, _, _ in SHEETS:
        ra = {"\t".join(r) for r in a[name][1:]}
        rb = {"\t".join(r) for r in b[name][1:]}
        if ra != rb or a[name][0] != b[name][0]:
            out[name] = dict(rows_v28d=len(a[name]) - 1, rows_v28f=len(b[name]) - 1,
                             rows_removed=len(ra - rb), rows_added=len(rb - ra))
    return out


def receipt_v28f(sheets, complete, changes, controls, receipt):
    """Write the v28f build receipt: per-sheet digests, the change summary and the controls."""
    import json
    dig = {}
    for name, _, _ in SHEETS:
        rows = sheets[name]
        h, b, _ = split_public(name, rows[0], rows[1:])
        dig[name] = dict(data_rows=len(rows) - 1, sha256_public=grid_digest([h] + b),
                         sha256_full=grid_digest(rows) if complete else None)
    rec = dict(unit="SData1", build="v28f", mode="complete" if complete else "public (withheld PREDICT cells absent)",
               positive_controls_on_331=controls, changes_from_v28d=changes, sheets=dig)
    Path(receipt).write_text(json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    return rec


#: English number words, for numbers the manuscript spells out.
WORDS = ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven",
         "twelve", "thirteen", "fourteen", "fifteen", "sixteen", "seventeen", "eighteen", "nineteen", "twenty"]


def number_word(n, capital=False):
    """Spell a count from 0 to 20 as the manuscript prints it ('Sixteen' with capital=True)."""
    w = WORDS[n]
    return w.capitalize() if capital else w


def numbers_table(numbers_tsv, n):
    """Recomputed value for every NUMBERS.tsv row, in file order.

    The ``code`` column names the key as ``panels.py:stat_numbers[<key>]`` and, for a
    spelled number, ends in ``|word`` or ``|Word``.

    Returns
    -------
    list of tuple
        ``(v28d_line, printed_value, recomputed_value, status)``.
    """
    out = []
    for r in as_dicts(read_tsv(numbers_tsv)):
        m = re.search(r"stat_numbers\[(\w+)\](?:\|(word|Word))?", r["code"])
        if not m:
            out.append((r["v28d_line"], r["printed_value"], "", r["status"]))
            continue
        v = n.get(m.group(1))
        if v is not None and m.group(2):
            v = number_word(v, capital=m.group(2) == "Word")
        out.append((r["v28d_line"], r["printed_value"], "" if v is None else str(v), r["status"]))
    return out


def discrepancy_census(audit):
    """Specimen-number discrepancy census from the RAW header audit.

    Port of the census in the Source Data package builder (derive_verified.py): every
    audited RAW file whose instrument record carries a different specimen number than
    its file name. Unit: acquisition (RAW file); n = 368 audited files.

    Parameters
    ----------
    audit : list of list of str
        ``CAPSCAN_raw_header_audit.tsv``.

    Returns
    -------
    tuple
        ``(rows, n_audited, n_participant15_audited_unchanged)``; rows in the published
        column layout with in_analysis_roster against the 331 roster.
    """
    num = lambda s: s.split("_T760_")[1].split("_")[0]
    a = as_dicts(audit)
    rows = [["acquisition", "in_analysis_roster", "roster_participant", "specimen_number_in_file_name",
             "specimen_number_in_instrument_record", "status"]]
    for r in a:
        if r["header_specimen_number_changed"] == "True":
            subj = r["subject"]
            subj = str(int(float(subj))) if subj else ""
            rows.append([r["sample"], "yes" if r["current_331"] == "True" else "no", subj, num(r["sample"]),
                         num(r["header_original_names"]), "unresolved: analysis labels follow the frozen roster"])
    p15 = [r for r in a if r["subject"] and int(float(r["subject"])) == 15 and r["current_331"] == "False"
           and r["header_status"] == "PARSED" and r["header_specimen_number_changed"] == "False"]
    return rows, len(a), len(p15)


# ---------------------------------------------------------------- v30: the five remaining CapScan rosters on 347
#: v30 re-indexing of the analyses left on 331 in v28f: v28d id -> (v30 id, population, source key, rule).
#: Population 'pm' = participants with protein and metabolite profiles in all five sampling groups;
#: 'dna' = the subset with DNA profiles in all five groups (``v30_capscan_crossomics_participants.tsv``).
V30_REINDEX = {
    "CAPSCAN_FUNCTIONAL_NONDNA14": ("CAPSCAN_FUNCTIONAL_NONDNA15", "pm", "specimens",
                                    "Valid microbial-OG specimen plus metabolite measurement, and all five participant regions"),
    "CAPSCAN_NATIVE12": ("CAPSCAN_NATIVE13", "dna", "participants",
                         "Same thirteen participants have all five eligible DNA regions; protein/metabolite specimens need not be identical to DNA specimens"),
    "CAPSCAN_FUNCTIONAL_DNA12": ("CAPSCAN_FUNCTIONAL_DNA13", "dna", "participants",
                                 "Positive OG/KO/pathway DNA profiles; eligible companion specimens; all five regions in thirteen participants"),
    "CAPSCAN_FUNCTIONAL_NONDNA12_SENS": ("CAPSCAN_FUNCTIONAL_NONDNA13_SENS", "dna", "participants",
                                         "Same thirteen DNA-complete participants; no DNA layer in the compared pair"),
}
#: Analyses whose membership also lists DNA runs (v28d producer: dnaruns= in enroll()).
V30_WITH_DNA_RUNS = ("CAPSCAN_NATIVE13", "CAPSCAN_FUNCTIONAL_DNA13")
#: Source files cited by the v30 rows (names as in v30_inputs_receipt.tsv).
V30_SOURCE_NAMES = {"specimens": "capscan_specimen_roster.tsv", "participants": "COMPLETE.json",
                    "histidine": "manifest_revision4.json"}
#: Reasons, worded as in the published workbook.
V30_SPECIMEN_REASON = "Specimen lacks a metabolite row or any positive microbial-OG acquisition; see specimen eligibility exclusions."
V30_NO_MET = "No metabolite row in the producing companion table."
V30_NO_OG = "No acquisition with positive microbial OG signal for this specimen."
V30_WHOLE_PARTICIPANT = ("Whole participant outside DNA-complete13 design: participant 5 lacks eligible colon DNA; "
                         "participant 7 lacks eligible duodenum DNA.")
V30_PARTICIPANT_EXCL = ("Excluded as a whole participant: no eligible DNA {region} profile after linking positive DNA runs "
                        "to valid companion specimens; complete five-region blocks required.")
V30_RUN_SALIVA = "Saliva sequencing run: outside the gut/stool proteomic cross-omics design; retained in the upstream completed DNA roster."
V30_RUN_SPECIMEN = "DNA specimen is not in the valid microbial-OG plus metabolite specimen set used by cross-omics."
V30_RUN_PARTICIPANT = "Participant excluded from complete five-region native/DNA design (5 lacks colon, 7 lacks duodenum)."
V30_H_SCOPE = "Input acquisition universe; endpoint-specific finite eligibility differs."
V30_H_QUANT_COLUMN = "in_quantitative_nonempty_roster"
V30_DESCRIPTIONS = dict(V28F_DESCRIPTIONS, **{
    "CAPSCAN H acquisitions": "347 histidine search acquisitions; quantitative zero is retained",
})
#: Analyses-row text changes for v30 (analysis id -> {column: text}); counts are computed.
V30_HISTIDINE_ROW = dict(
    inclusion_rule="Exact A20/A27/A28 fixed-347 search acquisition universe; do not apply quantitative exclusions",
    aggregation="Acquisition to specimen to participant-region; paired spatial potential universe has fifteen participants",
    provenance="manifest_revision4.json; complete_statistical_families.tsv")


def _sources(sd):
    """``{file name: sha256}`` of the files read by prepare_v30_inputs.py."""
    return {Path(r["file"]).name: r["sha256"] for r in as_dicts(read_tsv(Path(sd) / "v30_inputs_receipt.tsv"))}


def crossomics_membership(aid, roster, spec_ids, participants, runs, source, rule):
    """Membership rows of one cross-omics analysis in the producer's order (enroll() in build_rosters.py).

    Acquisitions in roster order, specimens by first appearance in the roster, participants
    numerically, participant-region cells participant by participant in region order
    (complete five-region blocks), then DNA runs in roster order.

    Parameters
    ----------
    aid : str
    roster : list of dict
        CAPSCAN acquisitions rows (unit: acquisition).
    spec_ids : set of str
        Included specimens (specimen_id_crossomics).
    participants : list of str
    runs : list of dict or None
        CAPSCAN DNA runs rows to list, or None.
    source : tuple
        ``(source_file, source_sha256)``.
    rule : str

    Returns
    -------
    list of dict
    """
    base = dict(analysis_id=aid, acquisition_id="", specimen_id_biology="", specimen_id_crossomics="",
                participant_id="", region="", run_accession="", source_file=source[0], source_sha256=source[1],
                inclusion_rule=rule)
    acq = [a for a in roster if a["specimen_id_crossomics"] in spec_ids]
    out = [dict(base, entity_type="acquisition", entity_id=a["acquisition_id"], acquisition_id=a["acquisition_id"],
                specimen_id_biology=a["specimen_id_biology"], specimen_id_crossomics=a["specimen_id_crossomics"],
                participant_id=a["participant_id"], region=a["region"]) for a in acq]
    first = {}
    for a in roster:
        first.setdefault(a["specimen_id_crossomics"], a)
    out += [dict(base, entity_type="specimen", entity_id=s, specimen_id_biology=a["specimen_id_biology"],
                 specimen_id_crossomics=s, participant_id=a["participant_id"], region=a["region"])
            for s, a in first.items() if s in spec_ids]
    ps = sorted(participants, key=int)
    out += [dict(base, entity_type="participant", entity_id=p, participant_id=p) for p in ps]
    out += [dict(base, entity_type="participant_region", entity_id=f"{p}|{r}", participant_id=p, region=r)
            for p in ps for r in REGION_ORDINAL]
    if runs is not None:
        bio = {a["specimen_id_crossomics"]: a["specimen_id_biology"] for a in roster}
        dna = [dict(base, entity_type="DNA_run", entity_id=r["run_accession"], run_accession=r["run_accession"],
                    specimen_id_crossomics=r["specimen_id_DNA"], specimen_id_biology=bio.get(r["specimen_id_DNA"], ""),
                    participant_id=r["participant_id_DNA"], region=r["region_DNA"]) for r in runs]
        out += dna
        seen = {}
        for d in dna:   # producer's reader-facing refinement: one DNA_specimen row per DNA specimen, first run's fields
            seen.setdefault(d["specimen_id_crossomics"], dict(d, entity_type="DNA_specimen", entity_id=d["specimen_id_crossomics"], run_accession=""))
        out += list(seen.values())
    return out


def reproduce_v28d_crossomics(v28d):
    """Positive control: the generator reproduces the published 331 rows of the four cross-omics analyses.

    The included specimens, participants and DNA runs are taken from the published rows;
    the control tests row construction and ordering, not the eligibility rules.

    Returns
    -------
    dict
        Per analysis: True if every generated row equals the published row.
    """
    roster = as_dicts(v28d["CAPSCAN acquisitions"])
    mem = as_dicts(v28d["CAPSCAN membership"])
    runs = as_dicts(v28d["CAPSCAN DNA runs"])
    res = {}
    for aid in V30_REINDEX:
        pub = [x for x in mem if x["analysis_id"] == aid]
        spec = {x["entity_id"] for x in pub if x["entity_type"] == "specimen"}
        parts = [x["entity_id"] for x in pub if x["entity_type"] == "participant"]
        rid = {x["entity_id"] for x in pub if x["entity_type"] == "DNA_run"}
        sel = [r for r in runs if r["run_accession"] in rid] if rid else None
        gen = crossomics_membership(aid, roster, spec, parts, sel, (pub[0]["source_file"], pub[0]["source_sha256"]),
                                    pub[0]["inclusion_rule"])
        key = lambda x: {k: x[k] for k in v28d["CAPSCAN membership"][0]}
        types = ("acquisition", "specimen", "participant", "participant_region", "DNA_run", "DNA_specimen")
        res[aid] = all([key(x) for x in gen if x["entity_type"] == t] == [key(x) for x in pub if x["entity_type"] == t] for t in types) \
            and len(gen) == len(pub)
    return res


def v30_populations(v28f, sd):
    """Cross-omics populations on the 347 roster.

    Valid specimens are recomputed (specimen with a positive microbial-OG acquisition in
    ``CAPSCAN_eligibility_audit_347.tsv`` and a metabolite row) and must equal the specimen
    roster of the 347 cross-omics run. DNA runs: gut/stool runs of the 301 whose specimen is
    valid and whose participant is DNA-complete; they must equal the native_profile runs of
    ``source_data/supplementary_methods_8_11/capscan_dna_runs.tsv`` when available.

    Returns
    -------
    dict
        ``valid`` (set of specimens), ``pm`` and ``dna`` (participant lists), ``runs`` (DNA run
        rows), ``met`` (specimens with a metabolite row), ``positive`` (specimens with a
        positive OG acquisition). Units: specimen, participant, DNA run.
    """
    sd = Path(sd)
    roster = as_dicts(v28f["CAPSCAN acquisitions"])
    elig = {r["acquisition"]: r["current_positive"] == "True" for r in as_dicts(read_tsv(sd / "CAPSCAN_eligibility_audit_347.tsv"))}
    positive = {a["specimen_id_crossomics"] for a in roster if elig[a["acquisition_id"]]}
    met = {r["specimen_id_crossomics"] for r in as_dicts(read_tsv(sd / "v30_capscan_metabolite_specimens.tsv"))}
    valid = positive & met
    saved = {r["specimen_id_crossomics"] for r in as_dicts(read_tsv(sd / "v30_capscan_crossomics_specimens.tsv"))}
    assert valid == saved, ("valid specimens differ from the 347 cross-omics roster", len(valid), len(saved))
    parts = as_dicts(read_tsv(sd / "v30_capscan_crossomics_participants.tsv"))
    pm = [p["participant_id"] for p in parts if p["protein_metabolite_complete"] == "True"]
    dna = [p["participant_id"] for p in parts if p["dna_complete"] == "True"]
    region = {a["specimen_id_crossomics"]: (a["participant_id"], a["region"]) for a in roster}
    for p in pm:   # the rule the participant list encodes: valid specimens in all five groups
        assert {region[s][1] for s in valid if region[s][0] == p} == set(REGION_ORDINAL), p
    runs = [r for r in as_dicts(v28f["CAPSCAN DNA runs"])
            if r["region_DNA"] != "saliva" and r["specimen_id_DNA"] in valid and r["participant_id_DNA"] in dna]
    for p in dna:
        assert {r["region_DNA"] for r in runs if r["participant_id_DNA"] == p} == set(REGION_ORDINAL), p
    native = Path(sd).parent / "supplementary_methods_8_11" / "capscan_dna_runs.tsv"
    if native.exists():
        nat = {r["run_accession"] for r in as_dicts(read_tsv(native)) if r["native_profile"] == "True"}
        assert nat == {r["run_accession"] for r in runs}, "DNA runs differ from the native_profile runs"
    return dict(valid=valid, pm=pm, dna=dna, runs=runs, met=met, positive=positive)


def build_v30(v28f, sd):
    """Apply the v30 CapScan re-indexing to the v28f data sheets.

    Re-indexes CAPSCAN_FUNCTIONAL_NONDNA14 -> NONDNA15, NATIVE12 -> NATIVE13,
    FUNCTIONAL_DNA12 -> DNA13, FUNCTIONAL_NONDNA12_SENS -> NONDNA13_SENS, and
    CAPSCAN_HISTIDINE_ROUTES with its three histidine sheets, plus CAPSCAN DNA eligibility,
    on the 347-acquisition roster. Exclusions and MS nonmembership of these analyses are
    regenerated with the published reason texts.

    Returns
    -------
    dict
        Data sheets of the v30 build.
    """
    sd = Path(sd)
    s = {k: [list(r) for r in v] for k, v in v28f.items()}
    roster = as_dicts(v28f["CAPSCAN acquisitions"])
    pop = v30_populations(v28f, sd)
    src = _sources(sd)
    source = {k: (n, src[n]) for k, n in V30_SOURCE_NAMES.items()}
    excluded_participants = sorted({a["participant_id"] for a in roster} - set(pop["dna"]), key=int)
    missing = {}
    for p in excluded_participants:
        have = {r["region_DNA"] for r in pop["runs"] if r["participant_id_DNA"] == p}
        allruns = {r["region_DNA"] for r in as_dicts(v28f["CAPSCAN DNA runs"])
                   if r["participant_id_DNA"] == p and r["specimen_id_DNA"] in pop["valid"]}
        missing[p] = [r for r in REGION_ORDINAL if r not in allruns]
    assert missing == {"5": ["colon"], "7": ["duodenum"]}, missing

    # membership
    mem_h = v28f["CAPSCAN membership"][0]
    new_mem, done, gen = [], set(), {}
    for aid, (nid, popkey, skey, rule) in V30_REINDEX.items():
        parts = pop[popkey]
        spec = {x for x in pop["valid"] if any(a["specimen_id_crossomics"] == x and a["participant_id"] in parts for a in roster)}
        runs = pop["runs"] if nid in V30_WITH_DNA_RUNS else None
        gen[aid] = crossomics_membership(nid, roster, spec, parts, runs, source[skey], rule)
    for r in as_dicts(v28f["CAPSCAN membership"]):
        aid = r["analysis_id"]
        if aid in V30_REINDEX:   # main block at its first row; DNA_specimen rows where the published ones stand
            block = "spec" if r["entity_type"] == "DNA_specimen" else "main"
            if (aid, block) not in done:
                new_mem += [x for x in gen[aid] if (x["entity_type"] == "DNA_specimen") == (block == "spec")]
                done.add((aid, block))
            continue
        new_mem.append(r)
    for aid in V30_REINDEX:
        assert (aid, "main") in done and ((aid, "spec") in done or not any(x["entity_type"] == "DNA_specimen" for x in gen[aid]))
    s["CAPSCAN membership"] = from_dicts(mem_h, new_mem)
    mem = as_dicts(s["CAPSCAN membership"])
    inc = lambda aid: {x["acquisition_id"] for x in mem if x["analysis_id"] == aid and x["entity_type"] == "acquisition"}

    # MS nonmembership
    non_h = v28f["CAPSCAN MS nonmembership"][0]
    new_non, done = [], set()
    for r in as_dicts(v28f["CAPSCAN MS nonmembership"]):
        aid = r["analysis_id"]
        if aid in V30_REINDEX:
            if aid not in done:
                nid = V30_REINDEX[aid][0]
                ids = inc(nid)
                parts = set(pop[V30_REINDEX[aid][1]])
                for a in roster:
                    if a["acquisition_id"] in ids:
                        continue
                    reason = V30_WHOLE_PARTICIPANT if a["participant_id"] not in parts else V30_SPECIMEN_REASON
                    new_non.append(dict(analysis_id=nid, acquisition_id=a["acquisition_id"], specimen_id_biology=a["specimen_id_biology"],
                                        specimen_id_crossomics=a["specimen_id_crossomics"], participant_id=a["participant_id"],
                                        region=a["region"], disposition_type="not_in_analysis", reason=reason))
                done.add(aid)
            continue
        new_non.append(r)
    s["CAPSCAN MS nonmembership"] = from_dicts(non_h, new_non)

    # exclusions
    ex_h = v28f["CAPSCAN exclusions"][0]
    blank = {k: "" for k in ex_h}
    new_ex, done = [], set()
    first = {}
    for a in roster:
        first.setdefault(a["specimen_id_crossomics"], a)
    for r in as_dicts(v28f["CAPSCAN exclusions"]):
        aid = r["analysis_id"]
        if aid in V30_REINDEX:
            if aid not in done:
                nid, popkey, skey, _ = V30_REINDEX[aid]
                if popkey == "dna":
                    for p in excluded_participants:
                        new_ex.append(dict(blank, analysis_id=nid, entity_type="participant", entity_id=p,
                                           reason=V30_PARTICIPANT_EXCL.format(region=missing[p][0]), source_file=source[skey][0],
                                           source_sha256=source[skey][1], missing_region=missing[p][0], disposition_type="analysis_eligibility"))
                else:
                    for sid, a in first.items():
                        if sid in pop["valid"]:
                            continue
                        new_ex.append(dict(blank, analysis_id=nid, entity_type="specimen", entity_id=sid,
                                           reason=V30_NO_MET if sid not in pop["met"] else V30_NO_OG,
                                           source_file=source[skey][0], source_sha256=source[skey][1],
                                           participant_id=a["participant_id"], region=a["region"], disposition_type="analysis_eligibility"))
                if nid == "CAPSCAN_NATIVE13":
                    keep = {x["run_accession"] for x in pop["runs"]}
                    for x in as_dicts(v28f["CAPSCAN DNA runs"]):
                        if x["run_accession"] in keep:
                            continue
                        reason = (V30_RUN_SALIVA if x["region_DNA"] == "saliva" else V30_RUN_SPECIMEN
                                  if x["specimen_id_DNA"] not in pop["valid"] else V30_RUN_PARTICIPANT)
                        new_ex.append(dict(blank, analysis_id=nid, entity_type="DNA_run", entity_id=x["run_accession"], reason=reason,
                                           source_file=source[skey][0], source_sha256=source[skey][1],
                                           participant_id=x["participant_id_DNA"], region=x["region_DNA"],
                                           specimen_id=x["specimen_id_DNA"], disposition_type="not_in_downstream_analysis"))
                done.add(aid)
            continue
        new_ex.append(r)
    s["CAPSCAN exclusions"] = from_dicts(ex_h, new_ex)

    # DNA eligibility (participant x region x representation), all fifteen participants
    el_h = [c if c != "in_complete12" else "in_complete13" for c in v28f["CAPSCAN DNA eligibility"][0]]
    gut = [r for r in as_dicts(v28f["CAPSCAN DNA runs"]) if r["region_DNA"] != "saliva" and r["specimen_id_DNA"] in pop["valid"]]
    rows = []
    for p in sorted({a["participant_id"] for a in roster}, key=int):
        for reg in REGION_ORDINAL:
            sel = [r for r in gut if r["participant_id_DNA"] == p and r["region_DNA"] == reg]
            for level in ("OG", "KO", "Pathway"):
                rows.append(dict(participant_id=p, region=reg, representation=level, eligible_DNA_runs=str(len(sel)),
                                 eligible_DNA_specimens=str(len({r["specimen_id_DNA"] for r in sel})),
                                 DNA_region_available=str(bool(sel)), in_complete13=str(p in pop["dna"])))
    s["CAPSCAN DNA eligibility"] = from_dicts(el_h, rows)

    # histidine sheets
    hd = {r["acquisition_id"]: r for r in as_dicts(read_tsv(sd / "v30_capscan_h_design.tsv"))}
    assert set(hd) == {a["acquisition_id"] for a in roster}
    nonempty = inc("CAPSCAN_QUANT_NONEMPTY")
    old_h = {r["acquisition_id"]: r for r in as_dicts(v28f["CAPSCAN H acquisitions"])}
    h_h = [c if c != "in_quantitative330_roster" else V30_H_QUANT_COLUMN for c in v28f["CAPSCAN H acquisitions"][0]]
    hrows = []
    for a in roster:
        d = hd[a["acquisition_id"]]
        assert (d["participant_id"], d["region"], d["histidine_specimen_id"]) == (a["participant_id"], a["region"], a["specimen_id_biology"]), a
        o = old_h.get(a["acquisition_id"])
        hrows.append({"acquisition_id": a["acquisition_id"], "histidine_specimen_id": d["histidine_specimen_id"],
                      "specimen_id_crossomics": a["specimen_id_crossomics"], "participant_id": a["participant_id"],
                      "region": a["region"], "aliquot": a["aliquot"], V30_H_QUANT_COLUMN: str(a["acquisition_id"] in nonempty),
                      "histidine_scope": V30_H_SCOPE, "source_file": o["source_file"] if o else source["histidine"][0]})
    s["CAPSCAN H acquisitions"] = from_dicts(h_h, hrows)
    s["CAPSCAN H finite n"] = read_tsv(sd / "v30_capscan_h_finite_n.tsv")
    s["CAPSCAN H participants"] = read_tsv(sd / "v30_capscan_h_participants.tsv")

    # Analyses
    an_h = v28f["Analyses"][0]
    count = lambda aid, et: str(sum(1 for x in mem if x["analysis_id"] == aid and x["entity_type"] == et))
    rows = []
    for r in as_dicts(v28f["Analyses"]):
        aid = r["analysis_id"]
        if aid in V30_REINDEX:
            nid, popkey, skey, rule = V30_REINDEX[aid]
            r = dict(r, analysis_id=nid, n_MS_acquisitions=count(nid, "acquisition"),
                     n_recorded_proteomic_specimens=count(nid, "specimen"), n_participants=count(nid, "participant"),
                     n_participant_region_cells=count(nid, "participant_region"),
                     roster_filter=r["roster_filter"].replace(aid, nid), exclusion_filter=r["exclusion_filter"].replace(aid, nid),
                     inclusion_rule=rule, provenance=source[skey][0],
                     limitations=r["limitations"].replace(V28F_NOT_REINDEXED, ""))
            if nid in V30_WITH_DNA_RUNS:
                r["n_DNA_runs"] = r["n_DNA_specimens"] = count(nid, "DNA_run")
            if nid == "CAPSCAN_NATIVE13":
                r["limitations"] = (f"{r['n_recorded_proteomic_specimens']} protein/metabolite specimens and {r['n_DNA_runs']} DNA "
                                    "specimens are separate layer inputs, combined at participant-region; do not assert "
                                    "specimen-level complete cases.")
        elif aid == "CAPSCAN_HISTIDINE_ROUTES":
            h = as_dicts(s["CAPSCAN H acquisitions"])
            r = dict(r, n_MS_acquisitions=str(len(h)), n_recorded_proteomic_specimens=str(len({x["histidine_specimen_id"] for x in h})),
                     n_participants=str(len({x["participant_id"] for x in h})),
                     n_participant_region_cells=str(len({(x["participant_id"], x["region"]) for x in h})),
                     limitations=r["limitations"].replace(V28F_NOT_REINDEXED, ""), **V30_HISTIDINE_ROW)
        rows.append(r)
    s["Analyses"] = from_dicts(an_h, rows)
    return s


def build_v30_all(v28f_all, sd):
    """All 34 sheets of the v30 workbook (every CapScan analysis on the 347-acquisition roster)."""
    data = {k: v for k, v in v28f_all.items() if k not in COMPUTED_SHEETS}
    sheets = build_v30(data, sd)
    sheets["Cohorts"] = cohorts_sheet(sheets, V28F_COHORT_DESIGN)
    sheets["Sheet directory"] = sheet_directory(SHEETS, sheets, V30_DESCRIPTIONS, [])
    return sheets


def stat_v30(sheets):
    """Roster counts of the v30 workbook.

    Returns
    -------
    dict
        Acquisitions per cohort and in total (unit: acquisition), and per re-indexed CapScan
        analysis its acquisitions, specimens, participants, cells and DNA runs.
    """
    coh = {r[0]: int(r[1]) for r in cohorts_sheet(sheets)[1:]}
    n = {f"acquisitions_{k}": v for k, v in coh.items()}
    n["acquisitions_total"] = sum(coh.values())
    for r in as_dicts(sheets["Analyses"]):
        if r["cohort"] == "CapScan":
            n[r["analysis_id"]] = "/".join(r[c] for c in ("n_MS_acquisitions", "n_recorded_proteomic_specimens", "n_participants",
                                                          "n_participant_region_cells", "n_DNA_runs"))
    return n


def receipt_v30(sheets, complete, changes, controls, receipt):
    """Write the v30 build receipt (per-sheet digests, change summary from v28f, controls)."""
    import json
    dig = {}
    for name, _, _ in SHEETS:
        rows = sheets[name]
        h, b, _ = split_public(name, rows[0], rows[1:])
        dig[name] = dict(data_rows=len(rows) - 1, sha256_public=grid_digest([h] + b),
                         sha256_full=grid_digest(rows) if complete else None)
    rec = dict(unit="SData1", build="v30", mode="complete" if complete else "public (withheld PREDICT cells absent)",
               positive_controls=controls, changes_from_v28f=changes, counts=stat_v30(sheets), sheets=dig)
    Path(receipt).write_text(json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    return rec
