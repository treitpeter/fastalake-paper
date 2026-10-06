#!/usr/bin/env python3
"""Copy the inputs of Supplementary Data 1 into ``source_data/SData1/`` (run once, by the authors).

This script is the only code in the unit that reads outside the repository. It reads the
published Source Data package copy of Supplementary Data 1 and the participant-15 roster
fix (upstream ``capscan_participant_roster`` and ``capscan_crossomics``)
from the lab project directory named by ``LAB_ROOT``, and writes

* ``source_data/SData1/v28d_<sheet>.tsv``: one public table per data sheet of the
  published workbook, cell values unchanged, with the PREDICT columns and rows listed in
  ``panels.WITHHELD_COLUMNS`` and ``panels.WITHHELD_ROWS`` removed;
* ``source_data/SData1/CAPSCAN_participant15_join.tsv``,
  ``CAPSCAN_eligibility_audit_347.tsv`` and ``CapScan_specimen_number_discrepancies.tsv``;
* ``source_data/SData1/CAPSCAN_raw_header_audit.tsv``: the RAW header audit of all 368
  CapScan RAW files (columns in ``RAW_AUDIT_COLUMNS``; local file paths and checksums dropped);
* ``tables/SData1/expected/``: the published ``ANALYSIS_INDEX.tsv`` (byte copy), the
  sha256 of the three published files, a digest of every published sheet, and the two
  sheets the code computes (Cohorts, Sheet directory) as published;
* ``<CONTROLLED_OUT>/SData1_PREDICT_withheld_cells.tsv``: the withheld PREDICT cells,
  one row per cell, for the controlled-access directory (never written into the repository).

Usage
-----
    LAB_ROOT=<lab project dir> python tables/SData1/prepare_source_data.py <CONTROLLED_OUT>
"""

import csv
import hashlib
import os
import shutil
import sys
from pathlib import Path

UNIT = Path(__file__).resolve().parent
REPO = UNIT.parents[1]
sys.path[:0] = [str(UNIT)]
import panels  # noqa: E402

LAB_ROOT = Path(os.environ["LAB_ROOT"])
FINAL = ("review/v15_integration/jess_integration_2026-09-29/analysis_followup/"
         "manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03")
PKG_SD1 = LAB_ROOT / FINAL / "final/figure_sd/source_data/package/FastaLake_Source_Data/Supplementary_Data_1"
C2 = LAB_ROOT / FINAL / "coauthor_comment_check_v25_2026-10-05/compute/C2_participant15"
RAW_AUDIT = (LAB_ROOT / "review/v15_integration/jess_integration_2026-09-29/cohort_reconciliation/"
             "raw_identity_audit/all_RAW_header_names.tsv")
#: Columns of the RAW header audit kept (file paths and checksums of local copies dropped).
RAW_AUDIT_COLUMNS = ["sample", "current_331", "subject", "header_original_names",
                     "header_specimen_number_changed", "header_status"]

#: Published files and their sha256 as recorded in the unit inventory and Codex census.
PUBLISHED = {
    "Sample_rosters.xlsx": "5ee76e28adb77d2edd9d424cf450538f72a8a6dc4f66bd07c091f524bafcff72",
    "ANALYSIS_INDEX.tsv": "0e10a9bc4d1dd724592a37bac503df5e156dfb67ecf8d1fe79d05b5c8b7be5cd",
    "CapScan_specimen_number_discrepancies.tsv": "2f4a14503dd71c1d5ce7440949f2406eefdc36efb1ce0910841de22ba83087e6",
}

#: Upstream files copied unchanged (source path relative to C2 -> source_data name).
C2_COPIES = {
    "roster/participant15_added.tsv": "CAPSCAN_participant15_join.tsv",
    "fig5g/crossomics_347/eligibility_audit.tsv": "CAPSCAN_eligibility_audit_347.tsv",
}


def sha(p):
    """sha256 of a file."""
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    """Write source_data, expected and the controlled extract."""
    controlled_out = Path(sys.argv[1])
    controlled_out.mkdir(parents=True, exist_ok=True)
    sd = REPO / "source_data" / "SData1"
    exp = UNIT / "expected"
    sd.mkdir(parents=True, exist_ok=True)
    exp.mkdir(exist_ok=True)
    for name, digest in PUBLISHED.items():
        assert sha(PKG_SD1 / name) == digest, name

    sheets = panels.read_xlsx_values(PKG_SD1 / "Sample_rosters.xlsx")
    withheld = []
    digests = []
    for name, rows in sheets.items():
        header, body = rows[0], rows[1:]
        assert all(len(r) == len(header) for r in body), name
        pub_header, pub_body, cells = panels.split_public(name, header, body)
        withheld += cells
        digests.append(dict(sheet=name, data_rows=len(body), columns=len(header),
                            sha256_full=panels.grid_digest(rows),
                            sha256_public=panels.grid_digest([pub_header] + pub_body)))
        if name in panels.COMPUTED_SHEETS:
            panels.write_tsv(exp / f"{panels.slug(name)}_v28d.tsv", rows)
            continue
        panels.write_tsv(sd / panels.sheet_file(name), [pub_header] + pub_body)

    panels.write_tsv(exp / "SHEET_DIGESTS_v28d.tsv",
                     [list(digests[0])] + [[str(d[k]) for k in digests[0]] for d in digests])
    panels.write_tsv(exp / "PUBLISHED_FILES_v28d.tsv",
                     [["file", "bytes", "sha256"]]
                     + [[n, str((PKG_SD1 / n).stat().st_size), d] for n, d in PUBLISHED.items()])
    shutil.copyfile(PKG_SD1 / "ANALYSIS_INDEX.tsv", exp / "ANALYSIS_INDEX_v28d.tsv")
    shutil.copyfile(PKG_SD1 / "CapScan_specimen_number_discrepancies.tsv",
                    sd / "CapScan_specimen_number_discrepancies.tsv")
    for src, dst in C2_COPIES.items():
        shutil.copyfile(C2 / src, sd / dst)
    audit = panels.read_tsv(RAW_AUDIT)
    keep = [audit[0].index(c) for c in RAW_AUDIT_COLUMNS]
    panels.write_tsv(sd / "CAPSCAN_raw_header_audit.tsv", [[r[j] for j in keep] for r in audit])

    panels.write_tsv(controlled_out / panels.CONTROLLED_FILE,
                     [["sheet", "row_index", "column_index", "column", "value"]]
                     + [[str(x) for x in c] for c in withheld])
    print(f"sheets {len(sheets)}; withheld cells {len(withheld)}; "
          f"controlled sha256 {sha(controlled_out / panels.CONTROLLED_FILE)}")


if __name__ == "__main__":
    main()
