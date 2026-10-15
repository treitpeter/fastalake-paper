#!/usr/bin/env python3
"""Scan the repository for content that must not be released.

Checks every git-tracked file (or every file under the repository if git lists none)
for:

* absolute lab or user paths (``/fs/``, ``/Volumes/``, ``/home/``, ``/Users/``);
* e-mail addresses;
* personal-genome identifiers (``Probant``, ``NBR`` followed by digits);
* clinical column names in tables (``.tsv``, ``.csv``, ``.txt``, ``.xlsx`` headers):
  ``idpatient``, ``Tube_ID``, ``meld``, ``crp``, ``antibiotic``, ``age``, ``sex``,
  ``visit`` (whole-word, case-insensitive, also as prefix, e.g. ``age_years``);
* files larger than 20 MB.

Lines of a notebook that hold only a stored base64 image payload are not text-scanned
(random base64 can spell a pattern such as ``nbr4``).

Files listed in ``ALLOW_PATH_FILES`` may contain the text patterns because they define
them (this script). ``PROVENANCE.md`` files may name original producer
locations only as ``LAB_ROOT/...``; absolute paths are flagged there too.

Usage
-----
    python tools/scan_release.py            # scan the repository, exit 1 on findings
    python tools/scan_release.py --self-test  # plant one error of each kind, expect each caught
"""

import argparse
import csv
import hashlib
import io
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

#: Size limit for any tracked file (bytes). Larger files go to Zenodo.
MAX_BYTES = 20 * 1024 * 1024

#: Text patterns that must not appear in any released file.
TEXT_PATTERNS = {
    "absolute_path": re.compile(r"(?<![A-Za-z0-9_.])/(?:fs|Volumes|home|Users)/"),
    "email": re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9-]+(?:\.[A-Za-z0-9-]+)*\.[A-Za-z]{2,}"),
    "personal_genome": re.compile(r"Probant|\bNBR\d|NBR66", re.IGNORECASE),
    "private_key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    "access_token": re.compile(r"\b(?:gh[pousr]_[A-Za-z0-9]{36,}|github_pat_[A-Za-z0-9_]{40,}|AKIA[A-Z0-9]{16})\b"),
}

#: Clinical covariate names that must not be table columns.
CLINICAL_COLUMNS = re.compile(
    r"^(?:idpatient|raw_patient_id|tube_?id|meld|crp|antibiotics?|age|sex|visit|clinical_group)(?:$|[_\s.\-])", re.IGNORECASE
)

#: Table formats whose header row is checked for clinical column names.
TABLE_SUFFIXES = {".tsv", ".csv", ".txt", ".tab"}

#: E-mail-like strings that are allowed (none of them identify a person).
EMAIL_ALLOW = {"noreply@github.com"}

#: Files allowed to contain the patterns because they define or test them.
ALLOW_PATH_FILES = {"tools/scan_release.py"}

# This published workbook contains eight cohort-mean composition rows (four
# exposure definitions x on/off), with no participant records. Permit only its
# reviewed aggregate header, and only while the entire workbook is unchanged.
AGGREGATE_WORKBOOK_HEADERS = {
    ('tables/STab3/expected/Supplementary_Table_3.xlsx', 'Composition', 'antibiotic_state'):
        '5d9815c575a8cc15f9f840994b2bf255fc95f1c0b3ddbb4aaf218c49a39a0f3d',
}

#: Notebook lines that hold only a base64 image payload (stored outputs); random base64 can
#: spell any pattern (e.g. "nbr4"), so these lines are not text-scanned.
NOTEBOOK_IMAGE_LINE = re.compile(r'^\s*"image/(?:png|jpeg|gif)": "[A-Za-z0-9+/=\\n]*",?\s*$')

#: Binary formats skipped for the text scan (size is still checked).
BINARY_SUFFIXES = {".png", ".pdf", ".jpg", ".jpeg", ".tif", ".tiff", ".gz", ".zip",
                   ".parquet", ".xlsx", ".otf", ".ttf", ".woff", ".woff2", ".svgz"}


def tracked_files(root):
    """Return repository files to scan, relative to ``root``.

    Parameters
    ----------
    root : pathlib.Path
        Repository root.

    Returns
    -------
    list of pathlib.Path
        Git-tracked plus untracked-but-not-ignored files; all files when not a git tree.
    """
    try:
        out = subprocess.run(
            ["git", "-C", str(root), "ls-files", "-co", "--exclude-standard", "-z"],
            check=True, capture_output=True,
        ).stdout.decode()
        files = [Path(p) for p in out.split("\0") if p]
        if files:
            return files
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass
    return [p.relative_to(root) for p in root.rglob("*") if p.is_file() and ".git" not in p.parts]


def table_header(path):
    """Return the header fields of a delimited text table, or an empty list."""
    try:
        with open(path, newline="", encoding="utf-8", errors="replace") as fh:
            first = fh.readline()
    except OSError:
        return []
    if first.startswith("#"):
        return []
    delim = "\t" if "\t" in first else ","
    return next(csv.reader(io.StringIO(first), delimiter=delim), [])


def scan(root):
    """Scan ``root`` and return a list of findings.

    Parameters
    ----------
    root : pathlib.Path
        Directory to scan.

    Returns
    -------
    list of tuple
        ``(relative_path, line_number, kind, matched_text)``; line 0 for whole-file findings.
    """
    findings = []
    for rel in tracked_files(root):
        path = root / rel
        if not path.is_file():
            continue
        size = path.stat().st_size
        if size > MAX_BYTES:
            findings.append((str(rel), 0, "file_over_20MB", f"{size} bytes"))
        if path.suffix.lower() == '.xlsx':
            from openpyxl import load_workbook
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            book = load_workbook(path, read_only=True, data_only=False)
            try:
                for sheet in book:
                    for rownum, row in enumerate(sheet.iter_rows(values_only=True), 1):
                        for value in row:
                            if not isinstance(value, str):
                                continue
                            if rownum == 1 and CLINICAL_COLUMNS.match(value.strip()):
                                reviewed = AGGREGATE_WORKBOOK_HEADERS.get((rel.as_posix(), sheet.title, value.strip()))
                                if reviewed != digest:
                                    findings.append((str(rel), rownum, 'clinical_column', sheet.title + ': ' + value))
                            for kind, pat in TEXT_PATTERNS.items():
                                for match in pat.finditer(value):
                                    if kind == 'email' and match.group(0).lower() in EMAIL_ALLOW:
                                        continue
                                    findings.append((str(rel), rownum, kind, sheet.title + ': ' + match.group(0)))
            finally:
                book.close()
        if path.suffix.lower() in TABLE_SUFFIXES:
            for col in table_header(path):
                if CLINICAL_COLUMNS.match(col.strip()):
                    findings.append((str(rel), 1, "clinical_column", col.strip()))
        if path.suffix.lower() in BINARY_SUFFIXES or rel.as_posix() in ALLOW_PATH_FILES:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        for lineno, line in enumerate(text.splitlines(), 1):
            if path.suffix == ".ipynb" and NOTEBOOK_IMAGE_LINE.match(line):
                continue
            for kind, pat in TEXT_PATTERNS.items():
                for m in pat.finditer(line):
                    if kind == "email" and m.group(0).lower() in EMAIL_ALLOW:
                        continue
                    findings.append((str(rel), lineno, kind, m.group(0)))
    return findings


def self_test():
    """Plant one violation of each kind in a temporary tree and require each be caught.

    Returns
    -------
    int
        0 if every planted error was found and the clean control produced no finding.
    """
    planted = {
        "absolute_path": ("notes.md", "data at /fs/pool/x.tsv\n"),
        "email": ("contact.md", "write to someone@example.org\n"),
        "personal_genome": ("genome.md", "sample Probant5 was used\n"),
        "absolute_path_in_notebook_text": ("nb.ipynb", '{\n "text/plain": ["/fs/pool/x"]\n}\n'),
        "clinical_column": ("source_data/X/t.tsv", "pid\tMELD\tvalue\n1\t2\t3\n"),
        "private_key": ("credential.key", "-----BEGIN PRIVATE KEY-----\n"),
        "access_token": ("credential.txt", "ghp_" + "x" * 36),
        "file_over_20MB": ("big.bin", None),
    }
    ok = True
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        (root / "clean.md").write_text("relative path figures/Fig1/README.md and LAB_ROOT/x\n")
        (root / "clean.ipynb").write_text('{\n "data": {\n  "image/png": "iVBORnbr4AAA/fs/xyz=",\n  "text/plain": ["x"]\n }\n}\n')
        (root / "source_data" / "Y").mkdir(parents=True)
        (root / "source_data" / "Y" / "ok.tsv").write_text("patient_pseudo_id\tvalue\tstage\n")
        clean = scan(root)
        if clean:
            print("SELF-TEST FAIL: clean control produced findings:", clean)
            ok = False
        for kind, (name, content) in planted.items():
            p = root / name
            p.parent.mkdir(parents=True, exist_ok=True)
            if content is None:
                with open(p, "wb") as fh:
                    fh.seek(MAX_BYTES + 1)
                    fh.write(b"\0")
            else:
                p.write_text(content)
        res = scan(root)
        found = {f[2] for f in res} | {"absolute_path_in_notebook_text" for f in res if f[0] == "nb.ipynb" and f[2] == "absolute_path"}
        for kind in planted:
            status = "caught" if kind in found else "MISSED"
            ok &= kind in found
            print(f"self-test {kind:16s} {status}")
        from openpyxl import Workbook
        book = Workbook()
        book.active.append(['clinical_group', 'value'])
        book.create_sheet('hidden').sheet_state = 'hidden'
        book['hidden'].append(['id', 'value'])
        book['hidden'].append(['x', '/home/private/data'])
        book.save(root / 'table.xlsx')
        workbook_findings = {r[2] for r in scan(root) if r[0] == 'table.xlsx'}
        workbook_ok = {'clinical_column', 'absolute_path'} <= workbook_findings
        ok &= workbook_ok
        print('self-test workbook_headers_and_hidden_cells', 'caught' if workbook_ok else 'MISSED')
        # An aggregate exception must stop applying as soon as workbook content changes.
        sample = root / 'aggregate.xlsx'
        book = Workbook()
        book.active.title = 'Summary'
        book.active.append(['antibiotic_state', 'mean'])
        book.active.append(['On', 0.3])
        book.save(sample)
        key = ('aggregate.xlsx', 'Summary', 'antibiotic_state')
        AGGREGATE_WORKBOOK_HEADERS[key] = hashlib.sha256(sample.read_bytes()).hexdigest()
        try:
            clean_aggregate = not any(r[0] == 'aggregate.xlsx' for r in scan(root))
            book.active.append(['Off', 0.2])
            book.save(sample)
            changed_aggregate = any(r[0] == 'aggregate.xlsx' and r[2] == 'clinical_column' for r in scan(root))
            aggregate_ok = clean_aggregate and changed_aggregate
            ok &= aggregate_ok
            print('self-test aggregate_checksum_boundary', 'caught' if aggregate_ok else 'MISSED')
        finally:
            del AGGREGATE_WORKBOOK_HEADERS[key]
    print("SELF-TEST PASS" if ok else "SELF-TEST FAIL")
    return 0 if ok else 1


def main(argv=None):
    """Command-line entry point."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--self-test", action="store_true", help="run the planted-error test")
    ap.add_argument("--root", default=str(REPO), help="directory to scan (default: repo)")
    args = ap.parse_args(argv)
    if args.self_test:
        return self_test()
    findings = scan(Path(args.root))
    for rel, line, kind, text in findings:
        print(f"{rel}:{line}: {kind}: {text}")
    print(f"{len(findings)} finding(s)")
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
