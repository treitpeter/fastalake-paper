#!/usr/bin/env python3
"""Check registered manuscript numbers against recorded recomputed values.

Reads ``figures/*/NUMBERS.tsv``, ``tables/*/NUMBERS.tsv`` and
``text_statistics/*/NUMBERS.tsv``. Each row has the columns (tab-separated, header
required, see CONVENTIONS.md)::

    v28d_line  quote  printed_value  recomputed_value  code  source_data  status

A row passes when

* ``status`` is ``MATCH`` and ``recomputed_value`` agrees with ``printed_value`` at the
  printed precision, or
* ``status`` is ``v28f_change:<new>`` and ``recomputed_value`` agrees with ``<new>`` at the
  precision of ``<new>`` (the manuscript value changes in v28f).

* ``status`` is ``SUPERSEDED_V33`` (optionally ``SUPERSEDED_V33:<previous status>``, e.g.
  ``SUPERSEDED_V33:v28f_change:1.40``): the v28d (or v28f) printed value was superseded by a v33 value
  computed on a changed roster (CapScan matched benchmark 282 -> 295 acquisitions); the row is
  then checked only through its v42 columns, which must be present.

``OPEN`` rows are listed and do not fail unless ``--strict``. Any other status, a row whose
status says MATCH but whose values disagree, or a malformed row, fails.

Agreement rule: the printed value is parsed as a number (thousands separators, ``%``,
``×``/``x``/``-fold`` suffixes, ``a × 10^b`` / ``a × 10−b`` / ``aeb`` notation and
Unicode minus are accepted). The recomputed value is rounded to the printed number of
decimals (or significant figures for scientific notation) and must then be equal.
Superscript exponents (``6.1 × 10⁻⁴``) are read as ordinary exponents. A printed integer
with trailing zeros (``250,000``) that fails the exact rule is accepted if the recomputed
value rounds to it at its significant figures (``250,197`` -> 2 s.f.); such rows are
reported as ROUNDED so the reader sees every one; the same rule and listing apply to
``v28f_change:<new>`` values and to ``v42_printed_value``. A printed bound ``< v``, ``> v``,
``more than v`` or ``fewer than v`` requires the recomputed value to satisfy it.
English number words zero to twenty are read as integers. Other non-numeric values are compared
as exact strings.

Current reference: v42. The four retained ``v36_*`` columns preserve the historical
anchors. Seven ``v42_*`` fields follow: ``para``, ``quote``, ``printed_value``,
``status``, ``recomputed_value``, ``code`` and ``source_data``. Current values and
input populations may differ from the historical columns. The current manuscript
has 494 paragraphs; paragraph indices are zero-based. MATCH checks recorded values
and, when supplied, exact manuscript quotations. OPEN records an unresolved claim;
NOT_IN_V42 records an unlocated quote. REMOVED_V30/V33/V36/V42 rows preserve an
explicit disposition for a claim removed from the corresponding manuscript.
The original tracked manuscript is never modified. Its accepted-text view and
identity are specified by ``docs/MANUSCRIPT.json``.

With ``--manuscript-text PATH`` (or env ``FASTALAKE_MANUSCRIPT_TXT``) the tool also checks
quotes. The format is detected from the file: a paragraph-indexed text (every non-empty
line ``<index>\t<paragraph>``, as ``visible_v42.txt``) checks ``v42_quote`` inside
paragraph ``v42_para`` and ``v42_printed_value`` inside ``v42_quote``; files without v42
columns are then reported as not checked against v42. Any other file is read as the
accepted v28d text with physical line numbers, and ``quote`` is checked on ``v28d_line``
and ``printed_value`` inside ``quote``. The manuscript text is not part of this repository.

Usage
-----
    python tools/check_numbers.py [--strict] [--unit Fig2] [--manuscript-text PATH]
"""

import argparse
import csv
import hashlib
import json
import math
import os
import re
import sys
from decimal import ROUND_HALF_UP, Decimal, InvalidOperation
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

#: Required NUMBERS.tsv columns, in order.
COLUMNS = ["v28d_line", "quote", "printed_value", "recomputed_value", "code", "source_data", "status"]

#: Optional reference-text (v42) columns, appended after COLUMNS.
REF = "v42"
REF_COLUMNS = [f"{REF}_{key}" for key in ("para", "quote", "printed_value", "status", "recomputed_value", "code", "source_data")]
HISTORICAL_REF_COLUMNS = [f"v36_{key}" for key in ("para", "quote", "printed_value", "status")]

#: Allowed v42_status values.
REF_STATUSES = ("MATCH", "OPEN", "NOT_IN_V42", "REMOVED_V30", "REMOVED_V33", "REMOVED_V36", "REMOVED_V42")

#: Glob patterns for NUMBERS.tsv files.
GLOBS = ["figures/*/NUMBERS.tsv", "tables/*/NUMBERS.tsv", "text_statistics/*/NUMBERS.tsv"]

_SCI = re.compile(
    r"^(?P<m>[-+]?\d*\.?\d+)\s*(?:[×x\*]\s*10\s*\^?\s*(?P<e1>[-+]?\d+)|[eE](?P<e2>[-+]?\d+))$"
)


_SUPER = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789-+")


def _clean(s):
    s = s.strip().translate(_SUPER).replace("−", "-").replace(" ", "").replace(" ", "").replace(" ", " ")
    s = s.replace("–", "-")
    return s


#: English number words read as integers ("nineteen" agrees with 19, "10" with "ten").
_WORDS = {w: i for i, w in enumerate("zero one two three four five six seven eight nine ten eleven twelve thirteen "
                                     "fourteen fifteen sixteen seventeen eighteen nineteen twenty".split())}


def parse(value):
    """Parse a printed value.

    Parameters
    ----------
    value : str
        Value as printed, e.g. ``"8,420"``, ``"0.318%"``, ``"2.5×"``, ``"3.2 × 10^-5"``,
        ``"P < 0.001"``.

    Returns
    -------
    tuple
        ``(kind, number, precision)`` where kind is ``"num"``, ``"lt"``, ``"le"``,
        ``"gt"``, ``"ge"``, ``"str"`` or ``"invalid"``. Non-finite numbers and
        numbers outside the floating-point range are invalid. Precision is
        ``("dp", n)`` decimals or ``("sf", n)`` significant figures.
    """
    s = _clean(value)
    if s.lower() == "two thousand":
        return "num", 2000.0, ("dp", 0)
    if s.lower() in _WORDS:
        return "num", float(_WORDS[s.lower()]), ("dp", 0)
    s = re.sub(r"^[A-Za-z][A-Za-z_ ]*\s*=\s*", "", s)  # drop a leading "P = ", "rho = "
    kind = "num"
    w = re.match(r"^(more|greater|fewer|less)(?:than|\s+than)\s*(.+)$", s, re.I)
    if w:
        s = (">" if w.group(1).lower() in ("more", "greater") else "<") + w.group(2)
    m = re.match(r"^(?:[A-Za-z][A-Za-z_ ]*)?\s*(<=|>=|[<>≤≥])\s*(.+)$", s)
    if m:
        kind = {"<": "lt", "≤": "le", "<=": "le", ">": "gt", "≥": "ge", ">=": "ge"}[m.group(1)]
        s = m.group(2).strip()
    s = re.sub(r"(%|×|-fold|fold|x)$", "", s).strip().replace(",", "")
    sci = _SCI.match(s.replace(" ", ""))
    try:
        if sci:
            mant = sci.group("m")
            exp = int(sci.group("e1") or sci.group("e2"))
            digits = len(mant.replace("-", "").replace("+", "").replace(".", "").lstrip("0")) or 1
            d = Decimal(f"{mant}e{exp}")
            precision = ("sf", digits)
        else:
            d = Decimal(s)
            precision = None
    except (InvalidOperation, ValueError):
        return "str", value.strip(), None
    if not d.is_finite():
        return "invalid", value.strip(), None
    number = float(d)
    if not math.isfinite(number) or (number == 0 and d != 0):
        return "invalid", value.strip(), None
    if precision is None:
        precision = ("dp", max(0, -d.as_tuple().exponent))
    return kind, number, precision


def rounded_integer_sf(printed):
    """Significant figures of a printed integer with trailing zeros (``"250,000"`` -> 2), else None."""
    kind, p, prec = parse(printed)
    if kind != "num" or prec != ("dp", 0) or p == 0:
        return None
    digits = str(int(abs(p)))
    sig = len(digits.rstrip("0"))
    return sig if sig < len(digits) else None


def agrees_rounded(printed, recomputed):
    """True if ``recomputed`` rounds to a trailing-zero integer ``printed`` at its significant figures."""
    n = rounded_integer_sf(printed)
    rk, r, _ = parse(recomputed)
    if n is None or rk != "num" or r == 0:
        return False
    _, p, _ = parse(printed)
    rr = round(r, -int(math.floor(math.log10(abs(r)))) + (n - 1))
    return math.isclose(rr, p, rel_tol=1e-12)


def agrees(printed, recomputed):
    """Return True if ``recomputed`` agrees with ``printed`` at the printed precision.

    Parameters
    ----------
    printed : str
        Value as printed in the manuscript (or the v28f replacement value).
    recomputed : str
        Full-precision value computed by the repository code.
    """
    kind, p, prec = parse(printed)
    if kind == "invalid":
        return False
    if kind == "str":
        return _clean(printed) == _clean(recomputed)
    rk, r, _ = parse(recomputed)
    if rk != "num":
        return False
    if kind == "lt":
        return r < p
    if kind == "le":
        return r <= p
    if kind == "gt":
        return r > p
    if kind == "ge":
        return r >= p
    mode, n = prec
    if mode == "dp":
        q = Decimal(1).scaleb(-n)
        return Decimal(repr(r)).quantize(q, ROUND_HALF_UP) == Decimal(repr(p)).quantize(q, ROUND_HALF_UP)
    if r == 0 or p == 0:
        return r == p
    rr = round(r, -int(math.floor(math.log10(abs(r)))) + (n - 1))
    return math.isclose(rr, p, rel_tol=1e-12)


def load_text(path):
    """Read a manuscript text and detect its indexing.

    Returns
    -------
    tuple
        ``("v42", {paragraph_index: text})`` for a paragraph-indexed text (every non-empty
        line ``<int>\t<text>``), else ``("v28d", {line_number: text})`` with physical,
        1-based line numbers of ``v28d_accepted.txt``.
    """
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    body = [ln for ln in lines if ln.strip()]
    if body and all(re.match(r"^\d+\t", ln) for ln in body):
        return "v42", {int(ln.split("\t", 1)[0]): ln.split("\t", 1)[1] for ln in body}
    return "v28d", dict(enumerate(lines, 1))


def agrees_any(printed, recomputed):
    """Exact rule or, for trailing-zero integers, the rounded rule."""
    return agrees(printed, recomputed) or agrees_rounded(printed, recomputed)


def check_ref(row, where, text, strict):
    """Check the v42 columns of one row; return (failures, open_rows)."""
    fails, opens = [], []
    row = dict(row, recomputed_value=row["v42_recomputed_value"])
    st = (row.get("v42_status") or "").strip()
    if st not in REF_STATUSES:
        return [f"{where}: unknown v42_status {st!r}"], []
    if st.startswith("REMOVED_"):
        # deliberately cut from the text; v42_quote holds the reason (edit id)
        print(f"{st} {where}: {row['quote'][:60]!r}; reason: {row['v42_quote']}")
        return fails, opens
    if st == "NOT_IN_V42":
        opens.append(f"{where}: NOT_IN_V42 {row['quote'][:60]!r}")
        if strict:
            fails.append(f"{where}: NOT_IN_V42 (strict)")
        return fails, opens
    if text is not None:
        try:
            para = text.get(int(row["v42_para"]), "")
        except ValueError:
            return [f"{where}: bad v42_para {row['v42_para']!r}"], []
        if _clean(row["v42_quote"]) not in _clean(para):
            fails.append(f"{where}: v42_quote not found in paragraph {row['v42_para']}")
        if _clean(row["v42_printed_value"]) not in _clean(row["v42_quote"]):
            fails.append(f"{where}: v42_printed_value not inside v42_quote")
    if st == "OPEN":
        opens.append(f"{where}: v42 OPEN printed {row['v42_printed_value']!r} recomputed {row['recomputed_value']!r}")
        if strict:
            fails.append(f"{where}: v42 OPEN (strict)")
    elif agrees_rounded(row["v42_printed_value"], row["recomputed_value"]) and not agrees(
        row["v42_printed_value"], row["recomputed_value"]
    ):
        print(f"ROUNDED {where}: v42 printed {row['v42_printed_value']!r} recomputed {row['recomputed_value']!r}")
    elif not agrees(row["v42_printed_value"], row["recomputed_value"]):
        fails.append(f"{where}: v42 MISMATCH printed {row['v42_printed_value']!r} recomputed {row['v42_recomputed_value']!r}")
    return fails, opens


def line_numbers(spec):
    """Expand ``"123"``, ``"L123"`` or ``"123-125"`` to a list of ints."""
    s = spec.strip().lstrip("Ll")
    if "-" in s:
        a, b = s.split("-", 1)
        return list(range(int(a), int(b.lstrip("Ll")) + 1))
    return [int(s)]


def check_file(path, text=None, strict=False, kind="v28d"):
    """Check one NUMBERS.tsv.

    ``kind`` is ``"v28d"`` or ``"v42"`` (see :func:`load_text`); with ``"v42"`` the
    ``text`` is paragraph-indexed and only the v42 columns are quote-checked.

    Returns
    -------
    tuple of (list, list)
        ``(failures, open_rows)``, each a list of human-readable strings.
    """
    failures, open_rows = [], []
    relp = path.relative_to(REPO).as_posix()
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames not in (COLUMNS, COLUMNS + HISTORICAL_REF_COLUMNS + REF_COLUMNS):
            return [f"{relp}: header must contain the base columns, retained v36 anchors and current v42 fields; found {reader.fieldnames}"], []
        has_ref = reader.fieldnames == COLUMNS + HISTORICAL_REF_COLUMNS + REF_COLUMNS
        if kind == "v42" and not has_ref:
            open_rows.append(f"{relp}: no v42 columns; not checked against the v42 text")
        for i, row in enumerate(reader, 2):
            where = f"{relp}:{i} (v28d line {row['v28d_line']})"
            status = (row["status"] or "").strip()
            if has_ref:
                fref, oref = check_ref(row, where + f" [v42 para {row['v42_para']}]", text if kind == "v42" else None, strict)
                failures += fref
                open_rows += oref
            if kind == "v42":
                pass
            elif text is not None and status != "OPEN":
                try:
                    lns = line_numbers(row["v28d_line"])
                    joined = " ".join(text.get(n, "") for n in lns)
                    if _clean(row["quote"]) not in _clean(joined):
                        failures.append(f"{where}: quote not found on that line")
                    if _clean(row["printed_value"]) not in _clean(row["quote"]):
                        failures.append(f"{where}: printed_value not inside quote")
                except ValueError:
                    failures.append(f"{where}: bad v28d_line {row['v28d_line']!r}")
            if status == "OPEN":
                open_rows.append(f"{where}: OPEN {row['quote'][:70]!r}")
                if strict:
                    failures.append(f"{where}: OPEN (strict)")
            elif status == "MATCH":
                if agrees_rounded(row["printed_value"], row["recomputed_value"]) and not agrees(
                    row["printed_value"], row["recomputed_value"]
                ):
                    print(f"ROUNDED {where}: printed {row['printed_value']!r} recomputed {row['recomputed_value']!r}")
                elif not agrees(row["printed_value"], row["recomputed_value"]):
                    failures.append(
                        f"{where}: MISMATCH printed {row['printed_value']!r} recomputed {row['recomputed_value']!r}"
                    )
            elif status.startswith("v28f_change:"):
                new = status.split(":", 1)[1].strip()
                if not new:
                    failures.append(f"{where}: v28f_change without a new value")
                elif agrees_rounded(new, row["recomputed_value"]) and not agrees(new, row["recomputed_value"]):
                    print(f"ROUNDED {where}: v28f value {new!r} recomputed {row['recomputed_value']!r}")
                elif not agrees(new, row["recomputed_value"]):
                    failures.append(f"{where}: v28f value {new!r} disagrees with recomputed {row['recomputed_value']!r}")
            elif status.split(":", 1)[0] == "SUPERSEDED_V33":
                if not has_ref:
                    failures.append(f"{where}: SUPERSEDED_V33 without v42 columns")
                else:
                    print(f"SUPERSEDED_V33 {where}: v28d printed {row['printed_value']!r}, v42 printed {row['v42_printed_value']!r}, recomputed {row['recomputed_value']!r}")
            else:
                failures.append(f"{where}: unknown status {status!r}")
    return failures, open_rows


def self_test():
    """Exercise the agreement rule on fixed cases; return 0 if all behave as expected."""
    cases = [
        ("8,420", "8420", True), ("0.318%", "0.31791", True), ("0.318%", "0.3186", False),
        ("2.5×", "2.46", True), ("3.2 × 10^-5", "3.24e-05", True), ("3.2 × 10−5", "3.3e-5", False),
        ("P < 0.001", "0.0004", True), ("P < 0.001", "0.002", False), ("P = 0.006", "0.0058", True),
        ("20.0%", "19.96", True), ("SIHUMIx", "SIHUMIx", True),
        ("6.1 × 10⁻⁴", "0.0006103515625", True), ("6.1 × 10⁻⁴", "0.00063", False),
        ("more than 10,000", "11428.6", True), ("more than 10,000", "9999", False),
        ("19", "nineteen", True), ("Nineteen", "19", True), ("10", "ten", True), ("10", "eleven", False),
    ]
    for printed, rec, want in [("250,000", "250197", True), ("29,700", "29663", True),
                               ("29,700", "29640", False), ("8,421", "8420", False)]:
        got = agrees_rounded(printed, rec)
        ok_r = got == want
        print(f"self-test rounded {printed!r:10} vs {rec!r:8} -> {got} {'ok' if ok_r else 'WRONG'}")
        if not ok_r:
            print("SELF-TEST FAIL")
            return 1
    ok = True
    for printed, rec, want in cases:
        got = agrees(printed, rec)
        ok &= got == want
        print(f"self-test {printed!r:16} vs {rec!r:10} -> {got} {'ok' if got == want else 'WRONG'}")
    print("SELF-TEST PASS" if ok else "SELF-TEST FAIL")
    return 0 if ok else 1


def main(argv=None):
    """Command-line entry point."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--strict", action="store_true", help="fail on OPEN rows")
    ap.add_argument("--unit", help="check only this unit directory name, e.g. Fig2")
    ap.add_argument("--manuscript-text", default=os.environ.get("FASTALAKE_MANUSCRIPT_TXT"))
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args(argv)
    if a.self_test:
        return self_test()
    kind, text = load_text(a.manuscript_text) if a.manuscript_text else ("v28d", None)
    if a.manuscript_text and kind == "v42":
        identity = json.loads((REPO / "docs/MANUSCRIPT.json").read_text())
        actual = hashlib.sha256(Path(a.manuscript_text).read_bytes()).hexdigest()
        if actual != identity["text_sha256"]:
            print("FAIL manuscript text differs from the pinned v42 accepted-text view")
            return 1
    files = sorted(p for g in GLOBS for p in REPO.glob(g) if not a.unit or p.parent.name == a.unit)
    if not files:
        print("FAIL no NUMBERS.tsv files matched" + (f" unit {a.unit!r}" if a.unit else ""))
        return 1
    failures, open_rows, nrows = [], [], 0
    for f in files:
        fl, op = check_file(f, text, a.strict, kind)
        failures += fl
        open_rows += op
        with open(f, encoding="utf-8") as fh:
            nrows += max(0, sum(1 for _ in fh) - 1)
    for o in open_rows:
        print(o)
    for x in failures:
        print("FAIL", x)
    print(f"{len(files)} file(s), {nrows} row(s), {len(open_rows)} OPEN, {len(failures)} failure(s)"
          + (f"; quotes checked against the {kind} text" if text is not None else "; quotes not checked (no --manuscript-text)"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
