"""The canonical peptide definition. Import this; do not re-implement it.

WHY THIS FILE EXISTS
--------------------
On 2026-08-07 an audit of every script that reads SAGE's `peptide` column found three
different definitions in simultaneous use, and most scripts used none of them
deliberately:

    24 scripts read the peptide column
     9 collapsed I/L
     8 stripped modifications
    15 did NEITHER

Measured consequence, PREDICT, 6 acquisitions: the raw definition counts a median
**+6.5% more peptides** than the canonical one. Modifications account for +6.2% of that
and I/L for +0.6% -- so the term everyone worried about is the small one, and the term
nobody mentioned is the large one.

The comparisons themselves survive: the FastaLake/metagenome ratio shifts by a median
+0.6% (max +3.1% over 10 paired acquisitions) between definitions, because both arms
inflate together. So ARM RATIOS are safe under either rule; ABSOLUTE COUNTS QUOTED FROM
DIFFERENT SCRIPTS ARE NOT COMPARABLE, and that is what produced four different peptide
counts for one acquisition in the working documents.

THE RULE
--------
A peptide is its BARE AMINO-ACID SEQUENCE with I and L collapsed.

    modifications  removed   -- M[+15.9949]PEPTLDEK and MPEPTLDEK are the same molecule
    I -> L         applied   -- isobaric; no mass spectrometer separates them
    gate           peptide_q <= 0.01 on TARGETS ONLY (label == "1")

⚠️ `spectrum_q` is NOT `peptide_q`. Gating on the wrong one is the 2026-06-14 FDR-gate
defect. `spectrum_q` controls PSMs; a peptide-level count needs `peptide_q`.

WHEN SOMETHING ELSE IS CORRECT
------------------------------
Two other forms are legitimate for specific purposes and must be NAMED when used:

    raw()        the modified form, as SAGE prints it. Correct ONLY for search-engine
                 index cost, because that is what the index actually holds.
    lfq-joined   canonical(), then intersected with lfq.tsv. Correct ONLY for the
                 quantitative matrix. It is ~5x smaller than the identified set.
"""
import csv
import re

_MOD = re.compile(r"\[[^\]]*\]")

#: The FDR field a peptide-level count must gate on. Never `spectrum_q`.
Q_FIELD = "peptide_q"
Q_THRESHOLD = 0.01


def canonical(peptide: str) -> str:
    """Bare sequence, I/L collapsed. THE definition. Use this everywhere."""
    return _MOD.sub("", peptide).replace("I", "L")


def raw(peptide: str) -> str:
    """The modified form as SAGE prints it. Search-engine index cost only."""
    return peptide


def bare(peptide: str) -> str:
    """Modifications stripped, I and L KEPT APART.

    Correct for two purposes, and it must be named when used:

    1. separating the two effects in an audit, which is why it was written; and
    2. **counting AMINO ACIDS rather than peptides.** A C-terminal residue census asks
       which residue a protease left behind, over a list that contains both I and L. The
       identity rule is wrong there: applying `canonical()` to a residue census makes
       isoleucine structurally 0.00% and relabels I+L as leucine, and the run still
       produces a full table and a verdict -- so it reads as a result rather than an
       error. A database search assigns the residue and DOES know which one it is.

    It is NOT a peptide-identity definition. Never report a peptide COUNT from it.
    """
    return _MOD.sub("", peptide)


def collapse_il(profile: dict) -> dict:
    """Fold I into L for a METHOD-AGREEMENT test only, and label the result honestly.

    Returns a profile over 17 keys where I and L are replaced by a single `I/L`. Use it
    when a residue profile is compared against de novo sequencing, which cannot support an
    I-vs-L distinction: scoring two methods on a distinction one of them cannot make
    manufactures disagreement. Report the 18-residue profile; TEST on this one, and apply
    it to BOTH sides -- collapsing one side alone is the defect this exists to prevent.
    """
    out = {r: v for r, v in profile.items() if r not in ("I", "L")}
    out["I/L"] = profile.get("I", 0.0) + profile.get("L", 0.0)
    return out


def iter_target_psms(path, q_field=Q_FIELD, q=Q_THRESHOLD):
    """Yield (canonical_peptide, row, header_index) for every PSM passing the gate.

    Centralises the gate as well as the definition, because the two have gone wrong
    independently: one script gated on `spectrum_q` while collapsing I/L correctly, and
    several gated correctly while counting modified forms as distinct peptides.
    """
    csv.field_size_limit(1 << 30)
    with open(path, newline="") as fh:
        r = csv.reader(fh, delimiter="\t")
        h = next(r)
        idx = {name: i for i, name in enumerate(h)}
        if q_field not in idx:
            raise KeyError(f"{path} has no {q_field!r} column; columns are {h[:12]}...")
        ip, il, iq = idx["peptide"], idx["label"], idx[q_field]
        for row in r:
            if row[il] != "1":
                continue
            try:
                if float(row[iq]) > q:
                    continue
            except (ValueError, IndexError):
                continue
            yield canonical(row[ip]), row, idx


def peptide_set(path, **kw):
    """The distinct canonical peptides identified in one search."""
    return {p for p, _, _ in iter_target_psms(path, **kw)}


def peptide_protein_graph(path, accumulate_into=None, **kw):
    """peptide -> set(protein accessions), under the canonical definition.

    ⚠️ ACCUMULATES. Pass the same dict across samples; never assign a fresh protein set
    per sample. Each sample searched a different database, so one peptide legitimately
    carries a different protein set in each sample it was identified in. Assignment kept
    only the last sample's view, dropped 21.3% of the pooled graph's edges, and inflated
    the cohort group count from 259,298 to 302,273 before it was caught on 2026-08-07.
    """
    g = {} if accumulate_into is None else accumulate_into
    for pep, row, idx in iter_target_psms(path, **kw):
        prots = [x for x in row[idx["proteins"]].split(";") if x]
        if prots:
            g.setdefault(pep, set()).update(prots)
    return g
