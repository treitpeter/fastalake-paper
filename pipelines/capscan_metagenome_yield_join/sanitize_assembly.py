#!/usr/bin/env python3
"""Prefix an assembly's headers with its sample, dropping damaged and duplicated records.

WHY THIS EXISTS
---------------
Two of the 335 CapScan assemblies are damaged at source (Dari's read-only inputs,
written Oct 2024). The damage is of two kinds and only the first was obvious:

  1. Malformed headers -- 15 per file, at the very end, where sequence text is
     spliced into header lines:
         >kLFGPYMSGNRTFAEGRVSRDIWYRypeQIAQRIR
         >Ne=EdYRQLKRHHDA959 # 1645 # 1YFPE*

  2. Large duplicated blocks in the MIDDLE of the file, far from the damaged tail.
     SAMN27256424 holds 114,644 records for 106,950 unique accessions; SAMN27256425
     holds 66,551 for 54,444. One accession (k141_207420_2) repeats hundreds of times
     with a byte-identical sequence.

Truncating at the first malformed line -- the obvious reading of "truncate" -- fixes
(1) and leaves (2) untouched, because the duplication is upstream of the damage. So
this drops malformed records AND collapses duplicate accessions, keeping the first
occurrence, and reports both counts separately so neither is silent.

Five accessions across the two files carry duplicates whose SEQUENCES DIFFER (3 in
...424, 2 in ...425). Keeping the first occurrence is a choice, not a deduction; it is
counted and printed under `dup_conflicting` so it appears in the log and can be revisited.

Clean assemblies pass through unchanged and report zero for every drop counter, so this
runs over all 335 without special-casing the two.

Originals are never modified. Output goes to a new directory.

Usage:  sanitize_assembly.py <in.faa> <out.faa> <COHORT> <SAMPLE>
Exit 0 on success; prints a one-line TSV summary to stdout.
"""

import hashlib
import re
import sys

VALID_HEADER = re.compile(r"^>k\d+_\d+_\d+ # ")


def main() -> int:
    src, dst, cohort, sample = sys.argv[1:5]

    seen: dict[str, str] = {}  # accession -> sha256 of its first sequence
    n_in = n_kept = 0
    n_malformed = n_dup_same = n_dup_conflicting = 0

    def flush(out, acc, header, seq):
        """Write one record if its accession has not been seen; return a drop reason."""
        nonlocal n_kept, n_dup_same, n_dup_conflicting
        body = "".join(seq)
        digest = hashlib.sha256(body.encode()).hexdigest()
        if acc in seen:
            if seen[acc] == digest:
                n_dup_same += 1
            else:
                n_dup_conflicting += 1
            return
        seen[acc] = digest
        n_kept += 1
        # Rewrite only the accession token. The Prodigal coordinate description is kept:
        # it is informative and costs nothing.
        out.write(f">{cohort}_{sample}_{header[1:]}\n")
        for i in range(0, len(body), 60):
            out.write(body[i : i + 60] + "\n")

    with open(src) as fh, open(dst, "w") as out:
        acc = header = None
        seq: list[str] = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if acc is not None:
                    flush(out, acc, header, seq)
                n_in += 1
                if not VALID_HEADER.match(line):
                    n_malformed += 1
                    acc = header = None
                    seq = []
                    continue
                header = line
                acc = line.split()[0]
                seq = []
            elif acc is not None:
                seq.append(line.strip())
        if acc is not None:
            flush(out, acc, header, seq)

    dropped = n_malformed + n_dup_same + n_dup_conflicting
    if n_in - dropped != n_kept:
        print(f"FATAL {sample}: {n_in} in - {dropped} dropped != {n_kept} kept", file=sys.stderr)
        return 1
    if n_kept == 0:
        print(f"FATAL {sample}: no records survived", file=sys.stderr)
        return 1

    print(f"{sample}\t{n_in}\t{n_kept}\t{n_malformed}\t{n_dup_same}\t{n_dup_conflicting}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
