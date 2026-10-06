#!/usr/bin/env python3
"""Shared tryptic digestion under the EXACT Sage 0.14.6 settings used in the searches.

Read from results.json of the final runs, not from documentation:

  enzyme          cleave_at KR, restrict P, missed_cleavages 2, min_len 5, max_len 50
  peptide mass    500.0 to 5000.0 Da
  static mods     C +57.0215
  variable mods   M +15.9949, max 2
  decoys          generated internally, decoy_tag rev_

Variable modifications are NOT expanded into distinct theoretical sequences here.
Sage expands them at scoring time; the unique theoretical *sequence* repertoire is the
quantity the peptide-space audit needs, and expanding M-oxidation would multiply the
same sequence rather than add new ones. Where a modification-expanded count is wanted
it is reported separately as a multiplier.

Serves Round 3D sections 1, 5, 8 and Round 3A section 4.
"""
import os
import sys

MONO = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "C": 103.00919 + 57.0215, "L": 113.08406, "I": 113.08406,
    "N": 114.04293, "D": 115.02694, "Q": 128.05858, "K": 128.09496,
    "E": 129.04259, "M": 131.04049, "H": 137.05891, "F": 147.06841,
    "R": 156.10111, "Y": 163.06333, "W": 186.07931,
}
H2O = 18.010565
MIN_LEN, MAX_LEN = 5, 50
MIN_MASS, MAX_MASS = 500.0, 5000.0
MISSED = 2


def mass(pep):
    m = H2O
    for a in pep:
        v = MONO.get(a)
        if v is None:
            return None
        m += v
    return m


def cleavage_sites(seq):
    """Trypsin: after K or R, not before P."""
    sites = [0]
    n = len(seq)
    for i, a in enumerate(seq):
        if a in "KR" and (i + 1 >= n or seq[i + 1] != "P"):
            sites.append(i + 1)
    if sites[-1] != n:
        sites.append(n)
    return sites


def digest(seq, missed=MISSED):
    """Yield peptides Sage would generate for this protein."""
    seq = seq.upper().replace("*", "")
    sites = cleavage_sites(seq)
    out = []
    for i in range(len(sites) - 1):
        for j in range(i + 1, min(i + 2 + missed, len(sites))):
            pep = seq[sites[i]:sites[j]]
            if MIN_LEN <= len(pep) <= MAX_LEN:
                m = mass(pep)
                if m is not None and MIN_MASS <= m <= MAX_MASS:
                    out.append(pep)
    return out


def read_fasta(path):
    acc, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line[0] == ">":
                if acc:
                    yield acc, "".join(buf)
                acc = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if acc:
        yield acc, "".join(buf)


def il(s):
    return s.replace("I", "L")


def digest_fasta(path, collapse_il=False):
    """Return (n_proteins, n_residues, total_peptides, pep2count) for one FASTA."""
    n_prot = n_res = n_pep = 0
    pep2n = {}
    for acc, seq in read_fasta(path):
        n_prot += 1
        n_res += len(seq)
        peps = digest(seq)
        n_pep += len(peps)
        for p in set(peps):
            k = il(p) if collapse_il else p
            pep2n[k] = pep2n.get(k, 0) + 1
    return n_prot, n_res, n_pep, pep2n


if __name__ == "__main__":
    p = sys.argv[1]
    a, b, c, d = digest_fasta(p, collapse_il=("--il" in sys.argv))
    print("%s\n  proteins %d  residues %d  peptide instances %d  unique sequences %d"
          % (os.path.basename(p), a, b, c, len(d)))
