"""Fragment-ion helpers for the spectrum evidence: parse() reads a peptide with modifications, match() matches theoretical b and y ions (charge 1 and 2) to observed peaks within a ppm tolerance. MASS and PROTON are the monoisotopic residue masses and proton mass."""
import re
import numpy as np
from numba import njit
MASS = np.array(
    [
        71.03711,
        0,
        103.00919,
        115.02694,
        129.04259,
        147.0684,
        57.02146,
        137.05891,
        113.08406,
        0,
        128.09496,
        113.08406,
        131.0405,
        114.04293,
        237.14774,
        97.05276,
        128.05858,
        156.1011,
        87.03203,
        101.04768,
        150.95363,
        99.06841,
        186.07932,
        0,
        163.06332,
        0,
    ],
    dtype=np.float32,
)

PROTON = np.float32(1.0072764)

TOK = re.compile(r"([A-Z])(?:\[([+-]?[0-9.]+)\])?")

def parse(s):
    t = TOK.findall(s)
    assert "".join(a + ("[" + m + "]" if m else "") for a, m in t) == s, s
    return "".join(a for a, m in t), np.array(
        [np.float32(MASS[ord(a) - 65] + np.float32(m or 0)) for a, m in t], dtype=np.float32
    )

def canon(s):
    return re.sub(r"\[[^\]]*\]", "", s).replace("I", "L")

@njit(cache=True)
def match(masses, ints, residues, total, charge):
    n = len(residues)
    b = np.float32(0)
    y = np.float32(total)
    bmask = np.uint64(0)
    ymask = np.uint64(0)
    matched = 0
    eligible = 0
    for i in range(n - 1):
        b = np.float32(b + residues[i])
        y = np.float32(y - residues[i])
        for series in range(2):
            neutral = b if series == 0 else y
            for z in range(1, max(2, min(charge, 3))):
                center = np.float32(neutral / np.float32(z))
                delta = np.float32(np.float32(center * np.float32(20)) / np.float32(1000000))
                lo = np.float32(center - delta)
                hi = np.float32(center + delta)
                if center + PROTON >= 200 and center + PROTON <= 2000:
                    eligible += 1
                index = np.searchsorted(masses, lo)
                found = False
                while index < len(masses) and masses[index] <= hi:
                    if ints[index] >= 0:
                        found = True
                    index += 1
                if found:
                    matched += 1
                    if series == 0:
                        bmask |= np.uint64(1) << np.uint64(i)
                    else:
                        ymask |= np.uint64(1) << np.uint64(i)
    return bmask, ymask, matched, eligible
