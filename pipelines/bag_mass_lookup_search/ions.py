"""Explicit fragment definitions; no universal confidence threshold."""
import re
import numpy as np
from pyteomics import mass

PROTON = 1.007276466621
WATER = mass.calculate_mass(formula="H2O")


def peptide_ions(spelling):
    sequence, masses, nterm = "", [], 0.0
    position = 0
    while position < len(spelling):
        char = spelling[position]
        if char == "[":
            end = spelling.index("]", position)
            delta = float(spelling[position + 1:end])
            if masses:
                masses[-1] += delta
            else:
                nterm += delta
            position = end + 1
        elif char == "-":
            position += 1
        else:
            if char not in "ACDEFGHIKLMNPQRSTVWY":
                raise ValueError(f"Unsupported residue: {char}")
            sequence += char
            masses.append(mass.std_aa_mass[char])
            position += 1
    assert len(sequence) > 1
    cumulative = np.cumsum(masses)
    neutral = float(cumulative[-1] + WATER + nterm)
    ions = {}
    for n in range(1, len(sequence)):
        b = float(cumulative[n - 1] + nterm)
        y = float(cumulative[-1] - cumulative[len(sequence) - n - 1] + WATER)
        for z in (1, 2):
            ions[("b", n, z)] = (b + z * PROTON) / z
            ions[("y", n, z)] = (y + z * PROTON) / z
    return sequence, neutral, ions


def matched(values, observed, ppm=20.0):
    observed = np.sort(np.asarray(observed, dtype=float))
    values = np.asarray(values, dtype=float)
    tolerance = values * ppm / 1e6
    return np.searchsorted(observed, values + tolerance, side="right") > np.searchsorted(observed, values - tolerance, side="left")


def legacy_ics(ions, observed, length):
    keys = [key for key in ions if key[2] == 1]
    hits = {key for key, good in zip(keys, matched([ions[k] for k in keys], observed)) if good}
    return len(hits) / (2 * (length - 1)), hits


def backbone(keys, length):
    return {n if series == "b" else length - n for series, n, z in keys
            if series in ("b", "y") and 0 < n < length}


def relation(left, right, other_label):
    if other_label != 1:
        return "decoy_comparator", set()
    if left == right:
        return "same_sequence", set()
    a, b = left.replace("I", "L"), right.replace("I", "L")
    if a == b:
        return "IL_equivalent", set()
    if len(a) == len(b):
        changed = {i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y}
        if sorted(a) == sorted(b):
            return "composition_preserving_rearrangement", changed
        if len(changed) == 1:
            return "single_substitution_candidate", changed
        if len(changed) == 2:
            return "two_substitution_candidate", changed
    return "other_sequence_no_counterpart_claim", set()


def site_metrics(keys, ions, other_ions, length, positions):
    spans = {key for key in keys if key in ions and any(
        (key[0] == "b" and p <= key[1]) or (key[0] == "y" and p > length - key[1]) for p in positions)}
    unique = {key for key in spans if not matched([ions[key]], list(other_ions.values()))[0]}
    return len({(s, n) for s, n, z in spans}), len({(s, n) for s, n, z in unique})


def self_test():
    seq, neutral, ions = peptide_ions("AC[+57.0215]DIK")
    assert abs(neutral - (mass.calculate_mass(sequence=seq) + 57.0215)) < 1e-7
    score, hits = legacy_ics(ions, list(ions.values()), len(seq))
    assert score == 1 and len(backbone(hits, len(seq))) == len(seq) - 1
    assert legacy_ics(ions, [], len(seq))[0] == 0
    assert backbone({("b", 2, 1), ("b", 2, 2), ("y", 3, 1)}, 5) == {2}
    for first, second in (("AILK", "ALLK"), ("AN[+0.984016]AK", "ADAK")):
        a, _, ai = peptide_ions(first)
        b, _, bi = peptide_ions(second)
        assert all(matched(list(ai.values()), list(bi.values())))
        assert site_metrics(set(ai), ai, bi, len(a), set(range(1, len(a)+1)))[1] == 0
    a, _, ai = peptide_ions("ACDEK")
    b, _, bi = peptide_ions("ADCEK")
    label, positions = relation(a, b, 1)
    assert label == "composition_preserving_rearrangement"
    common = {key for key in ai if matched([ai[key]], list(bi.values()))[0]}
    assert site_metrics(common, ai, bi, len(a), positions)[1] == 0
    assert site_metrics(set(ai), ai, bi, len(a), positions)[1] > 0
    assert relation(a, b, -1) == ("decoy_comparator", set())


if __name__ == "__main__":
    self_test()
    print("PASS: ion masses, modifications, legacy coverage, charge collapse, I/L, deamidation, site discrimination and decoy labelling")
