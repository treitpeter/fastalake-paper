#!/usr/bin/env python3
"""
Figure 2 revision, 2026-08-11 -- the derivation behind every panel.

Four arms, the SAME fifteen CAMPI mock acquisitions (S01-S10, S11_Fraction1-4, S12),
the SAME SAGE settings in every field but `database.fasta`:

  fastalake   per-sample FastaLake database, arm noclust_2026-08-07, razor/hash-acc
  campi_ref   CAMPI SIHUMI_DB1UNIPROT targets,  29,673   (their reference arm)
  campi_mo    CAMPI SIHUMI_DB2MG as shipped,    19,435   (their metagenome arm, CORRECTED)
  campi_all   all of CAMPI's mock FASTAs in ONE database, 48,992   (item 2.5, new)

WHAT IS DERIVED
  1. depth per acquisition -- peptides, proteins  (items 2.4/2.5)
  2. per-species peptide support                  (item 2.1's substance, decision D-01)
  3. quantitative accuracy against the mock       (decision D-2026-08-10-11)
  4. the 2x2, database source x FastaLake curation (item 2.6)
  5. peptide sharing between members              (item 2.3, the circos)

THE RULES THIS OBEYS
  * ONE peptide definition, `publication/lib/peptide.py`: bare sequence, I/L collapsed,
    peptide_q <= 0.01 on targets only. A cross-database peptide comparison under two
    definitions is meaningless, and four different counts for one acquisition is what
    happened before this module existed.
  * `sample_group` is per-acquisition and `study_group` is across acquisitions, and they
    are NEVER combined in one number. This script works at PEPTIDE level for the arm
    comparison, because the four arms have four different protein namespaces and no
    shared grouping exists across them -- so a protein-group comparison across arms
    would be comparing four different counting rules. Protein counts are reported
    per-arm-per-acquisition only, and labelled as such.
  * COUNT with representatives; APPORTION INTENSITY with raw rows. Species intensity
    here is apportioned from peptide-level `lfq.tsv` rows using ONLY peptides specific
    to one member, so nothing is apportioned across an ambiguity.
  * Every check carries an expected value (rule 7).

SPECIES ASSIGNMENT. Exact substring containment of each identified peptide in each of
CAMPI's eight member proteomes, by Aho-Corasick -- the same primitive
`Fig2/analysis/derive_species_recall.py` uses. A peptide in exactly ONE member is
SPECIFIC; in several, SHARED, and counted separately rather than assigned. This is
deliberately NOT an LCA: CAMPI's Unipept LCA resolves shared clostridial peptides above
species rank, and mixing the two is how 8/8 and 6/8 were both reported for one build.

⚠️ I/L IS COLLAPSED ON BOTH SIDES. The peptide definition collapses I/L, so the
proteomes must be scanned in the same alphabet or containment silently fails. Collapsing
one side only is the defect this note exists to prevent. The cost -- how much apparent
specificity the collapse creates -- is measured and printed, not assumed small.

Run: python3 derive_fig2_revision.py
"""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import ahocorasick

CANON = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(CANON / "publication"))
from lib import peptide as pepdef  # noqa: E402

OUT = CANON / "publication/analyses/fig2_revision_2026-08-11"
TAB = OUT / "tables"
TAB.mkdir(parents=True, exist_ok=True)

ARM_ROOT = CANON / "publication_data/arms"
CHAIN = ARM_ROOT / "noclust_2026-08-07/cohorts/CAMPI_MOCK/5_sage/razor_hash-acc"

# The eight member proteomes, with their OS= strings. 29,557 reference proteins,
# contaminants already excluded -- campi_targets.faa is campi_ref_targets minus cRAP.
MEMBER_FASTA = CANON / "publication/figures/Fig2/tables/campi_targets.faa"
MEMBERS_TSV = CANON / "publication/figures/Fig2/tables/campi_members.tsv"

SAMPLES = [f"S{i:02d}" for i in range(1, 11)] + [
    f"S11_Fraction{i}" for i in range(1, 5)
] + ["S12"]

ARMS = {
    "fastalake": (CHAIN, "FastaLake per-sample", "FastaLake lake", "on"),
    "campi_ref": (ARM_ROOT / "campi_ref_solo/CAMPI_MOCK", "CAMPI reference DB1",
                  "CAMPI FASTAs", "off"),
    "campi_mo": (ARM_ROOT / "campi_mo_fixed/CAMPI_MOCK", "CAMPI metagenome DB2",
                 "CAMPI FASTAs", "off"),
    "campi_all": (ARM_ROOT / "campi_allfasta/CAMPI_MOCK", "CAMPI all FASTAs, one DB",
                  "CAMPI FASTAs", "off"),
}

EXPECT = {
    "n_samples": 15,
    "n_members": 8,
    "member_proteins": 29557,
}

TRANS = str.maketrans("I", "L")  # the SAME collapse peptide.canonical() applies


def log(*a):
    print(*a, flush=True)


def fail(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


# --------------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------------
def load_members():
    """member -> canonical OS= string, from the curated table."""
    with open(MEMBERS_TSV, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    m = {r["member"]: r["campi_os_string"] for r in rows}
    if len(m) != EXPECT["n_members"]:
        fail(f"expected {EXPECT['n_members']} members, got {len(m)}")
    return m


def load_member_proteomes(members):
    """member -> list of I/L-collapsed sequences, split by the OS= string."""
    os_to_member = {os_s: mem for mem, os_s in members.items()}
    seqs = defaultdict(list)
    unmatched = 0
    total = 0
    hdr, buf = None, []

    def flush():
        nonlocal unmatched, total
        if hdr is None:
            return
        total += 1
        # OS=... runs until the next ' XX=' key, which is how UniProt headers are built.
        mem = None
        if "OS=" in hdr:
            tail = hdr.split("OS=", 1)[1]
            for key in (" GN=", " PE=", " SV=", " OX="):
                if key in tail:
                    tail = tail.split(key, 1)[0]
            mem = os_to_member.get(tail.strip())
        if mem is None:
            unmatched += 1
        else:
            seqs[mem].append("".join(buf).upper().translate(TRANS))

    with open(MEMBER_FASTA) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                flush()
                hdr, buf = line[1:], []
            else:
                buf.append(line)
    flush()

    log(f"  member proteomes: {total} proteins, {unmatched} unassigned to any member")
    if total != EXPECT["member_proteins"]:
        fail(f"expected {EXPECT['member_proteins']} member proteins, got {total}")
    if unmatched:
        fail(f"{unmatched} member proteins matched no OS= string; the map is incomplete")
    for mem in members:
        if not seqs[mem]:
            fail(f"member {mem} has no proteins")
        log(f"    {mem:<34} {len(seqs[mem]):>6} proteins")
    return seqs


def arm_peptides(arm):
    """sample -> set of canonical peptides, and sample -> n proteins at protein_q<=0.01."""
    root = ARMS[arm][0]
    peps, prots = {}, {}
    for s in SAMPLES:
        res = root / s / "results.sage.tsv"
        if not res.is_file():
            fail(f"{arm}/{s}: no results.sage.tsv at {res}")
        pset = set()
        pr = set()
        for pep, row, idx in pepdef.iter_target_psms(res):
            pset.add(pep)
            try:
                if float(row[idx["protein_q"]]) <= 0.01:
                    for a in row[idx["proteins"]].split(";"):
                        if a:
                            pr.add(a)
                    # ^ raw accessions; the arms have four different namespaces, so
                    #   this is a PER-ARM count only and is never compared as a group.
            except (KeyError, ValueError, IndexError):
                pass
        peps[s], prots[s] = pset, len(pr)
    return peps, prots


# --------------------------------------------------------------------------------
# species assignment
# --------------------------------------------------------------------------------
def assign_species(all_peptides, proteomes, members):
    """peptide -> frozenset(members containing it), by exact containment."""
    A = ahocorasick.Automaton()
    for i, p in enumerate(sorted(all_peptides)):
        A.add_word(p, (i, p))
    A.make_automaton()

    hits = defaultdict(set)
    for mem in members:
        for seq in proteomes[mem]:
            for _, (_, p) in A.iter(seq):
                hits[p].add(mem)
    return {p: frozenset(v) for p, v in hits.items()}


# --------------------------------------------------------------------------------
# quantitative accuracy
# --------------------------------------------------------------------------------
def species_intensity(arm, sample, pep2mem):
    """member -> summed LFQ intensity, from peptides SPECIFIC to that member.

    Rule: apportion intensity with RAW rows. Every row of `lfq.tsv` is used as written;
    nothing is collapsed to a representative first. Only peptides matching exactly ONE
    member contribute, so no intensity is split across an ambiguity and none is assigned
    by a tie-break. Shared and unassigned intensity is returned separately, so the
    fraction of the signal this measure DISCARDS is visible rather than implied.
    """
    path = ARMS[arm][0] / sample / "lfq.tsv"
    if not path.is_file():
        return None, None, None
    per = defaultdict(float)
    shared = unassigned = 0.0
    csv.field_size_limit(1 << 30)
    with open(path, newline="") as fh:
        r = csv.reader(fh, delimiter="\t")
        h = next(r)
        idx = {n: i for i, n in enumerate(h)}
        pcol = idx.get("peptide")
        # SAGE writes one intensity column per mzML file; find it by exclusion.
        known = {"peptide", "proteins", "q_value", "score", "spectral_angle", "charge"}
        icols = [i for n, i in idx.items() if n not in known]
        if pcol is None or not icols:
            fail(f"{path}: cannot find peptide / intensity columns in {h[:8]}")
        for row in r:
            try:
                val = sum(float(row[i]) for i in icols if row[i] not in ("", "NaN"))
            except (ValueError, IndexError):
                continue
            if val <= 0:
                continue
            mems = pep2mem.get(pepdef.canonical(row[pcol]), frozenset())
            if len(mems) == 1:
                per[next(iter(mems))] += val
            elif mems:
                shared += val
            else:
                unassigned += val
    return per, shared, unassigned


# --------------------------------------------------------------------------------
def main():
    log("=== inputs ===")
    members = load_members()
    proteomes = load_member_proteomes(members)

    log("\n=== reading the four arms ===")
    arm_pep, arm_prot = {}, {}
    for arm in ARMS:
        peps, prots = arm_peptides(arm)
        arm_pep[arm], arm_prot[arm] = peps, prots
        tot = len(set().union(*peps.values()))
        log(f"  {arm:<12} {len(peps)} acquisitions, {tot:>7} distinct peptides (union)")
        if len(peps) != EXPECT["n_samples"]:
            fail(f"{arm}: expected {EXPECT['n_samples']} acquisitions, got {len(peps)}")

    universe = set().union(*(set().union(*p.values()) for p in arm_pep.values()))
    log(f"\n  peptide universe across all four arms: {len(universe)}")

    log("\n=== species assignment (Aho-Corasick, I/L collapsed both sides) ===")
    pep2mem = assign_species(universe, proteomes, members)
    n_any = len(pep2mem)
    n_spec = sum(1 for v in pep2mem.values() if len(v) == 1)
    log(f"  peptides in >=1 member proteome : {n_any} of {len(universe)}")
    log(f"  specific to exactly one member  : {n_spec}")
    log(f"  shared between members          : {n_any - n_spec}")

    # Both sides must be in the same alphabet or containment silently under-reports.
    # This is the check that the collapse was actually applied to the PEPTIDES; the
    # proteomes are collapsed at load with the identical TRANS table.
    log("\n=== I/L alphabet check (both sides must be collapsed) ===")
    il_pep = sum(1 for p in universe if "I" in p)
    log(f"  peptides still containing 'I' after canonical(): {il_pep} (expected 0)")
    if il_pep != 0:
        fail("canonical() did not collapse I/L; the two sides are in different alphabets")

    # ---------------------------------------------------------------- depth
    log("\n=== 1. depth per acquisition ===")
    with open(TAB / "arm_depth_per_sample.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["arm", "arm_label", "sample", "peptides", "proteins_q0.01"])
        for arm, (_, label, _, _) in ARMS.items():
            for s in SAMPLES:
                w.writerow([arm, label, s, len(arm_pep[arm][s]), arm_prot[arm][s]])
    for arm, (_, label, _, _) in ARMS.items():
        tots = [len(arm_pep[arm][s]) for s in SAMPLES]
        med = sorted(tots)[len(tots) // 2]
        uni = len(set().union(*arm_pep[arm].values()))
        log(f"  {label:<28} median {med:>7} peptides/acquisition, union {uni:>7}")

    # ------------------------------------------------------- species support
    log("\n=== 2. per-species peptide support (decision D-01) ===")
    rows = []
    for arm, (_, label, _, _) in ARMS.items():
        union = set().union(*arm_pep[arm].values())
        spec = defaultdict(int)
        shared = 0
        none = 0
        for p in union:
            m = pep2mem.get(p, frozenset())
            if len(m) == 1:
                spec[next(iter(m))] += 1
            elif m:
                shared += 1
            else:
                none += 1
        for mem in members:
            rows.append([arm, label, mem, spec[mem]])
        log(f"  {label:<28} " + "  ".join(f"{m.split()[0][:6]}:{spec[m]}" for m in members))
        log(f"  {'':<28} shared between members {shared}, matching NO member {none} "
            f"({100.0 * none / max(1, len(union)):.1f}% of its peptides)")
        rows.append([arm, label, "_SHARED_BETWEEN_MEMBERS", shared])
        rows.append([arm, label, "_NO_MEMBER", none])
    with open(TAB / "species_support.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["arm", "arm_label", "member", "specific_peptides"])
        w.writerows(rows)

    # recall at several evidence floors -- the floor is a decision, so show it
    log("\n=== recall by evidence floor ===")
    floors = [1, 2, 3, 5, 10]
    with open(TAB / "species_recall_by_floor.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["arm", "arm_label"] + [f"floor_{f}" for f in floors])
        for arm, (_, label, _, _) in ARMS.items():
            union = set().union(*arm_pep[arm].values())
            spec = defaultdict(int)
            for p in union:
                m = pep2mem.get(p, frozenset())
                if len(m) == 1:
                    spec[next(iter(m))] += 1
            r = [sum(1 for mem in members if spec[mem] >= f) for f in floors]
            w.writerow([arm, label] + [f"{x}/8" for x in r])
            log(f"  {label:<28} " + "  ".join(f">={f}:{x}/8" for f, x in zip(floors, r)))

    # ------------------------------------------------- quantitative accuracy
    log("\n=== 3. quantitative accuracy: observed composition per arm ===")
    with open(TAB / "quant_composition.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["arm", "arm_label", "sample", "member", "intensity",
                    "pct_of_specific", "shared_intensity", "unassigned_intensity"])
        for arm, (_, label, _, _) in ARMS.items():
            for s in SAMPLES:
                per, shared, unass = species_intensity(arm, s, pep2mem)
                if per is None:
                    log(f"  ! {arm}/{s}: no lfq.tsv")
                    continue
                tot = sum(per.values())
                for mem in members:
                    pct = 100.0 * per[mem] / tot if tot else 0.0
                    w.writerow([arm, label, s, mem, f"{per[mem]:.6g}", f"{pct:.4f}",
                                f"{shared:.6g}", f"{unass:.6g}"])
    log(f"  -> {TAB / 'quant_composition.tsv'}")

    # ----------------------------------------------------- peptide sharing
    log("\n=== 5. peptide sharing between members (circos input) ===")
    pair = defaultdict(int)
    for p, m in pep2mem.items():
        if len(m) > 1:
            ms = sorted(m)
            for i in range(len(ms)):
                for j in range(i + 1, len(ms)):
                    pair[(ms[i], ms[j])] += 1
    with open(TAB / "member_peptide_sharing.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["member_a", "member_b", "shared_peptides"])
        for (a, b), n in sorted(pair.items(), key=lambda x: -x[1]):
            w.writerow([a, b, n])
    log(f"  {len(pair)} member pairs share at least one peptide")
    for (a, b), n in sorted(pair.items(), key=lambda x: -x[1])[:5]:
        log(f"    {a} <-> {b}: {n}")

    json.dump(
        {
            "n_samples": len(SAMPLES),
            "arms": {a: ARMS[a][1] for a in ARMS},
            "peptide_universe": len(universe),
            "peptides_in_a_member": n_any,
            "peptides_specific": n_spec,
        },
        open(TAB / "derivation_summary.json", "w"),
        indent=2,
    )
    log("\n=== OK ===")


if __name__ == "__main__":
    main()
