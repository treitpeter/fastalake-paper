"""Post hoc bounded-bag sweep on one finished acquisition of the native study.

Peter, 25 September 2026: "bound the expansion to peptides with small bags and few
substitutions ... we should test". The study's reservoir scan recorded, per protein and
original prediction, which route matched (exact 1, bag 2, mass 4). Bag size is a property of
the prediction, so the bag route can be masked per prediction after the fact, selection rerun
with the study's fixed-count rule, and the matched Sage searches repeated. No rescan.

The substitution half cannot be bounded post hoc: the native swap table pairs every
near-isobaric residue/dipeptide (W/AD/DA/EG/GE, R/GV/VG, N/GG, Q/AG/GA, ...), so every
lookup prediction has 10-40 one-substitution and 100+ two-substitution variants, and a hit
does not record which variant matched. That needs route tags per substitution level in the
query compiler and a rescan. Here mass_all and comb_all are kept only as reproduction checks.
"""
import gzip
import json
import math
import os
import sys
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
# PR 62 head 4215622: its search() accepts an arm list; the study revision's does not. The
# unbounded arms are asserted byte-identical to the study's FASTAs, which guards selection drift.
REPO = STUDY / "fdr_revision/repo"
sys.path.insert(0, str(REPO))
from fasta_lake.bags import dump, select_databases, sha  # noqa: E402
from tools.run_bags import search  # noqa: E402

EXACT = 1
B2, B4, B8, BALL = 1 << 4, 1 << 5, 1 << 6, 1 << 7        # bag orderings (incl. original) <= 2, 4, 8, any
MALL = 1 << 10                                          # mass route, unbounded (reproduction only)
ARMS_B = {
    "exact": EXACT,
    "bags_le2": EXACT | B2, "bags_le4": EXACT | B4, "bags_le8": EXACT | B8, "bags_all": EXACT | BALL,
    "mass_all": EXACT | MALL, "comb_all": EXACT | BALL | MALL,
}
REPRODUCE = {"exact": "exact", "bags_all": "bags", "mass_all": "mass", "comb_all": "combined"}


def permutations(elements):
    total = 1
    for element in elements:
        if element["unordered"]:
            letters = element["letters"]
            p = math.factorial(len(letters))
            for count in Counter(letters).values():
                p //= math.factorial(count)
            total *= p
    return total


def metrics(evidence_path, queries_path):
    """Per query unit: lookup flag; per real prediction: original, bag orderings, bagged residues, mass patterns."""
    lookup, original = {}, {}
    with gzip.open(evidence_path, "rt") as handle:
        for line in handle:
            d = json.loads(line)
            lookup[d["query"]] = d["lookup"]        # every query unit, including the two scrambled controls
            if d["control"] == 0:
                original[d["query"]] = d["original"]
    bag_perms, bag_residues, mass_patterns = {}, {}, Counter()
    with open(queries_path) as handle:
        for line in handle:
            d = json.loads(line)
            if d["route"] == 2:
                q = d["query"]
                bag_perms[q] = max(bag_perms.get(q, 1), permutations(d["elements"]))   # several recorded forms possible
                bag_residues[q] = max(bag_residues.get(q, 0),
                                      sum(len(e["letters"]) for e in d["elements"] if e["unordered"]))
            elif d["route"] == 4:
                mass_patterns[d["query"]] += 1
    return lookup, original, bag_perms, bag_residues, mass_patterns


def remap(q, routes, bag_perms):
    new = routes & EXACT
    if routes & 2:
        p = bag_perms.get(q, 1)
        new |= BALL
        if p <= 2:
            new |= B2
        if p <= 4:
            new |= B4
        if p <= 8:
            new |= B8
    if routes & 4:
        new |= MALL
    return new


def summary(task, index, records, lookup, original, bag_perms, bag_residues, mass_patterns, repro, seconds):
    eligible = [q for q in original if lookup[q]]   # real predictions only; control units carry no bound metrics
    return dict(
        sample=task["sample"], cohort=task["cohort"], index=index, reservoir_records=records,
        lookup_eligible_predictions=len(eligible),
        with_recorded_bag=sum(1 for q in eligible if q in bag_perms),
        bag_orderings_hist=dict(sorted(Counter(min(bag_perms.get(q, 1), 64) for q in eligible).items())),
        bag_residues_hist=dict(sorted(Counter(bag_residues.get(q, 0) for q in eligible).items())),
        mass_patterns_hist=dict(sorted(Counter(min(mass_patterns[q], 400) // 20 * 20 for q in eligible).items())),
        arms=ARMS_B, bits=dict(EXACT=EXACT, B2=B2, B4=B4, B8=B8, BALL=BALL, MALL=MALL),
        definitions=dict(
            bags_leN="bag route kept only for predictions whose recorded bag(s) admit at most N orderings including the original",
            mass_all="mass route unbounded; reproduction check only. A substitution bound needs compiler route tags and a rescan",
        ),
        reproduction=repro, selection_seconds=seconds,
    )


def main(index):
    manifest = json.loads((STUDY / "attempt_02/manifest.json").read_text())
    assert sha(STUDY / "attempt_02/manifest.json") == json.loads((HERE / "SOURCE.json").read_text())["study_manifest_sha256"]
    task = manifest["tasks"][index]
    src = STUDY / "attempt_02/outputs" / task["cohort"] / task["sample"]
    assert (src / "SEARCHES_COMPLETE.json").exists(), "source acquisition not finished"
    local = dict(manifest, output=str(HERE / "outputs"))
    out = HERE / "outputs" / task["cohort"] / task["sample"]
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    compile_receipt = json.loads((src / "compile.json").read_text())
    study_selection = json.loads((src / "selection.json").read_text())

    lookup, original, bag_perms, bag_residues, mass_patterns = metrics(src / "evidence.jsonl.gz", src / "queries.jsonl")

    bounded = out / "bounded_candidates.jsonl.gz"
    records = 0
    with gzip.open(src / "candidates.jsonl.gz", "rt") as handle, gzip.open(bounded, "wt", compresslevel=1) as dest:
        for line in handle:
            r = json.loads(line)
            records += 1
            expected = 0
            for q, routes in r["evidence"]:
                if lookup[q]:
                    expected |= routes
            assert expected == r["lookup_routes"], "lookup_routes rule violated in source"
            evidence = [[q, remap(q, routes, bag_perms)] for q, routes in r["evidence"]]
            lr = 0
            for q, routes in evidence:
                if lookup[q]:
                    lr |= routes
            dest.write(json.dumps({"header": r["header"], "sequence": r["sequence"],
                                   "lookup_routes": lr, "evidence": evidence}) + "\n")

    selection = select_databases(bounded, out, compile_receipt["query_units"], arms=ARMS_B)
    repro = {}
    for arm, study_arm in REPRODUCE.items():
        mine = selection["arms"][arm]["fasta"]["sha256"]
        theirs = study_selection["arms"][study_arm]["fasta"]["sha256"]
        repro[arm] = dict(study_arm=study_arm, identical=(mine == theirs), mine=mine, study=theirs)
    assert all(v["identical"] for v in repro.values()), "unbounded arms do not reproduce the study FASTAs: " + json.dumps(repro)
    dump(out / "CANDIDATES_COMPLETE.json", {"status": "CANDIDATES_COMPLETE", "selection": selection,
                                            "comparison": "bounded-bag arms versus fresh exact selection, same candidate stream"})
    dump(out / "BOUNDS.json", summary(task, index, records, lookup, original, bag_perms, bag_residues,
                                      mass_patterns, repro, time.monotonic() - started))

    search(local, index, arms=ARMS_B, baseline_arm="exact")
    bounded.unlink()   # reproducible from the source stream; FASTAs, assignments and BOUNDS.json are kept
    dump(out / "BOUNDED_COMPLETE.json", dict(status="BOUNDED_SWEEP_COMPLETE", sample=task["sample"], index=index,
                                             seconds=time.monotonic() - started, job_id=os.environ.get("SLURM_JOB_ID"),
                                             workflow_fdr="UNVERIFIED"))


if __name__ == "__main__":
    main(int(sys.argv[1]))
