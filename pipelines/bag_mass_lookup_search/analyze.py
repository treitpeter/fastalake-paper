"""Recount completed searches and locate changed peptides in selected databases."""
from pathlib import Path
from collections import Counter, defaultdict
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import json
import math
import os
import re
import statistics
import sys

import ahocorasick

HERE = Path(__file__).resolve().parent
ARMS = ("exact", "bags", "mass", "combined")
CONTRASTS = (("exact", "bags"), ("exact", "mass"), ("exact", "combined"),
             ("bags", "combined"), ("mass", "combined"))
MOD = re.compile(r"\[[^\]]*\]")
SOURCES = {}


def canonical(value):
    return MOD.sub("", value).replace("-", "").replace("I", "L")


def identity(path):
    path = Path(path).resolve()
    if str(path) not in SOURCES:
        before = path.stat()
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(block)
        after = path.stat()
        assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
        SOURCES[str(path)] = dict(path=str(path), bytes=after.st_size, mtime_ns=after.st_mtime_ns,
                                  sha256=digest.hexdigest())
    return SOURCES[str(path)]


def read(path):
    identity(path)
    return json.loads(Path(path).read_text())


def dump(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n")
    temporary.replace(path)


def tsv(path, rows):
    assert rows
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def search(path):
    identity(path)
    minimum, spectra, accepted, reported_targets = {}, {}, defaultdict(set), set()
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            assert int(row["rank"]) == 1
            peptide = canonical(row["peptide"])
            q, label = float(row["peptide_q"]), int(row["label"])
            key = (row["filename"], row["scannr"])
            assert key not in spectra, "Multiple rank-one rows for a spectrum"
            spectra[key] = (peptide, label, q)
            if label == 1:
                reported_targets.add(peptide)
            if label == 1 and math.isfinite(q):
                minimum[peptide] = min(q, minimum.get(peptide, math.inf))
                if q <= 0.01:
                    accepted[peptide].add(key)
    return dict(minimum=minimum, spectra=spectra, accepted=dict(accepted), reported_targets=reported_targets)


def fasta(path):
    header, parts = None, []
    with Path(path).open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:].split()[0], []
            else:
                assert header is not None
                parts.append(line.upper())
    if header is not None:
        yield header, "".join(parts)


def locate(path, peptides):
    identity(path)
    auto = ahocorasick.Automaton()
    for peptide in sorted(peptides):
        auto.add_word(peptide, peptide)
    auto.make_automaton()
    proteins, contexts, memberships = {}, defaultdict(set), defaultdict(set)
    for accession, original in fasta(path):
        assert accession not in proteins
        proteins[accession] = hashlib.sha256(original.encode()).hexdigest()
        sequence = original.replace("I", "L")
        if not peptides:
            continue
        for end, peptide in auto.iter(sequence):
            start = end + 1 - len(peptide)
            contexts[peptide].add((sequence[start - 1] if start else "^",
                                   sequence[end + 1] if end + 1 < len(sequence) else "$"))
            memberships[peptide].add(accession)
    return proteins, contexts, memberships


def category(change, other_present, other_q, other_reported=False):
    if not other_present:
        assert other_q is None and not other_reported, "Target PSM sequence absent from search FASTA"
        return "sequence_added" if change == "gained" else "sequence_removed"
    if other_q is not None:
        assert other_q > 0.01
        return "sequence_retained_above_q_cutoff"
    return "sequence_retained_nonfinite_q" if other_reported else "sequence_retained_no_reported_target"


def spectrum_states(keys, peptide, other):
    counts = Counter()
    for key in keys:
        row = other.get(key)
        if row is None:
            counts["no_reported_row"] += 1
        elif row[1] != 1:
            counts["decoy"] += 1
        elif row[0] == peptide:
            assert not math.isfinite(row[2]) or row[2] > 0.01
            counts["same_peptide_above_cutoff"] += 1
        elif math.isfinite(row[2]) and row[2] <= 0.01:
            counts["different_accepted_target"] += 1
        else:
            counts["different_unaccepted_target"] += 1
    return counts


def run_case(entry):
    case = Path(entry["directory"])
    task, candidate, receipt = read(case / "task.json"), read(case / "CANDIDATES_COMPLETE.json"), read(case / "SEARCHES_COMPLETE.json")
    assert candidate["baseline_reproduced"]
    assert read(case / "baseline_check.json")["identical_accessions_and_sequences"]
    selection, assignments = read(case / "selection.json"), read(case / "assignments.json")
    expected = {row["arm"]: row for row in receipt["rows"]}
    data, accepted, config_reference = {}, {}, None
    for arm in ARMS:
        path = case / "searches" / arm / "results.sage.tsv"
        data[arm] = search(path)
        execution = receipt["executions"][arm]
        if "results" in execution:
            assert identity(path)["sha256"] == execution["results"]["sha256"]
        accepted[arm] = set(data[arm]["accepted"])
        assert len(accepted[arm]) == expected[arm]["peptides"]
        assert len(set().union(*data[arm]["accepted"].values())) == expected[arm]["spectra"]
        config = read(case / "searches" / arm / "config.json")
        config.pop("output_directory", None)
        config["database"].pop("fasta")
        if config_reference is None:
            config_reference = config
        assert config == config_reference, "Search parameters differ beyond FASTA/output"
    identity(case / "peptide_changes.tsv.gz")
    saved = defaultdict(set)
    with gzip.open(case / "peptide_changes.tsv.gz", "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["arm"] in ARMS:
                saved[(row["arm"], row["change"])].add(row["peptide_IL"])
    for arm in ARMS:
        assert accepted[arm] - accepted["exact"] == saved[(arm, "gained")]
        assert accepted["exact"] - accepted[arm] == saved[(arm, "lost")]
        assert len(saved[(arm, "gained")]) == expected[arm]["gained"]
        assert len(saved[(arm, "lost")]) == expected[arm]["lost"]
    changed = set().union(*(accepted[a] ^ accepted[b] for a, b in CONTRASTS))
    locations = {}
    for arm in ARMS:
        path = case / (arm + ".fasta")
        assert identity(path)["sha256"] == selection["arms"][arm]["fasta"]["sha256"]
        locations[arm] = locate(path, changed)
        assert len(locations[arm][0]) == selection["arms"][arm]["selected_proteins"]
    null_sets, null_rows = {}, []
    null_case = Path(entry["null_directory"])
    null_report = read(null_case / "NULL_COMPARISON.json") if (null_case / "NULL_COMPARISON.json").exists() else None
    if null_report:
        for row in null_report["rows"]:
            if row["q"] == 0.01 and row["level"] == "canonical_peptide":
                method, seed = row["method"], row["seed"]
                path = (case / "searches" / f"shuffle_{seed}" / "results.sage.tsv" if method == "combined" else null_case / "searches" / f"shuffle_{seed}_{method}" / "results.sage.tsv")
                if path.exists():
                    values = set(search(path)["accepted"])
                else:
                    assert row["accepted"] == 0
                    values = set()
                assert len(values) == row["accepted"]
                null_sets[(method, seed)] = values
                null_rows.append(dict(cohort=task["cohort"], sample=task["sample"], **row))
    summaries, peptide_rows, turnover, thresholds = [], [], [], []
    meta = dict(cohort=task["cohort"], sample=task["sample"])
    for lower, upper in CONTRASTS:
        contrast = upper + "_vs_" + lower
        counts = Counter()
        for change, values, own, other in (("gained", accepted[upper] - accepted[lower], upper, lower),
                                          ("lost", accepted[lower] - accepted[upper], lower, upper)):
            for peptide in sorted(values):
                own_contexts, other_contexts = locations[own][1][peptide], locations[other][1][peptide]
                assert own_contexts, f"Accepted peptide absent from its FASTA: {peptide}"
                other_q = data[other]["minimum"].get(peptide)
                group = category(change, bool(other_contexts), other_q, peptide in data[other]["reported_targets"])
                counts[(change, group)] += 1
                own_q = data[own]["minimum"][peptide]
                psm = spectrum_states(data[own]["accepted"][peptide], peptide, data[other]["spectra"])
                row = dict(**meta, contrast=contrast, change=change, peptide_IL=peptide,
                           category=group, length=len(peptide), own_min_q=own_q,
                           other_min_q=other_q if other_q is not None else "",
                           own_accepted_spectra=len(data[own]["accepted"][peptide]),
                           own_protein_count=len(locations[own][2][peptide]),
                           other_protein_count=len(locations[other][2][peptide]),
                           shared_flanking_context=bool(own_contexts & other_contexts),
                           own_contexts=json.dumps(sorted(own_contexts)), other_contexts=json.dumps(sorted(other_contexts)),
                           own_null_seeds=sum(peptide in null_sets.get((own, seed), set()) for seed in (1, 2)),
                           null_controls_available=null_report is not None)
                row.update({"counter_spectra_" + key: psm[key] for key in ("same_peptide_above_cutoff", "different_accepted_target", "different_unaccepted_target", "decoy", "no_reported_row")})
                peptide_rows.append(row)
        gained, lost = len(accepted[upper] - accepted[lower]), len(accepted[lower] - accepted[upper])
        summary = dict(**meta, contrast=contrast, reference_peptides=len(accepted[lower]),
                       tested_peptides=len(accepted[upper]), retained=len(accepted[upper] & accepted[lower]),
                       gained=gained, lost=lost, net=gained-lost)
        for change in ("gained", "lost"):
            for group in ("sequence_added", "sequence_removed", "sequence_retained_above_q_cutoff", "sequence_retained_no_reported_target", "sequence_retained_nonfinite_q"):
                summary[change + "_" + group] = counts[(change, group)]
            assert sum(n for (c, _), n in counts.items() if c == change) == (gained if change == "gained" else lost)
        summaries.append(summary)
        p0, p1 = locations[lower][0], locations[upper][0]
        for accession in p0.keys() & p1.keys():
            assert p0[accession] == p1[accession]
        q0, q1 = assignments[lower], assignments[upper]
        turnover.append(dict(**meta, contrast=contrast, candidate_reference=selection["arms"][lower]["candidate_proteins"],
                             candidate_tested=selection["arms"][upper]["candidate_proteins"],
                             selected_reference=len(p0), selected_tested=len(p1), proteins_added=len(p1.keys()-p0.keys()),
                             proteins_removed=len(p0.keys()-p1.keys()), proteins_retained=len(p0.keys()&p1.keys()),
                             queries_newly_assigned=len(q1.keys()-q0.keys()), queries_no_longer_assigned=len(q0.keys()-q1.keys()),
                             shared_queries_reassigned=sum(q0[q] != q1[q] for q in q0.keys()&q1.keys())))
        for cutoff in (0.001, 0.005, 0.01):
            s0 = {p for p, q in data[lower]["minimum"].items() if q <= cutoff}
            s1 = {p for p, q in data[upper]["minimum"].items() if q <= cutoff}
            thresholds.append(dict(**meta, contrast=contrast, q=cutoff, reference=len(s0), tested=len(s1),
                                   gained=len(s1-s0), lost=len(s0-s1), net=len(s1)-len(s0)))
    print(json.dumps(dict(stage="CASE_COMPLETE", **meta, changed_peptides=len(changed))), flush=True)
    return summaries, peptide_rows, turnover, thresholds, null_rows


def main():
    manifest = read(HERE / "INPUTS.json")
    for record in manifest["frozen"]:
        assert identity(record["path"])["sha256"] == record["sha256"]
    output = HERE / "results"
    output.mkdir(exist_ok=False)
    tables = [[], [], [], [], []]
    for entry in manifest["cases"]:
        for target, rows in zip(tables, run_case(entry)):
            target.extend(rows)
    names = ("CONTRASTS.tsv", "PEPTIDE_MECHANISMS.tsv", "SELECTION_TURNOVER.tsv", "Q_SENSITIVITY.tsv", "NULL_COUNTS.tsv")
    for name, rows in zip(names, tables):
        tsv(output / name, rows)
    groups = defaultdict(list)
    for row in tables[1]:
        groups[(row["contrast"], row["change"], row["category"])].append(row)
    pooled = []
    for (contrast, change, group), rows in sorted(groups.items()):
        pooled.append(dict(contrast=contrast, change=change, category=group, peptide_acquisition_events=len(rows),
                           median_length=statistics.median(r["length"] for r in rows),
                           median_min_q=statistics.median(r["own_min_q"] for r in rows),
                           q_le_0001=sum(r["own_min_q"] <= 0.001 for r in rows),
                           multi_spectrum=sum(r["own_accepted_spectra"] >= 2 for r in rows),
                           shared_flanking_context=sum(r["shared_flanking_context"] for r in rows),
                           same_method_null_overlap=sum(r["own_null_seeds"] > 0 for r in rows)))
    tsv(output / "MECHANISM_TOTALS.tsv", pooled)
    for record in SOURCES.values():
        stat = Path(record["path"]).stat()
        assert (stat.st_size, stat.st_mtime_ns) == (record["bytes"], record["mtime_ns"]), "Input changed during analysis"
    dump(output / "SOURCE_IDENTITIES.json", list(SOURCES.values()))
    dump(output / "VALIDATION.json", dict(status="PASS", utc=datetime.now(timezone.utc).isoformat(),
         job_id=os.environ.get("SLURM_JOB_ID"), cases=len(manifest["cases"]), contrasts=len(tables[0]),
         changed_peptide_comparisons=len(tables[1]), sources=len(SOURCES),
         checks=["Raw canonical peptide and spectrum counts match saved receipts", "Exact comparisons match saved peptide identities",
                 "Search settings agree except FASTA/output", "Selected FASTA/result hashes match frozen receipts",
                 "Accepted changed peptides are present in their own FASTA", "Every changed peptide belongs to exactly one category",
                 "Matched null counts independently recounted", "Source file stats unchanged"],
         interpretation="Observed pipeline stages only; partial nonrepresentative acquisitions, no causal molecular correctness or FDR claim."))
    print(json.dumps(dict(status="PASS", output=str(output), cases=len(manifest["cases"]))))


if __name__ == "__main__":
    main()
