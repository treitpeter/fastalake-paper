#!/usr/bin/env python3
"""Regroup one manuscript cohort/arm, reaggregate annotations, and compare rules.

Only completed searches enter this workflow. Empty acquisitions remain explicit.
Original means the uploaded study rule replayed on exactly the same cohort.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import itertools
import json
import math
import re
import statistics
import subprocess
import sys
import traceback
from collections import Counter, defaultdict
from contextlib import contextmanager
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from fasta_lake import study  # noqa: E402

sys.path.insert(0, str(ROOT / "publication/analyses/simplified_grouping_2026-09-10"))
from hstool_full_comparison import coverage  # noqa: E402

METHODS = ("original", "simplified")
LEVELS = (
    "OG",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "COG_category",
    "CAZy",
    "Pfam",
    "GOs",
    "BRITE",
    "EC",
)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        yield from csv.DictReader(f, delimiter="\t")


@contextmanager
def table_writer(path, fields):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "xt", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        yield w


def write_table(path, records, fields=None):
    records = list(records)
    if fields is None:
        if not records:
            raise ValueError(f"No schema for empty table {path}")
        fields = list(records[0])
    with table_writer(path, fields) as w:
        w.writerows(records)


def median(values):
    values = [v for v in values if v is not None]
    return float(statistics.median(values)) if values else None


def corr(a, b):
    shared = sorted(a.keys() & b.keys())
    x = np.log2([a[k] for k in shared])
    y = np.log2([b[k] for k in shared])
    return (
        len(shared),
        float(np.corrcoef(x, y)[0, 1])
        if len(shared) >= 3 and np.std(x) > 0 and np.std(y) > 0
        else None,
    )


def annotation_terms(value, level):
    if not value or value in {"-", "NA", "N/A"}:
        return ()
    if level == "COG_category":
        return tuple(sorted(set(re.findall(r"[A-Z]", value))))
    values = {x.strip() for x in value.split(",") if x.strip() not in {"", "-"}}
    if level == "KEGG_ko":
        values = {x.removeprefix("ko:") for x in values}
    return tuple(sorted(values))


def fasta_records(path):
    header, parts = None, []
    with Path(path).open() as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:].strip(), []
            elif line.strip():
                parts.append("".join(c for c in line.upper() if "A" <= c <= "Z"))
    if header is not None:
        yield header, "".join(parts)


def load_evidence(task, names, registry):
    sequences, graph, fixed = {}, defaultdict(set), []
    for index, name in enumerate(names, 1):
        folder = Path(task["search_root"]) / name
        config = folder / task["config_name"]
        registry.append({"path": str(config), "sha256": sha(config)})
        fasta = Path(json.loads(config.read_text())["database"]["fasta"])
        members = study._read_fasta(fasta, sequences, registry)
        peptides, identified = set(), defaultdict(set)
        accepted_psms = 0
        for _, r, _ in study._rows(
            folder / "results.sage.tsv",
            {"peptide", "proteins", "label", "peptide_q", "protein_q"},
            registry=registry,
        ):
            if r["label"] != "1" or not study._passing(r["peptide_q"], 0.01):
                continue
            p = sys.intern(study.canonical_peptide(r["peptide"]))
            proteins = study._members(r["proteins"])
            if not proteins <= members:
                raise ValueError(f"Out-of-FASTA observation: {name}")
            peptides.add(p)
            accepted_psms += 1
            for accession in proteins:
                graph[sys.intern(accession)].add(p)
                if study._passing(r["protein_q"], 0.01):
                    identified[accession].add(p)
        fixed.append(
            dict(
                sample=name,
                database_sequences=len(members),
                database_residues=sum(len(sequences[p]) for p in members),
                accepted_canonical_peptides=len(peptides),
                accepted_target_psms=accepted_psms,
                identified_accessions_at_peptide_and_protein_q=len(identified),
                median_observed_peptides_per_identified_accession=median(
                    map(len, identified.values())
                ),
                mean_observed_peptides_per_identified_accession=(
                    statistics.mean(map(len, identified.values())) if identified else None
                ),
                median_observed_accession_coverage_pct=median(
                    coverage(sequences[p], peps)[0] for p, peps in identified.items()
                ),
            )
        )
        if index % 25 == 0 or index == len(names):
            print(f"Evidence reread: {index}/{len(names)} {task['id']}", flush=True)
    return sequences, graph, fixed


def representative_annotations(task, representatives, sequences, registry, out):
    source_fasta = Path(task["annotation_fasta"])
    source_sequences = {h.split()[0]: s for h, s in fasta_records(source_fasta)}
    source = {}
    for r in rows(task["annotations"]):
        p = r["accession"]
        if p in source:
            raise ValueError(f"Duplicate annotation accession: {p}")
        source[p] = {level: annotation_terms(r[level], level) for level in LEVELS}
    registry.extend({"path": str(p), "sha256": sha(p)} for p in [source_fasta, task["annotations"]])
    by_sequence = defaultdict(list)
    if task.get("project_annotations_by_exact_sequence"):
        wanted = {sequences[g] for g in representatives}
        for p, seq in source_sequences.items():
            if seq in wanted and p in source:
                by_sequence[seq].append(p)
    genomes = list(rows(task["genomes"]))
    fallback = {r["Species_rep"]: r["Lineage"] for r in genomes if r["Lineage"]}
    species = {}
    for r in genomes:
        lineage = r["Lineage"] or fallback.get(r["Species_rep"], "")
        terms = [t[3:] for t in lineage.split(";") if t.startswith("s__") and len(t) > 3]
        if terms:
            species[r["Genome"]] = terms[-1]
    registry.append({"path": task["genomes"], "sha256": sha(task["genomes"])})
    reference_species = defaultdict(set)
    if task.get("campi_reference"):
        member_names = {r["campi_os_string"]: r["member"] for r in rows(task["campi_members"])}
        for header, seq in fasta_records(task["campi_reference"]):
            name = header.split("OS=", 1)[1] if "OS=" in header else ""
            name = re.split(r" (?:GN|PE|SV|OX)=", name)[0].strip()
            if name in member_names:
                reference_species[seq].add(member_names[name])
        registry.extend(
            {"path": task[k], "sha256": sha(task[k])} for k in ["campi_reference", "campi_members"]
        )
    result, audit = {}, []
    for g in sorted(representatives):
        if task.get("project_annotations_by_exact_sequence"):
            aliases = sorted(by_sequence.get(sequences[g], []))
            origin = "exact_sequence_projection" if aliases else "unannotated"
        else:
            aliases = [g] if g in source else []
            origin = "same_accession" if aliases else "unannotated"
            if aliases and source_sequences.get(g) != sequences[g]:
                raise ValueError(f"Annotation sequence does not match representative {g}")
        annotation = {
            level: tuple(sorted({t for a in aliases for t in source[a][level]})) for level in LEVELS
        }
        reference = reference_species.get(sequences[g], set())
        genome = re.search(r"MGYG\d+", g)
        taxon = (
            next(iter(reference))
            if len(reference) == 1
            else "Shared_reference_members"
            if reference
            else species.get(genome.group(), "Unassigned")
            if genome
            else "Unassigned"
        )
        annotation["taxon"] = (taxon,) if taxon != "Unassigned" else ()
        result[g] = annotation
        audit.append(
            dict(
                study_group=g,
                sequence_sha256=hashlib.sha256(sequences[g].encode()).hexdigest(),
                annotation_origin=origin,
                annotation_accessions=";".join(aliases),
                representative_taxon=taxon,
                **{k: ";".join(v) for k, v in annotation.items() if k != "taxon"},
            )
        )
    write_table(out / "representative_annotations.tsv.gz", audit)
    return result


def write_matrix(path, names, matrix):
    by_group = defaultdict(dict)
    for sample, values in matrix.items():
        for group, value in values.items():
            by_group[group][sample] = value
    with table_writer(path, ["feature", *names]) as w:
        for group in sorted(by_group):
            w.writerow({"feature": group, **by_group[group]})


def group_metrics(sample, method, stats, sequences, features, total):
    support, covered = [], []
    incompatible_features, incompatible_intensity = 0, 0.0
    for g, (peps, _, _, bad_n, bad_v) in stats.items():
        support.append(len(peps))
        covered.append(coverage(sequences[g], peps)[0])
        incompatible_features += bad_n
        incompatible_intensity += bad_v
    two = sum(n >= 2 for n in support)
    return dict(
        sample=sample,
        method=method,
        quantified_groups=len(stats),
        lfq_features=features,
        intensity=total,
        median_assigned_peptides_per_group=median(support),
        mean_assigned_peptides_per_group=statistics.mean(support) if support else None,
        groups_with_at_least_two_assigned_peptides=two,
        groups_with_at_least_two_assigned_peptides_pct=100 * two / len(stats) if stats else None,
        median_compatible_assigned_coverage_pct=median(covered),
        incompatible_features=incompatible_features,
        incompatible_feature_pct=100 * incompatible_features / features if features else None,
        incompatible_intensity_pct=100 * incompatible_intensity / total if total else None,
    )


def compare_groups(names, direct, dictionary, sequences, out):
    matrices = {m: {s: {} for s in names} for m in METHODS}
    summaries, moves, profiles = [], [], []
    ledger = itertools.groupby(
        rows(direct / "feature_assignments.tsv.gz"), key=lambda r: r["sample"]
    )
    current = next(ledger, None)
    group_fields = [
        "sample",
        "method",
        "study_group",
        "assigned_peptides",
        "lfq_features",
        "intensity",
        "compatible_coverage_pct",
        "unsupported_peptides",
    ]
    diff_fields = [
        "sample",
        "study_group",
        "original_intensity",
        "simplified_intensity",
        "intensity_change",
        "original_assigned_peptides",
        "simplified_assigned_peptides",
    ]
    transfer_fields = [
        "sample",
        "input_line",
        "canonical_peptide",
        "original_group",
        "simplified_group",
        "intensity",
        "original_sequence_compatible",
    ]
    with (
        table_writer(out / "group_evidence.tsv.gz", group_fields) as evidence,
        table_writer(out / "group_diff.tsv.gz", diff_fields) as differences,
        table_writer(out / "changed_feature_assignments.tsv.gz", transfer_fields) as transfers,
        table_writer(
            out / "peptide_long.tsv.gz", ["sample", "canonical_peptide", "intensity"]
        ) as pw,
    ):
        for sample in names:
            stats = {m: defaultdict(lambda: [set(), 0, 0.0, 0, 0.0]) for m in METHODS}
            features, total, changed, changed_v, off_fasta, off_members = 0, 0.0, 0, 0.0, 0, 0
            peptide_values = defaultdict(float)
            present = current is not None and current[0] == sample
            for r in current[1] if present else ():
                pep, value = r["canonical_peptide"], float(r["intensity"])
                candidate = {
                    dictionary["protein_to_group"][p]
                    for p in r["member_set"].split(";")
                    if p in dictionary["protein_to_group"]
                }
                if not candidate:
                    raise ValueError(f"Original feature has no group: {sample}")
                original = min(candidate, key=lambda p: (dictionary["rank"][p], p))
                simple = dictionary["peptide_to_group"][pep]
                if simple != r["study_group"] or pep not in sequences[simple].replace("I", "L"):
                    raise ValueError("Direct dictionary/sequence agreement failed")
                features += 1
                total += value
                peptide_values[pep] += value
                off_fasta += r["representative_in_local_fasta"] == "0"
                off_members += r["representative_in_local_members"] == "0"
                if original != simple:
                    changed += 1
                    changed_v += value
                    transfers.writerow(
                        dict(
                            sample=sample,
                            input_line=r["input_line"],
                            canonical_peptide=pep,
                            original_group=original,
                            simplified_group=simple,
                            intensity=value,
                            original_sequence_compatible=(
                                pep in sequences[original].replace("I", "L")
                            ),
                        )
                    )
                for method, g in [("original", original), ("simplified", simple)]:
                    v = stats[method][g]
                    v[0].add(pep)
                    v[1] += 1
                    v[2] += value
                    bad = pep not in sequences[g].replace("I", "L")
                    v[3] += bad
                    v[4] += value * bad
            if present:
                current = next(ledger, None)
            pw.writerows(
                dict(sample=sample, canonical_peptide=p, intensity=v)
                for p, v in sorted(peptide_values.items())
            )
            moves.append(
                dict(
                    sample=sample,
                    lfq_features=features,
                    changed_features=changed,
                    changed_features_pct=100 * changed / features if features else None,
                    changed_intensity=changed_v,
                    changed_intensity_pct=100 * changed_v / total if total else None,
                    representative_outside_local_fasta_features=off_fasta,
                    representative_outside_local_members_features=off_members,
                )
            )
            for method in METHODS:
                matrix = matrices[method][sample]
                for g, (peps, n, value, _, _) in sorted(stats[method].items()):
                    matrix[g] = value
                    cov, bad_peps = coverage(sequences[g], peps)
                    evidence.writerow(
                        dict(
                            sample=sample,
                            method=method,
                            study_group=g,
                            assigned_peptides=len(peps),
                            lfq_features=n,
                            intensity=value,
                            compatible_coverage_pct=cov,
                            unsupported_peptides=bad_peps,
                        )
                    )
                if not math.isclose(sum(matrix.values()), total, rel_tol=1e-10, abs_tol=1e-6):
                    raise ValueError("Intensity conservation failed")
                summaries.append(
                    group_metrics(sample, method, stats[method], sequences, features, total)
                )
            for g in sorted(stats["original"].keys() | stats["simplified"].keys()):
                a, b = stats["original"].get(g), stats["simplified"].get(g)
                differences.writerow(
                    dict(
                        sample=sample,
                        study_group=g,
                        original_intensity=a[2] if a else 0.0,
                        simplified_intensity=b[2] if b else 0.0,
                        intensity_change=(b[2] if b else 0.0) - (a[2] if a else 0.0),
                        original_assigned_peptides=len(a[0]) if a else 0,
                        simplified_assigned_peptides=len(b[0]) if b else 0,
                    )
                )
            common, r = corr(matrices["original"][sample], matrices["simplified"][sample])
            profiles.append(dict(sample=sample, shared_groups=common, pearson_log2=r))
    if current is not None:
        raise ValueError("Unexpected or out-of-order acquisition in feature ledger")
    remaining = {(s, g) for s, values in matrices["simplified"].items() for g in values}
    for r in rows(direct / "study_group_long.tsv"):
        key = r["sample"], r["study_group"]
        if key not in remaining or not math.isclose(
            float(r["intensity"]), matrices["simplified"][key[0]][key[1]], rel_tol=1e-12
        ):
            raise ValueError("Direct table independent recount failed")
        remaining.remove(key)
    if remaining:
        raise ValueError("Direct table lacks quantified groups")
    for method in METHODS:
        write_matrix(out / f"{method}_study_group_matrix.tsv.gz", names, matrices[method])
    write_table(out / "per_acquisition.tsv", summaries)
    write_table(out / "assignment_changes.tsv", moves)
    write_table(out / "group_profile_agreement.tsv", profiles)
    return matrices, summaries, moves


def rollup_sample(values, annotations, level, mode):
    result = defaultdict(float)
    for g, value in values.items():
        terms = annotations[g][level] or ("Unassigned",)
        weight = value / len(terms) if mode == "split" else value
        for term in terms:
            result[term] += weight
    if mode == "split" and not math.isclose(
        sum(result.values()), sum(values.values()), rel_tol=1e-10, abs_tol=1e-6
    ):
        raise ValueError(f"Split annotation intensity conservation failed for {level}")
    return dict(result)


def reaggregate_functions(names, matrices, annotations, out):
    folder = out / "functional"
    folder.mkdir()
    agreements, resolution = [], []
    for level in (*LEVELS, "taxon"):
        for mode in ("full", "split"):
            rollups = {
                m: {s: rollup_sample(matrices[m][s], annotations, level, mode) for s in names}
                for m in METHODS
            }
            for method in METHODS:
                write_matrix(
                    folder / f"{method}_{level}_{mode}_matrix.tsv.gz", names, rollups[method]
                )
            with table_writer(
                folder / f"{level}_{mode}_diff.tsv.gz",
                ["sample", "term", "original", "simplified", "change"],
            ) as w:
                for sample in names:
                    a, b = rollups["original"][sample], rollups["simplified"][sample]
                    for term in sorted(a.keys() | b.keys()):
                        w.writerow(
                            dict(
                                sample=sample,
                                term=term,
                                original=a.get(term, 0.0),
                                simplified=b.get(term, 0.0),
                                change=b.get(term, 0.0) - a.get(term, 0.0),
                            )
                        )
                    common, r = corr(a, b)
                    agreements.append(
                        dict(
                            sample=sample,
                            level=level,
                            mode=mode,
                            shared_terms=common,
                            pearson_log2=r,
                            terms_original=len(a),
                            terms_simplified=len(b),
                        )
                    )
            if mode == "split":
                for method in METHODS:
                    for sample in names:
                        v = rollups[method][sample]
                        total = sum(matrices[method][sample].values())
                        missing = v.get("Unassigned", 0.0)
                        resolution.append(
                            dict(
                                method=method,
                                sample=sample,
                                level=level,
                                total_intensity=total,
                                unassigned_intensity=missing,
                                annotated_intensity_pct=100 * (total - missing) / total
                                if total
                                else None,
                            )
                        )
            del rollups
        print(f"Functional reaggregation complete: {level}", flush=True)
    write_table(out / "functional_profile_agreement.tsv", agreements)
    write_table(out / "annotation_resolution.tsv", resolution)
    leaders = {}
    with table_writer(
        out / "KO_taxon_long.tsv.gz",
        ["method", "sample", "KO", "representative_taxon", "intensity"],
    ) as w:
        for method in METHODS:
            for sample in names:
                values = defaultdict(lambda: defaultdict(float))
                for g, v in matrices[method][sample].items():
                    taxon = (annotations[g]["taxon"] or ("Unassigned",))[0]
                    for ko in annotations[g]["KEGG_ko"]:
                        values[ko][taxon] += v
                for ko, contribution in sorted(values.items()):
                    for taxon, v in sorted(contribution.items()):
                        w.writerow(
                            dict(
                                method=method,
                                sample=sample,
                                KO=ko,
                                representative_taxon=taxon,
                                intensity=v,
                            )
                        )
                    resolved = {
                        k: v
                        for k, v in contribution.items()
                        if k not in {"Unassigned", "Shared_reference_members"}
                    }
                    total, known = sum(contribution.values()), sum(resolved.values())
                    ranked = sorted(resolved, key=lambda k: (-resolved[k], k))
                    leaders[(method, sample, ko)] = dict(
                        leader=ranked[0] if ranked else "",
                        resolved_fraction=known / total,
                        leader_share=resolved[ranked[0]] / known if known else None,
                    )
    all_cells = sorted({(s, k) for _, s, k in leaders})
    changes = []
    for sample, ko in all_cells:
        a = leaders.get(("original", sample, ko), {})
        b = leaders.get(("simplified", sample, ko), {})
        aa, bb = a.get("leader", ""), b.get("leader", "")
        changes.append(
            dict(
                sample=sample,
                KO=ko,
                original_leader=aa,
                simplified_leader=bb,
                both_resolved=bool(aa and bb),
                leader_changed=(aa != bb) if aa and bb else None,
                resolution_status_changed=bool(aa) != bool(bb),
                original_resolved_fraction=a.get("resolved_fraction"),
                simplified_resolved_fraction=b.get("resolved_fraction"),
            )
        )
    write_table(
        out / "KO_leader_diff.tsv.gz",
        changes,
        fields=[
            "sample",
            "KO",
            "original_leader",
            "simplified_leader",
            "both_resolved",
            "leader_changed",
            "resolution_status_changed",
            "original_resolved_fraction",
            "simplified_resolved_fraction",
        ],
    )
    return agreements, changes


def detection(names, matrices, out):
    rows_out = []
    summaries = []
    for method in METHODS:
        counts = Counter(g for s in names for g in matrices[method][s])
        for group, n in sorted(counts.items()):
            rows_out.append(
                dict(
                    method=method,
                    study_group=group,
                    detected_acquisitions=n,
                    cohort_acquisitions=len(names),
                    detection_pct=100 * n / len(names),
                )
            )
        summaries.append(
            dict(
                method=method,
                detected_anywhere=len(counts),
                detected_at_least_half=sum(
                    n >= math.ceil(0.5 * len(names)) for n in counts.values()
                ),
                detected_at_least_80pct=sum(
                    n >= math.ceil(0.8 * len(names)) for n in counts.values()
                ),
                detected_in_all=sum(n == len(names) for n in counts.values()),
            )
        )
    write_table(
        out / "group_detection.tsv.gz",
        rows_out,
        fields=[
            "method",
            "study_group",
            "detected_acquisitions",
            "cohort_acquisitions",
            "detection_pct",
        ],
    )
    write_table(out / "detection_summary.tsv", summaries)
    return summaries


def report(
    task, direct_summary, summaries, fixed, changes, agreements, leaders, detection_rows, out
):
    metrics = [k for k in summaries[0] if k not in {"method", "sample"}]
    medians = {
        m: {k: median(r[k] for r in summaries if r["method"] == m) for k in metrics}
        for m in METHODS
    }
    diff = [
        dict(
            metric=k,
            original=medians["original"][k],
            simplified=medians["simplified"][k],
            difference=(medians["simplified"][k] - medians["original"][k])
            if None not in (medians["original"][k], medians["simplified"][k])
            else None,
        )
        for k in metrics
    ]
    write_table(out / "median_diff.tsv", diff)
    resolved = [r for r in leaders if r["both_resolved"]]
    features = sum(r["lfq_features"] for r in changes)
    summary = dict(
        status="PASS",
        task=task["id"],
        cohort=task["cohort"],
        arm=task["arm"],
        acquisitions=len(direct_summary["acquisition_names"]),
        acquisition_names=direct_summary["acquisition_names"],
        inferred_groups=direct_summary["inferred_study_groups"],
        pooled_canonical_peptides=direct_summary["pooled_canonical_peptides"],
        eligible_lfq_features=features,
        changed_feature_pct=100 * sum(r["changed_features"] for r in changes) / features
        if features
        else None,
        cohort_medians=medians,
        fixed_pipeline_medians={k: median(r[k] for r in fixed) for k in fixed[0] if k != "sample"},
        detection=detection_rows,
        median_KO_profile_correlation=median(
            r["pearson_log2"] for r in agreements if r["level"] == "KEGG_ko" and r["mode"] == "full"
        ),
        resolved_KO_leader_comparisons=len(resolved),
        changed_KO_leaders=sum(r["leader_changed"] for r in resolved),
        empty_quantified_acquisitions=[
            r["sample"] for r in summaries if r["method"] == "simplified" and not r["lfq_features"]
        ],
        checks=dict(
            exact_cohort=True,
            same_peptide_dictionary=True,
            direct_representative_sequence_compatible=True,
            direct_table_independently_recounted=True,
            group_intensity_conserved=True,
            split_annotation_intensity_conserved=True,
        ),
        scope="Post-search regrouping, complete annotation reaggregation, and original-rule diff",
    )
    (out / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = [
        f"# {task['cohort']} — {task['arm']}",
        "",
        f"Completed comparison on exactly {summary['acquisitions']} declared acquisitions.",
        "Existing search results were reused; no new de novo prediction or Sage search was run.",
        "",
        (
            "Original means the uploaded original rule replayed on this exact "
            "cohort; archived matrices"
        ),
        (
            "may have used a different population or first-accession partition and "
            "are not relabelled here."
        ),
        "",
        "| Per-acquisition median | Original | Simplified | Difference |",
        "|---|---:|---:|---:|",
    ]

    def fmt(value):
        return f"{value:,.4g}" if value is not None else "undefined"

    lines.extend(
        f"| {r['metric']} | {fmt(r['original'])} | "
        f"{fmt(r['simplified'])} | {fmt(r['difference'])} |"
        for r in diff
    )
    lines += [
        "",
        "## Reaggregated outputs",
        "",
        (
            "- `../direct/`: primary simplified group matrix, long table, "
            "dictionary, feature ledger and manifests."
        ),
        (
            "- `original_study_group_matrix.tsv.gz` and "
            "`simplified_study_group_matrix.tsv.gz`: paired matrices."
        ),
        (
            "- `group_diff.tsv.gz`: every acquisition/group in the union, both "
            "intensities and peptide counts."
        ),
        (
            "- `group_evidence.tsv.gz`: assigned peptide support and compatible "
            "coverage for both methods."
        ),
        (
            "- `changed_feature_assignments.tsv.gz`: every changed LFQ assignment "
            "and its original sequence compatibility."
        ),
        "- `peptide_long.tsv.gz`: canonical peptide intensity, fixed between assignment methods.",
        (
            "- `functional/`: all eleven annotation levels plus representative "
            "taxon; original, simplified and diffs."
        ),
        (
            "- `KO_taxon_long.tsv.gz`, `KO_leader_diff.tsv.gz`: contributors, "
            "resolved leaders and changes."
        ),
        (
            "- `annotation_resolution.tsv`, `representative_annotations.tsv.gz`: "
            "coverage and annotation provenance."
        ),
        "- `group_detection.tsv.gz`, `detection_summary.tsv`: cohort-level completeness.",
        (
            "- `group_profile_agreement.tsv`, `functional_profile_agreement.tsv`: "
            "paired log2 profile agreement."
        ),
        (
            "- `fixed_pipeline_metrics.tsv`: search FASTA sizes, accepted evidence,"
            " observed accession support/coverage."
        ),
        "",
        "## Definitions and interpretation",
        "",
        (
            "Peptide support in quantified groups counts distinct assigned "
            "canonical peptides with eligible LFQ signal."
        ),
        (
            "Accession support counts distinct observed peptides on target rows "
            "passing peptide and protein q <= 0.01."
        ),
        (
            "Coverage unions exact I/L-equivalent occurrences; incompatible "
            "peptides contribute no residues."
        ),
        (
            "Higher quantified group counts do not imply additional "
            "identifications. Blank matrix cells remain missing."
        ),
        (
            "Zero is used only for absent entries in diff arithmetic. Empty "
            "acquisitions remain explicit."
        ),
        "",
        (
            "Functional `full` mode contributes a group's entire intensity to each "
            "distinct annotation term; it is nonadditive."
        ),
        (
            "Functional `split` mode divides intensity equally between distinct "
            "terms within a level, retaining Unassigned."
        ),
        (
            "Both modes are applied identically to original and simplified group "
            "assignments; their differences isolate assignment."
        ),
        (
            "These are representative-based rollups at the retained peptide gate, "
            "not the old redundant threshold products."
        ),
        (
            "CAMPI standalone-control annotations, where requested, use exact "
            "full-sequence projection from the current"
        ),
        (
            "reference-plus-metagenome annotation source; identical-sequence alias "
            "terms are unioned and recorded."
        ),
        (
            "Taxon labels use the same representative lookup in both methods: exact"
            " CAMPI reference sequence when available,"
        ),
        (
            "otherwise UHGG genome metadata, otherwise Unassigned. These labels do "
            "not prove unique organism origin."
        ),
        (
            "Pooled detection and acquisition-profile correlations are descriptive;"
            " they do not replace patient/specimen-aware tests."
        ),
        "Study-group FDR and unique biological origin are not established by this rerun.",
        "",
        "## Tomorrow's figure review",
        "",
        (
            "Use these diffs to rebuild affected figure source tables and inspect "
            "their biological conclusions."
        ),
        (
            "Existing manuscript artwork is preserved. Identification-only "
            "endpoints retain their original units."
        ),
        (
            "Group-level quantities and functional contributors need reassessment; "
            "downstream captions must name the rule."
        ),
        "",
    ]
    (out / "REPORT.md").write_text("\n".join(lines))
    return summary


def run_task(task, output, baseline_path):
    output.mkdir(parents=True, exist_ok=False)
    names = Path(task["samples_file"]).read_text().splitlines()
    if len(names) != task["acquisitions"] or len(set(names)) != len(names):
        raise ValueError("Cohort roster changed")
    names = sorted(names)
    registry = [
        {"path": task["samples_file"], "sha256": sha(task["samples_file"])},
        {"path": str(baseline_path), "sha256": sha(baseline_path)},
    ]
    direct = output / "direct"
    print(f"Starting direct grouping: {task['id']} ({len(names)} acquisitions)", flush=True)
    direct_summary = study.group_searches(
        task["search_root"], direct, config_name=task["config_name"], acquisition_names=names
    )
    sequences, graph, fixed = load_evidence(task, names, registry)
    spec = importlib.util.spec_from_file_location("original_study_rule", baseline_path)
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    dictionary = original.study_dictionary(graph, set(sequences))
    peptide_map = {
        r["canonical_peptide"]: r["study_group"]
        for r in rows(direct / "peptide_to_study_group.tsv")
    }
    if peptide_map != dictionary["peptide_to_group"]:
        raise ValueError("Original and simplified inference dictionaries differ")
    del graph, peptide_map
    out = output / "reaggregated"
    out.mkdir()
    write_table(out / "fixed_pipeline_metrics.tsv", fixed)
    write_table(
        out / "original_protein_to_group.tsv.gz",
        [dict(protein=p, study_group=g) for p, g in sorted(dictionary["protein_to_group"].items())],
    )
    matrices, summaries, changes = compare_groups(names, direct, dictionary, sequences, out)
    annotations = representative_annotations(
        task, set(dictionary["rank"]), sequences, registry, out
    )
    del sequences, dictionary
    agreement, leaders = reaggregate_functions(names, matrices, annotations, out)
    detection_rows = detection(names, matrices, out)
    summary = report(
        task, direct_summary, summaries, fixed, changes, agreement, leaders, detection_rows, out
    )
    for path in [
        Path(__file__),
        Path(study.__file__),
        ROOT / "publication/analyses/simplified_grouping_2026-09-10/hstool_full_comparison.py",
    ]:
        registry.append({"path": str(path), "sha256": sha(path)})
    registry.append(
        {"path": str(direct / "input_manifest.json"), "sha256": sha(direct / "input_manifest.json")}
    )
    (out / "PROVENANCE.json").write_text(json.dumps(registry, indent=2) + "\n")
    (output / "COMPLETE.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--task-index", type=int, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    task = manifest["tasks"][args.task_index]
    output = Path(manifest["output_root"]) / "cohorts" / task["id"]
    if output.exists():
        raise FileExistsError(f"Refusing to reuse task output: {output}")
    commit = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True
    ).strip()
    if commit != manifest["code_commit"]:
        raise ValueError("Runtime code commit changed")
    for record in task["static_inputs"]:
        if sha(record["path"]) != record["sha256"]:
            raise ValueError(f"Static input changed: {record['path']}")
    for record in task["search_inputs"]:
        stat = Path(record["path"]).stat()
        if stat.st_size != record["bytes"] or stat.st_mtime_ns != record["mtime_ns"]:
            raise ValueError(f"Search input changed: {record['path']}")
    try:
        summary = run_task(task, output, Path(manifest["baseline_module"]))
        print(json.dumps(summary, indent=2), flush=True)
    except Exception as error:
        output.mkdir(parents=True, exist_ok=True)
        (output / "FAILED.json").write_text(
            json.dumps(
                dict(task=task["id"], error=str(error), traceback=traceback.format_exc()), indent=2
            )
            + "\n"
        )
        raise


if __name__ == "__main__":
    main()
