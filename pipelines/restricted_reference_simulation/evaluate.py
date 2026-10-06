"""Truth evaluation after all searches, with a separately implemented recount."""

import copy
import csv
import json
from pathlib import Path
import sys

import pandas as pd

from common import (
    PREVIOUS,
    THRESHOLDS,
    digest_intervals,
    execute_main,
    fasta,
    legacy,
    location,
    module,
    read_table,
    save,
    sha,
    table,
)

sys.path.insert(0, str(PREVIOUS))
truth_api = module("previous_truth_evaluator", PREVIOUS / "truth_evaluator.py")


def truth_for_search(truth, configured_inputs, reported_filenames):
    """Resolve Sage basenames only within an explicitly single-input search."""
    if len(configured_inputs) != 1 or len(reported_filenames) != 1:
        raise ValueError("Ambiguous input/filename mapping for independent truth")
    source = configured_inputs[0]
    reported = next(iter(reported_filenames))
    if {filename for filename, _ in truth} != {source}:
        raise ValueError("Truth does not describe the configured input acquisition")
    if reported not in {source, Path(source).name}:
        raise ValueError("Reported search filename does not identify the configured input")
    return {(reported, scan): sequences for (_, scan), sequences in truth.items()}


def independent_peptides(path, generating, threshold=0.01, low=5, high=50):
    """Pandas recount independent of the existing row-based truth evaluator."""
    frame = pd.read_csv(path, sep="\t")
    frame = frame.loc[frame.label == 1].copy()
    frame["sequence"] = frame.peptide.str.replace(r"\[[^\]]*\]", "", regex=True).str.replace(
        "I", "L"
    )
    best = frame.groupby("sequence").peptide_q.min()
    accepted = {p for p, q in best.items() if q <= threshold and low <= len(p) <= high}
    return accepted, accepted - generating


def partition_discoveries(first, second, generating):
    return {
        name: dict(
            discoveries=len(peptides),
            false_discoveries=len(peptides - generating),
            realised_FDP=len(peptides - generating) / len(peptides) if peptides else None,
        )
        for name, peptides in [
            ("first_only", first - second),
            ("second_only", second - first),
            ("shared", first & second),
        ]
    }


def diagnostic_with_truth(base, name, generating):
    path = base / "searches" / name / "sage/results.sage.tsv"
    bundle = base / "search_inputs" / name
    pairs = {r["target"]: r["entrapment"] for r in read_table(bundle / "peptide_pairs.tsv")}
    scores, classes = legacy.read_scores(path)
    paired = set(pairs) | set(pairs.values())
    control_digest, target_digest = set(), set()
    for accession, sequence in fasta(bundle / "search.fasta"):
        dest = control_digest if accession.startswith("E_") else target_digest
        dest.update(sequence[a:b].replace("I", "L") for a, b in digest_intervals(sequence))
    collisions = dict(
        control_digest=len(control_digest),
        target_digest=len(target_digest),
        control_sequences_generated=len(control_digest & generating),
        control_only_sequences_generated=len((control_digest - target_digest) & generating),
        paired_control_sequences_generated=len(set(pairs.values()) & generating),
    )
    rows = []
    for factor in (1, 0.1):
        modified_scores = {p: q * factor for p, q in scores.items()}
        for q in THRESHOLDS:
            for low in (5, 9):
                report = legacy.estimates(modified_scores, classes, pairs, q, low)
                accepted = {
                    p for p, score in modified_scores.items() if score <= q and low <= len(p) <= 50
                }
                unpaired = accepted - paired
                rows.append(
                    dict(
                        arm=name,
                        q_multiplier=factor,
                        **report,
                        actual_false=len(accepted - generating),
                        actual_fdp=len(accepted - generating) / len(accepted) if accepted else None,
                        unpaired_true=len(unpaired & generating),
                        unpaired_false=len(unpaired - generating),
                    )
                )
    score_rows = [
        dict(
            arm=name,
            peptide=p,
            peptide_q=q,
            membership=";".join(sorted(classes[p])),
            generated=p in generating,
            paired=p in paired,
        )
        for p, q in sorted(scores.items())
    ]
    return rows, score_rows, collisions


def main(phase, index):
    unit, base = location(phase, index)
    complete = json.loads((base / "ROUTES_COMPLETE.json").read_text())
    generated = json.loads((base / "GENERATED.json").read_text())
    for path, expected in complete["searches"].items():
        assert sha(base / path) == expected
    for acq in ("donor", "test"):
        name = f"truth/{acq}.tsv"
        assert sha(base / name) == generated["files"][name]
    truth_path = base / "truth/test.tsv"
    truth = truth_api.load_truth(truth_path)
    generating = set().union(*truth.values())
    assert len(truth) == len(generating) == 2048
    donor_truth = truth_api.load_truth(base / "truth/donor.tsv")
    assert set().union(*donor_truth.values()) == generating
    output = base / "evaluation"
    output.mkdir(exist_ok=False)
    errors, primary, sensitivity, scores_all, diagnostic, ledgers = [], [], [], [], [], []
    collisions, accepted_sets = {}, {}
    ordinary_gate, audited_files, searches_joined = True, {}, {}
    for path in sorted((base / "searches").glob("*/COMPLETE.json")):
        record = json.loads(path.read_text())
        name, search = record["name"], Path(record["results"])
        assert sha(search) == record["results_sha256"]
        assert sha(record["database"]) == record["database_sha256"]
        assert sha(path.parent / "config.json") == record["config_sha256"]
        config = json.loads((path.parent / "config.json").read_text())
        assert config["mzml_paths"] == [str(base / "spectra/test.mzML")]
        with search.open() as handle:
            all_rows = list(csv.DictReader(handle, delimiter="\t"))
        search_truth = truth_for_search(
            truth, config["mzml_paths"], {r["filename"] for r in all_rows}
        )
        augmented = name not in ("fixed", "self", "independent")
        rows, seen = truth_api.read_identifications(
            search, search_truth, "augmented" if augmented else "ordinary"
        )
        evaluated = truth_api.evaluate(search_truth, rows)
        errors.extend(dict(arm=name, **r) for r in evaluated)
        targets, false = independent_peptides(search, generating)
        primary_row = next(
            r
            for r in evaluated
            if r["unit"] == "canonical_peptide_at_peptide_q"
            and r["population"] == "augmented_or_ordinary_all"
            and r["min_residues"] == 5
            and r["max_residues"] == 50
            and r["q_threshold"] == 0.01
        )
        assert (len(targets), len(false)) == (
            primary_row["discoveries"],
            primary_row["false_discoveries"],
        )
        accepted_sets[name] = targets
        decoys = sum(r["label"] == "-1" for r in all_rows)
        primary.append(
            dict(
                arm=name,
                discoveries=len(targets),
                false_discoveries=len(false),
                realised_FDP=len(false) / max(len(targets), 1),
                decoy_winners=decoys,
                empty_discovery_set=not targets,
            )
        )
        if not augmented:
            ordinary_gate = ordinary_gate and bool(targets) and decoys > 0
        for row in rows:
            if row["peptide_q"] <= 0.05:
                ledgers.append(
                    dict(
                        arm=name,
                        filename=row["key"][0],
                        scan=row["key"][1],
                        peptide=row["peptide"],
                        peptide_q=row["peptide_q"],
                        psm_correct=row["psm_correct"],
                        peptide_correct=row["peptide"] in generating,
                        membership=";".join(sorted(row["classes"])),
                    )
                )
        perturbed = copy.deepcopy(rows)
        for row in perturbed:
            row["peptide_q"] *= 0.1
            row["spectrum_q"] *= 0.1
        sensitivity.extend(
            dict(arm=name, q_multiplier=0.1, **r)
            for r in truth_api.evaluate(search_truth, perturbed)
        )
        if augmented:
            diag, score_rows, collision = diagnostic_with_truth(base, name, generating)
            diagnostic.extend(diag)
            scores_all.extend(score_rows)
            collisions[name] = collision
        audited_files[str(search)] = sha(search)
        searches_joined[name] = dict(
            winners=len(seen),
            truth_scans_without_winner=len(set(search_truth) - seen),
            configured_input=config["mzml_paths"][0],
            reported_filename=next(iter(search_truth))[0],
        )
    # Verify donor exclusion in the actual recorded inputs and all three commands.
    selection_audit = {}
    for path in sorted((base / "selection").glob("*/SELECTION.json")):
        record = json.loads(path.read_text())
        prediction = Path(record["prediction"])
        expected_acquisition = "donor" if "independent" in path.parent.name else "test"
        assert (
            prediction
            == base
            / "predictions"
            / expected_acquisition
            / f"{expected_acquisition}_predictions.csv"
        )
        assert sha(prediction) == record["prediction_sha256"]
        for source, expected in (record["inputs"] | record["stage_fastas"]).items():
            assert sha(source) == expected
            audited_files[source] = expected
        pool = path.parent / "prediction_input"
        assert len(list(pool.iterdir())) == 1
        assert next(pool.iterdir()).resolve() == prediction.resolve()
        for stage in ("1_cohort", "3_acquisition"):
            command = json.loads((path.parent / stage / "COMMAND.json").read_text())["argv"]
            assert command[command.index("--peptides-dir") + 1] == str(pool)
        command = json.loads((path.parent / "4_razor_logs/COMMAND.json").read_text())["argv"]
        assert command[command.index("--peptides") + 1] == str(prediction)
        selection_audit[path.parent.name] = expected_acquisition
    incremental = []
    for comparator in ("fixed", "independent"):
        for q in THRESHOLDS:
            for low, high in ((5, 50), (9, 50), (5, 8)):
                first, _ = independent_peptides(
                    base / "searches/self/sage/results.sage.tsv", generating, q, low, high
                )
                second, _ = independent_peptides(
                    base / f"searches/{comparator}/sage/results.sage.tsv", generating, q, low, high
                )
                for name, counts in partition_discoveries(first, second, generating).items():
                    incremental.append(
                        dict(
                            comparison=f"self_vs_{comparator}",
                            subset=name,
                            q_threshold=q,
                            min_residues=low,
                            max_residues=high,
                            **counts,
                        )
                    )
    for name, data in [
        ("true_errors", errors),
        ("primary", primary),
        ("misleading_q_check", sensitivity),
        ("entrapment_vs_truth", diagnostic),
        ("observed_scores", scores_all),
        ("identification_ledger", ledgers),
        ("incremental", incremental),
    ]:
        table(output / f"{name}.tsv", data)
    save(output / "CONTROL_COLLISIONS.json", collisions)
    save(
        output / "AUDIT.json",
        dict(
            unit=unit,
            independent_recount_pass=True,
            matched_truth_pass=True,
            donor_exclusion_pass=True,
            selection_inputs=selection_audit,
            ordinary_outputs_assessable=ordinary_gate,
            searches_joined=searches_joined,
            input_hashes=audited_files,
            outputs={p.name: sha(p) for p in output.iterdir() if p.is_file()},
            whole_procedure_fdr_control_established=False,
        ),
    )
    print(unit, "truth evaluation and independent recount complete", flush=True)


if __name__ == "__main__":
    execute_main(main)
