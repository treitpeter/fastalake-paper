#!/usr/bin/env python3
"""Run a frozen-database peptide entrapment assay on a recorded HSTOOL transfer.

FDRBench generates disjoint target/entrapment/decoy quartets after the donor
database is frozen. Sage searches the explicit peptide FASTA without digestion.
This is a peptide-level assay variant, not an unchanged protein database search.
All inputs, commands and raw outputs are retained. No nominal threshold is tuned.
"""

import argparse
import collections
import csv
import hashlib
import json
import math
import os
import platform
import re
import subprocess
import time
from pathlib import Path

THRESHOLDS = (0.001, 0.005, 0.01, 0.02, 0.05)


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n")


def read_pairs(path):
    sequences, pairs = {}, collections.defaultdict(dict)
    with Path(path).open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            seq = row["sequence"]
            # Sage 0.14.6 includes selenocysteine U and pyrrolysine O.
            if seq != seq.replace("I", "L") or not re.fullmatch("[ACDEFGHKLMNOPQRSTUVWY]+", seq):
                raise ValueError(f"Noncanonical generated sequence: {seq}")
            kind, index = row["peptide_type"], row["peptide_pair_index"]
            if seq in sequences or kind in pairs[index]:
                raise ValueError("Duplicate generated sequence or pairing member")
            if (row["decoy"] == "Yes") != (kind in {"decoy", "p_decoy"}):
                raise ValueError("Incorrect decoy label")
            sequences[seq] = (kind, index)
            pairs[index][kind] = seq
    for pair in pairs.values():
        if set(pair) != {"target", "p_target", "decoy", "p_decoy"}:
            raise ValueError("Incomplete quartet")
        target = pair["target"]
        for seq in pair.values():
            if collections.Counter(seq) != collections.Counter(target) or seq[-1] != target[-1]:
                raise ValueError("Composition or terminal residue differs within quartet")
    if not pairs:
        raise ValueError("Empty peptide database")
    return sequences, pairs


def summarize(results, pair_file, out):
    sequences, pairs = read_pairs(pair_file)
    best = {}
    psm_rows = 0
    with Path(results).open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            psm_rows += 1
            peptide = re.sub(r"\[[^\]]*\]", "", row["peptide"]).replace("I", "L")
            if peptide not in sequences:
                raise ValueError(f"Search sequence outside audited database: {peptide}")
            kind, index = sequences[peptide]
            label = int(row["label"])
            if label not in {-1, 1} or (label == -1) != (kind in {"decoy", "p_decoy"}):
                raise ValueError("Search label disagrees with quartet")
            q, pep = float(row["peptide_q"]), float(row["posterior_error"])
            if not (math.isfinite(q) and 0 <= q <= 1 and math.isfinite(pep)):
                raise ValueError("Invalid search q value or posterior error")
            if label == -1:
                continue
            # Collapse modifications and charge states to the declared canonical
            # sequence population; rank by search q, then posterior error.
            score = (q, pep)
            if peptide not in best or score < best[peptide][0]:
                best[peptide] = (score, kind, index)
    with (out / "canonical_scores.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["peptide", "q_value", "posterior_error", "kind", "pair_index"])
        for peptide, (score, kind, index) in sorted(best.items(), key=lambda item: item[1][0]):
            writer.writerow([peptide, *score, kind, index])
    summaries = []
    for floor in (5, 9):
        for threshold in THRESHOLDS:
            accepted = {p: v for p, v in best.items() if len(p) >= floor and v[0][0] <= threshold}
            targets = sum(v[1] == "target" for v in accepted.values())
            entrapments = len(accepted) - targets
            missing = higher = ties = 0
            for peptide, (score, kind, index) in accepted.items():
                if kind != "p_target":
                    continue
                target = pairs[index]["target"]
                if target not in accepted:
                    missing += 1
                elif score < accepted[target][0]:
                    higher += 1
                elif score == accepted[target][0]:
                    ties += 1
            total = targets + entrapments
            numerator = entrapments + missing + 2 * higher
            summaries.append(
                dict(
                    min_length=floor,
                    q_threshold=threshold,
                    n_target=targets,
                    n_entrapment=entrapments,
                    partner_undiscovered=missing,
                    entrapment_ranks_higher=higher,
                    exact_score_ties=ties,
                    paired_numerator=numerator,
                    lower_fraction=entrapments / total if total else None,
                    combined_estimate=2 * entrapments / total if total else None,
                    paired_estimate=numerator / total if total else None,
                    paired_ties_upper=(numerator + 2 * ties) / total if total else None,
                )
            )
    with (out / "summary.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)
    return dict(
        psm_rows=psm_rows,
        canonical_target_labelled_sequences=len(best),
        paired_quartets=len(pairs),
        summaries=summaries,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical", type=Path, required=True)
    parser.add_argument("--sage", type=Path, required=True)
    parser.add_argument("--fdrbench", type=Path, required=True)
    parser.add_argument("--java", default="java", help="Java 11+ executable on the compute node")
    parser.add_argument("--target", required=True)
    parser.add_argument("--tier", choices=["T1_within_group", "T0_self"], default="T1_within_group")
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    manifest_path = (
        args.canonical.resolve()
        / "final_revision/thread_selection_math/manifests"
        / (f"stool_{args.target}__{args.tier}.json")
    )
    manifest = json.loads(manifest_path.read_text())
    if args.tier == "T1_within_group" and manifest["target"] == manifest["donor"]:
        raise ValueError("Transfer construction must use another acquisition")
    db = Path(manifest["database"])
    if sha(db) != manifest["database__sha256"]:
        raise ValueError("Frozen donor database changed")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    sage, jar = args.sage.resolve(), args.fdrbench.resolve()
    record = dict(
        status="running",
        target=manifest["target"],
        donor=manifest["donor"],
        group=manifest["target_group"],
        tier=args.tier,
        seed=args.seed,
        frozen_database=str(db),
        frozen_database_sha256=sha(db),
        sage_sha256=sha(sage),
        fdrbench_jar_sha256=sha(jar),
        manifest_sha256=sha(manifest_path),
        source_script_sha256=sha(__file__),
        source_commit=subprocess.check_output(
            ["git", "-C", str(Path(__file__).resolve().parent), "rev-parse", "HEAD"], text=True
        ).strip(),
        host=platform.node(),
        python=platform.python_version(),
        slurm_job_id=os.environ.get("SLURM_JOB_ID"),
        commands=[],
    )
    write_json(out / "RUN.json", record)

    def run(stage, argv):
        start = time.monotonic()
        with (out / (stage + ".log")).open("x") as log:
            code = subprocess.run(
                list(map(str, argv)), stdout=log, stderr=subprocess.STDOUT
            ).returncode
        record["commands"].append(
            dict(
                stage=stage,
                argv=list(map(str, argv)),
                returncode=code,
                seconds=time.monotonic() - start,
            )
        )
        write_json(out / "RUN.json", record)
        if code:
            raise RuntimeError(f"{stage} failed; see retained log")

    try:
        pair_file = out / "pairs.tsv"
        run(
            "build_entrapment",
            [
                args.java,
                "-Xmx8G",
                "-XX:ActiveProcessorCount=4",
                "-jar",
                jar,
                "-db",
                db,
                "-o",
                pair_file,
                "-I2L",
                "-level",
                "peptide",
                "-enzyme",
                "1",
                "-miss_c",
                "2",
                "-minLength",
                "5",
                "-maxLength",
                "50",
                "-fix_nc",
                "c",
                "-seed",
                args.seed,
                "-fix_seed",
                "-decoy",
            ],
        )
        sequences, pairs = read_pairs(pair_file)
        # Use the audited table to make labels explicit to Sage. Decoys are
        # generated by FDRBench alongside controls, with global nonoverlap.
        search_fasta = out / "search.fasta"
        with search_fasta.open("w") as stream:
            for peptide, (kind, index) in sequences.items():
                prefix = "rev_" if kind in {"decoy", "p_decoy"} else ""
                stream.write(f">{prefix}{kind}_{index}\n{peptide}\n")
        config_path = db.parent / "config.json"
        config = json.loads(config_path.read_text())
        config["database"].update(fasta=str(search_fasta), generate_decoys=False, decoy_tag="rev_")
        config["database"]["enzyme"].update(cleave_at="$", restrict=None, missed_cleavages=0)
        search = out / "search"
        search.mkdir()
        config["output_directory"] = str(search)
        config["mzml_paths"] = [manifest["mzml"]]
        write_json(out / "config.json", config)
        record.update(
            pair_quartets=len(pairs),
            pairs_sha256=sha(pair_file),
            search_fasta_sha256=sha(search_fasta),
            search_config_sha256=sha(out / "config.json"),
            original_config_sha256=sha(config_path),
        )
        write_json(out / "RUN.json", record)
        run("search", [sage, out / "config.json"])
        result = search / "results.sage.tsv"
        record.update(summarize(result, pair_file, out))
        record.update(
            status="COMPLETE",
            result_sha256=sha(result),
            mzml_sha256=sha(manifest["mzml"]),
            interpretation="Conditional paired-peptide assay; not known-truth FDR. "
            "Canonical sequence counts, technical dependencies, possible biological "
            "control overlap, null comparability and assay changes require assessment.",
        )
    except Exception as error:
        record.update(status="FAILED", error=str(error))
        write_json(out / "RESULT.json", record)
        raise
    write_json(out / "RESULT.json", record)
    print(json.dumps({k: record[k] for k in ["status", "target", "donor", "seed", "summaries"]}))


if __name__ == "__main__":
    main()
