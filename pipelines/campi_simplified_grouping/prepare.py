#!/usr/bin/env python3
"""Freeze explicit four-cohort/eight-arm inputs for smoke and scheduled runs."""

import argparse
import copy
import csv
import datetime
import json
import subprocess
from pathlib import Path

from run_cohort import sha


def read_manifest(path, key):
    with path.open() as f:
        return [r[key] for r in csv.DictReader(f, delimiter="\t")]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fastalake-root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    fl, out = args.fastalake_root.resolve(), args.out.resolve()
    code = Path(__file__).resolve().parents[3]
    commit = subprocess.check_output(
        ["git", "-C", str(code), "rev-parse", "HEAD"], text=True
    ).strip()
    if subprocess.check_output(
        ["git", "-C", str(code), "status", "--porcelain", "--untracked-files=no"], text=True
    ).strip():
        raise ValueError("Commit code before preparing execution manifests")
    can = fl / "canonical"
    arms = can / "publication_data/arms"
    baseline = fl / "release_prep/github_biorxiv_2026-09-07/repo/software/fasta_lake/study.py"
    sources = {
        "HSTOOL": code / "publication/analyses/simplified_grouping_2026-09-10/hstool_samples.txt",
        "PREDICT": code / "publication/cohorts/predict_sample_list.tsv",
        "CAPSCAN": code / "publication/cohorts/CAPSCAN/conservative_331.tsv",
        "CAMPI": code / "publication/cohorts/campi_mock_sample_list.tsv",
    }
    names = {
        "HSTOOL": sources["HSTOOL"].read_text().splitlines(),
        "PREDICT": read_manifest(sources["PREDICT"], "sample"),
        "CAPSCAN": read_manifest(sources["CAPSCAN"], "sample_id"),
        "CAMPI": read_manifest(sources["CAMPI"], "sample"),
    }
    expected = {"HSTOOL": 58, "PREDICT": 306, "CAPSCAN": 331, "CAMPI": 15}
    for c, values in names.items():
        if len(values) != expected[c] or len(set(values)) != expected[c]:
            raise ValueError(f"Incorrect manuscript roster: {c}")
        if any(Path(n).name != n or n in {"", ".", ".."} for n in values):
            raise ValueError(f"Invalid sample names: {c}")
    configs = [
        (c, "noclust", "noclust_2026-08-07", c, False) for c in ["HSTOOL", "PREDICT", "CAPSCAN"]
    ]
    configs += [
        ("CAMPI", label, arm, "CAMPI_MOCK", False)
        for label, arm in [
            ("FL_UNIPROT", "uniprotkb_2026_03_noclust_2026-09-09"),
            ("FL_REF", "uniprotkb_2026_03_noclust_campi_ref_append_2026-09-09"),
            ("FL_REF_MG", "uniprotkb_2026_03_noclust_campi_ref_metag_append_2026-09-09"),
        ]
    ]
    configs += [
        ("CAMPI", label, arm, "CAMPI_MOCK", True)
        for label, arm in [("REF", "campi_ref_solo"), ("MG", "campi_mo_fixed")]
    ]
    out.mkdir(parents=True, exist_ok=False)
    (out / "logs").mkdir()
    (out / "rosters").mkdir()
    (out / "smoke").mkdir()
    (out / "smoke/logs").mkdir()
    (out / "smoke/rosters").mkdir()
    for c, values in names.items():
        (out / "rosters" / f"{c}.txt").write_text("\n".join(sorted(values)) + "\n")
    tasks = []
    hash_cache = {}

    def hashed(path):
        path = str(path)
        if path not in hash_cache:
            hash_cache[path] = sha(path)
        return {"path": path, "sha256": hash_cache[path]}

    genomes = (
        fl.parent.parent
        / "Peter_Microbpred/metalab_magnify_databases/human_gut/genomes-all_metadata.tsv"
    )
    # fl is .../Peter/FastaLake; Peter_Microbpred is beside Peter.
    for index, (cohort, label, arm, folder_cohort, control) in enumerate(configs):
        search = (
            arms / arm / folder_cohort
            if control
            else arms / arm / "cohorts" / folder_cohort / "5_sage/razor_hash-acc"
        )
        annotation_arm = (
            "uniprotkb_2026_03_noclust_campi_ref_metag_append_2026-09-09" if control else arm
        )
        annotation_dir = arms / annotation_arm / "cohorts" / folder_cohort / "11_eggnog"
        annotation = annotation_dir / f"{folder_cohort}_eggnog_accession_keyed.tsv"
        annotation_fasta = annotation_dir / f"{folder_cohort}_selected.faa"
        inputs = {}
        for name in names[cohort]:
            folder = search / name
            config = folder / "results.json"
            fasta = Path(json.loads(config.read_text())["database"]["fasta"])
            if not fasta.is_absolute():
                raise ValueError(f"Nonabsolute searched FASTA: {config}")
            for path in [config, folder / "results.sage.tsv", folder / "lfq.tsv", fasta]:
                if not path.is_file() or not path.stat().st_size:
                    raise ValueError(f"Missing/empty required input: {path}")
                st = path.stat()
                inputs[str(path)] = dict(path=str(path), bytes=st.st_size, mtime_ns=st.st_mtime_ns)
        task = dict(
            index=index,
            id=f"{cohort}_{label}",
            cohort=cohort,
            arm=label,
            acquisitions=expected[cohort],
            samples_file=str(out / "rosters" / f"{cohort}.txt"),
            search_root=str(search),
            config_name="results.json",
            annotations=str(annotation),
            annotation_fasta=str(annotation_fasta),
            project_annotations_by_exact_sequence=control,
            genomes=str(genomes),
            search_inputs=list(inputs.values()),
        )
        static = [
            Path(task["samples_file"]),
            sources[cohort],
            annotation,
            annotation_fasta,
            genomes,
            baseline,
        ]
        if cohort == "CAMPI":
            task["campi_reference"] = str(can / "publication/figures/Fig2/tables/campi_targets.faa")
            task["campi_members"] = str(can / "publication/figures/Fig2/tables/campi_members.tsv")
            static.extend([Path(task["campi_reference"]), Path(task["campi_members"])])
        task["static_inputs"] = [hashed(p) for p in static]
        tasks.append(task)
        print(
            f"Preflight {task['id']}: {expected[cohort]} acquisitions, all files present",
            flush=True,
        )
    manifest = dict(
        schema_version=1,
        mode="full_manuscript_cohorts",
        code_commit=commit,
        baseline_module=str(baseline),
        output_root=str(out),
        cohort_counts=expected,
        start_local="2026-09-10T22:00:00+02:00",
        timezone="Europe/Berlin",
        partitions=["p.hpcl8", "p.hpcl93", "p.hpcl92"],
        max_concurrent_tasks=8,
        prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        tasks=tasks,
        code_source=str(code),
        preparation_script_sha256=sha(Path(__file__)),
    )
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    smoke = copy.deepcopy(manifest)
    smoke["mode"] = "smoke_subsets_not_full_cohorts"
    smoke["output_root"] = str(out / "smoke")
    smoke["start_local"] = "immediate"
    for task in smoke["tasks"]:
        values = names[task["cohort"]]
        by_size = sorted(
            values, key=lambda n: (Path(task["search_root"]) / n / "lfq.tsv").stat().st_size
        )
        selected = sorted({by_size[0], by_size[-1]})
        roster = out / "smoke/rosters" / f"{task['id']}.txt"
        roster.write_text("\n".join(selected) + "\n")
        old_roster = task["samples_file"]
        task["samples_file"] = str(roster)
        task["acquisitions"] = len(selected)
        task["static_inputs"] = [r for r in task["static_inputs"] if r["path"] != old_roster]
        task["static_inputs"].append(hashed(roster))
        # Search preflight covers the full input inventory; computation uses only this roster.
    (out / "smoke/manifest.json").write_text(json.dumps(smoke, indent=2) + "\n")
    print(out / "manifest.json", flush=True)


if __name__ == "__main__":
    main()
