"""Shared paths, provenance and small format helpers for this experiment."""

import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[1]
PREVIOUS = PROJECT / "analyses/fdr_validation_2026-09-08"
FUSION = PROJECT / "analyses/fusion_entrapment_2026-09-07"
RELEASE = PROJECT / "release_prep/github_biorxiv_2026-09-07"
RUNTIME = RELEASE / "logic_review_2026-09-10/campi_clean/_runtime/bin"
ALPHANOVO = Path("LAB_ROOT/an_window_af3ac89b")
CHECKPOINT = Path(
    "LAB_ROOT/checkpoints"
    "multislope_fragann_comp02_peptdeep_retrained_2pp/"
    "model-epoch12-step058020-peprecall0.5156.ckpt"
)
REFERENCE = PREVIOUS / "results/campi_reference_v2/REF_targets.fasta"
NOISE = PREVIOUS / (
    "vendor/FusionEntrapment/scripts/Simulation_Pipeline/3.AddNoise/mgf_noise_injector.py"
)
CONFIG = PROJECT / (
    "canonical/publication_data/arms/noclust_2026-08-07/cohorts/"
    "CAMPI_MOCK/5_sage/razor_hash-acc/S01/config.json"
)
THRESHOLDS = (0.001, 0.005, 0.01, 0.02, 0.05)
N_SPECTRA = 2048
PART_SIZE = 128
ACQUISITIONS = ("donor", "test")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def table(path, rows, fields=None):
    rows = list(rows)
    fields = fields or list(rows[0])
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_table(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


legacy = module("frozen_fusion_helpers", FUSION / "run_experiment.py")
fasta = legacy.fasta
digest_intervals = legacy.digest_intervals
canonical = legacy.canonical


def units(phase):
    if phase == "development":
        seeds = (2026091011, 2026091012)
    elif phase == "validation":
        seeds = range(2026090911, 2026090935)
    else:
        raise ValueError(phase)
    return [
        dict(id=f"u{i + 1:02d}", seed=seed, regime="baseline" if i % 2 == 0 else "stress")
        for i, seed in enumerate(seeds)
    ]


def location(phase, index):
    unit = units(phase)[int(index)]
    return unit, ROOT / "runs" / phase / unit["id"]


def guard():
    """Verify the frozen files, including binaries, before each execution stage."""
    manifest = json.loads((ROOT / "FROZEN.json").read_text())
    for path, expected in manifest["dependencies"].items():
        if sha(path) != expected:
            raise RuntimeError(f"Frozen dependency changed: {path}")
    for name, expected in manifest["source"].items():
        if sha(ROOT / name) != expected:
            raise RuntimeError(f"Frozen experiment source changed: {name}")
        committed = subprocess.check_output(["git", "show", f"HEAD:{name}"], cwd=ROOT)
        if hashlib.sha256(committed).hexdigest() != expected:
            raise RuntimeError(f"Source is not committed: {name}")
    return manifest


def binary(name):
    names = {
        "extractor": RUNTIME / "fasta_extractor",
        "parsimony": RUNTIME / "parsimony_engine",
        "sage": RUNTIME / "sage",
        "builder": FUSION / "bin/build_fused_v2",
    }
    return names[name]


def run(args, output, cwd=None, model=False):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    threads = str(int(os.environ.get("SLURM_CPUS_PER_TASK", "4")))
    env = dict(
        os.environ,
        WANDB_MODE="offline",
        OMP_NUM_THREADS=threads if model else "1",
        MKL_NUM_THREADS=threads if model else "1",
        OPENBLAS_NUM_THREADS="1",
        RAYON_NUM_THREADS=threads,
    )
    record = dict(
        argv=list(map(str, args)),
        start=time.time(),
        cwd=str(cwd) if cwd else None,
        slurm_job=os.getenv("SLURM_JOB_ID"),
    )
    save(output / "COMMAND.json", record)
    with (output / "stdout.log").open("w") as stdout, (output / "stderr.log").open("w") as stderr:
        result = subprocess.run(record["argv"], cwd=cwd, env=env, stdout=stdout, stderr=stderr)
    record.update(returncode=result.returncode, seconds=time.time() - record["start"])
    save(output / "COMMAND.json", record)
    if result.returncode:
        raise RuntimeError(f"Command failed; evidence retained in {output}")
    return record


def execute_main(function):
    """Keep stage command lines uniform: phase, zero-based unit, optional arguments."""
    guard()
    function(sys.argv[1], int(sys.argv[2]), *sys.argv[3:])
