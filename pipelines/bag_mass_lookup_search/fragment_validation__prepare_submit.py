"""Freeze and run a CAMPI pilot followed by a dependency-gated five-case array."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def record(path):
    return dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def main():
    assert not (HERE / "SUBMITTED.json").exists()
    queue = subprocess.check_output(["squeue", "-u", "treit", "-h", "-o", "%i|%j|%T|%R"], text=True)
    (HERE / "QUEUE_BEFORE.txt").write_text(queue)
    assert "FL_bag_fragments" not in queue
    result = subprocess.run([sys.executable, "-B", str(HERE / "ions.py")], capture_output=True, text=True)
    (HERE / "METHOD_CHECK.log").write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stderr
    cases = json.loads((HERE.parent / "INPUTS.json").read_text())["cases"]
    assert len(cases) == 6 and cases[0]["cohort"] == "CAMPI"
    previous = {r["path"]: r["sha256"] for r in json.loads((HERE.parent / "results/SOURCE_IDENTITIES.json").read_text())}
    for case in cases:
        root = Path(case["directory"])
        paths = [root / "task.json", root / "selection.json"]
        for arm in ("exact", "bags", "mass", "combined"):
            paths += [root / (arm + ".fasta"), root / "searches" / arm / "config.json", root / "searches" / arm / "results.sage.tsv"]
        case["frozen"] = [record(path) for path in paths]
        for entry in case["frozen"]:
            if entry["path"] in previous:
                assert entry["sha256"] == previous[entry["path"]], entry["path"]
    sage = json.loads((HERE.parent.parent / "attempt_02/manifest.json").read_text())["sage"]
    assert record(Path(sage["path"]))["sha256"] == sage["sha256"]
    code = [HERE / n for n in ("run.py", "ions.py", "prepare_submit.py", "PROTOCOL.md")]
    code += [HERE.parent / n for n in ("analyze.py", "INPUTS.json", "results/PEPTIDE_MECHANISMS.tsv")]
    dump(HERE / "MANIFEST.json", dict(cases=cases, sage=sage,
         code_sha256={str(path): record(path)["sha256"] for path in code}, utc=datetime.now(timezone.utc).isoformat()))
    script = HERE / "validate.sbatch"
    script.write_text("#!/bin/bash\nset -euo pipefail\nexport PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\n"
                      'exec LAB_ROOT/python -B -u run.py "${SLURM_ARRAY_TASK_ID:-0}"\n')
    (HERE / "logs").mkdir(exist_ok=True)
    commands, jobs = [], []
    for stage in ("pilot", "remaining"):
        command = ["sbatch", "--parsable", "--job-name=FL_bag_fragments_" + stage,
                   "--partition=p.hpcl8,p.hpcl92,p.hpcl93", "--nice=10000", "--cpus-per-task=4",
                   "--mem=16G", "--time=2:00:00", f"--chdir={HERE}",
                   f"--output={HERE}/logs/%A_%a.out", f"--error={HERE}/logs/%A_%a.err"]
        if stage == "remaining":
            command += ["--array=1-5%2", "--dependency=afterok:" + jobs[0], "--kill-on-invalid-dep=yes"]
        command.append(str(script))
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        job = result.stdout.strip().split(";")[0]
        assert job.isdigit()
        commands.append(command)
        jobs.append(job)
        dump(HERE / "SUBMITTED.json", dict(jobs=jobs, commands=commands, utc=datetime.now(timezone.utc).isoformat()))
    print(json.dumps(dict(jobs=jobs, status="SUBMITTED_DEPENDENCY_GATED")))


if __name__ == "__main__":
    main()
