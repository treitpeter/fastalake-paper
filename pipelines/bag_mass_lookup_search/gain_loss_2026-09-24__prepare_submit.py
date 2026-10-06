"""Freeze every currently completed acquisition, then submit one attribution job."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
BASE = HERE.parent


def record(path):
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def main():
    assert not (HERE / "SUBMITTED.json").exists()
    queue = subprocess.check_output(["squeue", "-u", "treit", "-h", "-o", "%i|%j|%T|%R"], text=True)
    (HERE / "QUEUE_BEFORE.txt").write_text(queue)
    assert "FL_bag_gainloss" not in queue
    cases, frozen = [], []
    for folder in (BASE / "attempt_02/outputs", BASE / "start_now_2026-09-23/environment_retry/retry_544_outputs"):
        for path in sorted(folder.glob("*/*/SEARCHES_COMPLETE.json")):
            case = path.parent
            task = json.loads((case / "task.json").read_text())
            candidate = json.loads((case / "CANDIDATES_COMPLETE.json").read_text())
            assert candidate["baseline_reproduced"]
            cases.append(dict(cohort=task["cohort"], sample=task["sample"], directory=str(case),
                              null_directory=str(BASE / "paired_nulls_2026-09-23/outputs" / task["cohort"] / task["sample"])))
            for name in ("task.json", "CANDIDATES_COMPLETE.json", "SEARCHES_COMPLETE.json", "selection.json", "baseline_check.json"):
                frozen.append(record(case / name))
    assert len({(c["cohort"], c["sample"]) for c in cases}) == len(cases)
    assert {c["cohort"] for c in cases} == {"CAMPI", "HSTOOL", "PREDICT", "CAPSCAN"}
    for name in ("analyze.py", "check_methods.py", "PROTOCOL.md", "prepare_submit.py"):
        frozen.append(record(HERE / name))
    root = BASE.parents[2]
    for name in ("V7_COORDINATION/MM_COMMENT_CROSSWALK.md", "MM_review_analysis/SECTION_EXECUTION.md", "MM_review_analysis/STATUS.md"):
        frozen.append(record(root / name))
    dump(HERE / "INPUTS.json", dict(utc=datetime.now(timezone.utc).isoformat(), cases=cases, frozen=frozen,
                                    scope="All completed validated native search suites at this timestamp; no outcome selection"))
    test = subprocess.run([sys.executable, "-B", str(HERE / "check_methods.py")], capture_output=True, text=True)
    (HERE / "METHOD_CHECK.log").write_text(test.stdout + test.stderr)
    assert test.returncode == 0, test.stderr
    script = HERE / "analysis.sbatch"
    script.write_text("#!/bin/bash\nset -euo pipefail\nexport PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\n"
                      "exec LAB_ROOT/python -B -u analyze.py\n")
    command = ["sbatch", "--parsable", "--job-name=FL_bag_gainloss", "--partition=p.hpcl8,p.hpcl92,p.hpcl93", "--nice=10000",
               "--cpus-per-task=2", "--mem=16G", "--time=1:00:00", f"--chdir={HERE}",
               f"--output={HERE}/%x_%j.out", f"--error={HERE}/%x_%j.err", str(script)]
    submitted = subprocess.run(command, check=True, capture_output=True, text=True)
    job = submitted.stdout.strip().split(";")[0]
    assert job.isdigit()
    dump(HERE / "SUBMITTED.json", dict(job_id=job, command=command, stdout=submitted.stdout, stderr=submitted.stderr,
                                      utc=datetime.now(timezone.utc).isoformat(), cases=len(cases)))
    print(json.dumps(dict(job_id=job, cases=len(cases), status="SUBMITTED")))


if __name__ == "__main__":
    main()
