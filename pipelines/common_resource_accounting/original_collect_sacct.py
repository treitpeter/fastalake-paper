#!/usr/bin/env python3
"""Collect SLURM accounting for the noclust_2026-08-07 chain, all six cohorts.

Reads each cohort's .jobids (HSTOOL: stages 0-3 from .jobids.pre_restart, stages 4+ from
.jobids.restart, because parsimony OOMed at 8G and the tail was resubmitted at 16G), pulls
sacct for every job and array task, and writes:

  data/compute_jobs.tsv     one row per job / array task
  data/compute_summary.tsv  one row per cohort x stage

MaxRSS comes from the .batch step. SLURM samples RSS every 30 s (JobAcctGatherFrequency=30),
so a task shorter than ~30 s carries only the batch shell (~2.5 MB); anything below 10 MB is written as empty, never as 0.
Elapsed and CPU time are job-level. cpu_alloc_h = NCPUS x elapsed; cpu_used_h = TotalCPU.
"""
import csv, re, subprocess, statistics as st
from collections import defaultdict
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
ARM = Path("LAB_ROOT/canonical/publication_data/arms/noclust_2026-08-07/cohorts")
COHORTS = ["SIHUMI", "CAMPI_MOCK", "SALIVA", "HSTOOL", "PREDICT", "CAPSCAN"]
STAGE_RE = re.compile(r"^(\d{2}_\S+)\s+(\S+)")

def read_jobids(path, keep=None):
    out = []
    for line in open(path):
        m = STAGE_RE.match(line)
        if not m or m.group(2) == "-":
            continue
        if keep and not keep(m.group(1)):
            continue
        out.append((m.group(1), m.group(2)))
    return out

jobs = {}
for c in COHORTS:
    if c == "HSTOOL":
        lst = read_jobids(ARM / c / ".jobids.pre_restart", keep=lambda s: s[:2] in ("01", "02", "05"))
        lst += read_jobids(ARM / c / ".jobids.restart")
    else:
        lst = read_jobids(ARM / c / ".jobids")
    for stage, jid in lst:
        jobs[jid] = (c, stage)

def gb(s):
    if not s:
        return None
    m = re.match(r"([\d.]+)([KMGT])?", s)
    if not m.group(2):
        return float(m.group(1)) / 1024**3
    return float(m.group(1)) * {"K": 1 / 1024**2, "M": 1 / 1024, "G": 1, "T": 1024}[m.group(2)]

def hms(s):
    d, _, rest = s.rpartition("-")  # [D-]HH:MM:SS or MM:SS.sss
    parts = [float(x) for x in rest.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0.0)
    h, m, sec = parts
    return (int(d) if d else 0) * 86400 + h * 3600 + m * 60 + sec

fields = "JobID,JobName%40,State,ElapsedRaw,ReqMem,MaxRSS,NCPUS,TotalCPU,CPUTimeRAW,Submit,Start,End,NodeList,Partition"
raw = subprocess.run(["sacct", "-j", ",".join(jobs), "--format", fields, "-P"], capture_output=True, text=True, check=True).stdout
(DATA / "sacct_raw.tsv").write_text(raw)

rows = {}
for r in csv.DictReader(raw.splitlines(), delimiter="|"):
    jid = r["JobID"]; base = jid.split(".")[0]; job = base.split("_")[0]
    if job not in jobs:
        continue
    d = rows.setdefault(base, {"cohort": jobs[job][0], "stage": jobs[job][1], "jobid": base, "maxrss_gb": None})
    if "." not in jid:
        d.update(name=r["JobName"], state=r["State"], elapsed_s=int(r["ElapsedRaw"]), reqmem_gb=gb(r["ReqMem"]),
                 ncpus=int(r["NCPUS"]), cpu_alloc_s=int(r["CPUTimeRAW"]), cpu_used_s=hms(r["TotalCPU"]),
                 submit=r["Submit"], start=r["Start"], end=r["End"], node=r["NodeList"], partition=r["Partition"])
    elif jid.endswith(".batch"):
        v = gb(r["MaxRSS"])
        # A task shorter than the 30 s sample interval reports only the batch shell (~2.5 MB);
        # the stage process was never sampled. Below 10 MB is "no sample", not "no memory".
        d["maxrss_gb"] = v if v and v >= 0.01 else None

def key(d):
    j = d["jobid"].split("_")
    return (COHORTS.index(d["cohort"]), d["stage"], int(j[0]), int(j[1]) if len(j) > 1 else -1)

out = sorted(rows.values(), key=key)
cols = ["cohort", "stage", "jobid", "name", "state", "elapsed_s", "cpu_alloc_s", "cpu_used_s", "ncpus", "reqmem_gb",
        "maxrss_gb", "submit", "start", "end", "node", "partition"]
with open(DATA / "compute_jobs.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(out)

agg = defaultdict(list)
for d in out:
    agg[(d["cohort"], d["stage"])].append(d)
scols = ["cohort", "stage", "n_tasks", "n_not_completed", "elapsed_med_min", "elapsed_max_min", "n_rss_sampled",
         "maxrss_med_gb", "maxrss_max_gb", "reqmem_gb", "ncpus", "cpu_alloc_h", "cpu_used_h", "stage_wall_h"]
with open(DATA / "compute_summary.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(scols)
    for (c, s), L in sorted(agg.items(), key=lambda kv: (COHORTS.index(kv[0][0]), kv[0][1])):
        e = [d["elapsed_s"] / 60 for d in L]
        m = [d["maxrss_gb"] for d in L if d["maxrss_gb"] is not None]
        wall = (max(datetime.fromisoformat(d["end"]) for d in L) - min(datetime.fromisoformat(d["start"]) for d in L)).total_seconds() / 3600
        w.writerow([c, s, len(L), sum(d["state"] != "COMPLETED" for d in L), f"{st.median(e):.2f}", f"{max(e):.2f}", len(m),
                    f"{st.median(m):.2f}" if m else "", f"{max(m):.2f}" if m else "", f"{L[0]['reqmem_gb']:.0f}", L[0]["ncpus"],
                    f"{sum(d['cpu_alloc_s'] for d in L) / 3600:.2f}", f"{sum(d['cpu_used_s'] for d in L) / 3600:.2f}", f"{wall:.2f}"])
print(f"{len(out)} jobs/tasks from {len(jobs)} job ids; states:", {s: sum(d['state'] == s for d in out) for s in {d['state'] for d in out}})
