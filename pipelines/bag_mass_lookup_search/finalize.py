from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import subprocess

here = Path(__file__).resolve().parent
validation = json.loads((here / "REPORT_VALIDATION.json").read_text())
assert validation["status"] == "PASS_ARITHMETIC_SOURCE_HASHES_PDF_STRUCTURE"
pdf = here / "Bags_Flips_Verification_Peter_Philipp_2026-09-24.pdf"
assert hashlib.sha256(pdf.read_bytes()).hexdigest() == validation["pdf_sha256"]
fragment = json.loads((here / "FRAGMENT_SUMMARY.json").read_text())
assert fragment["cases"] == 6 and fragment["replay_q_differences"] == 0
scheduler = subprocess.check_output(["sacct", "-j", "6000306,6000307,6000119", "--format=JobID,State,ExitCode,Elapsed,MaxRSS", "-P"], text=True)
(here / "FINAL_SACCT.txt").write_text(scheduler)
states = {line.split("|")[0]: line.split("|")[1:3] for line in scheduler.splitlines()[1:] if line}
for job in ["6000306"] + ["6000307_" + str(i) for i in range(1,6)]:
    assert states[job] == ["COMPLETED", "0:0"], (job, states[job])
result = dict(utc=datetime.now(timezone.utc).isoformat(), status="REPORT_AND_SIX_CASE_FRAGMENT_REVIEW_COMPLETE",
              pdf=str(pdf), pdf_sha256=validation["pdf_sha256"], pages=17, fragment=fragment,
              pending="Three full-reservoir Rust equivalence checks and ongoing 710-acquisition cohort work; downstream biological and FDR endpoints remain open.",
              github_actions="No commit, merge or push in this review; software worktree stays on bags at 1f0fc0e.")
(here / "DELIVERY_READY.json").write_text(json.dumps(result, indent=2)+"\n")
print(json.dumps(result))
