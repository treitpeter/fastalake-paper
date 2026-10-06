"""Wrapper that maps a global Slurm array index onto the PREDICT task list and runs search.py for that task (used by the array-index repair of job 5929446)."""
import json,runpy,sys
from pathlib import Path
r=Path(sys.argv[1]);sys.path.insert(0,str(r/"frozen"))
from lib import verify
m=json.loads((r/"inputs/predict_submission_manifest.json").read_text())
for x in m["code"]+m["inputs"]:verify(x)
i=int(sys.argv[2]);assert 332<=i<=636
sys.argv[2]=str(i+637)
runpy.run_path(str(r/"frozen/search.py"),run_name="__main__")
