"""Gate for the full candidate run: requires a PASS pilot of 2 acquisitions produced by the same worker.py (sha256), then runs worker.py full."""
from pathlib import Path
import json,hashlib,subprocess,sys
h=Path(__file__).resolve().parent
r=json.loads((h/'pilot/COMPLETE.json').read_text())
assert r['status']=='PASS' and r['mode']=='pilot' and r['acquisitions']==2
assert r['code_sha256']==hashlib.sha256((h/'worker.py').read_bytes()).hexdigest()
subprocess.run([sys.executable,str(h/'worker.py'),'full'],check=True)
