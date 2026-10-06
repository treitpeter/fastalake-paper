from pathlib import Path
import hashlib,json,subprocess,sys
HERE=Path(__file__).resolve().parent
p=HERE/'manifest.json'
assert hashlib.sha256(p.read_bytes()).hexdigest()=='6f02c2677bb6b0b91b58cf41542fc154174877b452de7d2cd928382618a5a0c8'
m=json.loads(p.read_text())
for r in m['code']:
 assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256'],r['path']
stage=sys.argv[1]
if stage in ['gate','collect']:
 command=[sys.executable,'-B',str(HERE/'finish.py'),stage]
else:
 index=int(sys.argv[2])
 if index not in m['paired_pilot_indices']:
  assert json.loads((HERE/'PILOT_GATE.json').read_text())['status']=='PASS'
 command=[sys.executable,'-B',str(HERE/'runtime/tools/run_bag_nulls.py'),str(p),stage,str(index)]
subprocess.run(command,check=True)
