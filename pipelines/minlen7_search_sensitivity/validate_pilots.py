from pathlib import Path
import json
O=Path(__file__).resolve().parent.parent
ids=json.loads((O/'pilots.json').read_text())
r=[json.loads((O/'receipts/pilot_'/('%03d.json'%i)).read_text()) for i in ids]
assert all(x['status']=='PASS' for x in r)
assert {x['kind'] for x in r}=={'common','fusion','joint'}
(O/'PILOT_VALIDATION.json').write_text(json.dumps(dict(status='PASS',indices=ids,scope='Search/config/length and corrected joint LFQ audit engineering gate; not a scientific result'),indent=2)+'\n')
print('PASS three wave2 pilot routes')
