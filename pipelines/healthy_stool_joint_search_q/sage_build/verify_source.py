"""Verify a source tree against the preserved producing source manifest."""
from pathlib import Path
import hashlib,json,sys
manifest=json.loads(Path(__file__).with_name('PRODUCING_SOURCE_MANIFEST.json').read_text())
root=Path(sys.argv[1])
bad=[]
for rel,want in manifest['source_files'].items():
    p=root/rel
    actual=hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None
    if actual!=want:bad.append(dict(path=rel,expected=want,actual=actual))
print(json.dumps(dict(status='FAIL'if bad else'PASS',files=len(manifest['source_files']),mismatches=bad),indent=2))
raise SystemExit(1 if bad else 0)
