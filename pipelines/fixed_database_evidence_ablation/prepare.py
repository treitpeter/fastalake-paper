"""fixed_database_evidence_ablation: prepare. Moved from the original analysis as it ran, apart from path scrubbing; see README.md."""
from pathlib import Path
import json,hashlib,gzip,csv
b=Path(__file__).resolve().parents[1];j=json.loads((b/'inputs/PREPARE.json').read_text());p=b/'inputs/MG_accepted_peptides.tsv.gz'
with gzip.open(p,'wt') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['sample','canonical_peptide'])
 for t in j['tasks']:
  r=t['sets']['MG'];data=Path(r['path']).read_bytes();assert hashlib.sha256(data).hexdigest()==r['sha256'];w.writerows((t['sample'],s) for s in sorted(set(data.decode().splitlines())))
sources=[Path(j[k]) for k in ['features','original_long','annotations']]+[p];t={k:j[k] for k in ['features','original_long','annotations']};t.update(cohort='PREDICT',samples=[x['sample'] for x in j['tasks']],mg_peptides=str(p),sources=[dict(path=str(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in sources]);(b/'manifest.json').write_text(json.dumps(dict(tasks=[t]),indent=2)+'\n');print('PREPARE PASS',len(t['samples']),flush=True)
