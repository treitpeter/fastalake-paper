from pathlib import Path
import json,hashlib,sys,os,subprocess,time
from full_route_entrapment import summarize_search
H=Path(__file__).resolve().parent;M=json.loads((H/'MANIFEST.json').read_text());r=M['tasks'][int(sys.argv[1])];sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(M['sage'])==M['sage_sha256'] and sha(r['config'])==r['config_sha256'] and sha(r['database'])==r['database_sha256']
code=json.loads((H/'CODE_SHA256.json').read_text());assert all(sha(H/k)==v for k,v in code.items())
o=Path(r['config']).parent;assert not (o/'COMPLETE.json').exists();t=time.time()
with (o/'stdout.log').open('w') as so,(o/'stderr.log').open('w') as se:subprocess.run([M['sage'],r['config'],'--write-pin'],stdout=so,stderr=se,check=True)
summary=summarize_search(o/'results.sage.tsv',json.loads((H/'control_classes.json').read_text()),o)
(o/'COMPLETE.json').write_text(json.dumps({'task':r,'seconds':time.time()-t,'job':os.getenv('SLURM_JOB_ID'),'array_job':os.getenv('SLURM_ARRAY_JOB_ID'),'summary':summary,'code_sha256':code},indent=2));print('Complete',r['task'],r['arm'],r['sample'])
