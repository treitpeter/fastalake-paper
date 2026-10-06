#!/usr/bin/env python3
"""Matched searches after clustering, with an archived-baseline reproduction gate."""
import csv,gzip,json,os,re,subprocess,sys,time
from pathlib import Path
from h01_eligibility import sha
def peptides(path):
    ans=set()
    with Path(path).open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['label']=='1' and 0<=float(r['peptide_q'])<=.01:
                p=re.sub(r'\[[^\]]*\]','',r['peptide']).replace('I','L')
                if 5<=len(p)<=50:ans.add(p)
    return ans
base=Path(sys.argv[1]);index=int(sys.argv[2]);m=json.loads((base/'manifest.json').read_text());t=m['tasks'][index]
out=base/'results'/t['sample'];out.mkdir(parents=True,exist_ok=False)
assert sha(m['H04_manifest']['path'])==m['H04_manifest']['sha256']
for r in m['binaries']+t['sources']:
    p=Path(r['path']);st=p.stat();assert st.st_size==r['bytes'] and st.st_mtime_ns==r['mtime_ns'],str(p)
    if 'sha256' in r:assert sha(p)==r['sha256'],str(p)
producer=Path(t['H04'])/'COMPLETE.json';d=json.loads(producer.read_text())
assert d['status']=='PASS' and d['baseline_sequence_set_reproduced']
for r in d['outputs']:assert sha(r['path'])==r['sha256']
sources=[dict(path=r['path'],sha256=sha(r['path'])) for r in t['sources']]
sources.append(dict(path=str(producer),sha256=sha(producer)))
cfg0=json.loads(Path(t['config']).read_text());assert cfg0['version']=='0.14.6'
assert subprocess.check_output([m['sage'],'--version'],text=True).strip()=='sage 0.14.6'
expected=peptides(t['archived_search']);reference=set(Path(t['reference_peptides']).read_text().splitlines())
membership={}
with Path(t['members']).open() as f:
    for r in csv.DictReader(f,delimiter='\t'):
        if r['n_members']=='1':membership[r['peptide']]=r['members']
arms={'unclustered':Path(cfg0['database']['fasta']),'cluster_097':Path(t['H04'])/'cluster_097_parsimony'/(t['sample']+'_razor.fasta'),'cluster_099':Path(t['H04'])/'cluster_099_parsimony'/(t['sample']+'_razor.fasta')}
sets={};counts=[];by_member=[];commands=[]
for arm,fasta in arms.items():
    folder=out/arm;folder.mkdir();cfg=json.loads(json.dumps(cfg0));cfg.pop('output_paths',None);cfg['database']['fasta']=str(fasta)
    config=folder/'config.json';config.write_text(json.dumps(cfg,indent=2)+'\n')
    sources.append(dict(path=str(fasta),sha256=sha(fasta)))
    argv=[m['sage'],str(config),'-o',str(folder),'--write-pin','--disable-telemetry-i-dont-want-to-improve-sage']
    started=time.time()
    with (folder/'search.log').open('w') as f:r=subprocess.run(argv,env=dict(os.environ,RAYON_NUM_THREADS='2'),stdout=f,stderr=subprocess.STDOUT)
    commands.append(dict(arm=arm,argv=argv,returncode=r.returncode,seconds=time.time()-started));(out/'COMMANDS.json').write_text(json.dumps(commands,indent=2)+'\n')
    assert r.returncode==0,arm
    p=peptides(folder/'results.sage.tsv');sets[arm]=p
    (folder/'accepted_canonical_5_50.txt').write_text(''.join(x+'\n' for x in sorted(p)))
    if arm=='unclustered':
        gate=dict(status='PASS' if p==expected else 'FAIL',reproduced=len(p),archived=len(expected),gained=len(p-expected),lost=len(expected-p),jaccard=len(p&expected)/max(1,len(p|expected)))
        (out/'BASELINE_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
        assert p==expected,'Baseline identity drift: do not interpret or expand the intervention'
    counts.append(dict(sample=t['sample'],arm=arm,accepted=len(p),reference_denominator=len(reference),reference_recovered=len(p&reference),lost_vs_unclustered=len(sets['unclustered']-p),gained_vs_unclustered=len(p-sets['unclustered'])))
    for member in sorted(set(membership.values())):
        specific={p for p,v in membership.items() if v==member}&reference
        by_member.append(dict(sample=t['sample'],arm=arm,member=member,member_specific_REF_denominator=len(specific),recovered=len(specific&p),lost_vs_unclustered=len((specific&sets['unclustered'])-p),gained_vs_unclustered=len((specific&p)-sets['unclustered'])))
for name,rows in [('summary.tsv',counts),('member_specific_recovery.tsv',by_member)]:
    with (out/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
with gzip.open(out/'peptide_fates.tsv.gz','wt') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['peptide','REF','member_specific',*arms])
    for p in sorted(set.union(reference,*sets.values())):w.writerow([p,int(p in reference),membership.get(p,''),*[int(p in sets[a]) for a in arms]])
(out/'SOURCE_HASHES.json').write_text(json.dumps(sources,indent=2)+'\n')
(out/'COMPLETE.json').write_text(json.dumps(dict(status='PASS',sample=t['sample'],searches=3,baseline_exact=True,scope='Same-spectrum/source-config searches after candidate clustering; reference member-specific evidence is not strain truth or equal-error sensitivity. No end-to-end MetaNovo comparison.',outputs=[dict(path=str(p),sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file()]),indent=2)+'\n')
