#!/usr/bin/env python3
"""Partition measured LFQ features, fixing the existing reporting dictionary."""
import csv,gzip,json,math,re,sys
from collections import defaultdict,Counter
from pathlib import Path
import pandas as pd
from hashlib import sha256
def sha(p):
    h=sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(8<<20),b""):h.update(b)
    return h.hexdigest()
def rows(p):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:yield from csv.DictReader(f,delimiter='\t')
base=Path(sys.argv[1]);stage=sys.argv[2];index=int(sys.argv[3]);m=json.loads((base/'manifest.json').read_text());t=m['tasks'][index]
out=base/(stage+'_'+t['cohort']);out.mkdir(exist_ok=False)
hashes=[]
for r in t['sources']:
    p=Path(r['path']);st=p.stat();assert st.st_size==r['bytes'] and st.st_mtime_ns==r['mtime_ns']
    hashes.append(dict(path=str(p),sha256=sha(p)))
(out/'SOURCE_HASHES.json').write_text(json.dumps(hashes,indent=2)+'\n')
samples=set(t['samples'] if stage=='full' else t['samples'][:1]);mg=defaultdict(set)
for r in rows(t['mg_peptides']):
    if r['sample'] in samples:mg[r['sample']].add(r['canonical_peptide'])
annotations={r['study_group']:set(re.findall(r'K\d{5}',r['KEGG_ko'])) for r in rows(t['annotations'])}
groups=defaultdict(lambda:[0.0,0.0]);features=defaultdict(Counter)
for r in rows(t['features']):
    s=r['sample']
    if s not in samples:continue
    value=float(r['intensity']);assert math.isfinite(value) and value>0
    kind=0 if r['canonical_peptide'] in mg[s] else 1
    groups[s,r['study_group']][kind]+=value
    features[s]['shared' if kind==0 else 'FL_only']+=1
    features[s]['outside_local_representative']+=r['representative_in_local_fasta']=='0'
assert set(features)==samples
# Independently reconstruct every original group, including mixed groups.
seen=set()
for r in rows(t['original_long']):
    if r['sample'] not in samples:continue
    key=r['sample'],r['study_group'];assert key not in seen;seen.add(key)
    total=sum(groups[key]);original=float(r['intensity'])
    assert abs(total-original)<=1e-8*max(original,1),(key,total,original)
assert seen==set(groups)
ko=defaultdict(lambda:[0.0,0.0]);totals=defaultdict(lambda:[0.0,0.0]);missing=defaultdict(lambda:[0.0,0.0])
with gzip.open(out/'group_ablation.tsv.gz','wt') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['sample','study_group','shared_feature_intensity','FL_only_feature_intensity','all_feature_intensity','KO_terms'])
    for (s,g),values in sorted(groups.items()):
        assert g in annotations,g
        terms=annotations[g];w.writerow([s,g,*values,sum(values),';'.join(sorted(terms))])
        for i,v in enumerate(values):totals[s][i]+=v
        if not terms:
            for i,v in enumerate(values):missing[s][i]+=v
        for k in terms:
            for i,v in enumerate(values):ko[s,k][i]+=v/len(terms)
with gzip.open(out/'KO_ablation_equal_split.tsv.gz','wt') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['sample','KO','shared_feature_intensity','FL_only_feature_intensity','all_feature_intensity'])
    for (s,k),v in sorted(ko.items()):w.writerow([s,k,*v,sum(v)])
summary=[]
for s in sorted(samples):
    v=[x for (sample,k),x in ko.items() if sample==s];a=sum(x[0] for x in v);b=sum(x[1] for x in v)
    for i,total in enumerate([a,b]):assert abs(total-(totals[s][i]-missing[s][i]))<=1e-8*max(totals[s][i],1)
    # Bray–Curtis between normalized nonnegative KO profiles; no pseudo-counts.
    distance=sum(abs((x[0]/a if a else 0)-(sum(x)/(a+b) if a+b else 0)) for x in v)/2 if a>0 and a+b>0 else None
    summary.append(dict(sample=s,all_LFQ_intensity=sum(totals[s]),shared_LFQ_intensity=totals[s][0],FL_only_LFQ_intensity=totals[s][1],fraction_intensity_removed=totals[s][1]/sum(totals[s]),all_KOs=len(v),shared_KOs=sum(x[0]>0 for x in v),FL_only_supported_KOs=sum(x[0]==0 and x[1]>0 for x in v),normalized_KO_profile_BrayCurtis_all_vs_shared=distance,**features[s]))
pd.DataFrame(summary).to_csv(out/'acquisition_summary.tsv',sep='\t',index=False)
(out/'REPORT.md').write_text('# Peptide-evidence ablation\n\nEvery original study-group intensity was reconstructed before partitioning. Shared versus FL-only classes use canonical identities accepted in the same acquisition\'s MG comparison. Mixed groups are now divided using their measured LFQ features, not an arbitrary split of group totals. KO projection divides each group equally among its fixed annotations to conserve annotated intensity. The dictionary, features and annotations are held fixed: this is a reporting/evidence ablation, not a new identification search, re-inference of group representatives or a causal test of the clinical conclusions. The all-versus-shared KO-profile distance quantifies profile dependence on added features. Inferential clinical/own-stool analyses require the same roster and original models; they are not replaced by this descriptive distance. Preserve unmatched annotation intensity and representative-locality flags.\n')
(out/'COMPLETE.json').write_text(json.dumps(dict(status='PASS',cohort=t['cohort'],stage=stage,samples=len(samples),all_original_group_intensities_reproduced=True,outputs=[dict(path=str(p),sha256=sha(p)) for p in sorted(out.iterdir())]),indent=2)+'\n')
