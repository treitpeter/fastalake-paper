#!/usr/bin/env python3
"""Predeclared measured-spectrum examples for newly recovered H boundaries."""
import csv,itertools,json,os,re,sys
from collections import defaultdict
from pathlib import Path
import ahocorasick
import numpy as np
import pandas as pd
from pyteomics import mzml
RUN=Path(sys.argv[1]).resolve();sys.path.insert(0,str(RUN/'frozen'))
from lib import meta,verify,dump
from peptide import canonical,bare
from context_helpers import fasta
os.environ['NUMBA_CACHE_DIR']=str(RUN/'results/numba_cache')
from fragment_helpers import parse,match,MASS,PROTON
m=json.loads((RUN/'inputs/manifest_revision2.json').read_text());index=int(sys.argv[2]);pilot=index==0
cohort=['PILOT','CAPSCAN','PREDICT'][index];out=RUN/'results_revision2'/('pilot' if pilot else cohort);out.mkdir(exist_ok=True)
assert not (out/'COMPLETE.json').exists()
for x in m['code']+m['inputs']:verify(x)
tasks=m['tasks'][:4] if pilot else [t for t in m['tasks'] if t['cohort']==cohort]
byname={t['sample']:t for t in tasks};candidates={};counts=[];sources=[]
for t in tasks:
    p=Path(t['receipt']);c=json.loads(p.read_text());assert c['status']=='PASS';sources.append(meta(p))
    x=next(x for x in c['outputs'] if Path(x['path']).name=='site_evidence.tsv.gz');verify(x)
    try:f=pd.read_csv(x['path'],sep='\t',keep_default_na=False)
    except pd.errors.EmptyDataError:f=pd.DataFrame()
    verify(t['baseline_peptides']);baseline=set(Path(t['baseline_peptides']['path']).read_text().splitlines())
    compartment='stool' if t['cohort']=='PREDICT' or t['design']['region']=='stool' else 'capsule'
    for origin,side in itertools.product(['Human','Microbial'],['N','C']):
        key=(t['cohort'],compartment,origin,side);pp=set()
        if len(f):
            z=f[f.route.eq('sage_unspecific_bothq')&f.policy.eq('unambiguous')&f.origin.eq(origin)]
            z=z[z.motifs.map(lambda s:bool(s) and all(w[3]=='H' for w in s.split(';')))]
            for cell in z[side+'_peptides']:pp.update(p for p in str(cell).split(';') if p)
        pp-=baseline;candidates[key,t['sample']]=pp
        counts.append(dict(cohort=t['cohort'],compartment=compartment,origin=origin,side=side,sample=t['sample'],candidates=len(pp)))
pd.DataFrame(counts).to_csv(out/'complete_candidate_counts.tsv',sep='\t',index=False)
psmcache={};selected=[];strata=[];eligible_rows=[]
def psms(t):
    if t['sample'] not in psmcache:
        receipt=json.loads(Path(t['unspecific_receipt']).read_text());assert receipt['status']=='PASS';sources.append(meta(t['unspecific_receipt']))
        oo={Path(x['path']).name:x for x in receipt['outputs']}
        verify(oo['results.sage.tsv']);verify(oo['results.json']);sources.extend([oo['results.sage.tsv'],oo['results.json']])
        f=pd.read_csv(oo['results.sage.tsv']['path'],sep='\t');f=f[(f['rank']==1)&(f.label==1)&(f.peptide_q<=.01)&(f.spectrum_q<=.01)].copy()
        f['canonical_peptide']=f.peptide.map(canonical);f=f[f.canonical_peptide.str.len().between(9,50)]
        psmcache[t['sample']]=(f,oo)
    return psmcache[t['sample']]
for key in sorted({k for k,s in candidates}):
    choices=sorted([(len(v),sample) for (k,sample),v in candidates.items() if k==key and v])
    if not choices:strata.append(dict(zip(['cohort','compartment','origin','side'],key),status='NO_ELIGIBLE_SITE'));continue
    midpoint=len(choices)//2;order=sorted(range(len(choices)),key=lambda i:(abs(i-midpoint),i));chosen=None
    for rank in order:
        n,sample=choices[rank];t=byname[sample];f,oo=psms(t);g=f[f.canonical_peptide.isin(candidates[key,sample])].copy()
        g=g.sort_values(['peptide_q','hyperscore','scannr','peptide'],ascending=[True,False,True,True]).drop_duplicates('canonical_peptide')
        if not len(g):continue
        g=g.sort_values(['matched_intensity_pct','canonical_peptide','scannr']).reset_index(drop=True);selected_rank=len(g)//2
        row=g.iloc[selected_rank].to_dict();case='HIS_'+str(len(selected)+1).zfill(3)
        row.update(case_id=case,cohort=key[0],compartment=key[1],origin=key[2],side=key[3],sample=sample,
                   acquisition_rank=rank,acquisition_median_rank=midpoint,eligible_peptides=len(g),selected_rank=selected_rank)
        selected.append(row);chosen=case
        for i,record in g.iterrows():eligible_rows.append(dict(case_id=case,sample=sample,selection_rank=i,selected=i==selected_rank,**record.to_dict()))
        break
    strata.append(dict(zip(['cohort','compartment','origin','side'],key),status='SELECTED' if chosen else 'NO_BOTHQ_PSM',case_id=chosen))
pd.DataFrame(strata).to_csv(out/'all_strata.tsv',sep='\t',index=False)
pd.DataFrame(selected).to_csv(out/'selected_cases.tsv',sep='\t',index=False)
pd.DataFrame(eligible_rows).to_csv(out/'eligible_case_selection_psms.tsv.gz',sep='\t',index=False)
summary=[];comparators=[]
for row in selected:
    t=byname[row['sample']];f,oo=psms(t);case=out/row['case_id'];case.mkdir(exist_ok=True)
    verify(t['raw']);verify(t['fasta']);cfg=json.loads(Path(oo['results.json']['path']).read_text())
    assert cfg['max_peaks']==150 and cfg['max_fragment_charge']==2 and cfg['fragment_tol']=={'ppm':[-20,20]} and not cfg['deisotope']
    assert cfg['database']['ion_kinds']==['b','y'] and cfg['database']['fragment_min_mz']==200 and cfg['database']['fragment_max_mz']==2000
    assert cfg['database']['enzyme']['cleave_at']=='' and cfg['database']['static_mods']=={'C':57.0215} and cfg['database']['variable_mods']=={'M':[15.9949]}
    assert cfg['database']['max_variable_mods']==2 and cfg['precursor_tol']=={'ppm':[-20,20]}
    with mzml.MzML(t['raw']['path'],use_index=True,decode_binary=True) as reader:sp=reader.get_by_id(row['scannr'])
    assert sp['ms level']==2 and sp['id']==row['scannr']
    mz=np.asarray(sp['m/z array'],np.float32);intensity=np.asarray(sp['intensity array'],np.float32)
    assert np.isfinite(mz).all() and np.isfinite(intensity).all() and (intensity>=0).all()
    pd.DataFrame({'mz':mz,'intensity':intensity}).to_csv(case/'observed_peaks.tsv.gz',sep='\t',index=False)
    k=(mz>=200)&(mz<=2000);ma=mz[k];it=intensity[k];ix=np.argsort(it,kind='stable')[-150:];ma=ma[ix];it=it[ix];ix=np.argsort(ma);ma=ma[ix];it=it[ix];masses=np.asarray(ma-PROTON,np.float32)
    pd.DataFrame({'mz':ma,'intensity':it}).to_csv(case/'processed_peaks.tsv',sep='\t',index=False)
    seq,res=parse(row['peptide']);total=np.float32(row['calcmass']);charge=int(row['charge']);bm,ym,n,eligible=match(masses,it,res,total,charge)
    assert n==int(row['matched_peaks']),(row['case_id'],n,row['matched_peaks'])
    b=np.float32(0);y=total;annotations=[]
    for i in range(len(res)-1):
        b=np.float32(b+res[i]);y=np.float32(y-res[i])
        for series,neutral,ordinal in [('b',b,i+1),('y',y,len(res)-i-1)]:
            for z in range(1,max(2,min(charge,3))):
                center=np.float32(neutral/np.float32(z));delta=np.float32(center*np.float32(20)/np.float32(1000000));ix=np.flatnonzero((masses>=np.float32(center-delta))&(masses<=np.float32(center+delta)))
                if len(ix):
                    j=ix[np.argmax(it[ix])];annotations.append(dict(series=series,ordinal=ordinal,charge=z,observed_mz=float(ma[j]),calculated_mz=float(center+PROTON),intensity=float(it[j]),ppm=float((masses[j]-center)/center*1e6)))
    assert len(annotations)==n;pd.DataFrame(annotations).to_csv(case/'matched_fragments.tsv',sep='\t',index=False)
    # Bounded alternatives present in the ACTUAL local FASTA: same sequence,
    # single replacements and adjacent swaps, with all Sage C/M modification forms.
    original=canonical(row['peptide']);variants={original:'same_canonical'};alphabet='ACDEFGHKLMNPQRSTVWY'
    for i,a in enumerate(original):
        for b in alphabet:
            if a!=b:variants.setdefault(original[:i]+b+original[i+1:],'single_substitution')
    for i in range(len(original)-1):
        if original[i]!=original[i+1]:variants.setdefault(original[:i]+original[i+1]+original[i]+original[i+2:],'adjacent_swap')
    machine=ahocorasick.Automaton()
    for p in variants:machine.add_word(p,p)
    machine.make_automaton();found=defaultdict(set);audit={}
    for header,protein in fasta(t['fasta']['path'],audit):
        for end,p in machine.iter(protein.replace('I','L')):found[p].add(header.split()[0])
    assert original in found and audit['sha256']==t['fasta']['sha256']
    alt=[]
    for p,parents in sorted(found.items()):
        meth=[i for i,a in enumerate(p) if a=='M']
        for nm in range(min(2,len(meth))+1):
            for oxidized in itertools.combinations(meth,nm):
                form=''.join(a+('[+57.0215]' if a=='C' else '[+15.9949]' if i in oxidized else '') for i,a in enumerate(p))
                _,residues=parse(form);mass=np.float32(18.010565)
                for v in residues:mass=np.float32(mass+v)
                if not 500<=mass<=5000:continue
                isotope=[i for i in range(-1,4) if abs(float(row['expmass'])-i*1.0033548378-float(mass))/float(mass)*1e6<=20]
                if not isotope:continue
                is_original=form.replace('I','L')==row['peptide'].replace('I','L')
                if is_original:mass=total
                ab,ay,an,ae=match(masses,it,residues,mass,charge)
                alt.append(dict(peptide=form,canonical_peptide=p,change=variants[p],parents=';'.join(sorted(parents)),
                    isotope_errors=';'.join(map(str,isotope)),calculated_mass=float(mass),matched_ions=an,
                    original=is_original,at_least_original_matches=an>=n))
    assert any(a['original'] and a['matched_ions']==n for a in alt),'Original hypothesis missing from bounded-alternative reconstruction'
    pd.DataFrame(alt).to_csv(case/'bounded_sequence_alternatives.tsv',sep='\t',index=False)
    for route,source in dict(t['routes'],unspecific=dict(search=oo['results.sage.tsv'])).items():
        verify(source['search'])
        with open(source['search']['path']) as fh:
            for r in csv.DictReader(fh,delimiter='\t'):
                if r['scannr']==row['scannr']:comparators.append(dict(case_id=row['case_id'],route=route,**r))
    raw_identity=meta(t['raw']['path']);verify(t['raw'])
    dump(case/'PROVENANCE.json',dict(raw=raw_identity,fasta=t['fasta'],unspecific_config=oo['results.json'],selected_psm=row,
        interpretation='Measured spectrum; N-side H is a flanking protein-context assignment. Alternative match counts are bounded diagnostics, not recalibrated Sage scores/q values.'))
    summary.append(dict(case_id=row['case_id'],cohort=row['cohort'],compartment=row['compartment'],origin=row['origin'],side=row['side'],
        matched_ions=n,backbone_coverage=int(bm|ym).bit_count()/(len(seq)-1),alternatives=len(alt),
        nonoriginal_at_least_original=sum(not a['original'] and a['at_least_original_matches'] for a in alt)))
pd.DataFrame(summary).to_csv(out/'spectrum_audit.tsv',sep='\t',index=False)
pd.DataFrame(comparators).to_csv(out/'same_scan_all_search_calls.tsv',sep='\t',index=False)
(out/'REPORT.md').write_text('''# New unspecific H-boundary measured spectra

All requested cohort/compartment/origin/side strata are retained, including empties.
Cases are absent from same-acquisition accepted tryptic identities and have an
unambiguous internal H boundary, target rank1, peptide-q<=.01 and spectrum-q<=.01.
Choose the acquisition at median positive candidate-count rank; deterministic nearest
rank fallback only if its stricter PSM set is empty. Within it take the median
matched-intensity representative of canonical-peptide best PSMs. Full candidate
counts and selection tables are saved; cases are not selected for visual quality.

Observed peak lists and b/y annotations independently reproduce Sage matched-ion
counts. All same-scan search calls, including rejected/decoy calls, are retained.
Bounded alternatives cover local-FASTA single substitutions, adjacent swaps and
Sage-supported C/M modification placements with permitted isotope/precursor ranges.
These are fragment-count diagnostics, not exhaustive alternatives, new q values or
a replacement for search scoring. Other PTMs, multiple substitutions, chimericity,
I/L molecular discrimination and cleavage timing are not resolved by these cases.
No polished figure or protease assignment is produced here. Scientific review pending.
''')
dump(out/'COMPLETE.json',dict(status='PASS',pilot=pilot,cohort=cohort,acquisitions=len(tasks),cases=len(selected),
    scientific_review='PENDING',source_receipts=sources,job_id=os.environ.get('SLURM_JOB_ID'),
    outputs=[meta(p) for p in out.rglob('*') if p.is_file() and p.name!='COMPLETE.json']))
print('HISTIDINE SPECTRUM AUDIT PASS',cohort,len(selected),flush=True)
