#!/usr/bin/env python3
import csv,json,math,os,pathlib,re,shutil,sys
from collections import Counter,defaultdict
RUN=pathlib.Path(sys.argv[1]).resolve()
sys.path.insert(0,str(RUN/'frozen'))
from lib import meta,sha,verify,dump,table,command
from peptide import canonical
def fa(path):
    out={};h=None;seq=[]
    for l in open(path):
        if l.startswith('>'):
            if h is not None:assert h not in out;out[h]=''.join(seq)
            h=l[1:].split()[0];seq=[]
        else:seq.append(l.strip().upper())
    if h is not None:assert h not in out;out[h]=''.join(seq)
    return out
def peptide(s):
    mt=re.fullmatch(r'[A-Z*-]\.(.*)\.[A-Z*-]',s)
    if mt:s=mt.group(1)
    s=canonical(s)
    assert re.fullmatch('[A-Z]+',s),s
    return s
def run(index):
    m=json.loads((RUN/'inputs/manifest.json').read_text());t=m['tasks'][index]
    for x in m['frozen_code']+[m['jar'],m['percolator'],t['source_config'],t['fasta'],t['raw']]:verify(x)
    case=RUN/'results'/t['id'];case.mkdir(parents=True,exist_ok=True)
    assert not (case/'COMPLETE.json').exists()
    targets=fa(t['fasta']['path']);assert targets and all(not p.startswith('rev_') for p in targets)
    db=case/'target_decoy.fasta'
    with db.open('w') as f:
        for p,s in targets.items():f.write('>'+p+'\n'+s+'\n>rev_'+p+'\n'+s[::-1]+'\n')
    check=fa(db);assert len(check)==2*len(targets)
    assert all(check['rev_'+p]==s[::-1] for p,s in targets.items())
    dump(case/'database.json',dict(source=t['fasta'],target_entries=len(targets),decoy_entries=len(targets),decoy_policy='One full-sequence reverse per original target; explicit rev_ prefix; no selection on search outcomes',target_decoy_identical_sequences=len(set(targets.values())&{s[::-1] for s in targets.values()}),fasta=meta(db)))
    source=json.loads(pathlib.Path(t['source_config']['path']).read_text())
    cfg={}
    for line in (RUN/'frozen/closed_fragger.params').read_text().splitlines():
        line=line.split('#',1)[0].strip()
        if '=' in line:k,v=line.split('=',1);cfg[k.strip()]=v.strip()
    mods=source['database']
    assert mods['static_mods']=={'C':57.0215} and mods['variable_mods']=={'M':[15.9949]}
    e=mods['enzyme'];assert e['cleave_at']=='KR' and e['restrict']=='P' and not e.get('semi_enzymatic')
    cfg.update(database_name=str(db),num_threads='8',precursor_mass_lower=str(source['precursor_tol']['ppm'][0]),precursor_mass_upper=str(source['precursor_tol']['ppm'][1]),precursor_true_tolerance='20',fragment_mass_tolerance='20',calibrate_mass='0',deisotope='0',deneutralloss='0',isotope_error='/'.join(str(v) for v in range(source['isotope_errors'][0],source['isotope_errors'][1]+1)),search_enzyme_nocut_1='P',allowed_missed_cleavage_1=str(e['missed_cleavages']),num_enzyme_termini='2',variable_mod_01='15.9949 M 2',variable_mod_02='0.0 site_02 0',max_variable_mods_per_peptide='2',digest_min_length=str(e['min_len']),digest_max_length=str(e['max_len']),digest_mass_range='500.0 5000.0',add_C_cysteine='57.0215',precursor_charge='2 4',output_format='tsv_pepxml_pin',output_report_topN='1',output_max_expect='50',report_alternative_proteins='1',minimum_peaks='15',use_topN_peaks='150',min_matched_fragments='4',minimum_ratio='0',remove_precursor_peak='0')
    # Fully disable the additional default acetylation rather than retaining an
    # undeclared modification. Zero mass is the documented inactive slot.
    cfg['variable_mod_02']='0.0 site_02 1'
    params=case/'fragger.params';params.write_text('\n'.join(k+' = '+v for k,v in cfg.items())+'\n')
    raw=case/(t['sample']+'.mzML')
    if not (case/'NATIVE_SEARCH_CACHE.json').exists():shutil.copy2(t['raw']['path'],raw)
    rawmeta=meta(raw);assert rawmeta['sha256']==sha(t['raw']['path'])
    if (case/'NATIVE_SEARCH_CACHE.json').exists():
        cache=json.loads((case/'NATIVE_SEARCH_CACHE.json').read_text())
        for x in cache['files']:assert sha(x['path'])==x['sha256'],'Native search cache changed'
    else:
        command(['java','-Xmx24G','-jar',m['jar']['path'],params,raw],case/'msfragger')
    pins=[case/(t['sample']+'.pin')];assert pins[0].is_file()
    labels=Counter()
    with pins[0].open() as f:
        rd=csv.DictReader(f,delimiter='\t')
        assert 'Label' in rd.fieldnames or 'label' in rd.fieldnames
        for r in rd:
            label=r.get('Label',r.get('label'))
            if label in ['1','-1']:labels[label]+=1
    assert labels['1']>=100 and labels['-1']>=100,labels
    # MSFragger's native PIN appends precursor charge to the peptide string.
    # Remove this encoding BEFORE Percolator's peptide-level competition, retaining
    # the independent charge features. Merely stripping it in final display would
    # incorrectly report peptide-ion q-values as peptide q-values.
    normalized=case/'peptide_level.pin';converted=Counter()
    with pins[0].open() as fi,normalized.open('w') as fo:
        rd=csv.reader(fi,delimiter='\t');header=next(rd)
        w=csv.writer(fo,delimiter='\t',lineterminator='\n');w.writerow(header)
        pi=header.index('Peptide')
        for row in rd:
            match=re.fullmatch(r'([A-Z*-])\.(.+?)(\d+)\.([A-Z*-])',row[pi])
            assert match,('Unexpected native PIN peptide notation',row[pi])
            left,body,charge,right=match.groups();charge=int(charge)
            feature='charge_'+str(charge) if charge<7 else 'charge_7_or_more'
            assert float(row[header.index(feature)])==1,('PIN charge suffix disagrees with charge feature',row[pi])
            assert re.fullmatch('[A-Z]+',canonical(body))
            row[pi]=left+'.'+body+'.'+right;w.writerow(row);converted[str(charge)]+=1
    dump(case/'pin_peptide_identity_audit.json',dict(status='PASS',original=meta(pins[0]),normalized=meta(normalized),charge_counts=dict(converted),rule='Remove validated charge suffix before peptide-level Percolator competition; preserve every row, label and feature'))
    pc=[m['percolator']['path'],'--post-processing-tdc','--seed','913','--num-threads','8','--results-peptides',case/'target_peptides.tsv','--decoy-results-peptides',case/'decoy_peptides.tsv','--results-psms',case/'target_psms.tsv','--decoy-results-psms',case/'decoy_psms.tsv',normalized]
    command(pc,case/'percolator')
    accepted={};allq=Counter()
    with (case/'target_peptides.tsv').open() as f:
        rd=csv.reader(f,delimiter='\t');h=next(rd);qi=h.index('q-value');pi=h.index('peptide')
        for r in rd:
            q=float(r[qi]);assert math.isfinite(q) and 0<=q<=1
            allq['reported_target_peptide_forms']+=1
            if q<=.01:
                p=peptide(r[pi]);accepted[p]=min(q,accepted.get(p,1))
    assert accepted,'No accepted pilot/production peptides; inspect validation'
    import ahocorasick
    a=ahocorasick.Automaton()
    for p in accepted:a.add_word(p,p)
    a.make_automaton();parents=defaultdict(set)
    for pid,s in targets.items():
        for p in {p for _,p in a.iter(s.replace('I','L'))}:parents[p].add(pid)
    assert set(parents)==set(accepted),'An accepted target peptide is absent from this target FASTA'
    rows=[dict(peptide=p,length=len(p),peptide_q=accepted[p],parents=';'.join(sorted(parents[p]))) for p in sorted(accepted)]
    table(case/'accepted_canonical_peptides.tsv',rows)
    for lo in [5,9]:
        (case/('accepted_'+str(lo)+'_50.txt')).write_text(''.join(p+'\n' for p in sorted(accepted) if lo<=len(p)<=50))
    dump(case/'COMPLETE.json',dict(status='PASS',sample=t['sample'],arm=t['arm'],target_entries=len(targets),labels=dict(labels),accepted_canonical_peptides=len(accepted),accepted_length9_50=sum(9<=len(p)<=50 for p in accepted),source_raw=dict(t['raw'],sha256=rawmeta['sha256']),source_fasta=meta(t['fasta']['path']),confidence='MSFragger 4.1; independent per-acquisition Percolator 3.06.5 target-decoy competition, peptide q <=0.01',job_id=os.environ.get('SLURM_JOB_ID'),outputs=[meta(p) for p in sorted(case.iterdir()) if p.is_file() and p.name!='COMPLETE.json']))
    print('VALIDATED',t['id'],len(accepted),flush=True)
if __name__=='__main__':run(int(sys.argv[2]))
