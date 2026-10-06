"""Reproject the frozen complete parent graph after source annotation repair.

The quantitative emission below is copied unchanged from the validated initial
projection. Parent edges, peptides, acquisition roster and LFQ inputs are fixed.
"""
from pathlib import Path
from collections import defaultdict, Counter
import argparse, csv, gzip, hashlib, json, os, sys, traceback
import numpy as np
from scipy import sparse
HERE=Path(__file__).resolve().parent
N=HERE.parent
E=N.parent/'evidence_all_cohorts'
sys.path.insert(0,str(E))
from run_projection import records,digest,dump
from evidence_rules import LEVELS,PeptideEvidence,canonical,common_prefix,lineage,signature,valid_intensity
csv.field_size_limit(1 << 30)

def domain_qc(tax,domains):
    if not tax:return 'no_resolved_lineage'
    if not domains:return 'no_eggnog_domain_annotation'
    if len(domains)>1:return 'orthology_domains_disagree_flagged'
    return 'domain_compatible' if tax[0][3:] in domains else 'orthology_domain_differs_flagged'

def merge_annotations(old,extra):
    result={k:{tuple(v) for v in values} for k,values in old.items()}
    for level,values in extra.items():result.setdefault(level,set()).update(tuple(v) for v in values)
    return result

def load_adapter():
    root=N/'adapter';receipt=json.loads((root/'COMPLETE.json').read_text());assert receipt['status']=='PASS'
    for name,sha in receipt['outputs'].items():assert digest(root/name)==sha,name
    taxa={};tax_status=defaultdict(set);ann={};dom={}
    for r in records(root/'taxonomy.tsv.gz'):
        key=(r['accession'],r['sequence_sha256_IL']);value=lineage(r['lineage'])
        tax_status[key].add(r['status'])
        # Missing annotation in an alternate copy of the SAME source protein is
        # not an additional unknown biological parent. Conflicting positive
        # source lineages are intersected; genuine crosswalk failures are kept.
        if value:
            taxa[key]=common_prefix(taxa[key],value) if key in taxa else value
        elif r['status']=='source_crosswalk_unresolved_or_conflicting':taxa[key]=()
    with gzip.open(root/'functions.jsonl.gz','rt') as f:
        for r in map(json.loads,f):
            sha=r['sequence_sha256_IL'];ann[sha]=r['annotations'];dom[sha]=set(r['domains'])
    return taxa,tax_status,ann,dom

def run(task,mode,out):
    manifest=json.loads((HERE/'INPUTS.json').read_text())
    frozen=json.loads((HERE/'CODE_MANIFEST.json').read_text())
    for name,sha in frozen.items():assert digest(HERE/name)==sha,name
    source_manifest=json.loads((N/'REPROJECTION_INPUTS.json').read_text())
    for path,sha in source_manifest.items():assert digest(path)==sha,path
    cfg=manifest['tasks'][task];assert cfg['index']==task
    out.mkdir(parents=True,exist_ok=False)
    if mode=='full':assert json.loads((HERE/'pilot_gate/COMPLETE.json').read_text())['status']=='PASS'
    original=E/'full'/cfg['id'];receipt=json.loads((original/'COMPLETE.json').read_text())
    assert receipt['status']=='PASS' and receipt['mode']=='full'
    checked={}
    for name in ['parent_evidence.tsv.gz','peptide_parent_edges.tsv.gz','peptide_evidence.tsv.gz']:
        path=original/name;got=digest(path);assert got==receipt['outputs'][name]['sha256'];checked[str(path)]=got
    for path in [cfg['feature_ledger'],cfg['acquisition_summary']]:
        got=digest(path);assert got==manifest['files'][path]['sha256'];checked[path]=got
    summary=list(records(cfg['acquisition_summary']));sample_names=[r['sample'] for r in summary]
    assert len(set(sample_names))==len(sample_names)==cfg['acquisitions']
    sample_ix={s:i for i,s in enumerate(sample_names)}
    originals=[{k:r[k] for k in ('peptide_id','canonical_peptide','compatible_parents')}
               for r in records(original/'peptide_evidence.tsv.gz')]
    assert all(int(r['peptide_id'])==i for i,r in enumerate(originals))
    ordered=[r['canonical_peptide'] for r in originals]
    if mode=='pilot' and cfg['id']!='CAMPI_REF':
        ordered=sorted(sorted(ordered,key=lambda p:hashlib.sha256(p.encode()).digest())[:256])
    ids={p:i for i,p in enumerate(ordered)}
    old_to_new={int(r['peptide_id']):ids[r['canonical_peptide']] for r in originals if r['canonical_peptide'] in ids}
    expected_parents={int(r['peptide_id']):int(r['compatible_parents']) for r in originals if int(r['peptide_id']) in old_to_new}
    del originals
    evidence=[PeptideEvidence() for _ in ordered];needed=set();edges=0
    qc_flags=np.zeros((len(ordered),2),dtype=bool)
    with gzip.open(out/'peptide_parent_edges.tsv.gz','wt',compresslevel=3) as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['peptide_id','parent_id'])
        for row in records(original/'peptide_parent_edges.tsv.gz'):
            old=int(row['peptide_id'])
            if old in old_to_new:
                w.writerow([old_to_new[old],row['parent_id']]);needed.add(row['parent_id']);edges+=1
    taxa,statuses,extra,extra_domains=load_adapter()
    parent_info={};domain_checks=Counter();native_checks=Counter();repair_checks=Counter()
    with gzip.open(out/'parent_evidence.tsv.gz','wt',compresslevel=3) as f:
        w=None
        for row in records(original/'parent_evidence.tsv.gz'):
            if row['parent_id'] not in needed:continue
            key=(row['accession'],row['sequence_sha256_IL']);sha=key[1]
            old_tax=lineage(row['source_lineage']);old_accepted=row['accepted_lineage'];old_check=row['eggnog_domain_check']
            tax=taxa.get(key,old_tax)
            annotations=merge_annotations(json.loads(row['annotation_sets_json']),extra.get(sha,{}))
            domains=set(filter(None,row['eggnog_domains'].split(',')))|extra_domains.get(sha,set())
            check=domain_qc(tax,domains)
            if row['native_eggnog_status']=='gene_sequence_unverified_or_conflicting':
                tax=();check='genome_sequence_conflict_withheld'
            source_status=','.join(sorted(statuses.get(key,()))) or 'original_source_genome_taxonomy'
            repair_checks[source_status]+=1;domain_checks[check]+=1;native_checks[row['native_eggnog_status']]+=1
            row.update(source_lineage=';'.join(tax),accepted_lineage=';'.join(tax),eggnog_domain_check=check,
                eggnog_domains=','.join(sorted(domains)),annotation_sets_json=json.dumps({k:sorted(v) for k,v in annotations.items()},sort_keys=True),
                annotation_repair_status=source_status,original_source_lineage=';'.join(old_tax),
                original_accepted_lineage=old_accepted,
                original_eggnog_domain_check=old_check,
                common_function_annotation_index_available=int(sha in extra))
            if w is None:w=csv.DictWriter(f,fieldnames=list(row),delimiter='\t');w.writeheader()
            w.writerow(row);parent_info[row['parent_id']]=(tax,annotations,
                check.startswith('orthology_'),bool(old_tax and not old_accepted and tax))
    assert set(parent_info)==needed
    parent_count=len(parent_info)
    for edge in records(out/'peptide_parent_edges.tsv.gz'):
        tax,annotations,flag,restored=parent_info[edge['parent_id']];pi=int(edge['peptide_id'])
        evidence[pi].add(tax,annotations);qc_flags[pi]|=(flag,restored)
    assert all(e.parents for e in evidence)
    assert all(evidence[new].parents==expected_parents[old] for old,new in old_to_new.items())
    if mode=='full':assert edges==receipt['parent_edges'] and parent_count==receipt['compatible_parents']
    del parent_info,needed,taxa,statuses,extra,extra_domains
    reference_records=0
    qc_signal=np.zeros((2,len(sample_names)));qc_species_signal=np.zeros_like(qc_signal)
    species_supported=np.array([e.taxonomy()['rank']=='species' for e in evidence])
    catalog={};lookup={level:[] for level in LEVELS};tax_catalog={};tax_ids=[]
    with gzip.open(out/'peptide_evidence.tsv.gz','wt',compresslevel=3) as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['peptide_id','canonical_peptide','compatible_parents','taxonomy_json','function_evidence_json'])
        for i,(pep,e) in enumerate(zip(ordered,evidence)):
            tax=e.taxonomy();fun={level:e.functions[level].result() for level in LEVELS}
            w.writerow([i,pep,e.parents,json.dumps(tax,sort_keys=True),json.dumps(fun,sort_keys=True)])
            taxkey=json.dumps({k:tax[k] for k in ('state','rank','supported_lineage')},sort_keys=True)
            taxid=tax_catalog.setdefault(taxkey,len(tax_catalog));tax_ids.append(taxid)
            for level in LEVELS:
                key,record=signature(level,fun[level],tax)
                if key not in catalog:catalog[key]=record
                lookup[level].append(key)
    totals=np.zeros(len(sample_names));counts=np.zeros(len(sample_names),dtype=np.int64)
    source_totals=np.zeros(len(sample_names));source_counts=np.zeros(len(sample_names),dtype=np.int64)
    matrices={level:defaultdict(float) for level in LEVELS};taxmatrix=defaultdict(float)
    last_line={};seen_sample_blocks=set();current_sample=None
    for r in records(cfg['feature_ledger']):
        s=r['sample'];si=sample_ix[s];line=int(r['input_line']);v=valid_intensity(r['intensity'])
        # The accepted ledger is source-order sorted. Fail rather than risk double counting.
        if current_sample!=s:
            assert s not in seen_sample_blocks,('Non-contiguous sample',s)
            seen_sample_blocks.add(s);current_sample=s
        assert line>last_line.get(s,-1),('Duplicated or unsorted LFQ feature',s,line)
        last_line[s]=line
        p=canonical(r['canonical_peptide'])
        source_totals[si]+=v;source_counts[si]+=1
        if mode=='pilot' and p not in ids:continue
        assert p in ids,('Quantified peptide missing from dictionary',p)
        pi=ids[p];totals[si]+=v;counts[si]+=1
        qc_signal[:,si]+=qc_flags[pi]*v
        if species_supported[pi]:qc_species_signal[:,si]+=qc_flags[pi]*v
        for level in LEVELS:matrices[level][lookup[level][pi],si]+=v
        taxmatrix[tax_ids[pi],si]+=v
    expected=np.array([float(r['intensity']) for r in summary])
    expected_counts=np.array([int(r['lfq_features']) for r in summary])
    assert np.allclose(source_totals,expected,rtol=1e-9,atol=1e-6)
    assert np.array_equal(source_counts,expected_counts)
    np.savez_compressed(out/'sample_totals.npz',samples=np.array(sample_names),intensity=totals,features=counts,source_intensity=source_totals,source_features=source_counts)
    with gzip.open(out/'quantity_definitions.jsonl.gz','wt',compresslevel=3) as f:
        for key in sorted(catalog):f.write(json.dumps(dict(quantity_id=key,**catalog[key]),sort_keys=True)+'\n')
    category_rows=[]
    for level,d in matrices.items():
        names=sorted({key for key,si in d});ix={k:i for i,k in enumerate(names)}
        mat=sparse.coo_matrix(([v for v in d.values()],([ix[k] for k,si in d],[si for k,si in d])),shape=(len(names),len(sample_names))).tocsr()
        assert np.allclose(np.asarray(mat.sum(axis=0)).ravel(),totals,rtol=1e-9,atol=1e-6)
        sparse.save_npz(out/f'{level}_evidence_signal.npz',mat)
        (out/f'{level}_quantity_ids.json').write_text(json.dumps(names)+'\n')
        # Single-function rows remain annotated with their supported taxonomic rank.
        single=np.array([catalog[k]['function']['state']=='single_consistent' for k in names])
        sparse.save_npz(out/f'{level}_single_consistent_signal.npz',mat[single])
        (out/f'{level}_single_consistent_quantity_ids.json').write_text(json.dumps([k for k,b in zip(names,single) if b])+'\n')
        categories=defaultdict(lambda:np.zeros(len(sample_names)))
        for ri,k in enumerate(names):
            rec=catalog[k];category=(rec['function']['state'],rec['taxonomy']['state'],rec['taxonomy']['rank'])
            categories[category]+=mat.getrow(ri).toarray().ravel()
        for (fs,ts,rank),vals in categories.items():
            for si,v in enumerate(vals):category_rows.append([level,sample_names[si],fs,ts,rank,v,totals[si],v/totals[si] if totals[si] else ''])
    with open(out/'sample_coverage.tsv','w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['level','sample','function_state','taxonomy_state','taxonomic_rank','intensity','total_intensity','fraction']);w.writerows(category_rows)
    tx=sorted(tax_catalog,key=tax_catalog.get)
    tm=sparse.coo_matrix(([v for v in taxmatrix.values()],([k for k,s in taxmatrix],[s for k,s in taxmatrix])),shape=(len(tx),len(sample_names))).tocsr()
    assert np.allclose(np.asarray(tm.sum(axis=0)).ravel(),totals,rtol=1e-9,atol=1e-6)
    sparse.save_npz(out/'taxonomy_evidence_signal.npz',tm);dump(out/'taxonomy_rows.json',[json.loads(x) for x in tx])
    with open(out/'orthology_taxonomy_QC_signal.tsv','w') as f:
        w=csv.writer(f,delimiter='\t')
        w.writerow(['sample','total_LFQ','any_parent_orthology_domain_flag_signal','species_supported_with_any_parent_domain_flag_signal','any_parent_previous_domain_veto_restored_signal','species_supported_with_previous_domain_veto_restored_signal'])
        for si,s in enumerate(sample_names):w.writerow([s,totals[si],qc_signal[0,si],qc_species_signal[0,si],qc_signal[1,si],qc_species_signal[1,si]])
    outputs={p.name:dict(bytes=p.stat().st_size,sha256=digest(p)) for p in out.iterdir() if p.is_file()}
    result=dict(status='PASS',arm=cfg['id'],mode=mode,job_id=os.environ.get('SLURM_JOB_ID'),
        acquisitions=len(sample_names),positive_acquisitions=int((source_totals>0).sum()),
        query_peptides=len(ordered),compatible_parents=parent_count,parent_edges=edges,
        reference_records_scanned=reference_records,quantified_features=int(counts.sum()),
        source_features=int(source_counts.sum()),total_intensity=float(totals.sum()),
        conservation_each_level=True,source_summary_reproduced=True,
        pilot_independent_containment=False, frozen_complete_parent_graph_preserved=True, annotation_repair_status_counts=dict(repair_checks),eggnog_domain_checks=dict(domain_checks),native_eggnog_checks=dict(native_checks),checked_inputs=checked,outputs=outputs,
        scope='Evidence signal relative to the declared searched-protein union plus frozen full UHGG reference. Other upstream reservoirs are not exhaustively rescanned. No species biomass, new FDR or biological regional conclusion.',
        taxonomy_limit='Verified UHGG/GMGC/GMSC/UniProt native taxonomy; GTDB r202 backbone; missing assembly parents retained. Orthology-domain differences are QC flags; affected signal is separately reported.',
        grouping='Existing frozen cohort protein groups remain reporting labels. They are not used to choose taxonomic or functional parents.')
    dump(out/'COMPLETE.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('checked_inputs','outputs')},indent=2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('task',type=int);parser.add_argument('mode',choices=['pilot','full']);args=parser.parse_args()
    cfg=json.loads((HERE/'INPUTS.json').read_text())['tasks'][args.task]
    out=HERE/args.mode/cfg['id']
    try:run(args.task,args.mode,out)
    except Exception:
        if out.exists() and not (out/'COMPLETE.json').exists():dump(out/'FAILURE.json',dict(status='FAIL',traceback=traceback.format_exc(),job_id=os.environ.get('SLURM_JOB_ID')))
        raise
