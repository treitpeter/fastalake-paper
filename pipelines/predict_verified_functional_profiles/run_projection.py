"""Frozen-reference peptide evidence reconstruction, with no forced species split."""
import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import traceback
from functools import lru_cache
import re

import ahocorasick
import numpy as np
from scipy import sparse
from evidence_rules import LEVELS, PeptideEvidence, canonical, common_prefix, lineage, parse_terms, signature, valid_intensity, eggnog_domains, taxonomy_check

HERE=Path(__file__).resolve().parent
csv.field_size_limit(1 << 30)

def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()

def fasta(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:
        header=None;buf=[]
        for line in f:
            if line.startswith('>'):
                if header is not None:yield header,''.join(buf).upper()
                header=line[1:].strip();buf=[]
            else:buf.append(line.strip())
        if header is not None:yield header,''.join(buf).upper()

def dump(path, obj):
    path.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')

def records(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:yield from csv.DictReader(f,delimiter='\t')

def genome_lineages(path):
    rows=list(records(path));rep={}
    for r in rows:
        lin=lineage(r['Lineage'])
        if lin:
            old=rep.get(r['Species_rep'])
            if old is not None and old != lin:
                raise ValueError('Conflicting species-representative lineages: '+r['Species_rep'])
            rep[r['Species_rep']]=lin
    return {r['Genome']:lineage(r['Lineage']) or rep.get(r['Species_rep'],()) for r in rows}

def taxonomy(acc, genomes):
    import re
    ids=set(re.findall(r'MGYG\d+',acc))
    if ids:
        vals=[genomes.get(x,()) for x in sorted(ids)]
        if not all(vals):return ()
        result=vals[0]
        for v in vals[1:]:result=common_prefix(result,v)
        return result
    # Only the unambiguous UniProt species suffix is recognized here.
    # Other source namespaces remain explicit unknowns, never parsed as species guesses.
    if re.search(r'_HUMAN(?:\||$)',acc):
        return lineage('d__Eukaryota;p__Chordata;c__Mammalia;o__Primates;f__Hominidae;g__Homo;s__Homo sapiens')
    return ()

def annotation_index(sources):
    byseq={};domains={}
    for src in sources:
        anns={}
        for r in records(src['table']):
            key=r['accession']
            vals={level:parse_terms(r.get({'KO':'KEGG_ko','OG':'OG','Pathway':'KEGG_Pathway'}[level],''),level) for level in LEVELS}
            vals['_domains']=eggnog_domains(r.get('OG',''))
            anns.setdefault(key,[]).append(vals)
        for header,seq in fasta(src['fasta']):
            key=header.split()[0]
            if key not in anns:continue
            sid=hashlib.sha256(canonical(seq).encode()).digest()
            dest=byseq.setdefault(sid,{level:set() for level in LEVELS})
            for row in anns[key]:
                for level in LEVELS:dest[level].add(row[level])
                domains.setdefault(sid,set()).update(row['_domains'])
    return byseq,domains

class NativeEggnog:
    def __init__(self,files):
        self.files=files

    @lru_cache(maxsize=32)
    def genome(self,genome_id):
        p=self.files.get(genome_id)
        if p is None:return {}
        result={}
        for row in records(p):
            key=row['#query']
            vals={level:parse_terms(row.get({'KO':'KEGG_ko','OG':'eggNOG_OGs','Pathway':'KEGG_Pathway'}[level],''),level) for level in LEVELS}
            vals['_domains']=eggnog_domains(row.get('eggNOG_OGs',''))
            if key in result:assert result[key]==vals,('Conflicting native eggNOG query',key)
            result[key]=vals
        return result

    def get(self,acc):
        if not re.fullmatch(r'MGYG\d+_\d+',acc):return None
        return self.genome(acc.split('_')[0]).get(acc)

def run(task,mode,out):
    manifest=json.loads((HERE/'INPUTS.json').read_text())
    expected_code=json.loads((HERE/'CODE_MANIFEST.json').read_text())
    for name,sha in expected_code.items():
        assert digest(HERE/name)==sha,('Code changed after submission',name)
    cfg=manifest['tasks'][task]
    assert cfg['index']==task
    out.mkdir(parents=True,exist_ok=False)
    if mode=='full':
        gate=json.loads((HERE/'pilot_gate/COMPLETE.json').read_text())
        assert gate['status']=='PASS'
    checked={}
    for p in cfg['verify_paths']:
        spec=manifest['files'][p]
        assert Path(p).stat().st_size==spec['bytes'],p
        got=digest(p);assert got==spec['sha256'],p
        checked[p]=got
    print('INPUTS VERIFIED',cfg['id'],len(checked),flush=True)
    summary=list(records(cfg['acquisition_summary']))
    sample_names=[r['sample'] for r in summary]
    assert len(set(sample_names))==len(sample_names)==cfg['acquisitions']
    sample_ix={s:i for i,s in enumerate(sample_names)}
    queries={canonical(r['canonical_peptide']) for r in records(cfg['peptide_table'])}
    if mode=='pilot':
        queries=set(sorted(queries,key=lambda p:hashlib.sha256(p.encode()).digest())[:256])
    ordered=sorted(queries);ids={p:i for i,p in enumerate(ordered)}
    automaton=ahocorasick.Automaton()
    for p,i in ids.items():automaton.add_word(p,i)
    automaton.make_automaton()
    evidence=[PeptideEvidence() for _ in ordered]
    anns,ann_domains=annotation_index(cfg['annotations'])
    native=NativeEggnog(manifest.get('native_eggnog',{}))
    genomes=genome_lineages(manifest['genomes'])
    print('QUERIES',len(ordered),'ANNOTATED FULL SEQUENCES',len(anns),flush=True)
    parent_info={};parent_count=0;edges=0;reference_records=0;domain_checks=defaultdict(int)
    native_checks=defaultdict(int);uhgg_gene_sha={}
    native_sequences=set()
    # Pilot independently verifies all-parent containment for eight query peptides.
    probes=ordered[:8] if mode=='pilot' else []
    naive={p:set() for p in probes};fast={p:set() for p in probes}
    with gzip.open(out/'peptide_parent_edges.tsv.gz','wt',compresslevel=3) as ef, gzip.open(out/'parent_evidence.tsv.gz','wt',compresslevel=3) as pf:
        ew=csv.writer(ef,delimiter='\t');pw=csv.writer(pf,delimiter='\t')
        ew.writerow(['peptide_id','parent_id']);pw.writerow(['parent_id','accession','sequence_sha256_IL','first_reference_file','source_lineage','accepted_lineage','eggnog_domain_check','eggnog_domains','native_eggnog_status','annotation_sets_json','original_header'])
        for reference in cfg['reference_fastas']:
            for header,seq in fasta(reference):
                reference_records+=1
                acc=header.split()[0]
                assert not acc.startswith(('rev_','DECOY_')),('Decoy in target reference',acc)
                seq=canonical(seq)
                sh=None
                if reference==manifest.get('uhgg'):
                    sh=hashlib.sha256(seq.encode()).digest()
                    if acc in uhgg_gene_sha:assert uhgg_gene_sha[acc]==sh,('UHGG identifier sequence conflict',acc)
                    uhgg_gene_sha[acc]=sh
                hits={i for _,i in automaton.iter(seq)}
                if mode=='pilot':
                    slow=[p for p in probes if p in seq]
                    assert {ids[p] for p in slow}==hits.intersection(ids[p] for p in probes)
                if not hits:continue
                sh=sh or hashlib.sha256(seq.encode()).digest()
                parent_id=hashlib.sha256(acc.encode()+b'\0'+sh).hexdigest()
                if parent_id in parent_info:continue
                parent_count+=1
                source_tax=taxonomy(acc,genomes)
                parent_info[parent_id]=(acc,sh,reference,source_tax,header)
                for i in sorted(hits):
                    ew.writerow([i,parent_id]);edges+=1
                    if ordered[i] in fast:fast[ordered[i]].add(parent_id)
                if mode=='pilot':
                    for p in slow:naive[p].add(parent_id)
            print('REFERENCE DONE',Path(reference).name,'parents',parent_count,'edges',edges,flush=True)
        # The full reference gene-sequence index is now available regardless
        # of where the UHGG FASTA occurred in the reference list.
        for parent_id,info in parent_info.items():
                acc,sh,reference,source_tax,header=info
                native_status='not_available'
                ext=native.get(acc)
                if ext is not None:
                    if uhgg_gene_sha.get(acc)!=sh:
                        native_status='gene_sequence_unverified_or_conflicting'
                    else:
                        native_status='gene_sequence_verified'
                        native_sequences.add(sh)
                        for level in LEVELS:anns.setdefault(sh,{}).setdefault(level,set()).add(ext[level])
                        ann_domains.setdefault(sh,set()).update(ext['_domains'])
                parent_info[parent_id]=(acc,sh,reference,source_tax,header,native_status)
        # Finalize only after all reference records have supplied annotations.
        # Thus full-sequence transfer and conflict handling do not depend on
        # FASTA order, a group's representative, or which sample saw a parent.
        for parent_id,info in parent_info.items():
            acc,sh,reference,source_tax,header,native_status=info
            annotations=anns.get(sh,{})
            domains=ann_domains.get(sh,set())
            if native_status=='not_available' and sh in native_sequences:
                native_status='full_sequence_annotation_transfer'
            native_checks[native_status]+=1
            check=taxonomy_check(source_tax,domains)
            if native_status=='gene_sequence_unverified_or_conflicting':check='genome_sequence_conflict_withheld'
            domain_checks[check]+=1
            tax=() if check.endswith('_withheld') else source_tax
            pw.writerow([parent_id,acc,sh.hex(),reference,';'.join(source_tax),';'.join(tax),check,','.join(sorted(domains)),native_status,json.dumps({k:sorted(v) for k,v in annotations.items()},sort_keys=True),header])
            parent_info[parent_id]=(tax,annotations)
    for edge in records(out/'peptide_parent_edges.tsv.gz'):
        tax,annotations=parent_info[edge['parent_id']]
        evidence[int(edge['peptide_id'])].add(tax,annotations)
    del parent_info,uhgg_gene_sha,anns,ann_domains
    if mode=='pilot':assert naive==fast
    assert all(e.parents for e in evidence),[ordered[i] for i,e in enumerate(evidence) if not e.parents][:10]
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
    outputs={p.name:dict(bytes=p.stat().st_size,sha256=digest(p)) for p in out.iterdir() if p.is_file()}
    result=dict(status='PASS',arm=cfg['id'],mode=mode,job_id=os.environ.get('SLURM_JOB_ID'),
        acquisitions=len(sample_names),positive_acquisitions=int((source_totals>0).sum()),
        query_peptides=len(ordered),compatible_parents=parent_count,parent_edges=edges,
        reference_records_scanned=reference_records,quantified_features=int(counts.sum()),
        source_features=int(source_counts.sum()),total_intensity=float(totals.sum()),
        conservation_each_level=True,source_summary_reproduced=True,
        pilot_independent_containment=bool(mode=='pilot'),eggnog_domain_checks=dict(domain_checks),native_eggnog_checks=dict(native_checks),checked_inputs=checked,outputs=outputs,
        scope='Evidence signal relative to the declared searched-protein union plus frozen full UHGG reference. Other upstream reservoirs are not exhaustively rescanned. No species biomass, new FDR or biological regional conclusion.',
        taxonomy_limit='UHGG lineage and explicit human identifiers only; other namespaces stay unresolved. No inferred crosswalk, OS-name parsing or exclusion of unknown parents.',
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
