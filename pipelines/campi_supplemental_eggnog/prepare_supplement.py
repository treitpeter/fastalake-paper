"""Find never-queried reference sequences, preserving their actual I/L residues."""
from pathlib import Path
from collections import defaultdict,Counter
import hashlib,json,csv,gzip,os,sys
N=Path(__file__).resolve().parent
E=N.parent/'evidence_all_cohorts'
sys.path.insert(0,str(E))
from run_projection import fasta,records,digest,dump

def sh(seq):return hashlib.sha256(seq.encode()).hexdigest()

def main():
    spec=json.loads((E/'INPUTS.json').read_text());frozen=json.loads((N/'SUPPLEMENT_PREP_MANIFEST.json').read_text())
    for p,sha in frozen['code'].items():assert digest(N/p)==sha,p
    assert digest(E/'INPUTS.json')==frozen['original_input_manifest_sha256']
    out=N/'supplemental_functions';out.mkdir(exist_ok=False)
    candidates={};origins=defaultdict(set);input_hashes={}
    for p in (N/'inventory').glob('*/non_uhgg_reference_proteins.tsv.gz'):
        assert digest(p)==json.loads((p.parent/'COMPLETE.json').read_text())['output_sha256']
        for r in records(p):
            # GMGC/GMSC native paired query FASTAs are being exhaustively reused.
            if r['namespace'] in ('GMGC','GMSC'):continue
            candidates.setdefault(r['sequence_sha256_IL'],set()).add(r['namespace'])
            origins[r['sequence_sha256_IL']].add(p.parent.name)
    initial=len(candidates);queried=set();selected_files=sorted({a['fasta'] for t in spec['tasks'] for a in t['annotations']})
    for p in selected_files:
        assert digest(p)==spec['files'][p]['sha256'];input_hashes[p]=spec['files'][p]['sha256']
        for h,seq in fasta(p):
            key=sh(seq.upper().replace('I','L'))
            if key in candidates:queried.add(key)
    for key in queried:candidates.pop(key)
    # UHGG's native genome annotations have already been reused by complete
    # protein identity. Do not spend a fresh search on another name for it.
    p=spec['uhgg'];assert digest(p)==spec['files'][p]['sha256'];input_hashes[p]=spec['files'][p]['sha256'];uhgg=set()
    for h,seq in fasta(p):
        key=sh(seq.upper().replace('I','L'))
        if key in candidates:uhgg.add(key)
    for key in uhgg:candidates.pop(key)
    queries={};provenance=defaultdict(set)
    refs=sorted({p for t in spec['tasks'] for p in t['reference_fastas'] if p!=spec['uhgg']})
    for p in refs:
        assert digest(p)==spec['files'][p]['sha256'];input_hashes[p]=spec['files'][p]['sha256']
        for h,seq in fasta(p):
            key=sh(seq.upper().replace('I','L'))
            if key not in candidates:continue
            # Annotate real amino acid sequences, NOT artificially all-L proteins.
            raw=seq.upper()
            assert '*' not in raw[:-1], ('Internal stop in searched target',h)
            q='FLSEQ_'+sh(raw);queries[q]=raw;provenance[q].add((h.split()[0],p,key))
    assert {sh(s.replace('I','L')) for s in queries.values()}==set(candidates)
    ordered=sorted(queries)
    # Deterministic hash ordering; not chosen for expected annotation success.
    pilot=ordered[:64];chunks=[ordered[i:i+4000] for i in range(0,len(ordered),4000)]
    for name,ids in [('pilot',pilot)]+[(f'chunk_{i:03d}',ids) for i,ids in enumerate(chunks)]:
        with open(out/(name+'.faa'),'w') as f:
            for q in ids:f.write('>'+q+'\n'+queries[q]+'\n')
    with gzip.open(out/'query_provenance.tsv.gz','wt',compresslevel=3) as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['query','raw_sequence_sha256','sequence_sha256_IL','source_accession','source_fasta','cohorts'])
        for q in ordered:
            for acc,p,key in sorted(provenance[q]):w.writerow([q,q[6:],key,acc,p,','.join(sorted(origins[key]))])
    counts=Counter()
    for key in candidates:
        for cohort in origins[key]:counts[cohort]+=1
    report=dict(status='PASS',initial_non_catalog_IL_sequences=initial,already_queried_selected=len(queried),already_native_UHGG=len(uhgg),
        missing_IL_sequences=len(candidates),raw_annotation_queries=len(queries),cohort_missing_sequences=dict(counts),chunks=len(chunks),pilot_queries=len(pilot),
        input_hashes=input_hashes,outputs={p.name:digest(p) for p in out.iterdir() if p.is_file()},job_id=os.environ.get('SLURM_JOB_ID'),
        method='Same manuscript eggNOG 5.0.2 / mapper2.1.12. Only never-queried reference sequences; no retries selected for missing KO or adverse results. Raw I/L variants are annotated separately and combined only by the existing evidence rules.')
    dump(out/'PREPARE_COMPLETE.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('input_hashes','outputs')},indent=2),flush=True)

if __name__=='__main__':main()
