#!/usr/bin/env python3
import csv,gzip,json,math,pathlib,re,sys
from collections import Counter,defaultdict
import ahocorasick
RUN=pathlib.Path(sys.argv[1]).resolve();sys.path.insert(0,str(RUN/'frozen'))
from lib import meta,sha,verify,dump,table,load_queries,low,clean
from peptide import canonical
from search_helpers import read_fasta
def run(index):
    m=json.loads((RUN/'inputs/manifest.json').read_text());t=m['tasks'][index]
    for x in m['frozen_code']+[t['prediction'],t['fasta'],t['search']]:verify(x)
    out=RUN/'results'/t['cohort']/t['sample'];out.mkdir(parents=True,exist_ok=True)
    assert not (out/'COMPLETE.json').exists()
    counts,filtering,retained=load_queries(t['prediction']['path']);top=set(counts)
    rows=list(csv.DictReader(open(t['prediction']['path'])));qual=defaultdict(set);all_predictions=set();beams=Counter();schema=Counter()
    for i,r in enumerate(rows,1):
        s=(r.get('peptide_prediction_detokenized_unmodified') or '').strip().upper()
        if re.fullmatch('[A-Z]+',s):all_predictions.add(s.replace('I','L'))
        beams[r.get('beam_rank','unrecorded')]+=1
        try:score=float(r.get('score','0.5'))
        except ValueError:score=.5;schema['score_parse_fallback']+=1
        assert math.isfinite(score),'Non-finite score requires an explicit parsimony evidence audit'
        if not s or not 9<=len(s)<=50 or score<=0 or low(s):continue
        if not re.fullmatch('[A-Z]+',s):schema['literal_unmatchable_parsimony_syntax']+=1
        qual[s.replace('I','L')].add(r.get('spec_idx') or str(i))
    schema['lookup_queries_different_from_parsimony_literals']=len(top-set(qual))
    top_spectra=defaultdict(set)
    for i,q,_ in retained:top_spectra[q].add(rows[i-1].get('spec_idx') or str(i))
    accepted=defaultdict(set);identified=defaultdict(set);native_confident=set()
    for r in csv.DictReader(open(t['search']['path']),delimiter='\t'):
        if r['label']!='1' or r['rank']!='1' or not float(r['peptide_q'])<=.01:continue
        p=canonical(r['peptide']);assert re.fullmatch('[A-Z]+',p)
        for prot in r['proteins'].split(';'):
            if not prot.startswith('rev_'):accepted[p].add(prot);identified[prot].add(p)
        if float(r['protein_q'])<=.01:native_confident.update(p for p in r['proteins'].split(';') if not p.startswith('rev_'))
    correct_top=set();correct_all=set();refpaired=None
    if t.get('reference_fate'):
        verify(t['reference_fate']);n=0;z=Counter();strata=defaultdict(Counter)
        with gzip.open(t['reference_fate']['path'],'rt') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                n+=1
                if r['reference_accepted']!='1':continue
                z['paired_predictions']+=1;ok=r['category']=='exact_IL';z['exact_predictions']+=ok
                q=r['prediction_IL']
                if ok:correct_all.add(q)
                if ok and r['production_row_retained']=='1':correct_top.add(q)
                k=(len(q),r['production_row_retained'],r['category']);strata[k]['n']+=1
        assert n==len(rows)
        table(out/'reference_agreement_by_length.tsv',[dict(length=k[0],production_retained=k[1],category=k[2],n=v['n']) for k,v in sorted(strata.items())])
        refpaired=dict(z)
    qkeys=set(qual);aseqs=set(accepted)
    fasta=read_fasta(t['fasta']['path']);assert set(identified)<=set(fasta)
    a=ahocorasick.Automaton()
    for q in sorted(qkeys|aseqs|top):a.add_word(q,q)
    a.make_automaton();matches={};multiplicity=Counter()
    for p,s in fasta.items():
        qs={q for _,q in a.iter(s)};matches[p]=qs;multiplicity.update(qs&qkeys)
    records=[];edgefile=gzip.open(out/'prediction_protein_edges.tsv.gz','wt')
    w=csv.writer(edgefile,delimiter='\t',lineterminator='\n');w.writerow(['query','protein','production_lookup','parsimony_qualifying','same_spectrum_REF_correct','selected_protein_multiplicity','qualifying_spectra'])
    for p,qs in matches.items():
        tq=qs&top;pq=qs&qkeys;ids=identified.get(p,set())
        assert tq and pq,('Selected protein has no reconstructed selection evidence',p)
        for q in sorted(pq):w.writerow([q,p,int(q in top),1,int(q in correct_all) if t.get('reference_fate') else '',multiplicity[q],len(qual[q])])
        records.append(dict(protein=p,identified_accession=int(bool(ids)),native_protein_q01=int(p in native_confident),production_prediction_sequences=len(tq),production_prediction_spectra=len(set().union(*(top_spectra[q] for q in tq))),parsimony_prediction_sequences=len(pq),parsimony_prediction_spectra=len(set().union(*(qual[q] for q in pq))),unique_among_selected_proteins=sum(multiplicity[q]==1 for q in pq),shared_among_selected_proteins=sum(multiplicity[q]>1 for q in pq),same_spectrum_REF_correct_top=len(tq&correct_top) if t.get('reference_fate') else '',same_spectrum_REF_correct_all=len(pq&correct_all) if t.get('reference_fate') else '',predictions_corroborated_by_FL_peptides=len(pq&aseqs),final_assigned_peptides=len(ids),final_peptides_not_predicted_any=len(ids-all_predictions),final_peptides_not_in_top30=len(ids-top),final_peptides_not_in_parsimony=len(ids-qkeys)))
    edgefile.close();table(out/'protein_support.tsv',records)
    populations={}
    for name,rs in [('all_selected',records),('identified_accessions',[r for r in records if r['identified_accession']]),('native_protein_q01',[r for r in records if r['native_protein_q01']])]:
        z=dict(n=len(rs),one_top_prediction=sum(r['production_prediction_sequences']==1 for r in rs),one_parsimony_prediction=sum(r['parsimony_prediction_sequences']==1 for r in rs),with_FL_corroborated_prediction=sum(r['predictions_corroborated_by_FL_peptides']>0 for r in rs))
        if t.get('reference_fate'):z.update(with_REF_correct_top=sum(r['same_spectrum_REF_correct_top']>0 for r in rs),with_REF_correct_all=sum(r['same_spectrum_REF_correct_all']>0 for r in rs))
        populations[name]=z
    dump(out/'COMPLETE.json',dict(status='PASS',sample=t['sample'],cohort=t['cohort'],filtering=filtering,beam_rank_counts=dict(beams),schema_audit=dict(schema),populations=populations,reference_agreement=refpaired,accepted_peptides=len(accepted),accepted_not_predicted_any=len(aseqs-all_predictions),accepted_not_in_top30=len(aseqs-top),accepted_not_in_parsimony=len(aseqs-qkeys),inputs=[meta(t[k]['path']) for k in ['prediction','fasta','search']],outputs=[meta(p) for p in out.iterdir() if p.is_file()]))
    print('SUPPORT VALIDATED',t['cohort'],t['sample'],len(records),flush=True)
if __name__=='__main__':run(int(sys.argv[2]))
