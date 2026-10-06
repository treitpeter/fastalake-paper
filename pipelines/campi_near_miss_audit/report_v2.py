import csv,gzip,json,pathlib
from collections import Counter,defaultdict
from lib import dump,meta,sha,table
def build(run,m):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    summaries=[];checks=[];audit=[];pooled=defaultdict(Counter);total=paired=0
    metrics=['reference_present_final','reference_identified_fl_anywhere','reference_identified_fl_same_spectrum','same_candidate_contains_reference','same_selected_contains_reference','family_any_compatible','family_exclusively_compatible','selected_family_any_compatible','selected_family_exclusively_compatible']
    for sm in m['samples']:
        out=run/'results'/sm['sample'];check=json.loads((out/'COMPLETE.json').read_text());assert check['status']=='PASS'
        for r in check['outputs']:assert sha(r['path'])==r['sha256'],r['path']
        rows=list(csv.DictReader((out/'summary.tsv').open(),delimiter='\t'));summaries+=rows;checks.append(check)
        observed=defaultdict(Counter);n=nr=0
        with gzip.open(out/'spectrum_prediction_fate.tsv.gz','rt') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                n+=1
                if r['reference_accepted']!='1':continue
                nr+=1;pop=['all_reference_paired']
                if r['production_row_retained']=='1':
                    pop.append('production_retained_reference_paired')
                    if r['also_spectrum_q_01']=='1':pop.append('production_retained_reference_paired_both_q01')
                for p in pop:
                    z=observed[p,r['category']];z['n']+=1
                    for k in metrics:
                        if r[k]!='':z[k+'_evaluable']+=1;z[k]+=int(r[k])
                    if r['candidate_proteins']!='':
                        z['any_lookup']+=int(int(r['candidate_proteins'])>0)
                        z['all_candidates_annotated']+=int(int(r['candidate_proteins'])>0 and r['candidate_proteins']==r['annotated_candidates'])
        assert n==check['prediction_rows'] and nr==check['paired_prediction_rows']
        for r in rows:
            z=observed[r['population'],r['category']]
            for k in r:
                if k not in ['sample','population','category']:assert int(r[k] or 0)==z[k],(sm['sample'],r['population'],r['category'],k)
            pooled[r['population'],r['category']].update(z)
        total+=n;paired+=nr
        audit.append(dict(sample=sm['sample'],prediction_rows=n,reference_paired_rows=nr,retained_rows=check['retained_rows'],rows_with_bags=check['rows_with_bags'],reference_peptides_with_family=check['reference_peptides_with_family'],reference_unique_peptides=check['reference_unique_peptides'],candidate_fasta_entries=check['candidate_fasta']['entries'],final_fasta_entries=check['final_fasta']['entries'],scan_offset=check['scan_offset'],annotation_conflicts=check['annotation_conflicts'],output_hashes='PASS',row_reconciliation='PASS'))
    assert len(checks)==15
    results=run/'results'
    table(results/'acquisition_summary.tsv',summaries,sorted(set().union(*(set(r) for r in summaries))))
    table(results/'acquisition_validation.tsv',audit)
    prows=[dict(population=p,category=c,**dict(v)) for (p,c),v in sorted(pooled.items())]
    table(results/'pooled_summary.tsv',prows,sorted(set().union(*(set(r) for r in prows))))
    cats=['exact_IL','adjacent_swap','rearrangement','substitution','length_difference','other_disagreement']
    labels=['Exact I/L','Adjacent swap','Rearrangement','Substitution','Length difference','Other discrepancy']
    plt.rcParams.update({'font.size':9,'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(16.5,5.2))
    y=np.arange(len(cats));alln=[pooled['all_reference_paired',c]['n'] for c in cats];kept=[pooled['production_retained_reference_paired',c]['n'] for c in cats]
    ax=axes[0];ax.barh(y,alln,color='#cdd6da',label='All reference-paired predictions');ax.barh(y,kept,color='#177b8f',label='Production-retained predictions');ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlabel('Prediction rows');ax.set_title('a  What disagreements occur?',loc='left');ax.legend(loc='lower right',fontsize=7);ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    def vals(key):return [100*pooled['production_retained_reference_paired',c][key]/n if n else np.nan for c,n in zip(cats,kept)]
    ax=axes[1]
    for delta,key,color,label in [(-.22,metrics[0],'#177b8f','REF peptide in final FASTA'),(0,metrics[1],'#df983a','FL identifies it in the acquisition'),(.22,metrics[2],'#704d8c','FL identifies it on the same spectrum')]:ax.barh(y+delta,vals(key),height=.22,color=color,label=label)
    ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlim(0,100);ax.set_xlabel('% of retained, REF-paired predictions');ax.set_title('b  Does reference evidence survive?',loc='left');ax.legend(loc='lower right',fontsize=7)
    ax=axes[2];ax.barh(y-.15,vals(metrics[3]),height=.3,color='#177b8f',label='Same candidate contains REF peptide');ax.barh(y+.15,vals(metrics[5]),height=.3,color='#95b499',label='Observed REF-compatible family');ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlim(0,100);ax.set_xlabel('% of retained, REF-paired predictions');ax.set_title('c  What did exact lookup retrieve?',loc='left');ax.legend(loc='lower right',fontsize=7)
    fig.suptitle('CAMPI: from de novo disagreement to searchable reference evidence',x=.06,ha='left',fontsize=13)
    fig.text(.06,.015,'15 acquisitions · target rank-1 peptide q ≤ 0.01 · I/L-normalized bare sequences · incomplete family annotation; coverage and exclusive-compatibility counts in source tables',fontsize=8)
    fig.tight_layout(rect=(0,.07,1,.91),w_pad=2)
    for ext in ['png','pdf','svg']:fig.savefig(results/('A03_near_misses.'+ext),dpi=220)
    plt.close(fig)
    nonexact=Counter()
    for c in cats[1:]:nonexact.update(pooled['production_retained_reference_paired',c])
    def ratio(z,k):return f"{z[k]:,}/{z['n']:,} ({100*z[k]/z['n']:.2f}%)" if z['n'] else 'Undefined'
    lines=['# Analysis 3 — realistic prediction disagreements','',
    'All 15 CAMPI acquisitions were evaluated against the archived REF search and the September 9 blind UniProt unclustered FL arm. Every original prediction row is retained in the evidence tables. Spectrum joins, production query sets, original hit counts, final FASTA scans and row-level summaries passed validation.','',
    f'Available predictions: {total:,}; paired to an accepted REF spectrum: {paired:,}; unpaired: {total-paired:,}.','',
    '| Class | All REF-paired | Production-retained | REF peptide in final FASTA | FL identification anywhere | FL identification on same spectrum |','|---|---:|---:|---:|---:|---:|']
    for c,l in zip(cats,labels):
        z=pooled['production_retained_reference_paired',c];lines.append('| '+l+' | '+f"{pooled['all_reference_paired',c]['n']:,} | {z['n']:,} | "+' | '.join(ratio(z,k) for k in metrics[:3])+' |')
    lines+=['','## Interpretation against Matthias’s request','',
    f"Among non-exact retained REF-paired predictions, {ratio(nonexact,metrics[3])} have at least one exact-lookup candidate containing the reference peptide as well. After parsimony, {ratio(nonexact,metrics[4])} retain that co-occurrence in one selected protein. Separately, {ratio(nonexact,metrics[0])} retain the reference peptide somewhere in the final FASTA.",'',
    'These endpoints separate useful retrieval by the particular prediction from retention supplied by any prediction in the acquisition. Other spectra can retrieve the protein containing the reference peptide; that does not demonstrate rescue by the discrepant prediction itself.','',
    f"Observed family compatibility: {ratio(nonexact,metrics[5])}. REF families are known for {nonexact['family_any_compatible_evaluable']:,}/{nonexact['n']:,} non-exact retained paired rows. All retrieved candidates are annotated in {nonexact['all_candidates_annotated']:,}/{nonexact['n']:,} rows. Exclusive compatibility requires complete candidate annotation and every observed family to be REF-compatible.",'',
    'Exact peptide co-occurrence was examined across the complete acquisition candidate FASTA. Pre-existing family annotations cover a subset of this space. Unknown or conflicting assignments are never automatically called unrelated. Broad orthologous-group agreement does not identify a strain or establish the generating protein.','',
    'The REF search uses a community reference independent of FL selection, but contains surrogate strains. Agreement is an operational comparator conditional on REF identification. It does not establish whole-model recall or calibrated identification error. An additional spectrum-q sensitivity population is provided separately from the primary peptide-q rule.','',
    'Beam ranks, scores, bag coverage, source hashes and validation counts are retained. Missing beams cannot be reconstructed from rank-zero-only exports. Repeated rows and acquisitions are descriptive observations, not independent biological replicates.','',
    '[Figure](A03_near_misses.pdf) · [Acquisition summaries](acquisition_summary.tsv) · [Pooled summaries](pooled_summary.tsv) · [Validation and coverage](acquisition_validation.tsv). Each acquisition directory contains all prediction rows, query/reference-pair traces, lookup counts, reference retention and a checksum manifest.','',
    'Matthias’s k-mer and ambiguity expansion idea is an additional intervention. The archived CSVs contain sequence_with_bags; the FL extraction evaluated here remains exact/I-L. Retrieval gain, family compatibility, candidate inflation and identification must be measured separately before any benefit is claimed.']
    (results/'REPORT.md').write_text('\n'.join(lines)+'\n')
    dump(results/'COMPLETE.json',dict(status='PASS',n_acquisitions=15,prediction_rows=total,reference_paired_rows=paired,nonexact_retained=dict(nonexact),all_row_tables_reconciled=True,manifest_sha256=sha(run/'inputs/manifest.json'),outputs=[meta(p) for p in results.iterdir() if p.is_file() and p.name!='COMPLETE.json']))
