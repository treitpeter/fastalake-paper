import csv,json,math,os,pathlib,time
from collections import Counter
import numpy as np
from lib import *

def rows(path):
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))
def probability(query,frequencies,windows):
    logp=0.0
    for aa in query:
        f=frequencies.get(aa,0.0)
        if f<=0:return 0.0
        logp+=math.log(f)
    return -math.expm1(-windows*math.exp(logp))
def render(out,data,summary,seeds):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
    fig,(ax,bx)=plt.subplots(2,1,figsize=(10.6,7.2),sharex=True,gridspec_kw={'height_ratios':[3.2,1],'hspace':.12})
    lens=np.arange(9,51);real=[r for r in data if r['arm']=='REAL']
    rr=np.array([np.nan if r['distinct_queries']==0 else r['fraction_matched'] for r in real])
    ax.plot(lens,rr,'o-',ms=3.5,lw=1.6,color='#147c77',label='Authentic predictions')
    for seed in seeds:
        rs=[r for r in data if r['seed']==seed]
        ax.plot(lens,[np.nan if r['distinct_queries']==0 else r['fraction_matched'] for r in rs],color='#cd9a6b',alpha=.4,lw=.75,zorder=1)
    ax.fill_between(lens,[np.nan if r['scrambled_min'] is None else r['scrambled_min'] for r in summary],[np.nan if r['scrambled_max'] is None else r['scrambled_max'] for r in summary],color='#ce8841',alpha=.15)
    ax.plot(lens,[np.nan if r['scrambled_median'] is None else r['scrambled_median'] for r in summary],'o-',color='#be7329',ms=3,lw=1.5,label='Scrambled: median and range of five seeds')
    ax.plot(lens,[r['uniform20_reference'] for r in summary],':',color='#737b86',lw=1.4,label='Uniform 20 residues; no I/L equivalence')
    ax.plot(lens,[r['uniform20_il_median'] if r['uniform20_il_median'] is not None else np.nan for r in summary],'--',color='#608bb7',lw=1.4,label='Uniform residues; actual query composition + I/L')
    ax.plot(lens,[r['composition_median'] if r['composition_median'] is not None else np.nan for r in summary],'--',color='#8061a3',lw=1.6,label='Measured reservoir composition + I/L')
    ax.set_yscale('symlog',linthresh=1e-5,linscale=.45);ax.set_ylim(-1e-6,1.2)
    ax.set_yticks([0,1e-5,1e-4,1e-3,1e-2,1e-1,1]);ax.set_yticklabels(['0','0.00001','0.0001','0.001','0.01','0.1','1'])
    ax.set_ylabel('Fraction of distinct queries with ≥1 reservoir match')
    ax.set_title('Accidental lookup as a function of peptide length',loc='left',fontweight='bold',pad=15)
    ax.legend(frameon=False,fontsize=8,loc='upper right')
    ax.grid(axis='y',color='#e2e5e8',lw=.5)
    ax.text(.99,.025,'Log scale above 10⁻⁵; linear near zero',transform=ax.transAxes,ha='right',fontsize=8,color='#626b74')
    bx.plot(lens,[r['distinct_queries'] for r in real],color='#147c77',lw=1.5)
    ns=np.array([[r['distinct_queries'] for r in data if r['seed']==seed] for seed in seeds])
    bx.fill_between(lens,np.min(ns,axis=0),np.max(ns,axis=0),color='#ce8841',alpha=.2)
    bx.plot(lens,np.median(ns,axis=0),color='#be7329',lw=1.5)
    bx.set_yscale('symlog',linthresh=1);bx.set_ylabel('Queries tested');bx.set_xlabel('Peptide length (residues)');bx.set_xlim(8.5,50.5)
    bx.set_xticks([9,12,15,20,25,30,35,40,45,50])
    fig.text(.10,.02,'Worked HSTOOL acquisition · full reference reservoir · production top-30% filtering · 9–50 residues · I/L equivalence',fontsize=8,color='#626b74')
    fig.subplots_adjust(left=.12,right=.98,top=.92,bottom=.12)
    for ext in ['png','svg','pdf']:fig.savefig(out/('A02_lookup_by_length.'+ext),dpi=220,bbox_inches='tight',facecolor='white')
    plt.close(fig)
def build(run,m):
    out=run/'results';out.mkdir(exist_ok=True)
    assert not (out/'COMPLETE.json').exists(),'Completed report already exists'
    p=json.loads((run/'derived/PREPARED.json').read_text());lk=json.loads((run/'derived/LOOKUP.json').read_text())
    comp=json.loads((run/'derived/reservoir_composition.json').read_text())
    verify(lk['hits']);verify(comp['source_fasta'])
    assert comp['entries']==lk['stats']['proteins_searched']==233497702
    occurrences,blocks,_=hit_counts(lk['hits']['path'])
    lengths={int(k):int(v) for k,v in comp['length_histogram'].items()}
    windows={l:sum(max(k-l+1,0)*n for k,n in lengths.items()) for l in range(9,51)}
    f={k:v/comp['residues'] for k,v in comp['aa_counts'].items()};f['L']+=f['I'];f['I']=0
    uniform={c:1/20 for c in 'ACDEFGHIKLMNPQRSTVWY'};uniform['L']+=uniform['I'];uniform['I']=0
    data=[];queries=[];archive=[]
    union=set()
    for a in p['arms']:
        verify(a['queries']);qs=rows(a['queries']['path']);bylen={l:[] for l in range(9,51)};aset=set()
        for row in qs:
            q=row['query'];n=int(row['retained_rows']);length=int(row['length']);aset.add(q);union.add(q)
            rec=dict(arm=a['arm'],seed=a['seed'],query=q,length=length,retained_rows=n,matched=q in occurrences,match_occurrences=occurrences.get(q,0),uniform20_il_probability=probability(q,uniform,windows[length]),measured_composition_probability=probability(q,f,windows[length]))
            queries.append(rec);bylen[length].append(rec)
        if a['arm'] in ['REAL','SIM_20260901']:
            oldarm='REAL' if a['arm']=='REAL' else 'SIM';oldpath=run/'inputs'/(oldarm+'_sample_hits.tsv')
            original=next(x for x in m['sources'] if x['frozen']['path']==str(oldpath.resolve()));verify(original['frozen'])
            old_occ,_,_=hit_counts(oldpath)
            assert set(old_occ)<=aset,oldarm+' historical hit not present in audited query set'
            current=Counter({q:occurrences[q] for q in aset if q in occurrences})
            assert current==old_occ,oldarm+' full-reservoir matching differs from archived query/occurrence results'
            archive.append(dict(arm=a['arm'],matched_query_count=len(current),occurrence_counts_identical=True,archived_hits=meta(oldpath)))
        for length in range(9,51):
            rs=bylen[length];n=len(rs);matched=sum(r['matched'] for r in rs);nr=sum(r['retained_rows'] for r in rs);mr=sum(r['retained_rows'] for r in rs if r['matched'])
            data.append(dict(arm=a['arm'],seed=a['seed'],length=length,distinct_queries=n,matching_queries=matched,fraction_matched=matched/n if n else None,retained_rows=nr,matching_rows=mr,row_weighted_fraction=mr/nr if nr else None,within_protein_windows=windows[length],uniform20_reference=-math.expm1(-windows[length]*20.0**(-length)),uniform20_il_mean=float(np.mean([r['uniform20_il_probability'] for r in rs])) if n else None,measured_composition_mean=float(np.mean([r['measured_composition_probability'] for r in rs])) if n else None))
    assert set(occurrences)<=union
    table(out/'query_source_data.tsv',queries);table(out/'length_by_arm.tsv',data)
    summary=[]
    for l in range(9,51):
        rs=[r for r in data if r['length']==l and r['seed'] is not None and r['distinct_queries']>0]
        observed=[r['fraction_matched'] for r in rs]
        summary.append(dict(length=l,seeds_with_queries=len(rs),scrambled_min=min(observed) if rs else None,scrambled_median=float(np.median(observed)) if rs else None,scrambled_max=max(observed) if rs else None,uniform20_reference=-math.expm1(-windows[l]*20.0**(-l)),uniform20_il_median=float(np.median([r['uniform20_il_mean'] for r in rs])) if rs else None,composition_median=float(np.median([r['measured_composition_mean'] for r in rs])) if rs else None))
    table(out/'five_seed_summary.tsv',summary);dump(out/'ARCHIVED_REPRODUCTION.json',dict(status='PASS',comparisons=archive))
    dump(out/'theory_definition.json',dict(window_formula='sum(max(protein_length-L+1,0))',probability_approximation='1-exp(-W_L*product(residue_frequencies_in_query))',frequency_denominator='all sequence symbols, including any stop symbols; all counted explicitly',merged_il_frequencies=f,protein_entries=comp['entries'],sequence_residues=comp['residues'],windows_by_length=windows,assumptions=['residues independent under the approximation','Poisson occurrence approximation; overlapping windows and homologous proteins are dependent in the actual reservoir','scrambling is conditional on prediction composition; it is not a model of all real prediction errors']))
    render(out,data,summary,m['seeds'])
    def percent(x):return 'not estimable' if x is None else f'{100*x:.4g}%'
    lines=['# Analysis 2 — §4.2: lookup specificity by peptide length','',f"All six query populations were validated against the production Rust extractor. The historical seed 20260901 was reproduced byte-for-byte, and the fresh full-reservoir run reproduces the archived authentic and original scrambled per-query occurrence counts.",'',f"The freshly counted reservoir contains {comp['entries']:,} FASTA entries and {comp['residues']:,} sequence symbols. The full-query denominator is retained for every length and arm, including sequences with no hit.",'','| Arm | Length | Distinct queries | Matching queries | Match fraction |','|---|---:|---:|---:|---:|']
    for length in [9,12]:
        for r in data:
            if r['length']==length:lines.append(f"| {r['arm']} | {length} | {r['distinct_queries']:,} | {r['matching_queries']:,} | {percent(r['fraction_matched'])} |")
    lines+=['','## What this establishes','', 'This is a direct measurement of how often the tested authentic and composition-preserving scrambled sequences occur in the reservoir under the production lookup rules. It separates the event of matching at least one protein from the number of sequence occurrences. The complete curve and denominators determine how strongly length supports selectivity; no monotonic decline or particular threshold is assumed as an acceptance criterion.','', 'Five-seed envelopes describe variation between these fixed scramble realizations. They are not confidence intervals over acquisitions, biological samples or all possible prediction errors. Zero observed matches at a length are reported as zero events among a stated number of queries. A length with no queries has an undefined rate and is omitted from connected curves rather than represented as zero.','', 'The uniform 20-residue curve without I/L is the illustrative reference. The I/L-aware uniform and measured-composition curves are averaged over each arm’s actual query compositions. They use valid within-protein windows. Their independence and Poisson assumptions do not capture all effects of homologous proteins, repetitive windows or sequence correlations. Departures from those approximations remain visible.','', 'The scope is one fixed primary HSTOOL acquisition and five controls. It establishes neither a cohort-wide error rate nor the fraction of authentic AlphaNovo predictions that are incorrect. Realistic near-misses and final identification error are the subsequent §4.2 analyses.','', '## Figure caption','',f"**Length dependence of authentic and scrambled lookup.** Distinct I/L-normalized queries from the worked HSTOOL acquisition were matched by the production Rust extractor to the full {comp['entries']:,}-entry reservoir. All arms used the production positive-score, low-complexity, 9–50-residue and top-30% rules, retaining score ties at the cutoff. The authentic arm is shown alongside five composition-preserving scrambles (seeds {', '.join(map(str,m['seeds']))}). Thin control lines and shading show the five realizations and their minimum–maximum range; the central line is their median. Reference curves show uniform-residue and measured-composition approximations as specified in the source data. The probability axis is logarithmic above 10⁻⁵ and linear near zero, so observed zeros remain visible. The lower panel reports query denominators; a missing denominator yields an undefined fraction. No shared evidence prefilter or pooled DC predictions enter this experiment.",'', '[Figure](A02_lookup_by_length.pdf) · [All query outcomes](query_source_data.tsv) · [Per-length counts](length_by_arm.tsv) · [Five-seed summary](five_seed_summary.tsv) · [Archived reproduction](ARCHIVED_REPRODUCTION.json)']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    dump(out/'COMPLETE.json',dict(status='PASS',job_id=os.environ.get('SLURM_JOB_ID'),utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),manifest_sha256=sha(run/'inputs/manifest.json'),script_sha256=sha(__file__),n_arms=6,n_seeds=5,outputs=[meta(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='COMPLETE.json']))

