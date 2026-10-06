#!/usr/bin/env python3
"""Audit actual FASTA entries and render the approved four-cohort comparison."""
import csv, hashlib, json, os, pathlib, re, sys, time
import numpy as np

def sha(p):
    h=hashlib.sha256()
    with pathlib.Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def metadata(p):
    p=pathlib.Path(p);s=p.stat()
    return dict(path=str(p.resolve()),bytes=s.st_size,mtime_ns=s.st_mtime_ns)

def verify(r):
    assert metadata(r['path'])=={k:r[k] for k in ['path','bytes','mtime_ns']}, 'Changed input '+r['path']
    if 'sha256' in r:assert sha(r['path'])==r['sha256'], 'Changed content '+r['path']

def scan(path, unique=False, human=False):
    """Record counts are primary; optional exact sequence identities are a separate audit."""
    path=pathlib.Path(path);before=metadata(path);raw=hashlib.sha256()
    n=residues=length=0;header=None;parts=[];seqs=set();ids=set();duplicates=isoforms=nonhuman=nonreviewed=0
    def finish():
        if header is None:return
        if length==0:raise ValueError('Empty FASTA entry: '+str(path))
        if unique:seqs.add(hashlib.sha256(b''.join(parts).upper()).digest())
    with path.open('rb',buffering=8*1024*1024) as f:
        for line in f:
            raw.update(line)
            if line.startswith(b'>'):
                finish();header=line[1:].strip();length=0;parts=[];n+=1
                if not header:raise ValueError('Empty header: '+str(path))
                if unique:
                    ident=header.split()[0];duplicates+=int(ident in ids);ids.add(ident)
                    if ident.startswith((b'ENTRAP',b'DECOY',b'REV_')):raise ValueError('Control/decoy in selected target FASTA: '+str(path))
                if human:
                    fields=header.split(b'|');nonreviewed+=int(not header.startswith(b'sp|'))
                    isoforms+=int(len(fields)>1 and b'-' in fields[1]);nonhuman+=int(b'OX=9606' not in header)
            else:
                s=line.strip()
                if not s:continue
                if header is None:raise ValueError('Sequence before header: '+str(path))
                if b' ' in s or b'\t' in s:raise ValueError('Internal sequence whitespace: '+str(path))
                length+=len(s);residues+=len(s)
                if unique:parts.append(s)
    finish()
    assert metadata(path)==before, 'Changing FASTA: '+str(path)
    d=dict(**before,sha256=raw.hexdigest(),entries=n,residues=residues)
    if unique:d.update(unique_exact_sequences=len(seqs),duplicate_accession_count=duplicates,sequence_set_sha256=hashlib.sha256(b''.join(sorted(seqs))).hexdigest())
    if human:d.update(isoform_accessions=isoforms,nonhuman_headers=nonhuman,nonreviewed_headers=nonreviewed)
    return d

def write_json(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def table(path,rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

def render(out,rows,summaries,human,worked):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
    fig,(ax,bx)=plt.subplots(1,2,figsize=(12.8,5.5),gridspec_kw={'width_ratios':[1,1.8]})
    values=[worked[k]['entries'] for k in ['reservoir','candidates','selected']]
    ax.plot(range(3),values,color='#294f66',lw=1.5,zorder=1)
    ax.scatter(range(3),values,s=[100,75,55],c=['#879ba8','#407c91','#158d80'],zorder=2)
    for x,y in enumerate(values):ax.annotate(f'{y:,}',(x,y),xytext=(0,12),textcoords='offset points',ha='center',fontweight='bold')
    ax.set_yscale('log');ax.set_ylim(min(values)*0.25,max(values)*6)
    ax.set_xlim(-0.45,2.4);ax.set_xticks(range(3),['Reference\nreservoir','Exact-match\ncandidates','After\nparsimony'])
    ax.set_ylabel('FASTA entries (log scale)');ax.set_title('a   Worked stool acquisition',loc='left',fontweight='bold',pad=24)
    ax.text(0.02,-0.27,f"58-acquisition evidence filter: {worked['shared']['entries']:,} entries\nReservoir to selected database: {values[0]/values[2]:,.0f}-fold reduction",transform=ax.transAxes,fontsize=9,va='top')
    colors=['#7664a0','#168b80','#4587b8','#c68b3e'];order=['CAMPI','HSTOOL','PREDICT','CapScan']
    rng=np.random.default_rng(20260913)
    for x,(name,color) in enumerate(zip(order,colors)):
        ys=np.array([r['entries'] for r in rows if r['cohort']==name])
        bx.scatter(x+rng.uniform(-0.24,0.24,len(ys)),ys,s=12,alpha=0.58,color=color,edgecolors='none',zorder=2)
        med=float(np.median(ys));bx.plot([x-0.30,x+0.30],[med,med],color='#18212b',lw=2.2,zorder=3)
    bx.scatter([1],[worked['selected']['entries']],s=75,facecolors='none',edgecolors='#18212b',linewidths=1.2,zorder=4,label='Worked acquisition')
    bx.legend(loc='upper left',frameon=False,fontsize=8)
    bx.axhline(human['entries'],ls=(0,(4,3)),color='#8c4651',lw=1.3,zorder=1)
    bx.text(0.99,0.98,f"Human Swiss-Prot canonical: {human['entries']:,}",color='#8c4651',ha='right',va='top',transform=bx.transAxes,fontsize=9)
    bx.set_xticks(range(4),[f"{s['cohort']}\nn = {s['n']}\nR = {s['reservoir_entries']/1e6:.1f}M" for s in summaries])
    ymax=max(human['entries'],max(r['entries'] for r in rows));bx.set_ylim(0,ymax*1.20);bx.set_xlim(-0.6,3.6)
    bx.yaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:,.0f}'))
    bx.set_ylabel('Final selected FASTA entries');bx.set_title('b   All primary acquisitions',loc='left',fontweight='bold',pad=24)
    bx.text(0,-0.27,'Each point is an acquisition; black bars mark medians.\nR denotes the cohort’s own reservoir size.',transform=bx.transAxes,fontsize=9,va='top')
    fig.subplots_adjust(left=0.08,right=0.98,top=0.86,bottom=0.29,wspace=0.33)
    for ext in ['svg','pdf','png']:fig.savefig(out/('cohort_database_sizes.'+ext),dpi=220,facecolor='white',bbox_inches='tight')
    plt.close(fig)

def main(run):
    run=pathlib.Path(run).resolve();root=run.parents[1];out=run/'results'
    if (out/'COMPLETE.json').exists():raise RuntimeError('Completed analysis exists; do not overwrite')
    m=json.loads((run/'inputs/manifest.json').read_text());extra=json.loads((run/'inputs/audit_sources.json').read_text())
    for r in extra['sources']:verify(r['frozen'])
    for stage in ['gate','evidence']:
        assert json.loads((run/'status'/(stage+'.json')).read_text())['status']=='PASS'
    for stage in ['extract','select']:
        for i in range(58):assert json.loads((run/'status'/(stage+'_'+str(i)+'.json')).read_text())['status']=='PASS'
    verify(m['human_reference']['frozen'])
    human=scan(m['human_reference']['frozen']['path'],True,True)
    assert human['entries']>0 and human['isoform_accessions']==human['nonhuman_headers']==human['nonreviewed_headers']==0
    write_json(out/'human_reference.json',dict(label=m['human_reference_label'],**human))
    reservoirs=extra['reservoirs'];records=[];comparisons=[]
    labels={'CAMPI_MOCK':'CAMPI','HSTOOL':'HSTOOL','PREDICT':'PREDICT','CAPSCAN':'CapScan'}
    hroot=run/'rebuild/cohorts/HSTOOL'
    for a in m['acquisitions']:
        verify(a['fasta']);oldpath=pathlib.Path(a['fasta']['path']);p=oldpath
        if a['cohort']=='HSTOOL':p=hroot/'4_parsimony/razor_hash-acc'/a['sample']/'parsimony.fasta'
        new=scan(p,True)
        if a['cohort']=='HSTOOL':
            old=scan(oldpath,True)
            ns=json.loads((hroot/'3_persample'/a['sample']/'sample_lake.stats.json').read_text())
            os_=json.loads(pathlib.Path(extra['hstool_candidate_stats'][a['sample']]).read_text())
            comparisons.append(dict(sample=a['sample'],old_shared_pool_n=97,new_shared_pool_n=58,old_candidate_entries=os_['total_output_proteins'],new_candidate_entries=ns['total_output_proteins'],candidate_counts_equal=os_['total_output_proteins']==ns['total_output_proteins'],old_selected_entries=old['entries'],new_selected_entries=new['entries'],selected_counts_equal=old['entries']==new['entries'],selected_bytes_identical=old['sha256']==new['sha256'],selected_sequence_sets_identical=old['sequence_set_sha256']==new['sequence_set_sha256'],old_sha256=old['sha256'],new_sha256=new['sha256']))
        records.append(dict(cohort=labels[a['cohort']],sample=a['sample'],reservoir_entries=reservoirs[a['cohort']]['entries'],below_human_reference=new['entries']<human['entries'],**new))
    assert len(records)==710 and len(comparisons)==58
    table(out/'acquisition_sizes.tsv',records);table(out/'hstool_58_vs_97.tsv',comparisons)
    print('Counted all 710 primary acquisition FASTAs',flush=True)
    summaries=[]
    for name in ['CAMPI','HSTOOL','PREDICT','CapScan']:
        rs=[r for r in records if r['cohort']==name];ys=np.array([r['entries'] for r in rs]);n=len(ys)
        summaries.append(dict(cohort=name,n=n,reservoir_entries=rs[0]['reservoir_entries'],min=int(ys.min()),q25=float(np.quantile(ys,.25)),median=float(np.median(ys)),q75=float(np.quantile(ys,.75)),max=int(ys.max()),n_below=int(sum(ys<human['entries'])),fraction_below=float(np.mean(ys<human['entries'])),n_equal=int(sum(ys==human['entries'])),n_above=int(sum(ys>human['entries'])),human_reference_entries=human['entries']))
    table(out/'cohort_summary.tsv',summaries)
    worked={};sample=m['worked_sample']
    paths={'reservoir':pathlib.Path(m['reservoir']['path']),'shared':hroot/'1_evidence/evidence_lake.fasta','candidates':hroot/'3_persample'/sample/'sample_lake.fasta','selected':hroot/'4_parsimony/razor_hash-acc'/sample/'parsimony.fasta'}
    verify(m['reservoir'])
    for stage,p in paths.items():
        print('Recounting worked stage',stage,str(p),flush=True)
        worked[stage]=scan(p,unique=(stage=='selected'))
        write_json(out/('worked_'+stage+'.json'),worked[stage])
    shared_stats=json.loads((hroot/'1_evidence/evidence_lake.stats.json').read_text());candidate_stats=json.loads((hroot/'3_persample'/sample/'sample_lake.stats.json').read_text())
    assert worked['reservoir']['entries']==reservoirs['HSTOOL']['entries']==shared_stats['proteins_searched']
    assert worked['shared']['entries']==shared_stats['total_output_proteins']==candidate_stats['proteins_searched']
    assert worked['candidates']['entries']==candidate_stats['total_output_proteins']
    assert worked['selected']['entries']==next(r['entries'] for r in records if r['cohort']=='HSTOOL' and r['sample']==sample)
    write_json(out/'worked_example.json',dict(sample=sample,shared_pool_n=58,stages=worked,shared_extraction=shared_stats,acquisition_extraction=candidate_stats))
    render(out,records,summaries,human,worked)
    allbelow=sum(r['below_human_reference'] for r in records);allmin=min(r['entries'] for r in records);allmax=max(r['entries'] for r in records)
    lines=['# Analysis 1 — results for §4.1','',f"The worked stool acquisition reduces {worked['reservoir']['entries']:,} reservoir FASTA entries to {worked['candidates']['entries']:,} exact-match candidates and {worked['selected']['entries']:,} selected entries. These counts were read afresh from the files. Its shared intermediate contains {worked['shared']['entries']:,} entries and was built from the 58 primary HSTOOL acquisitions, excluding all 39 DC acquisitions.",'',f"Across the four primary cohorts, {allbelow}/710 acquisitions ({100*allbelow/710:.1f}%) have fewer selected entries than the {human['entries']:,}-entry archived human reference. The complete acquisition range is {allmin:,}–{allmax:,}.",'','| Cohort | Acquisitions | Median entries | Full range | Below human reference |','|---|---:|---:|---:|---:|']
    for s in summaries:lines.append(f"| {s['cohort']} | {s['n']} | {s['median']:,.1f} | {s['min']:,}–{s['max']:,} | {s['n_below']}/{s['n']} ({100*s['fraction_below']:.1f}%) |")
    nbytes=sum(c['selected_bytes_identical'] for c in comparisons);nsets=sum(c['selected_sequence_sets_identical'] for c in comparisons);ncandidates=sum(c['candidate_counts_equal'] for c in comparisons)
    lines+=['','## Effect of excluding DC evidence','',f"Candidate counts agree with the archived 97-pool run in {ncandidates}/58 acquisitions. Final selected sequence sets agree in {nsets}/58, and complete FASTA bytes agree in {nbytes}/58. Equal candidate counts alone are not proof of candidate-set identity; the final-sequence comparison uses exact uppercase sequence hashes, without I/L collapsing.",'','## Interpretation and counting definitions','', 'The result measures the size of the database made available to a conventional search. It does not count identified proteins or infer how many distinct protein species the sample expresses. The reference is the archived UniProtKB/Swiss-Prot reviewed human canonical FASTA, with no isoform-suffixed accessions; its exact file and checksum are recorded in human_reference.json. A formal UniProt release identifier has not been inferred from the local archive filename.','',f"The human file contains {human['entries']:,} entries and {human['unique_exact_sequences']:,} distinct exact sequences. The plotted endpoint counts entries consistently on both sides. Distinct sequence counts for each selected FASTA are provided separately. The large reservoir recount verifies FASTA records and residues; its sequence-deduplicated construction is documented by the original builder, rather than re-established through a new 233-million-sequence uniqueness computation.",'', 'The strongest defensible heading must follow the complete distributions. If any acquisition exceeds the human reference, a universal claim that every final database is smaller would be false. The worked-example claim and the cohort-specific fractions remain directly assessable. Acquisition-weighted summaries are descriptive; repeated acquisitions are not treated as independent biological samples.','', 'CAMPI, PREDICT and CapScan retain their archived production databases and their own reservoir sizes. Only HSTOOL was rebuilt under the newly agreed 58-acquisition restriction. The study-wide merged search databases are outside this analysis.','', '## Figure caption','', f"**A de novo lookup produces compact acquisition databases from broad sequence reservoirs.** **a**, Worked HSTOOL acquisition ({sample}), with freshly counted reservoir, acquisition-candidate and final selected FASTA entries. The shared evidence filter was rebuilt using only the 58 primary HSTOOL acquisitions; its intermediate size is stated below the panel. The vertical axis is logarithmic. **b**, Final selected FASTA entries for CAMPI (15), HSTOOL (58), PREDICT (306) and CapScan (331), with every acquisition shown, cohort medians marked in black, the worked acquisition circled, and the {human['entries']:,}-entry reviewed canonical human Swiss-Prot reference as a dashed line. R denotes the actual reservoir size for each cohort. Point jitter is deterministic and has no analytical meaning. All acquisitions are retained, including those above the reference.",'','[Figure](cohort_database_sizes.pdf) · [All acquisition counts](acquisition_sizes.tsv) · [Cohort summaries](cohort_summary.tsv) · [58-versus-97 comparison](hstool_58_vs_97.tsv)']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    outputs=[dict(**metadata(p),sha256=sha(p)) for p in sorted(out.iterdir()) if p.is_file() and p.name!='COMPLETE.json']
    write_json(out/'COMPLETE.json',dict(status='PASS',job_id=os.environ.get('SLURM_JOB_ID'),utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),n_acquisitions=710,n_hstool=58,n_excluded_dc=39,manifest_sha256=sha(run/'inputs/manifest.json'),audit_sources_sha256=sha(run/'inputs/audit_sources.json'),script_sha256=sha(__file__),outputs=outputs))
    print('AUDIT COMPLETE',flush=True)

if __name__=='__main__':main(sys.argv[1])
