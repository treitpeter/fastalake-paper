"""Require a complete cohort, collect control traces and scan target overlap."""
from pathlib import Path
import argparse,collections,csv,json,statistics,subprocess,sys
import pandas as pd
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'src'))
from full_route_entrapment import canonical,write_json,sha

def write_tsv(path,rows,fields=None):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def baseline_set(path):
    peps=set()
    for d in pd.read_csv(path,sep='\t',usecols=['peptide','label','peptide_q'],chunksize=100000):
        peps.update(canonical(p) for p in d.loc[(d.label==1)&(d.peptide_q<=.01),'peptide'])
    return peps

def main():
    ap=argparse.ArgumentParser();ap.add_argument('cohort');a=ap.parse_args()
    cfg=json.loads((R/'experiment.json').read_text());c=cfg['cohorts'][a.cohort];root=Path(cfg['output_root'])/a.cohort
    samples=list(csv.DictReader(open(c['sample_list']),delimiter='\t'));out=R/'tables'/a.cohort;out.mkdir(exist_ok=True)
    ev=json.loads((root/'1_evidence/COMPLETE.json').read_text());assert ev['n_pooled']==c['n']
    traces=[];summary=[];discoveries=[];queries=set();comparisons=[]
    for r in samples:
        s=r['sample_id'];p=root/'samples'/s
        done=json.loads((p/'COMPLETE.json').read_text());assert done['sample']==s
        for stage,v in [('0_control_bank',{'counts':cfg['control_class_counts']}),('1_cohort',ev),('3_acquisition',done['selection']),('4_razor',done['razor'])]:
            for cls in ['shuffle','thermophilic','plant']:
                traces.append(dict(cohort=a.cohort,sample=s,stage=stage,control_class=cls,n_proteins=v['counts'].get(cls,0)))
        for row in csv.DictReader((p/'5_sage/summary.tsv').open(),delimiter='\t'):
            summary.append(dict(cohort=a.cohort,sample=s,**row))
        current=set()
        for row in csv.DictReader((p/'5_sage/discoveries.tsv').open(),delimiter='\t'):
            if float(row['q_threshold'])==.01 and 'target' in row['classes'].split(';'):current.add(row['peptide'])
            if row['assignment']=='target':continue
            discoveries.append(dict(cohort=a.cohort,sample=s,**row));queries.add(row['peptide'])
        old=baseline_set(r['baseline_results'])
        comparisons.append(dict(cohort=a.cohort,sample=s,baseline_canonical=len(old),full_route_target_supported=len(current),intersection=len(old&current),lost=len(old-current),gained=len(current-old),count_ratio=len(current)/len(old) if old else None,own_set_recovery=len(old&current)/len(old) if old else None))
    assert len(comparisons)==c['n']
    qpath=out/'control_queries.txt';qpath.write_text(''.join(p+'\n' for p in sorted(queries)))
    result=out/'original_reservoir_overlap.tsv'
    subprocess.run([str(R/'src/scan_overlap'),str(qpath),c['reservoir'],str(result),'16'],check=True)
    overlap={r['peptide']:int(r['present_in_original_reservoir']) for r in csv.DictReader(result.open(),delimiter='\t')}
    assert set(overlap)==queries
    for r in discoveries:r['present_in_original_reservoir']=overlap[r['peptide']]
    by=collections.defaultdict(list)
    for r in discoveries:by[r['sample'],float(r['q_threshold'])].append(r)
    for r in summary:
        n=int(r['n_canonical']);floor=int(r['min_residues']);rr=[x for x in by[r['sample'],float(r['q_threshold'])] if len(x['peptide'])>=floor]
        for cls in ['shuffle','plant','thermophilic']:
            only=[x for x in rr if x['assignment']==cls]
            novel=sum(not x['present_in_original_reservoir'] for x in only)
            r[cls+'_only_not_in_original_reservoir']=novel
            r[cls+'_only_not_in_original_reservoir_fraction']=novel/n if n else None
    write_tsv(out/'stage_control_counts.tsv',traces)
    write_tsv(out/'per_acquisition.tsv',summary)
    write_tsv(out/'control_discoveries.tsv',discoveries,['cohort','sample','peptide','q_threshold','assignment','classes','present_in_original_reservoir'])
    write_tsv(out/'baseline_comparison.tsv',comparisons)
    groups=collections.defaultdict(list)
    for r in summary:groups[float(r['q_threshold']),int(r['min_residues'])].append(r)
    medians=[]
    for (q,floor),rr in sorted(groups.items()):
        for cls in ['shuffle','thermophilic','plant']:
            for metric in [cls+'_only_fraction',cls+'_only_not_in_original_reservoir_fraction']:
                v=[float(r[metric]) for r in rr if r[metric] not in ['',None]]
                qu=statistics.quantiles(v,method='exclusive') if len(v)>1 else [None]*3
                medians.append(dict(cohort=a.cohort,q_threshold=q,min_residues=floor,control_class=cls,metric=metric,n_attempted=len(rr),n_defined=len(v),median=statistics.median(v) if v else None,q1=qu[0],q3=qu[2],max=max(v) if v else None))
    write_tsv(out/'summary.tsv',medians)
    write_json(out/'COMPLETE.json',{'cohort':a.cohort,'n':len(comparisons),'all_stages_and_searches_complete':True,'original_reservoir_queries_scanned':len(queries),'observed_control_fraction_only':True,'paired_estimate_reported':False,'whole_procedure_calibration_established':False,'sources':{p.name:sha(p) for p in out.glob('*.tsv')}})
    print(json.dumps([r for r in medians if r['q_threshold']==.01 and r['min_residues']==5],indent=2))

if __name__=='__main__':main()
