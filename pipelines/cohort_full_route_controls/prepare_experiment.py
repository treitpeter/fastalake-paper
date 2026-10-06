"""Freeze a full-route experiment after successful control generation."""
from pathlib import Path
import csv,hashlib,json,shutil,subprocess,sys
R=Path(__file__).resolve().parents[1];CAN=R.parents[2]
sys.path.insert(0,str(CAN))
from fasta_lake.full_route_entrapment import sha,fasta,write_json
F=CAN.parent
ARCH=F/'FIGURATOR/Supplementary/ENTRAPMENT/rerun_2026-09-03'
stats=json.loads((R/'controls/shuffle_stats.json').read_text())
assert stats['all_written_cleavage_coordinates_verified'] and stats['all_written_block_compositions_verified'] and stats['written']>100000
controls=R/'controls/all_controls.fasta';mapping=R/'controls/all_controls.tsv'
if controls.exists():raise FileExistsError('Frozen controls already assembled')
shutil.copyfile(R/'controls/shuffle.fasta',controls);shutil.copyfile(R/'controls/shuffle_map.tsv',mapping)
counts={'shuffle':stats['written']}
foreign={'thermophilic':Path('LAB_ROOT/pyrococcus_furiosus_swissprot.fasta'),'plant':Path('LAB_ROOT/arabidopsis_thaliana_swissprot.fasta')}
foreign_manifest=[]
with controls.open('a') as out,mapping.open('a') as mf:
    writer=csv.writer(mf,delimiter='\t',lineterminator='\n')
    for cls,src in foreign.items():
        counts[cls]=0
        frozen=R/'inputs'/(cls+'.fasta');shutil.copy2(src,frozen)
        foreign_manifest.append({'class':cls,'origin':str(src),'sha256':sha(frozen)})
        for acc,s in fasta(frozen):
            name='FLX'+hashlib.sha256((cls+'\0'+acc+'\0'+s).encode()).hexdigest()[:32]
            out.write('>'+name+'\n'+s+'\n');writer.writerow([name,cls,acc,len(s),'reference']);counts[cls]+=1
assert counts['thermophilic']==512 and counts['plant']==16418
config={'experiment':'full_route_entrapment_2026-09-06','output_root':str(R/'runs'),'controls_fasta':str(controls),'controls_map':str(mapping),'controls_sha256':sha(controls),'control_sequences':sum(counts.values()),'control_class_counts':counts,'control_source_stats':stats,'foreign_sources':foreign_manifest,'binaries':{},'cohorts':{}}
binarypaths={'extractor':CAN/'rust/fasta_extractor_v2/target/release/fasta_extractor','parsimony':CAN/'rust/parsimony_engine/target/release/parsimony_engine','sage':Path('LAB_ROOT/sage')}
(R/'bin').mkdir(exist_ok=True)
for key,p in binarypaths.items():
    if key!='sage':
        crate=p.parents[2]
        state=subprocess.check_output(['git','-C',str(CAN),'status','--porcelain','--',str(crate)],text=True)
        assert not state,'Production crate has changed: '+state
    dest=R/'bin'/p.name;shutil.copy2(p,dest)
    config['binaries'][key]={'path':str(dest),'origin':str(p),'sha256':sha(dest),'version':subprocess.check_output([str(dest),'--version'],text=True).strip()}
for cohort,n in [('CAMPI_MOCK',15),('HSTOOL',58),('PREDICT',306),('CAPSCAN',331)]:
    rows=list(csv.DictReader((ARCH/f'lists/{cohort}_ruled.tsv').open(),delimiter='\t'));assert len(rows)==n
    for row in rows:
        baseline=CAN/f"publication_data/arms/noclust_2026-08-07/cohorts/{cohort}/5_sage/razor_hash-acc/{row['sample_id']}/config.json"
        assert baseline.is_file() and Path(row['denovo_csv']).is_file() and Path(row['mzml']).is_file()
        frozen_cfg=R/'inputs/configs'/cohort/(row['sample_id']+'.json');frozen_cfg.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(baseline,frozen_cfg)
        row['sage_config']=str(frozen_cfg)
        row['baseline_results']=str(baseline.parent/'results.sage.tsv')
        st=Path(row['denovo_csv']).stat();row['denovo_bytes']=st.st_size;row['denovo_mtime_ns']=st.st_mtime_ns
    path=R/'inputs'/(cohort+'_samples.tsv')
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    ev=json.loads((CAN/f'publication_data/arms/noclust_2026-08-07/cohorts/{cohort}/1_evidence/evidence_lake.stats.json').read_text())
    reservoir=CAN/f'publication_data/cohorts/{cohort}/0_lake/{cohort}_lake.fasta'
    assert Path(ev['database']).resolve()==reservoir.resolve()
    config['cohorts'][cohort]={'sample_list':str(path),'sample_list_sha256':sha(path),'n':n,'reservoir':str(reservoir),'reservoir_sequences':ev['proteins_searched']}
config['canonical_commit']=subprocess.check_output(['git','-C',str(CAN),'rev-parse','HEAD'],text=True).strip()
module=CAN/'fasta_lake/full_route_entrapment.py';shutil.copy2(module,R/'src/full_route_entrapment.py')
config['runner_sha256']=sha(module)
write_json(R/'experiment.json',config)
print(json.dumps({'control_counts':counts,'cohorts':list(config['cohorts']),'frozen_binaries':config['binaries']},indent=2))
