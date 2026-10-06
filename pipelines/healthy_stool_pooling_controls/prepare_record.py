from pathlib import Path
import shutil,json,csv,copy,hashlib
P=Path('LAB_ROOT/FastaLake');B=P.parent/'BiorXiv';H=P/'FIGURATOR/Figure3/review_2026-09-07';L=B/'05_references/introduction_planning_2026-09-07';L.mkdir(parents=True,exist_ok=True)
shutil.copytree(H/'literature',L/'verified_figure3_literature',dirs_exist_ok=True)
shutil.copy2(H/'FIGURE3_REVIEW.md',L/'FIGURE3_REVIEW_at_conception.md')
(L/'INTRODUCTION_EVIDENCE.md').write_text('''# Introduction evidence saved 7 September 2026

Preserved for the later Introduction revision; the Introduction has not been rewritten.

The opening should acknowledge established stool metaproteomic reproducibility and longitudinal profiling (Gonzalez et al., 2020), conserved and variable functions in healthy people with taxon-specific contributions (Tanca et al., 2024), and functional persistence despite taxonomic/protein variation (Blakeley-Ruiz et al., 2019; Crohn’s remission, not healthy controls).

Our contribution is an evidence-guided search workflow accompanied by a reusable benchmark: transfer and pooling of selected sequence evidence on unchanged target spectra, five distinct measurement/sampling designs, and functional-contributor variation interpreted against technical variation. Do not claim that healthy stool variability or functional redundancy is newly discovered, that 58 acquisitions are 58 people, or that unmatched metrics establish superiority. The repeated-sampling dates still need confirmation.

The technical other-eleven pool raises median canonical peptide counts by 9.1% and retains 94.9% of own-search peptide identities. These descriptive results support final-database reuse. Shared upstream cohort evidence prevents a whole-procedure independence claim. Entrapment outcomes assess control behaviour; they do not establish nominal 1% whole-procedure FDR calibration.

## Verified primary sources

- Gonzalez et al. (2020), https://doi.org/10.1128/mSystems.00200-20 — 145 specimens, 29 participants, five time points, duplicate processing; technical/preparation reproducibility precedent. Reported R² and our LFQ Pearson r are not interchangeable.
- Tanca et al. (2024), https://doi.org/10.1038/s41522-024-00526-4 — ten datasets, 134 healthy individuals; 182 functions and 83 pathways detected in all individuals. Direct precedent for the healthy community functional view.
- Blakeley-Ruiz et al. (2019), https://doi.org/10.1186/s40168-019-0631-8 — 25 metaproteomes from five Crohn’s-remission participants over one year; functional persistence alongside taxonomic/protein variability.

Full primary texts, short verified excerpts, and comparison limits are retained in verified_figure3_literature/. These are selected verified precedents, not an exhaustive priority search. Figure3 review details and benchmark provenance are in the accompanying saved review and the FIGURATOR/Figure3/review_2026-09-07 directory.
''')
E=H/'entrapment_pooling';E.mkdir(exist_ok=True);(E/'logs').mkdir(exist_ok=True);(E/'databases').mkdir(exist_ok=True)
R=P/'canonical/publication/analyses/full_route_entrapment_2026-09-06';old=json.loads((H/'pooling/MANIFEST.json').read_text());samples=[r['sample'] for r in old['source_fastas']];sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
# Import definitions only, not the original experiment preparation.
ns={};exec((H/'pooling/prepare.py').read_text().split('rows=')[0].replace('H=Path(__file__).resolve().parent','H=Path("'+str(H/'pooling')+'")'),ns);read=ns['read_fasta']
classes={r['accession']:r['class'] for r in csv.DictReader((R/'runs/HSTOOL/1_evidence/controls.tsv').open(),delimiter='\t')};(E/'control_classes.json').write_text(json.dumps(classes,indent=2))
dbs={};sources=[];members={}
for s in samples:
 f=R/'runs/HSTOOL/samples'/s/'4_razor'/f'{s}_razor.fasta';dbs[s]=read(f);sources.append({'sample':s,'fasta':str(f),'sha256':sha(f),'n':len(dbs[s])})
 for a,v in dbs[s].items():
  if a not in members:members[a]={'record':v,'donors':[]}
  assert members[a]['record']==v
  members[a]['donors'].append(s)
def write(p,keys):
 with p.open('x') as f:
  for a in sorted(keys):h,q=members[a]['record'];f.write('>'+h+'\n'+q+'\n')
tasks=[]
for i,s in enumerate(samples):
 for arm,donors in [('pool_all12',samples),('pool_other11',[x for x in samples if x!=s])]:
  keys={a for a,v in members.items() if set(v['donors'])&set(donors)};fp=E/'databases'/('all12.fasta' if arm=='pool_all12' else f'other11_{i:02d}.fasta')
  if not fp.exists():write(fp,keys)
  out=E/'search'/arm/s;out.mkdir(parents=True,exist_ok=True)
  orig=R/'runs/HSTOOL/samples'/s/'5_sage/config.json';cfg=json.loads(orig.read_text());cfg['database']['fasta']=str(fp);cfg['output_directory']=str(out);cp=out/'config.json';cp.write_text(json.dumps(cfg,indent=2))
  tasks.append({'task':len(tasks),'sample':s,'arm':arm,'config':str(cp),'config_sha256':sha(cp),'database':str(fp),'database_sha256':sha(fp),'donors':donors,'target_in_final_pool':s in donors,'n_accessions':len(keys),'control_accessions':{c:sum(classes.get(a)==c for a in keys) for c in ['shuffle','plant','thermophilic']},'original_config':str(orig),'original_config_sha256':sha(orig)})
shutil.copy2(R/'src/full_route_entrapment.py',E/'full_route_entrapment.py')
(E/'MANIFEST.json').write_text(json.dumps({'sage':old['sage'],'sage_sha256':old['sage_sha256'],'sources':sources,'tasks':tasks,'n_tasks':24,'full_route_parent':str(R),'parent_experiment_sha256':sha(R/'experiment.json'),'scope':'Controls entered before cohort and acquisition evidence selection and razor in the parent experiment. Pooling reuses those final databases. Other11 excludes target final FASTA but shares upstream cohort evidence. No calibrated whole-procedure FDR estimate.'},indent=2))
print('Introduction evidence:',L);print('Control pooling tasks:',len(tasks));print(tasks[0]['control_accessions'])
