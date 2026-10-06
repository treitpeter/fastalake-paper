#!/usr/bin/env python3
import ast,json,pathlib,shutil,subprocess,time
from lib import *
RUN=pathlib.Path(__file__).resolve().parent
ROOT=RUN.parents[1]
def freeze(src,dst):
    dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(str(src),str(dst))
    a,b=meta(src),meta(dst);assert a['sha256']==b['sha256']
    return dict(source=a,frozen=b)
def main():
    assert not (RUN/'inputs/manifest.json').exists(),'Manifest already frozen'
    for d in ['inputs','frozen','bin','jobs','logs','status','derived','results']:(RUN/d).mkdir(parents=True,exist_ok=True)
    a1=ROOT/'MM_review_analysis/A01_4_1_database_sizes'
    m1=json.loads((a1/'inputs/manifest.json').read_text());sample=m1['worked_sample']
    old=ROOT/'canonical/publication_data/arms/fig1b_shuffle_control_2026-09-01'
    real=freeze(old/'REAL_pep'/(sample+'_predictions.csv'),RUN/'inputs/authentic.csv')
    archived=freeze(old/'SIM_pep'/(sample+'_predictions.csv'),RUN/'inputs/archived_seed20260901.csv')
    pred=next(x for x in m1['predictions'] if x['sample']==sample)
    assert real['frozen']['sha256']==pred['frozen']['sha256'],'Different worked acquisition predictions'
    b=next(x for x in m1['binaries'] if x['source']['path'].endswith('/fasta_extractor'))
    verify(b['source']);verify(b['frozen'])
    binary=freeze(b['frozen']['path'],RUN/'bin/fasta_extractor')
    sources=[freeze(a1/'inputs/manifest.json',RUN/'inputs/analysis1_manifest.json')]
    for p,name in [(old/'acquisition.json','archived_acquisition.json'),(ROOT/'FIGURATOR/Figure1/Fig1B/shuffle_one.py','archived_shuffle_one.py'),(ROOT/'FIGURATOR/Figure1/Fig1B/METHODS_fig1b.md','archived_methods.md'),(ROOT/'FIGURATOR/Figure1/Fig1B/arm.sbatch','archived_arm.sbatch')]:
        sources.append(freeze(p,RUN/'inputs'/name))
    for arm in ['REAL','SIM']:
        for name in ['sample_lake.stats.json','sample_hits.tsv']:
            sources.append(freeze(old/arm/name,RUN/'inputs'/(arm+'_'+name)))
    for p in (a1/'frozen/source/fasta_extractor_v2').rglob('*'):
        if p.is_file():sources.append(freeze(p,RUN/'frozen/source'/p.relative_to(a1/'frozen/source')))
    code=[]
    for p in sorted(RUN.glob('*.py')):
        ast.parse(p.read_text());code.append(freeze(p,RUN/'frozen'/p.name))
    code.append(freeze(RUN/'composition.cpp',RUN/'frozen/composition.cpp'))
    compiler=shutil.which('g++');assert compiler,'g++ missing'
    argv=[compiler,'-O3','-std=c++11',str(RUN/'frozen/composition.cpp'),'-o',str(RUN/'bin/composition')]
    subprocess.run(argv,check=True)
    # The FASTA census must respect multiline sequences and protein boundaries.
    toy=RUN/'inputs/census_fixture.fasta';toy.write_text('>one\nACIL\nGG\n>two\nK\n')
    subprocess.run([str(RUN/'bin/composition'),str(toy),str(RUN/'status/census_fixture')],check=True)
    st=json.loads((RUN/'status/census_fixture.json').read_text())
    assert st['entries']==2 and st['residues']==7 and st['aa_counts']['I']==1 and st['aa_counts']['L']==1 and st['length_histogram']=={'1':1,'6':1},st
    manifest=dict(analysis='Analysis 2 / section 4.2 / register A04',approved_utc_date='2026-09-13',sample=sample,seeds=[20260901,20260913,20260914,20260915,20260916],authentic=real,archived_scramble=archived,extractor=binary,extractor_original_provenance=b,canonical_commit=m1['canonical_commit'],sources=sources,code=code,composition_binary=meta(RUN/'bin/composition'),composition_compile_command=argv,reservoir=m1['reservoir'],settings=dict(min_length=9,max_length=50,top_fraction=.30,normalize_il=True,score_ties='all at cutoff',low_complexity='production'),prepared_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
    dump(RUN/'inputs/manifest.json',manifest)
    dump(RUN/'status/PREFLIGHT.json',dict(status='PASS',manifest=meta(RUN/'inputs/manifest.json'),checks=['same worked predictions as A01','production binary matches frozen committed-source provenance','Python syntax','multiline reservoir census fixture']))
    print('FROZEN',RUN/'inputs/manifest.json',flush=True)
if __name__=='__main__':main()

