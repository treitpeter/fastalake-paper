#!/usr/bin/env python3
import ast,csv,json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];MM=ROOT/'MM_review_analysis'
sys.path.insert(0,str(MM/'code'))
from setup_review_followups import submit,meta,sha,dump
run=MM/'A15_4_4_transfer_pool'
assert not (run/'inputs/manifest.json').exists(),'Already prepared; inspect journal'
for d in ['inputs','frozen','results','jobs','logs','status']:(run/d).mkdir(parents=True,exist_ok=True)
(run/'README.md').write_text('''# §4.4 — transfer and other-eleven pooling audit

R16/register A10. Independently recount all 58 original target acquisitions with
their own, within-design-group and interindividual-donor databases (174 saved
searches). Recount all 12 technical targets with own, all-twelve and other-eleven
final-FASTA unions (36 saved searches). No new identification search is needed.

Verify matched target mzML identity and all search parameters; verify each
transferred FASTA equals its named donor FASTA and each pool equals its explicit
source membership. Compare reconstructed I/L-canonical target peptide q≤0.01
sets to archived identity files and the corrected 58-row manuscript recount.
Keep 5–50 primary and nested 9–50 populations, distinct gained/lost identities,
counts, and two different retention measures. Independently scan each searched
FASTA for every own-arm peptide's sequence and compatible Sage KR/P digest.
Distinguish sequence absence, incompatible cleavage, and available-but-unidentified.

The other-eleven control excludes the target final FASTA; the upstream cohort
evidence lake remains shared. It does not establish whole-workflow independence.
T1 in the interindividual group is already between donors. T4 draws from the
interindividual group; donor relationships outside that group are not asserted.
Technical injections and preparations are not independent biological participants.
Counts are conditional identification results, not a new FDR calibration.

A transfer and pooling pilot gate the remaining 68 acquisition/experiment audits.
Two concurrent 1-CPU/8G jobs follow the public-comparison jobs, preserving bounded
concurrency. Freeze all scripts and source identities before execution. Preserve
all original data; scientific/visual review and DOCX integration follow validation.
''')
sources=[(MM/'code/transfer_pool_audit.py','search.py'),(MM/'code/public_comparison.py','parameters.py'),
    (MM/'A14_4_3_campi_gap/frozen/search.py','digest.py'),(MM/'A02_4_2_lookup_specificity/lib.py','lib.py'),
    (ROOT/'canonical/publication/lib/peptide.py','peptide.py')]
for src,name in sources:shutil.copy2(src,run/'frozen'/name);ast.parse((run/'frozen'/name).read_text())
stool=ROOT/'canonical/final_revision/thread_selection_math/stool_v2'
pairs=list(csv.DictReader((stool/'PAIRS.tsv').open(),delimiter='\t'));assert len(pairs)==174
old=ROOT/'FIGURATOR/Figure3/edit3_pages_2026-09-09/source_data/03_source_data/review_repairs_2026-09-06/transfer/per_acquisition.tsv'
assert sha(old)=='2431b30bafb08fefe3f53709c2d49b7515b9422417f14fd8eeaa17f21968a691'
pool=ROOT/'FIGURATOR/Figure3/review_2026-09-07/pooling';pm=json.loads((pool/'MANIFEST.json').read_text())
records=[]
def arm(p,historical,extra):
    c=json.loads((p/'results.json').read_text());f=Path(c['database']['fasta']).resolve()
    return dict(config=meta(p/'results.json'),search=meta(p/'results.sage.tsv',False),fasta=meta(f),historical_set=meta(historical),**extra),c
sourcepath=ROOT/'canonical/publication_data/arms/noclust_2026-08-07/cohorts/HSTOOL/4_parsimony/razor_hash-acc'
for sample in sorted({r['target'] for r in pairs}):
    pp=[r for r in pairs if r['target']==sample];assert len(pp)==3
    arms={};sourcefastas={};raw=None
    for r in pp:
        tier=r['tier'];a,c=arm(stool/tier/sample,stool/'sets'/(sample+'__'+tier+'.pep'),dict(donor=r['donor']))
        arms['own' if tier=='T0_self' else tier]=a
        sourcefastas[r['donor']]=meta(sourcepath/r['donor']/'parsimony.fasta')
        mz=Path(c['mzml_paths'][0]).resolve();assert raw is None or raw==mz;raw=mz
    records.append(dict(id='transfer__'+sample,kind='transfer',sample=sample,group=pp[0]['target_group'],raw=str(raw),arms=arms,source_fastas=sourcefastas))
sourcefastas={r['sample']:meta(r['fasta']) for r in pm['source_fastas']}
for r in pm['source_fastas']:assert sourcefastas[r['sample']]['sha256']==r['fasta_sha256']
poolrecords=[]
for sample in sorted(sourcefastas):
    arms={};raw=None
    for r in [r for r in pm['tasks'] if r['sample']==sample]:
        p=Path(r['config']).parent;a,c=arm(p,p/'accepted_peptides.txt',dict(donors=r['donor_samples']))
        assert sha(r['config'])==r['config_sha256'] and a['fasta']['sha256']==r['database_sha256']
        arms[r['arm']]=a;mz=Path(c['mzml_paths'][0]).resolve();assert raw is None or raw==mz;raw=mz
    assert len(arms)==3
    poolrecords.append(dict(id='pool__'+sample,kind='pool',sample=sample,group='TechnicalReplicates',raw=str(raw),arms=arms,source_fastas=sourcefastas))
tasks=[records[0],poolrecords[0]]+records[1:]+poolrecords[1:];assert len(tasks)==70
dump(run/'inputs/manifest.json',dict(analysis=run.name,register=['A10'],tasks=tasks,
    historical_transfer=meta(old),historical_pool=meta(pool/'per_acquisition.tsv'),source_pairs=meta(stool/'PAIRS.tsv'),source_pool=meta(pool/'MANIFEST.json'),
    frozen_code=[meta(p) for p in (run/'frozen').glob('*.py')],protocol=meta(run/'README.md')))
dump(run/'status/PREFLIGHT.json',dict(status='PASS',transfer_targets=58,pool_targets=12,searches=210,
    source_rosters_and_configs_frozen=True,raw_recount_and_digest_validation='pilot-gated'))
prior=json.loads((MM/'A16_4_5_public_comparison/status/SUBMITTED.json').read_text())['summary']
pilot=submit(run,'pilot','search.py',indices='0-1',dependency=prior,cap=2)
full=submit(run,'audit','search.py',indices='2-69',dependency=pilot,cap=2)
dump(run/'status/SUBMITTED.json',dict(status='SUBMITTED_NOT_COMPLETED',pilot=pilot,audit=full,tasks=70,searches=210))
