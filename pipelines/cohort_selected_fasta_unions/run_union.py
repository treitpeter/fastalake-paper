"""Build the exact-sequence union with the frozen workflow.build_union (byte-identical copy in ../code).
Modes: baseline331 = the production CAPSCAN_noclust task FASTAs; fixed347 = those + the 16 participant-15 parsimony.fasta;
negctl = fixed347 with one planted wrong byte count (verify() must raise -> exit 1)."""
import sys,json,csv,shutil
from pathlib import Path
O=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(O/'code'))
import workflow as W
mode=sys.argv[1]
SRC=Path('LAB_ROOT/sample_specific_merge_all_raw/2026-09-11')
m=json.loads((SRC/'manifest.json').read_text())
t=[x for x in m['tasks'] if x['id']=='CAPSCAN_noclust'][0]
fastas=list(t['fastas']);assert len(fastas)==331
if mode in('fixed347','negctl'):
    added=[r['sample_id'] for r in csv.DictReader(open('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/roster/participant15_added.tsv'),delimiter='\t')]
    assert len(added)==16
    base=Path('LAB_ROOT/canonical/publication_data/arms/noclust_2026-08-07/cohorts/CAPSCAN/4_parsimony/razor_hash-acc')
    have={f['path'] for f in fastas}
    for s in added:
        st=W.stamp(base/s/'parsimony.fasta');assert st['path'] not in have;fastas.append(st)
    assert len(fastas)==347
    if mode=='negctl':fastas[-1]=dict(fastas[-1],bytes=fastas[-1]['bytes']+1)
folder=O/mode;folder.mkdir(exist_ok=True)
for p in ['union.fasta','union_summary.json','sequence_provenance.tsv.gz']:
    if (folder/p).exists():raise SystemExit(f'refusing: {folder/p} exists')
shutil.copy2(SRC/'cohorts/CAPSCAN_noclust/sage.json',folder/'sage.json')  # build_union reads only database.decoy_tag from it
W.build_union(dict(folder=str(folder),fastas=fastas))
