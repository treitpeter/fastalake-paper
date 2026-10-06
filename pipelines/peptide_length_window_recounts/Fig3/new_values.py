"""Write NEW_VALUES.tsv for Fig3 and SFig07 (old 5-50 vs new 9-50) from the same tables the builds used."""
import os  # release: paths come from environment variables
import sys
from pathlib import Path
import numpy as np,pandas as pd
D=Path(__file__).resolve().parent.parent;REPO=Path((os.environ.get('FASTALAKE_PAPER_REPO') or str(Path(__file__).resolve().parents[3])))
sys.path.insert(0,str(D/'Fig3'));import panels_len950 as P
c=pd.read_csv(REPO/'source_data/Fig3/transfer_pool_comparisons.tsv',sep='\t')
SRC3='repo source_data/Fig3/transfer_pool_comparisons.tsv (A15 all_comparisons.tsv), rows min_length 5 vs 9 (9-50 = accepted 5-50 set filtered to >=9 residues by A15 search.py)'
rows=[]
def med(L,kind,arm,field,g=None):
    x=c[c.kind.eq(kind)&c.min_length.eq(L)&c.arm.eq(arm)]
    if g:x=x[x.group.eq(g)]
    return float(np.median(x[field]))
lab={'T1_within_group':'same series (T1)','T4_different_donor':'different donor (T4)'}
for pan,field,q in [('b','count_retention_percent','median peptide count, % of own'),('c','own_set_recovery_percent','median own identities recovered, %')]:
    for arm in P.TRANSFER_ARMS:
        for g in P.DESIGNS:
            rows.append((pan,f'{q}; {g}; {lab[arm]}',med(5,'transfer',arm,field,g),med(9,'transfer',arm,field,g),SRC3))
for pan,field,q in [('d','count_retention_percent','median peptide count, % of own'),('e','own_set_recovery_percent','median own identities recovered, %')]:
    for arm in P.POOL_ARMS[1:]:
        rows.append((pan,f'{q}; {arm}',med(5,'pool',arm,field),med(9,'pool',arm,field),SRC3))
for f in ['gained','lost']:
    rows.append(('d (text L103)',f'median peptides {f}, other-11 union vs own',med(5,'pool','pool_other11',f),med(9,'pool','pool_other11',f),SRC3))
for arm in P.POOL_ARMS[1:]:
    x5=c[c.kind.eq('pool')&c.min_length.eq(5)&c.arm.eq(arm)];x9=c[c.kind.eq('pool')&c.min_length.eq(9)&c.arm.eq(arm)]
    rows.append(('d',f'acquisitions with more peptides than own; {arm}',int((x5.peptides>x5.own).sum()),int((x9.peptides>x9.own).sum()),SRC3))
pd.DataFrame(rows,columns=['panel','quantity','old_5_50','new_9_50','source_file_and_filter']).to_csv(D/'Fig3/NEW_VALUES.tsv',sep='\t',index=False)
# SFig07
rows=[]
ic=pd.read_csv(REPO/'source_data/SFig07/identity_comparison.tsv',sep='\t')
H=pd.read_csv(D/'SFig07/work/hstool/identity_per_acquisition_len.tsv',sep='\t');S=pd.read_csv(D/'SFig07/work/hstool/identity_summary_len.tsv',sep='\t');C=pd.read_csv(D/'SFig07/work/hstool/completeness_len.tsv',sep='\t')
SA='repo source_data/SFig07/identity_comparison.tsv (A21 compare_revision2.py, job 5928104), rows min_length 5 vs 9'
SH='len950/SFig07/work/hstool/*_len.tsv (hstool_recount.py, srun job 6072769) from A21 results/<arm>/accepted_identities.tsv.gz and A22 joint_search/out/joint_fileq (job 6048637); 9-50 = accepted 5-50 canonical peptides with >=9 residues'
SC='len950/SFig07/work/hstool/completeness_len.tsv (hstool_recount.py, srun 6072769) from same_file_id.features.npz; features kept if canonical peptide >=9 residues'
for col,arm in [('joint corrected study q','joint_corrected'),('joint original study q','joint_original'),('repeat individual','repeat_individual'),('union individual','union_individual')]:
    f=lambda L,v:float(np.median(ic[ic.min_length.eq(L)&ic.test.eq(arm)][v]))
    rows.append(('a',f'median original identities recovered, %; {col}',100*f(5,'original_identity_recovery'),100*f(9,'original_identity_recovery'),SA))
g=lambda L,arm,v:float(np.median(H[H.min_length.eq(L)&H.arm.eq(arm)][v]))
rows.append(('a',f'median original identities recovered, %; joint corrected acquisition q',100*g(5,'joint_fileq','original_identity_recovery'),100*g(9,'joint_fileq','original_identity_recovery'),SH))
for arm,lab in [('union_individual','separate common-FASTA searches'),('joint_fileq','joint, per-acquisition q'),('joint_corrected','joint, study-level q')]:
    s=lambda L,v:S[S.min_length.eq(L)&S.arm.eq(arm)][v].iloc[0]
    rows.append(('text L115',f'median peptide-count gain, %; {lab}',100*(s(5,'median_count_ratio')-1),100*(s(9,'median_count_ratio')-1),SH))
    rows.append(('text L115 / a',f'median original identities retained, %; {lab}',100*s(5,'median_identity_recovery'),100*s(9,'median_identity_recovery'),SH))
    rows.append(('text L115',f'acquisitions with more peptides (of 58); {lab}',s(5,'n_higher'),s(9,'n_higher'),SH))
for arm,lab in [('original','acquisition-specific'),('union_individual','separate common-FASTA searches'),('joint_fileq','joint, per-acquisition q'),('joint_corrected','joint, study-level q')]:
    s=lambda L:S[S.min_length.eq(L)&S.arm.eq(arm)].study_distinct.iloc[0]
    rows.append(('text L115',f'distinct peptides across 58 acquisitions; {lab}',s(5),s(9),SH))
    k=lambda L:C[C.min_length.eq(L)&C.arm.eq(arm)&C.tier.eq('TechnicalReplicates')].rows_complete_all.iloc[0]
    rows.append(('text L115',f'peptide features complete in all 12 technical replicates; {lab}',k(5),k(9),SC))
rows.append(('text L115','median gain with study-level q plus PSM PEP filter, %',5.9,'NOT RECOMPUTED','untraced in repo (SFig07 OPEN_ITEMS; NUMBERS.tsv OPEN); no producer found, so no 9-50 value'))
q5=pd.read_csv(REPO/'source_data/SFig07/common_search_participants.tsv',sep='\t');q9=pd.read_csv(D/'SFig07/build/source_data_len950/common_search_participants.tsv',sep='\t')
SB='len950/SFig07/work/clin/part_*.tsv (clin_recount.py, Slurm array 6072755) from each acquisition own and common-FASTA results.sage.tsv (rank 1, target, peptide q<=0.01, canonical; 9-50 = >=9 residues), aggregated as summarize_v2.py'
for co in ['PREDICT','CAPSCAN']:
    for stat in ['median','min','max']:
        rows.append(('b',f'participant accepted-peptide change, common vs own FASTA, % ({stat}); {co}',getattr(q5[q5.cohort.eq(co)].peptides_percent_difference,stat)(),getattr(q9[q9.cohort.eq(co)].peptides_percent_difference,stat)(),SB))
    rows.append(('b',f'participants with more peptides on common FASTA; {co}',int((q5[q5.cohort.eq(co)].peptides_percent_difference>0).sum()),int((q9[q9.cohort.eq(co)].peptides_percent_difference>0).sum()),SB))
pd.DataFrame(rows,columns=['panel','quantity','old_5_50','new_9_50','source_file_and_filter']).to_csv(D/'SFig07/NEW_VALUES.tsv',sep='\t',index=False)
