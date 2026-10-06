"""RES-17 Source Data table: every accepted peptide of the worked acquisition, flagged if absent from the archived AlphaNovo predictions.
Same reading rules as recheck_4258.py (A09 analysis.py): target, rank 1, peptide q<=0.01, modifications stripped, I->L; predictions = unmodified detokenized sequences, I->L.
usage: python res17_table.py LAB_ROOT OUT.tsv"""
import csv,hashlib,re,sys
L,out=sys.argv[1:3]
pred=L+'/MM_review_analysis/A01_4_1_database_sizes/inputs/predictions/20181231_QX7_PhGe_LC95_SA_WorkflowReplicates_2_predictions.csv'
sage=L+'/canonical/publication_data/arms/noclust_2026-08-07/cohorts/HSTOOL/5_sage/razor_hash-acc/20181231_QX7_PhGe_LC95_SA_WorkflowReplicates_2/results.sage.tsv'
sha=lambda p:hashlib.sha256(open(p,'rb').read()).hexdigest()
assert sha(pred)=='8826920da56b6ed2fa443f519270c830f02294e40f36a7c497741f64b13cb15d'
assert sha(sage)=='92339fc9c159bd00ba24c0edd07e3955120962a98b1d0951c3a67dc012a391d4'
MOD=re.compile(r'\[[^\]]*\]')
allp=set()
for r in csv.DictReader(open(pred)):
    s=(r.get('peptide_prediction_detokenized_unmodified') or '').strip().upper()
    if re.fullmatch('[A-Z]+',s):allp.add(s.replace('I','L'))
best={}
for r in csv.DictReader(open(sage),delimiter='\t'):
    if r['label']!='1' or r['rank']!='1' or not float(r['peptide_q'])<=.01:continue
    p=MOD.sub('',r['peptide']).replace('I','L'); q=float(r['peptide_q'])
    best[p]=min(q,best.get(p,1.0))
rows=sorted(best)
nmiss=sum(p not in allp for p in rows)
assert (len(rows),nmiss)==(29663,4258),(len(rows),nmiss)
with open(out,'w') as f:
    f.write('acquisition\tpeptide_IL\tlength\tbest_peptide_q\tin_alphanovo_predictions\n')
    for p in rows: f.write(f'20181231_QX7_PhGe_LC95_SA_WorkflowReplicates_2\t{p}\t{len(p)}\t{best[p]:.6g}\t{"yes" if p in allp else "no"}\n')
print(len(rows),nmiss,f'{100*nmiss/len(rows):.1f}%')
