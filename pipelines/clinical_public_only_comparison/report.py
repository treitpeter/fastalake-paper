#!/usr/bin/env python3
"""Collect validated §4.4/§4.5 results without changing their frozen outputs."""
import csv,json,os,sys
from pathlib import Path
import pandas as pd
RUN=Path(sys.argv[1]).resolve();sys.path.insert(0,str(RUN/'frozen'))
from lib import dump,meta,verify
m=json.loads((RUN/'inputs/report_manifest.json').read_text())
for x in m['code']+m['receipts']+m.get('additional_inputs',[]):verify(x)
out=RUN/'results';assert not (out/'REPORT_VALIDATION.json').exists()
for p in m['receipts']:
    c=json.loads(Path(p['path']).read_text());assert c['status']=='PASS'
    for x in c['outputs']:verify(x)
if RUN.name.startswith('A15'):
    f=pd.concat([pd.read_csv(Path(p['path']).with_name('comparisons.tsv'),sep='\t') for p in m['receipts']],ignore_index=True)
    assert len(m['receipts'])==70 and not f.duplicated(['sample','kind','arm','min_length']).any()
    f.to_csv(out/'all_comparisons.tsv',sep='\t',index=False)
    summary=f.groupby(['kind','group','arm','min_length']).agg(acquisitions=('sample','nunique'),median_count_retention_percent=('count_retention_percent','median'),median_own_identity_recovery_percent=('own_set_recovery_percent','median'),median_gained=('gained','median'),median_lost=('lost','median'),median_own_digest_compatible_unidentified=('own_digest_compatible_unidentified','median'),median_database_accessions=('database_accessions','median')).reset_index()
    summary.to_csv(out/'comparison_summary.tsv',sep='\t',index=False)
    lines=['# §4.4 — transfer and other-eleven pooling','',
      'All 70 target/experiment audits passed: 58 transfer targets and 12 pooling targets, encompassing 210 saved raw searches. Actual donor FASTAs, exact pool unions, raw acquisitions and search settings were checked. Counts and identities were independently recounted.', '',
      '| Design | Arm | Acquisitions | Median count retained | Median own identities retained |','|---|---|---:|---:|---:|']
    for r in summary[summary.min_length.eq(5)].itertuples():lines.append(f'| {r.group} | {r.arm} | {r.acquisitions} | {r.median_count_retention_percent:.2f}% | {r.median_own_identity_recovery_percent:.2f}% |')
    lines+=['', 'Count retention and identity recovery have different denominators. A transferred FASTA can lose original peptides and gain other ones. Values are medians of paired acquisition ratios; all gains/losses and digest-compatible-but-unidentified peptides remain in source data.', '',
      'T1 means another acquisition in the same design; T4 is the declared other-design donor. The biological relationship is not inferred from tier labels. In the interindividual design, within-design transfer is not a same-person control. The donor labels are retained explicitly.', '',
      'Other-eleven excludes the target final FASTA from the exact union; the upstream cohort evidence remains shared. This control therefore tests final database pooling, not a wholly independent end-to-end workflow. Twelve technical acquisitions are not twelve biological participants. Scientific/figure review and manuscript integration remain pending.']
else:
    assert RUN.name.startswith('A16') and len(m['receipts'])==2
    rows=[];frames=[]
    for p in m['receipts']:
        c=json.loads(Path(p['path']).read_text());co=c['cohort'];f=pd.read_csv(Path(p['path']).with_name('per_acquisition.tsv'),sep='\t');f['cohort']=co;frames.append(f)
        for (a,b,lo),g in f.groupby(['test_arm','reference_arm','min_length']):
            rows.append(dict(cohort=co,test=a,reference=b,min_length=int(lo),acquisitions=len(g),defined_pairs=int(g.count_ratio.notna().sum()),test_more=int(g.n_test.gt(g.n_reference).sum()),same_count=int(g.n_test.eq(g.n_reference).sum()),test_fewer=int(g.n_test.lt(g.n_reference).sum()),median_paired_count_ratio=float(g.count_ratio.median()),ratio_of_median_counts=float(g.n_test.median()/g.n_reference.median()),median_reference_identity_recovery=float(g.reference_recovery.median()),median_test_identity_recovery=float(g.test_recovery.median())))
    pd.concat(frames).to_csv(out/'all_acquisitions.tsv',sep='\t',index=False);pd.DataFrame(rows).to_csv(out/'count_definition_reconciliation.tsv',sep='\t',index=False)
    lines=['# §4.5 — public-only versus cohort-MG and metagenome searches','',
      'All 306 PREDICT and 331 CapScan acquisitions are retained. The metagenome comparison uses 296 and 282 matched acquisitions respectively. The public reservoir/source proof passed independently for both cohorts (A08 and A17). Shared cohort de novo evidence remains part of selection; no participant assembly is in the public-only reservoir.', '',
      '| Cohort | Test / reference | Paired acquisitions | More / equal / fewer | Median paired count ratio | Ratio of medians | Median reference identities recovered |','|---|---|---:|---:|---:|---:|---:|']
    for r in rows:
        if r['min_length']==5:lines.append(f"| {r['cohort']} | {r['test']} / {r['reference']} | {r['acquisitions']} | {r['test_more']} / {r['same_count']} / {r['test_fewer']} | {r['median_paired_count_ratio']:.4f} | {r['ratio_of_median_counts']:.4f} | {r['median_reference_identity_recovery']:.2%} |")
    lines+=['', 'FL_public excludes participant metagenome assemblies from the reservoir. FL_cohortMG allows them. The no-sequencing claim must use FL_public; the cohort-MG arm is a separate comparison. Recovery denotes the reference search\'s own canonical target peptide identities, not agreement with molecular ground truth.', '',
      'The archived PREDICT 1.23-fold result is the ratio of medians (1.233879) for FL_cohortMG versus MG. The median paired ratio is 1.207318. Public-only PREDICT gives 1.190153 paired and 1.217943 ratio-of-medians. Do not combine an estimator or arm from one cohort with a differently defined result in another.', '',
      'Canonical 5–50 peptide identities at target rank-1 peptide-q ≤0.01 are primary; nested 9–50 results remain in every table. Nominal search confidence is not an independent false-discovery calibration. Count gains and losses are both retained; the public-only arm does not contain every metagenome-identified peptide.', '',
      'Participant summaries are saved separately with equal acquisitions within specimen, then equal specimens within participant. Four PREDICT acquisitions without participant IDs remain acquisition-only. The CapScan grouping was independently checked against the canonical manuscript design: the same 292 participant/region/specimen groups in both directions, with no split or merge. Zero denominators remain undefined. Mean ratios can be driven by small reference counts; paired medians and denominator counts accompany them.', '',
      'Scientific/figure review, cost/turnaround sourcing and manuscript integration remain pending. Host/microbial clinical robustness is separately evaluated in A19 and terminal-source compatibility in A20; no public-arm clinical biology equivalence is asserted here.']
(out/'REPORT.md').write_text('\n'.join(lines)+'\n')
outputs=[meta(p) for p in out.iterdir() if p.is_file() and p.name not in ['REPORT_VALIDATION.json','AGGREGATE_VALIDATION.json']]
record=dict(status='PASS',source_receipts=m['receipts'],job_id=os.environ.get('SLURM_JOB_ID'),scientific_review='PENDING',outputs=outputs)
dump(out/'REPORT_VALIDATION.json',record);dump(out/'AGGREGATE_VALIDATION.json',record)
print('TRANSFER/PUBLIC REPORT PASS',RUN.name,flush=True)
