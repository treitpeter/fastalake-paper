"""Extend the A33/FIGURATOR CapScan design.tsv (331 rows) with the 16 join-fixed participant-15 capsules.

Covariates come from Capscan_Main_PeTr.xlsx through the unchanged canonical builder's own
load_xlsx_index()/acq_number() (canonical/publication/metadata/build_capscan_metadata.py, imported
read-only). Participant and region come from the join-fixed roster. specimen_id = 'main:<Main>' as for
every other capsule. Control: the same procedure regenerates every covariate of the 331 existing rows.
usage: build_design_347.py <out.tsv>
"""
import importlib.util, json, sys
from pathlib import Path
import pandas as pd
OUT = Path(__file__).resolve().parents[1]; C2 = OUT.parent
A33D = Path('LAB_ROOT/MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/CAPSCAN__simplified/source_data/Figure5/design.tsv')
FIG = Path('LAB_ROOT/FIGURATOR/Figure5/v5_2026-09-07/source_data/design.tsv')
B = Path('LAB_ROOT/canonical/publication/metadata/build_capscan_metadata.py')
spec = importlib.util.spec_from_file_location('capmeta', B); cm = importlib.util.module_from_spec(spec); spec.loader.exec_module(cm)
d = pd.read_csv(A33D, sep='\t', dtype=str, keep_default_na=False)
assert d.equals(pd.read_csv(FIG, sep='\t', dtype=str, keep_default_na=False)), 'A33 design differs from FIGURATOR design'
assert len(d) == 331
idx = cm.load_xlsx_index()
cov = [o for o, _ in cm.COVAR]
# control: regenerate covariates for the 331
bad = []
for _, r in d.iterrows():
    which, rec = idx[cm.acq_number(r.sample_id)]
    if which.lower() != r.aliquot or any(rec[c] != r[c] for c in cov): bad.append(r.sample_id)
assert not bad, ('covariate regeneration control failed', bad[:5])
add = pd.read_csv(C2/'roster/participant15_added.tsv', sep='\t')
rows = []
for _, a in add.iterrows():
    n = cm.acq_number(a.sample_id); which, rec = idx[n]
    assert which == 'Pellet' and rec['_subject'] == '15', (a.sample_id, which, rec['_subject'])
    assert rec['aliquot_main_id'] == str(a.workbook_Main), a.sample_id
    row = {c: '' for c in d.columns}
    row.update(sample_id=a.sample_id, cohort='CAPSCAN', participant='15', region=a.region, group_ordinal=str(a.region_ordinal),
               aliquot='pellet', **{c: rec[c] for c in cov},
               subject_source='C2 join-fixed roster (fixed_roster_347.tsv)', group_source='C2 join-fixed roster (fixed_roster_347.tsv)',
               covariate_source='Capscan_Main_PeTr.xlsx:Sheet1', subject_xlsx_crosscheck='agree', metadata_built='2026-10-05',
               specimen_id='main:' + rec['aliquot_main_id'], specimen_mapping_valid='True')
    p15 = d[d.participant.eq('15')]
    row['prior_antibiotic_months'] = p15.prior_antibiotic_months.iloc[0] if len(p15) else 'not reported'
    rows.append(row)
x = pd.concat([d, pd.DataFrame(rows)[d.columns]], ignore_index=True)
assert len(x) == 347 and x.sample_id.is_unique and x.specimen_id[x.participant.eq('15') & x.region.ne('stool')].is_unique
x.to_csv(sys.argv[1], sep='\t', index=False)
print(json.dumps(dict(status='PASS', rows=len(x), added=len(rows), control_331_covariates_regenerated=True,
                      p15_by_region=x[x.participant.eq('15')].region.value_counts().to_dict()), indent=1))
