"""Write A4 configs (control331, fixed347, negative) with input sha256. Run under Slurm."""
import json, hashlib
from pathlib import Path
A = Path(__file__).resolve().parents[1]
F = Path(__import__('os').environ['FASTALAKE_DATA_DIR'])
B = F/'review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit'
C2 = B/'claude_review_2026-10-03/comment_check_v25_2026-10-05/compute/C2_participant15'
MC = F/'00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23'
REC = F/'review/v15_integration/integration_2026-09-29/cohort_reconciliation/crossomics'
def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''): h.update(b)
    return h.hexdigest()
def files(c):
    fs = [c[k] for k in ['design', 'specimen_roster', 'metabolites', 'dna_cells', 'previous_rows'] + [k for k in ['metabolites_npy', 'dna_cached', 'recon_native'] if k in c]]
    for g in ['other_unresolved_groups', 'supported_human_groups']: fs += [c['groups'] + '/' + g + '.npz', c['groups'] + '/' + g + '_axes.json']
    return {p: sha(p) for p in fs}
x347 = json.loads((C2/'fig5g/crossomics_347/COMPLETE.json').read_text())
control = dict(mode='control331', out=str(A/'control331'), groups=str(B/'groups/CAPSCAN'), design=str(F/'analyses/enterotypes_2026-09-07/CAPSCAN_specimen_design.tsv'),
    specimen_roster=str(B/'cross_refit/capscan_specimen_roster.tsv'), donors=[1, 2, 3, 4, 6, 8, 9, 10, 11, 12, 13, 14],
    metabolites=str(F/'analyses/enterotypes_2026-09-07/CAPSCAN_metabolite_intensity.tsv'), metabolites_npy=str(MC/'method_comparison_mofa_2026-09-23/prepared/CAPSCAN/Metabolites_log_profiles.npy'),
    dna_cells=str(MC/'results/CAPSCAN/DNA_native_cells.npz'), dna_cached=str(MC/'method_comparison_mofa_2026-09-23/distances/CAPSCAN/distance_matrices.npz'),
    previous_rows=str(B/'distance_refit/comparisons.tsv'))
fixed = dict(mode='fixed347', out=str(A/'run347'), groups=str(C2/'fig5g/groups_347'), design=str(C2/'fig5g/CAPSCAN_specimen_design_347.tsv'),
    specimen_roster=str(C2/'fig5g/crossomics_347/capscan_specimen_roster.tsv'), donors=x347['dna_participants'],
    metabolites=str(REC/'metabolomics/CAPSCAN_metabolite_intensity.tsv'), dna_cells=str(REC/'native/current331_plus16_all/DNA_native_cells.npz'),
    recon_native=str(REC/'native/current331_plus16_all/distance_matrices.npz'), previous_rows=str(B/'distance_refit/comparisons.tsv'))
assert len(fixed['donors']) == 13 and 15 in fixed['donors']
for c in (control, fixed): c['sha256'] = files(c)
neg = dict(fixed, out=str(A/'negative_should_not_exist'), donors=sorted(fixed['donors'] + [5]))  # participant 5 lacks DNA in one region
(A/'control331.json').write_text(json.dumps(control, indent=1) + '\n'); (A/'fixed347.json').write_text(json.dumps(fixed, indent=1) + '\n')
(A/'negative_participant5.json').write_text(json.dumps(neg, indent=1) + '\n'); print('configs written', fixed['donors'])
