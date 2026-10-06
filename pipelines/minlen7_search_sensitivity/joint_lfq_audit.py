from pathlib import Path
import numpy as np
import pandas as pd
def audit_joint(out,cfg):
    ranks=pd.read_csv(out/'ranked_lfq_scores.tsv',sep='\t');assert np.array_equal(ranks['rank'],np.arange(len(ranks)));assert not ranks.duplicated(['precursor_id','charge','decoy']).any();assert set(ranks.decoy)=={0,1}
    scores=ranks.score_f32.to_numpy(dtype=np.float32);assert np.isfinite(scores).all() and (np.diff(scores)<=0).all()
    d=ranks.decoy.to_numpy(dtype=np.float32);den=np.cumsum(1-d,dtype=np.float32)
    with np.errstate(divide='ignore'):raw=(1+np.cumsum(d,dtype=np.float32))/den
    expected=np.minimum.accumulate(np.minimum(raw,1)[::-1])[::-1];np.testing.assert_allclose(ranks.q_value.to_numpy(),expected,atol=1e-7,rtol=1e-6)
    search=out;a=pd.read_csv(search/'lfq_target_decoy_audit.tsv',sep='\t');assert len(a)==len(ranks)
    joined=a.merge(ranks,on=['precursor_id','charge','decoy'],validate='one_to_one',suffixes=('_export','_rank'));assert len(joined)==len(a);np.testing.assert_array_equal(joined.score.to_numpy(dtype=np.float32),joined.score_f32.to_numpy(dtype=np.float32));np.testing.assert_allclose(joined.q_value_export,joined.q_value_rank,atol=1e-7,rtol=1e-6)
    primary=pd.read_csv(search/'lfq.tsv',sep='\t');targets=a[a.decoy==0].drop(columns=['precursor_id','decoy']);assert list(primary)==list(targets)
    # Exported target measurements and q-values must be exactly those in the complete audit.
    sort=['peptide','charge','proteins','score'];pd.testing.assert_frame_equal(primary.sort_values(sort).reset_index(drop=True),targets.sort_values(sort).reset_index(drop=True),check_exact=True)
    rawcols=[Path(p).name for p in cfg['mzml_paths']];v=targets[rawcols].to_numpy();finite=np.isfinite(v)&(v>0)
    audit=dict(status='PASS',ranked_peaks=len(ranks),target_peaks=len(targets),decoy_peaks=int((a.decoy==1).sum()),targets_q05=int((targets.q_value<=.05).sum()),targets_q01=int((targets.q_value<=.01).sum()),positive_finite_target_cells=int(finite.sum()),invalid_or_negative_cells=int((~np.isfinite(v)|(v<0)).sum()),actual_tie_order_retained=True,q_values_independently_reconstructed=True,all_target_decoy_assignments_verified=True,target_output_matches_complete_audit=True,confidence_unit='Pooled peptide-feature confidence, not per-acquisition transferred-cell confidence')

    return audit
