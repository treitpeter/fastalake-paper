#!/usr/bin/env python3
"""CapScan rows of the SFig08 b-d table, with the exact transform of
MM_review_analysis/end_of_day_2026-09-15/code/plot_corrected_ablation.py (line 13):
 intensity_removed_percent=100*fraction_intensity_removed; KO_support_lost_percent=100*(1-shared_KOs/all_KOs);
 BrayCurtis=normalized_KO_profile_BrayCurtis_all_vs_shared; cohort='CapScan'. PREDICT rows are taken unchanged from the
 committed Source Data table. Writes clinical_fixed_dictionary_ablation_<n>.tsv and medians_<n>.tsv (stat_bcd_medians rule)."""
import sys
from pathlib import Path
import pandas as pd
summary, n, ref, outdir = Path(sys.argv[1]), int(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
d = pd.read_csv(summary, sep='\t'); assert len(d) == n and d['sample'].is_unique
d['cohort'] = 'CapScan'; d['intensity_removed_percent'] = 100 * d.fraction_intensity_removed
d['KO_support_lost_percent'] = 100 * (1 - d.shared_KOs / d.all_KOs); d['BrayCurtis'] = d.normalized_KO_profile_BrayCurtis_all_vs_shared
r = pd.read_csv(ref, sep='\t'); pred = r[r.cohort.eq('PREDICT')]
assert list(pred.columns) == list(d.columns), (list(pred.columns), list(d.columns))
out = pd.concat([pred, d]); out.to_csv(outdir / f'clinical_fixed_dictionary_ablation_{n}.tsv', sep='\t', index=False)
cols = ['intensity_removed_percent', 'KO_support_lost_percent', 'BrayCurtis']
m = out.groupby('cohort')[cols].median(); m.insert(0, 'n', out.groupby('cohort')['sample'].nunique())
m.to_csv(outdir / f'medians_{n}.tsv', sep='\t'); print(m.to_string())
