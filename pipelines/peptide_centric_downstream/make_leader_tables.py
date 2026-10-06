#!/usr/bin/env python3
"""Rebuild the two participant leader tables Figure 5e reads, from the revised collector.

Same shape as the published MM_A19_CAPSCAN_profiles tables: one row per participant, one
column per grouping rule, agreement averaged over that participant's regions.
"""
import sys
from pathlib import Path
import pandas as pd
PACKAGE = Path(__file__).resolve().parents[1]
src = pd.read_csv(PACKAGE / 'results/matched_capsule_stool_donors.tsv', sep='\t')
out = PACKAGE / 'render/tables'
out.mkdir(parents=True, exist_ok=True)
for threshold, name in [('resolved', 'participant_leaders_resolved.tsv'),
                        ('strict_both_methods', 'participant_leaders_strict_both_methods.tsv')]:
    s = src[src.scope.eq('common_functions') & src.threshold.eq(threshold)]
    t = (s.groupby('participant')[['agreement_original', 'agreement_simplified']].mean()
         .reset_index())
    assert len(t) == 14, (threshold, len(t))
    t.to_csv(out / name, sep='\t', index=False)
    print(f'{name}: {len(t)} participants, original {t.agreement_original.mean():.4f}, '
          f'simplified {t.agreement_simplified.mean():.4f}')
