#!/usr/bin/env python3
"""Regenerate the SFig09 derived tables (per specimen, intestinal acquisitions with tertile, tertile summaries 5-50 and
9-50) from the 295-row acquisition tables with the unit's own derive_tables() (the moved-unchanged producer step), and
print the new EXPECTED_BIN_N. Run with cwd = phase2/repo/figures/SFig09."""
import sys
from pathlib import Path
U = Path.cwd(); sys.path[:0] = [str(U.parents[1]), str(U)]
import pandas as pd, panels
SD = U.parents[1] / 'source_data/SFig09'
acq = pd.read_csv(SD / 'capscan_assembly_dna_yield_gain_by_acquisition.tsv', sep='\t')
assert len(acq) == 295
spec, caps, summary = panels.derive_tables(acq)
spec.to_csv(SD / 'capsule_assembly_gain_per_specimen.tsv', sep='\t', index=False)
caps.to_csv(SD / 'capsule_assembly_gain_acquisitions.tsv', sep='\t', index=False)
summary.to_csv(SD / 'capsule_assembly_gain_tertiles.tsv', sep='\t', index=False)
a9 = pd.read_csv(SD / 'capscan_acquisitions_9_50.tsv', sep='\t')
_, _, s9 = panels.derive_tables(a9)
s9.to_csv(SD / 'capsule_assembly_gain_tertiles_9_50.tsv', sep='\t', index=False)
print('specimens', len(spec), 'intestinal acquisitions', len(caps))
print(summary.to_string()); print(s9.to_string())
print('EXPECTED_BIN_N', list(summary.n))
