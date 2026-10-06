"""Benjamini-Hochberg q values for the four direct-assignment EDS prediction tests.

The v30 manuscript reports EDS held-out prediction for direct peptide assignment only, so
the multiple-testing family shrinks from eight grouping-by-model tests to four model
tests. This script (1) recomputes the stored eight-test q values as a control
(``permutation_q_8``, written by ``analyze_score.py``), and (2) recomputes BH q values over
the four ``grouping == 'simplified'`` rows (protein, baseline, combined, distance).
Statistical unit: patient, n = 85. Permutation P values are plus-one with 199 label
permutations, first split seed (20260915).

Input: the directory holding ``model_metrics.tsv`` and ``phenotype_associations.tsv``
(argument 1; LAB_ROOT/.../iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit).
Output (current directory): ``eds_bh_direct4.tsv``.

Code as it ran (LAB_ROOT/.../final/v30/grouping_recount/eds_bh.py, Slurm job 6074462);
only this docstring was added.
"""
import sys, pandas as pd, numpy as np
from statsmodels.stats.multitest import multipletests
D=sys.argv[1]
m=pd.read_csv(D+'/model_metrics.tsv',sep='\t')
ok=m.permutation_p.notna()
q8=multipletests(m.loc[ok,'permutation_p'],method='fdr_bh')[1]
print('CONTROL 8-test max|diff| vs stored q_8:',np.abs(q8-m.loc[ok,'permutation_q_8']).max())
d=m[ok & m.grouping.eq('simplified')].copy()
d['q_4']=multipletests(d.permutation_p,method='fdr_bh')[1]
print(d[['grouping','seed','model','rho','permutation_p','permutation_q_8','q_4']].to_string(index=False))
d.to_csv('eds_bh_direct4.tsv',sep='\t',index=False)
a=pd.read_csv(D+'/phenotype_associations.tsv',sep='\t'); print(a.to_string(index=False))
