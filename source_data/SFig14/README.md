# Source data: Supplementary Figure 14

All files are copies of the Source Data package folder `Supplementary_Figure_14/` (identical
bytes), plus `spectra/spectrum_audit.tsv` and `spectra/same_scan_all_search_calls.tsv` from the
upstream analysis (A25). The package
`PROVENANCE.json` files carry file names only; the upstream originals held absolute paths and
are otherwise identical. Raw-file names (for example
`20191207_QX7_PeTr_LC02_SA_MicrobPred_224`) identify acquisitions; they carry no clinical
information. No clinical column is present.

## PREDICT_spectrum_case_review.tsv

One row per case (4 rows).

| column | meaning |
|---|---|
| case_id | HIS_001 to HIS_004 (panels a-d) |
| sample | acquisition (raw file name without extension) |
| compartment | stool |
| origin | Human or Microbial (source protein) |
| side | boundary side of the histidine: N (histidine just before the peptide in the source protein) or C (peptide ends in histidine) |
| peptide | peptide as reported by Sage (modifications in brackets, Da) |
| matched_ions | matched b/y fragment ions (Sage matched_peaks) |
| backbone_coverage | fraction of peptide bonds covered by at least one matched b or y ion |
| peptide_q, spectrum_q | Sage q values |
| accepted_in_semi | True if a semi-tryptic PSM on the same scan is a target with peptide and spectrum q <= 0.01 and the same canonical peptide (case audit, line 158). PREDICT had no semi-tryptic search in A25, so it is False for all four cases |
| nonoriginal_alternatives | bounded sequence alternatives other than the reported peptide |
| nonoriginal_matching_at_least_original | of those, how many match at least as many ions |
| N_flanking_H_inferred | side == N: the histidine is inferred from the source protein, not observed in the peptide (case audit, line 160) |
| parents | source protein accession(s) |

## spectra/HIS_00k/

| file | columns | meaning |
|---|---|---|
| observed_peaks.tsv.gz | mz, intensity | every peak of the MS2 scan, read from the mzML (float32) |
| processed_peaks.tsv | mz, intensity | peaks in 200-2000 m/z, the 150 most intense, sorted by m/z; the gray peaks |
| matched_fragments.tsv | series (b/y), ordinal, charge, observed_mz, calculated_mz, intensity, ppm | matched fragment ions; ppm = (observed - calculated)/(calculated - proton) x 1e6 |
| bounded_sequence_alternatives.tsv | peptide, canonical_peptide, change, parents, isotope_errors, calculated_mass, matched_ions, original, at_least_original_matches | alternatives present in the local FASTA within precursor tolerance; `original` marks the reported peptide |
| PROVENANCE.json | raw, fasta, unspecific_config (file name, bytes, mtime, sha256); selected_psm (the Sage PSM row and selection ranks); interpretation | identity of inputs and the selected PSM |

## spectra/ selection tables

| file | rows | meaning |
|---|---|---|
| complete_candidate_counts.tsv | 1,224 | per acquisition (306 PREDICT stool acquisitions) and stratum (origin x side): number of candidate peptides (columns cohort, compartment, origin, side, sample, candidates) |
| all_strata.tsv | 4 | strata and status (SELECTED) with case id |
| selected_cases.tsv | 4 | the selected Sage PSM rows plus case_id, acquisition_rank, acquisition_median_rank, eligible_peptides, selected_rank |
| eligible_case_selection_psms.tsv.gz | 173 | the best PSM per canonical peptide in each chosen acquisition, with selection_rank and selected flag |
| spectrum_audit.tsv | 4 | upstream summary: matched_ions, backbone_coverage, alternatives, nonoriginal_at_least_original |
| same_scan_all_search_calls.tsv | 6 | every Sage PSM on the four selected scans, with case_id and route (tryptic or unspecific) prepended to the Sage columns |

Sage PSM columns follow the Sage 0.14.6 `results.sage.tsv` format; `canonical_peptide` is the
unmodified sequence with I replaced by L.
