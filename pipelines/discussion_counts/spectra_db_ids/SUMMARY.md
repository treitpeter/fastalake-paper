# Spectra, database size, candidate competition and identifications per acquisition

Written 2026-10-06. Existing search outputs only, no new search. Slurm jobs: 6073318 (tasks.py), 6073323-6 (pilots),
6073330 (worker.py, 100-task array, peak RSS 1.14 GB), 6073437 (collect.py). Rows read: 332,406,424 Sage PSM rows;
701 mzML files (111,649,751 MS2 spectra in total).

Files: per_acquisition.tsv (701 rows: PREDICT 296, CapScan 347, healthy stool 58), summary_medians.tsv,
summary_paired_ratios.tsv, CHECKS.json; scripts tasks.py, worker.py, worker.sbatch, collect.py.

## Arms

| arm | what | source |
|---|---|---|
| FLpub | FL own, public-sequence FastaLake (Fig. 4 / Supp. Fig. 6d,e same-spectrum benchmark) | PREDICT abl_nometag_noclust; CapScan A08 (331 production; 16 participant-15 from public347/A08_347) |
| MG | matched metagenome + human | metag_hashed_human |
| FLcmg | FL own, cohort FastaLake (noclust, with cohort metagenomes): the own-FASTA side of the Supp. Fig. 7b comparison | noclust_2026-08-07 |
| FLunion | FL study-union FASTA, searched per acquisition (Supp. Fig. 7b; healthy stool: A06 study_union) | single_raw_common_fasta_2026-09-14; A06_7_hstool_second_pass |
| FLown | healthy stool FL own (noclust) | noclust_2026-08-07 HSTOOL |

FLcmg is included because the own-FASTA baseline of the Supp. Fig. 7b union searches is the cohort FastaLake search, not FLpub.

## Checks (all PASS, CHECKS.json)

- FLpub and MG accepted sets equal the A10/A08 set files at 5-50 and 9-50 (Supp. Fig. 6d reproduced: PREDICT 21.5k/18.8k, 272/296; CapScan 22.2k/13.7k, 281/282).
- FLcmg/FLunion counts equal len950/SFig07 clin_recount (627 acquisitions); healthy-stool counts equal hstool_recount.
- Paired scored_candidates equal task E1 for all 578 pairs.
- All arms of an acquisition point to the same mzML (realpath); MS1+MS2 equals the spectrumList count in all 701 files.
- Every accepted 9-50 peptide was located in every protein Sage assigned it to (0 unplaced).

## Cohort medians per arm (median over acquisitions)

| stratum | arm | n | MS2 | DB proteins | DB residues (M) | cand. paired | cand. all | pep 5-50 | pep 9-50 | proteins w/ pep | pep/protein mean | pep/protein median | coverage | IDs/1000 MS2 (9-50) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PREDICT | FLpub | 296 | 188,494 | 8,129 | 3.51 | 228 | 228 | 27,255 | 21,484 | 7,802 | 6.22 | 4 | 17.5% | 113.5 |
| PREDICT | MG | 296 | 188,494 | 241,596 | 60.9 | 5,128 | 4,414 | 22,285 | 18,847 | 12,247 | 2.83 | 2 | 8.2% | 99.4 |
| PREDICT | FLcmg | 296 | 188,494 | 8,122 | 3.49 | 226 | 226 | 27,433 | 21,576 | 7,793 | 6.26 | 4 | 17.8% | 114.2 |
| PREDICT | FLunion | 296 | 188,494 | 555,987 | 237.6 | 7,474 | 6,171 | 25,561 | 22,024 | 228,491 | 3.51 | 2 | 7.9% | 115.9 |
| CapScan matched | FLpub | 282 | 133,466 | 8,585 | 3.65 | 242 | 241 | 27,840 | 22,237 | 8,317 | 6.26 | 5 | 18.9% | 167.4 |
| CapScan matched | MG | 282 | 133,466 | 149,418 | 36.9 | 2,925 | 2,506 | 16,078 | 13,746 | 8,224 | 2.59 | 2 | 9.8% | 102.9 |
| CapScan matched | FLcmg | 282 | 133,466 | 8,613 | 3.61 | 238 | 238 | 27,846 | 22,331 | 8,347 | 6.30 | 5 | 19.3% | 168.5 |
| CapScan matched | FLunion | 282 | 133,466 | 452,255 | 189.7 | 6,496 | 5,689 | 27,099 | 23,278 | 199,663 | 3.45 | 2 | 8.5% | 175.0 |
| CapScan all | FLpub | 347 | 133,235 | 7,917 | 3.33 | (282 only) | 221 | 27,304 | 21,744 | 7,678 | 6.25 | 5 | 18.8% | 163.0 |
| CapScan | FLcmg | 331 | 133,190 | 7,941 | 3.35 | 218 | 218 | 26,777 | 21,318 | 7,648 | 6.24 | 5 | 19.0% | 159.4 |
| CapScan | FLunion | 331 | 133,190 | 452,255 | 189.7 | 6,434 | 5,561 | 26,171 | 22,144 | 196,362 | 3.39 | 2 | 8.3% | 166.8 |
| Healthy stool | FLown | 58 | 190,238 | 12,213 | 5.37 | 237 | 237 | 31,609 | 26,099 | 11,705 | 7.14 | 5 | 20.6% | 137.4 |
| Healthy stool | FLunion | 58 | 190,238 | 153,267 | 64.7 | 1,746 | 1,494 | 33,916 | 29,188 | 112,789 | 5.17 | 3 | 13.3% | 153.7 |

cand. paired = median Sage scored_candidates on rank-one spectra reported by both arms of the pair (FLpub-MG, FLcmg-FLunion, FLown-FLunion).
cand. all = all rank-one spectra reported by that arm. IDs/1000 MS2 at 5-50 is in the TSVs.

## Paired ratios (median of per-acquisition ratio; higher/lower/equal)

| metric | PREDICT FLpub/MG (296) | CapScan FLpub/MG (282) | PREDICT FLunion/FLcmg (296) | CapScan FLunion/FLcmg (330-331) | Healthy stool FLunion/FLown (58) |
|---|---|---|---|---|---|
| DB proteins | 0.034 (0/296/0) | 0.052 (0/282/0) | 68.5 (296/0/0) | 57.0 (331/0/0) | 12.6 (58/0/0) |
| DB residues | 0.055 (0/296/0) | 0.086 (0/282/0) | 68.2 (296/0/0) | 56.6 (331/0/0) | 12.1 (58/0/0) |
| scored_candidates, paired | 0.046 (0/296/0) | 0.069 (0/282/0) | 31.5 (296/0/0) | 27.6 (331/0/0) | 6.93 (58/0/0) |
| peptides 5-50 | 1.19 (294/2/0) | 1.52 (282/0/0) | 0.93 (30/266/0) | 0.96 (104/226/0) | 1.07 (58/0/0) |
| peptides 9-50 | 1.09 (272/24/0) | 1.42 (281/1/0) | 1.01 (166/130/0) | 1.04 (208/122/0) | 1.11 (58/0/0) |
| proteins with a 9-50 peptide | 0.61 (21/275/0) | 0.83 (94/188/0) | 26.9 (296/0/0) | 24.0 (330/0/0) | 9.29 (58/0/0) |
| peptides per protein, mean | 2.22 (296/0/0) | 2.41 (282/0/0) | 0.56 (0/296/0) | 0.54 (0/330/0) | 0.75 (0/58/0) |
| peptides per protein, median | 2.5 (290/0/6) | 3.0 (281/0/1) | 0.5 (0/296/0) | 0.5 (0/330/0) | 0.67 (0/57/1) |
| median coverage | 2.15 (296/0/0) | 1.92 (282/0/0) | 0.45 (0/296/0) | 0.44 (0/330/0) | 0.68 (0/58/0) |
| IDs per 1000 MS2, 9-50 | 1.09 (272/24/0) | 1.42 (281/1/0) | 1.01 (166/130/0) | 1.04 (208/122/0) | 1.11 (58/0/0) |

One CapScan acquisition (20210513_EXPL0_PeTr_SA_T760_1341_20210528032734) has 0 FLcmg peptides and is excluded from those ratios (n = 330).

## Is sequence coverage of identified proteins more complete with FL than MG?

Yes. Median per-protein coverage by accepted 9-50 peptides: PREDICT 17.5% (FL) vs 8.2% (MG), median paired ratio 2.15, FL higher in 296 of 296;
CapScan 18.9% vs 9.8%, ratio 1.92, FL higher in 282 of 282. Median peptides per identified protein: 4 vs 2 (PREDICT), 5 vs 2 (CapScan).

## What is measured, and its limits

Each row reads the existing Sage output of each arm on the same mzML. Accepted peptides: rank 1, target, peptide q <= 0.01, modifications
stripped, I/L-equivalent, distinct; 9-50 residues (5-50 alongside). scored_candidates is Sage 0.14.6's per-spectrum count of database
peptides within the precursor tolerance that were scored against the spectrum (targets and decoys). It exists only for spectra Sage
wrote out (report_psms = 1, rank one); paired medians use spectra reported by both arms. FLpub is a same-spectrum selection: its database
is drawn from the acquisition's own top-30% de novo peptides, so its gain is the uncorrected one (cross-fit and split-half in Supp. Fig. 6).
"Proteins with a 9-50 peptide" counts every database entry Sage lists for an accepted peptide. It is not a protein-group count, and redundant
entries (MG, study union) inflate it. Protein groups are omitted because no grouping was applied consistently to both the FL and MG arms
(the FastaLake post-search grouping exists only for FL databases). Coverage is per database entry in each arm's own searched FASTA, so the
two arms measure it on different sequence sets. Healthy-stool FLunion is a separate per-acquisition search against the 58-file union,
not the joint search.
