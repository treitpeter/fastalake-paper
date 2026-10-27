# Open items: Fig6

## Manuscript sentences that must change for v28f (values recomputed here, NUMBERS.tsv)

- MISMATCH: line 165 "5.0–5.5%" becomes "5.1–5.6%" (non-specific region means 5.092 to 5.639%), "all fourteen paired participants" becomes "all fifteen", "0.33%" becomes "0.34%" (0.3365%).
- MISMATCH: line 169 legend "(n=14)" becomes "(n=15)".
- MISMATCH: line 173 "40.3% and 40.8%" becomes "39.9% and 40.3%"; "all fourteen participants" becomes "all fifteen"; "18.0 percentage points" becomes "17.6" (15.2 unchanged).
- MISMATCH: lines 285 and 855 "all 331 CapScan" becomes "347" (347 acquisitions in `capscan_acquisition_design.tsv`).
- MISMATCH: lines 287, 863 and 869 "fourteen participants" becomes "fifteen"; line 869 "whereas descriptive stool summaries may include fifteen" then no longer describes a difference and could be removed.

## Statistical details the manuscript does not state (filled from the code in STATS.tsv)

- TEXT-GAP: exact P and q of the Fig. 6a,b tests are not printed: P = 6.1035e-05 (2 of 32,768 sign vectors, exact), q = 0.000290 for all five legend tests (STATS.tsv rows 1-5); the figure prints "q < 0.001".
- TEXT-GAP: paired effect sizes with intervals for Fig. 6a,b are not stated (for example non-specific route 2.22 pp, 95% CI in STATS.tsv row 3; host 15.20 pp, microbial 17.55 pp with CIs in rows 4-5).
- TEXT-GAP: seeds are not stated. A24: test exact (no seed); bootstrap weights from numpy default_rng seeded with the first 4 bytes (little-endian) of sha256("bootstrap|n|draws|20260914"), drawn once per n. A28: rng(20260914 + n) for random signs (n > 16), rng(20261014 + n) for the mean bootstrap, rng(20261114 + n) for the ratio bootstrap. Fig. 6c and point jitter: one default_rng(6) shared across panels a-c.
- TEXT-GAP: implementation is not named: sign-flip tests, bootstrap and Benjamini-Hochberg are custom NumPy code (moved to `panels.py`), not a library routine; recomputed with CPython 3.12.4, numpy 2.2.6, pandas 2.3.3 (originally CPython 3.10.9, numpy 1.26.4).
- TEXT-GAP: all tests are two-sided; the sign-flip statistic is the mean paired difference with |null mean| >= |observed mean| - 1e-12 counted as extreme; no degrees of freedom apply.
- TEXT-GAP: tests with fewer than eight paired participants remain in their family with p = 1 (stated in the A24/A28 reports, not in v28d).
- TEXT-GAP: per-stratum n for d-f: CapScan n = 15 in every displayed region (fluid and stool), PREDICT n = 93.
- TEXT-GAP: the legend says "Bars (c) and symbols (d–f) show participant or patient means with 95% bootstrap intervals"; panel c uses 2,000 resamples and d-f 20,000 (stated only in Supplementary Methods line 869).

## Not traced or not recomputable here

- UNTRACED: "9–50" upper bound 50 (lines 169, 285, 855) is applied in the upstream peptide filters (A20/A26 code), not visible in `source_data/` (NUMBERS.tsv row 7).
- UNTRACED: "306 PREDICT acquisitions" (lines 285, 855) needs the PREDICT acquisition roster, which is not part of this unit's source data.
- UNTRACED: "n≤16" (line 869) is a code constant (`EXACT_MAX_N`). The rule is supported on both sides: the PREDICT family contains exact-enumeration tests at n = 8 and n = 11 (for example P = 0.603515625) and the CapScan tests at n = 15 are exact, while tests with n >= 18 use random signs. No test has n = 16 or 17, so only the cutoff value itself is not exercised.
- MISMATCH: two of 18,144 recomputed A24 tests differ by one in the positive/negative sign count (tryptic route, Human origin: `H_fraction_nonKR` proximal minus stool, and `residue_P_fraction_nonKR` duodenum minus colon in the composition family). In each, one participant's difference is zero up to floating-point rounding, and summation order decides its sign. P, q, means and intervals agree to 1e-15. Neither test is shown or cited.
- UNTRACED: PREDICT panel f values (two rows) are read from the shipped table; patient-level complementary-product differences are not in this repository, so their means, intervals and P were not recomputed (STATS.tsv last two rows).
- UNTRACED: Slurm job ids of the figure build (`run_all.sh`) and of the PREDICT A28 run are not recorded.
- AUTHOR: line 165 "0.33%" becomes "0.34%" (0.3365%) under the definition the figure uses: the mean over all participant-region cells of regions 1-4 and stool (75 cells, n = 15 each), which is the dashed tryptic line in panel a. The sentence places it next to the fluid-region and stool values; the mean over regions 1-4 only is 0.3215% (0.32%). The v28d text does not say which definition 0.33% used. The author chooses the wording and the definition.
- AUTHOR: Methods line 855 states Sage v0.14.7, while the participant-15 non-specific search preflight (C2 fig6/A05) reports Sage 0.14.6; confirm the binary (inventory gap, not checked by this unit).
- AUTHOR: the Methods (lines 287, 861) say mixed and unknown origins are kept in Source Data; here they are rows of `capscan_terminal_participant_region.tsv` (origins Mixed, Unknown, NotInSelectedFASTA), not a dedicated table.
- AUTHOR: in the v28f composite, the overlaid labels "Region 3" and "Region 4" in panel d are wider than the words they replace and run into the y-axis line. Reproduced as is (no restyling); an author decision.

## v30 sync (2026-10-06)

- TEXT-GAP: NUMBERS.tsv keeps the seven v28d columns and appends the v30 columns `v30_para`, `v30_quote`, `v30_printed_value`, `v30_status` that `tools/check_numbers.py` accepts; every v30 quote is checked against FINAL `v30/visible_v30f.txt` (0 failures; 9 OPEN = 4 v28d OPEN, 4 v30 OPEN on the same values, 1 NOT_IN_V30). Every `v28f_change` value equals the v30 printed value.
- TEXT-GAP: v30 prints the paired counts as digits ("all 15 paired participants", paragraph 85; "in all 15 participants", paragraph 89); rows 165 and 173 now carry `v28f_change:15` and the code returns 15. Paragraphs 150, 441 and 444 keep the word "fifteen".
- TEXT-GAP: the v28d line 869 phrase "descriptive stool summaries may include fifteen" is not in v30 paragraph 444; its NUMBERS.tsv row is kept (MATCH against v28d) and no longer has a v30 counterpart.
- TEXT-GAP: v30 paragraph 151 adds "exact enumeration for fifteen participants; 200,000 random sign assignments for 93 PREDICT patients"; these equal the existing rows (15 from `stat_terminal_families`, 200,000 and 93 from `stat_predict_patient`), no separate row added.
- PIXEL: the panel b key "Human-assigned" is now written by `panel_b` (6 Oct 2026), not covered and retyped on the composite as in the v30 producer. Against the v30 figure: 5,884 differing pixels at 300 dpi, all inside the key (x 1266-1882, y 535-589 px): glyph rasterisation and a ~1 px lower "Microbial" row. Wording, face, size and position agree. Panel b, builder and v28f checks now differ from their "Host" references in the key only; panels a, c-f identical. If pixel identity with the manuscript file is wanted, re-export Figure 6 from this notebook.
- No new source tables: the v30 figure uses the same n = 15 data as v28f; `source_data/Fig6/` is unchanged (PREDICT patients already P001-P093).

## Pixel checks

- PIXEL: none open. All six panel crops, the builder PDF/PNG, the v28f PDF/PNG and the v30 PDF/PNG (`expected/Fig6_v30.*`) are IDENTICAL to the producer outputs with Söhne fonts. The repository style module differs from the producer's (1 pt vs 0.6 pt lines), so `panels.producer_style` reproduces the producer's settings instead of calling only `fl.init()`.

## Shared requests

- ENV: the post-processing (`panels.postprocess`, `panels.relabel_regions`) embeds `soehne-buch.ttf` and `soehne-halbfett.ttf` from `FASTALAKE_FONT_DIR`; `style/README.md` lists only the `.otf` files. Please add the two `.ttf` files to the font instructions.
- ZENODO: the complete A24 participant-region table (all lengths, `per_participant_region_endpoints.tsv.gz`, 4.65 MB gzip, sha256 `6d998e903173e85cff0166be599d05fc0005d1b462bf29101a8a93fd3f746dac`) and the A24 acquisition and specimen tables (`per_acquisition_endpoints.tsv.gz`, `per_specimen_endpoints.tsv.gz`) for Source Data completeness; the repository keeps only `min_length == 9`.
- CONTROLLED: predict_pseudonyms.tsv | raw_patient_id, pseudonym | cbc1ae4a900e541c3422f0c7ce91742e19a34679cc5b7fb24094573c51e57f27 | figures/Fig6/prepare_source_data.py (writes it; no notebook cell reads it) | map between PREDICT raw patient ids and the pseudonyms P001-P093 in `source_data/Fig6/predict_histidine_patient_enrichment.tsv`; other units that pseudonymise PREDICT patients (for example SFig13) should use the same map.
- PIPELINE: terminal_participant_inference, terminal_histidine_inference, histidine_boundary_localisation, terminal_route_counts, denovo_threshold_routes and capscan_nonspecific_search have no `pipelines/` directory yet; their statistics parts used by Figure 6 are moved into `panels.py`, the upstream counting steps are not in the repository.
- Source Data package: the submission package `Figure6/` CapScan sheets are still the n = 14 build and `SOURCE_DATA_MAP_reanchored.tsv` rows RES-89 to RES-92, FIG6-f and FIG6-c quote n = 14 values; replace with the tables in `source_data/Fig6/`.

## v33 sync (2026-10-07)

- RESOLVED: the v33 Figure 6 (producer `final/v30/figures/Fig6_v31`: shared d-f route key raised 11 pt, the 6 pt blank band cut at the page bottom, 21 white rows cropped after flattening; 2162 x 1884 px) is `expected/Fig6_v33.{png,pdf}`. The notebook now builds it (`panels.raise_route_key`, `postprocess(..., crop_bottom_pt=6)`, `panels.crop_bottom`) after the v28f/v30 steps, which are unchanged. Pixel check against the v33 figure: 5,884 differing pixels (PNG), all in the panel b key (x 1266-1882, y 535-589 px), the same known difference as the v30 check (key drawn as "Human-assigned" by the code instead of cover-and-retype).
