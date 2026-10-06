# Open items: SFig01

- MISMATCH: the v28d legend says "all 710 acquisition-specific FASTAs"; the figure plots 726; the v30 legend (v30f) reads 726 (v28d line 405; NUMBERS.tsv row 11).
- MISMATCH: the v28d legend says "306 PREDICT and 331 CapScan acquisitions"; the v30 legend (v30f) reads 347 (v28d line 405; NUMBERS.tsv row 17).
- MISMATCH: Methods say "we used the same 710 selected FASTAs" for the full-reservoir context; v30 Methods (v30f) read "the same 726 selected FASTAs" (v28d line 675; NUMBERS.tsv row 20).
- TEXT-GAP: the text does not state that horizontal point positions are jittered (uniform, half-width 0.14, seed 20260928); display only, no number depends on it (panels.py JITTER_SEED).
- ZENODO: the Source Data package still holds only the 710-row `Review_v16/Database_sizes/database_sizes_710.tsv` and has no `Supplementary_Figure_01/` folder; the 726-row table (`source_data/Fig1/acquisition_sizes_726.tsv`, sha256 fd9ade3b...) should be exported as Supplementary Figure 1 Source Data and map rows SF1, SF1-n and RES-08a-d re-pointed.
- UNTRACED: no Slurm job id is recorded for the producer runs (`SF01/run.sh`, the 28 Sep `build_evidence.py` run and `reexport.py`); outputs are identified by sha256 only (PROVENANCE.md).
- PIPELINE: cohort_database_sizes: the upstream size audit (audit.py, frozen/audit_v2.py scan(), and the participant-15 extension sizes.py/assemble.py) is not in `pipelines/`; the reservoir entry counts it read are not in this repository (UPSTREAM.json).
- AUTHOR: the `sample` column of the size table carries raw acquisition file names (operator initials in the names, no patient identifiers or clinical data); confirm these names are the ones deposited with the raw data.
- shared requests: `tools/check_numbers.py` fails rows where the text rounds to the nearest hundred: "3,100" (median 3057.0) and "12,200" (median 12212.5) on v28d line 55 are both the recomputed medians rounded to the nearest hundred, but the tool reads trailing zeros as exact integers. Please accept a printed integer with trailing zeros at the precision of its last non-zero digit, or add a status for stated rounding (NUMBERS.tsv rows 4 and 5).
- shared requests: this unit reads `source_data/Fig1/acquisition_sizes_726.tsv` (owned by Fig1) and checks its sha256; the Fig1 unit's `source_data/Fig1/README.md` must carry its data dictionary, and the file must not be renamed or changed.
- UNTRACED: no Slurm job id is recorded for the v30 house-style build (`build_sfig01_style.py`); outputs are identified by sha256 (PROVENANCE.md, v30).
- AUTHOR: the v30 PNG is not yet among the media of the v30f docx; placement in the next build is pending (ASSEMBLY.md, v30).
