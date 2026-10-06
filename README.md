# FastaLake paper repository

Code, notebooks and source data that regenerate every figure panel and every number of the manuscript "FastaLake selects compact metaproteomic databases from large sequence reservoirs by de novo peptide lookup" (reference text v30). Every printed number is listed with the code and the source data table that produce it (`docs/REVIEWER_GUIDE.md`), every statistical test with its unit, n, correction and seed (`docs/STATISTICS.md`), and every choice (cutoffs, seeds, filters) in a Choices table in the README of the unit that makes it.

The manuscript text is not part of the repository. The reference text is v30 (`visible_v30f.txt`, one paragraph per line as `<index><TAB><text>`). Units in `tables/` and `text_statistics/` give for every number the v30 paragraph, wording, printed value and status (`NUMBERS.tsv` columns `v30_*`) and check them with `python tools/check_numbers.py --manuscript-text visible_v30f.txt`. `v28d_line` is kept as the physical line number of the earlier accepted text v28d (`CONVENTIONS.md`); figure units still refer to it.

## Install

Tested on Linux x86_64 (SUSE Linux Enterprise 15, kernel 6.4, glibc 2.38) with CPython 3.12.4. Measured install time: 66 s (4 cores, 8 GB RAM, packages from PyPI, 2026-10-06). No special hardware is needed for the notebooks.

    python3.12 -m venv .venv
    .venv/bin/python -m pip install --upgrade pip==25.2
    .venv/bin/python -m pip install -r environment/requirements.lock.txt
    .venv/bin/python -m ipykernel install --prefix .venv --name python3

Details, optional inputs and the dependency list are in `environment/README.md`.

## Rebuild

    .venv/bin/python paper.py list              # every unit with number and open-item counts
    .venv/bin/python paper.py build Fig2        # rebuild one figure (executes its notebook)
    .venv/bin/python paper.py build-all         # every figure, table and text-statistics notebook
    .venv/bin/python paper.py check             # release scan and number check

The notebooks read only `source_data/` (343 MB in the repository) and, where noted, controlled inputs. Figures are drawn with the same code as the submitted panels; schematics are included as supplied artwork and are never redrawn (`figures/<unit>/artwork/`, with their typed numbers checked in `ARTWORK_NUMBERS.md`). A rebuilt figure is compared with the reference render by `tools/check_pixels.py`; fonts are not included (`style/README.md`), so text can differ at pixel level without the fonts.

## Demo

The smallest notebook is `tables/STab2/STab2.ipynb` (regenerates Supplementary Table 2 from `source_data/STab2/`). Measured run time with `python tools/run_notebooks.py --no-write tables/STab2`: 93 s on 4 cores, 8 GB RAM, including kernel start (2026-10-06, Slurm). Large figure notebooks are given a limit of 15 minutes each.

## Layout

| path | content |
|---|---|
| `figures/<Fig1..Fig6, SFig01..SFig15>/` | notebook, plotting and statistics functions, `NUMBERS.tsv`, `STATS.tsv`, provenance, open items, reference renders |
| `tables/<STab1..STab6, SData1>/` | supplementary tables and data, regenerated from source data |
| `text_statistics/<section>/` | numbers stated in abstract, discussion, Methods and Supplementary Methods that no figure shows |
| `pipelines/<id>/` | the 94 upstream analyses (two of them placeholders for running analyses) that produced the source data (code as it ran, Slurm launchers, receipts with hashes and job states); they are not rerun by `paper.py` |
| `source_data/<unit>/` | the tables the notebooks read, each with a data dictionary |
| `controlled/` | list of controlled clinical inputs and how to request them; the data are not here |
| `style/`, `tools/`, `environment/` | figure style, checks, pinned dependencies |
| `docs/` | reviewer guide and statistics (v30 order), reporting summary draft, code and data availability, changes for v28f |

Pipelines added for v30: `clinical_second_engine` (MSFragger 4.1 and Percolator 3.06.5 on PREDICT, CapScan and healthy stool), `peptide_length_window_recounts` (9–50-residue recounts), `discussion_counts`, `cohort_database_accumulation`, `cohort_candidate_competition` and `metagenome_header_audit` (Discussion and Methods counts); `capscan_mg_benchmark_295` and `clinical_group_ko_homogeneity` are placeholders for analyses still running. `pipelines/INDEX.md` lists all of them with their manuscript section.

Start with `docs/REVIEWER_GUIDE.md`. Items not traced or not matching are in each unit `OPEN_ITEMS.md`; the repository does not hide them.

## Relation to the FastaLake software

The FastaLake software (database selection by de novo lookup, Rust and Python) is released separately at https://github.com/MannLabs/fasta-lake under the MIT licence. Methods and Supplementary Methods cite its version 1.1.0rc2 (commit 9ad09c9); `docs/CODE_AVAILABILITY.md` lists the cited paths and functions. This repository holds only the code and data of the paper: the analyses, figure notebooks and statistics built on the software outputs. AlphaNovo code is not included; its predictions are inputs (`docs/DATA_AVAILABILITY.md`).

## Controlled data

PREDICT patient-level clinical covariates are not in the repository. Public tables use random pseudonyms. Code that needs the controlled files reads them from the directory in `FASTALAKE_CONTROLLED_DIR`, skips the step with a message if a file is missing, and never substitutes other data. See `controlled/CONTROLLED_INPUTS.md`.

## Licence

Code: MIT (`LICENSE`). Source data: CC BY 4.0 (`LICENSE-DATA.md`). Both choices await confirmation by the authors. Citation: `CITATION.cff`.
