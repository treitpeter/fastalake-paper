# FastaLake paper resources

[**Read the paper on bioRxiv**](https://www.biorxiv.org/content/10.64898/2026.10.08.757395v1) · [DOI: 10.64898/2026.10.08.757395](https://doi.org/10.64898/2026.10.08.757395)

[![Checks](https://github.com/treitpeter/fastalake-paper/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/treitpeter/fastalake-paper/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.12-1f4e5f.svg)](environment/README.md)
[![Code licence: MIT](https://img.shields.io/badge/code-MIT-1f4e5f.svg)](LICENSE)

**[Quickstart](#quickstart) · [Figures and tables](#figures-and-tables) · [Reviewer guide](docs/REVIEWER_GUIDE.md) · [Data access](docs/DATA_AVAILABILITY.md) · [FastaLake software](https://github.com/MannLabs/fasta-lake)**

Analysis code, notebooks and source tables for *FastaLake selects compact
metaproteomic databases from large sequence reservoirs by de novo peptide lookup*.
Use this repository to inspect the evidence behind a result or rebuild figures
and tables. To run FastaLake on your own spectra, start with the
[software repository](https://github.com/MannLabs/fasta-lake).

This repository includes improvements since the posted preprint: stricter held-out
testing, additional error controls, corrected functional annotations and expanded
validation. **The bioRxiv manuscript update is pending.** See the short
[change notes](CHANGELOG.md) for the results that changed, and cite the commit used
alongside the article.

## Quickstart

Use Git and Python 3.12 on Linux. The pinned environment was tested with
CPython 3.12.4; the notebooks do not require a GPU.

```bash
git clone https://github.com/treitpeter/fastalake-paper.git
cd fastalake-paper
python3.12 -m venv .venv
.venv/bin/python -m pip install --upgrade pip==25.2
.venv/bin/python -m pip install -r environment/requirements.lock.txt
.venv/bin/python -m ipykernel install --prefix .venv --name python3
```

Rebuild Supplementary Table 2 from the included source data:

```bash
.venv/bin/python paper.py build STab2
```

The workbook is written to `build/STab2/Supplementary_Table_2.xlsx`. Its notebook
compares the rebuilt cells with the reference workbook. The recorded run took
93 seconds on four CPU cores with 8 GB RAM; timing varies by machine.
See [environment setup](environment/README.md) for dependencies and optional inputs.

## Figures and tables

| Resource | Start here |
|---|---|
| Main figures | [1](figures/Fig1/) · [2](figures/Fig2/) · [3](figures/Fig3/) · [4](figures/Fig4/) · [5](figures/Fig5/) · [6](figures/Fig6/) |
| Supplementary figures | [Browse all figure units](figures/) (`SFig01`–`SFig15`) |
| Supplementary tables and data | [Browse tables](tables/) (`STab1`–`STab6`, `SData1`) |
| Numbers stated in the text | [Text statistics](text_statistics/) |
| Input tables and definitions | [Source data](source_data/README.md) |
| Upstream searches and analyses | [Pipeline index](pipelines/INDEX.md) |
| Registered claims and statistical methods | [Reviewer guide](docs/REVIEWER_GUIDE.md) · [Statistics](docs/STATISTICS.md) |

Each analysis unit includes a README, a notebook, source references and its known
limitations in `OPEN_ITEMS.md`. Quantitative figure panels are rebuilt from source
tables; schematic artwork is supplied with the unit.

## Reproduce and check

Run these commands from the repository root after installation:

```bash
.venv/bin/python paper.py list          # List units and outstanding entries
.venv/bin/python paper.py build Fig2    # Rebuild one figure
.venv/bin/python paper.py build-all     # Execute all 33 notebooks
.venv/bin/python paper.py check         # Scan files and check registered values
```

`paper.py build` and `build-all` update notebook outputs and write generated files
under `build/`. To execute notebooks without rewriting their stored outputs:

```bash
.venv/bin/python tools/run_notebooks.py --no-write tables/STab2
```

The number checker compares registered values. To check quotations too, supply
the current reference text, `visible_current.txt`, with one paragraph per line
in `<index><TAB><text>` format:

```bash
.venv/bin/python tools/check_numbers.py --manuscript-text visible_current.txt
.venv/bin/python tools/check_updated_analyses.py
```

The reference manuscript text is supplied separately; its checksum is in
[MANUSCRIPT.json](docs/MANUSCRIPT.json). The `v42_*` fields describe the current
reference, with earlier references retained for reproducibility. Unresolved
entries are reported as `OPEN`; add `--strict` to make them fail the check.
The [updated-analysis tables](source_data/v42_updates/) include the additional
held-out and error-control results. Add `--bootstrap` to their checker to repeat
the participant-bootstrap intervals.

## Data and figure reproduction

The public notebooks use tables in `source_data/`. Large inputs, raw acquisitions
and controlled PREDICT clinical data have separate access requirements;
see [Data availability](docs/DATA_AVAILABILITY.md) and the
[controlled input list](controlled/CONTROLLED_INPUTS.md).
Set `FASTALAKE_CONTROLLED_DIR` to an authorized input directory for controlled
calculations. Missing inputs are reported as skipped, so public CI does not
validate those calculations. Upstream search pipelines are provided for reference
and are not run by the notebook commands.

The original Söhne fonts are licensed and are not redistributed. With the default
fonts, numeric checks still run and exact pixel comparisons are reported as
skipped. See [figure style](style/README.md) to reproduce the original typography.

## Citation and licence

Treit et al. (2026). *FastaLake selects compact metaproteomic databases from large
sequence reservoirs by de novo peptide lookup*.
[bioRxiv](https://doi.org/10.64898/2026.10.08.757395).
See [CITATION.cff](CITATION.cff) for the article metadata, and record the analysis
commit used. The FastaLake software cited in the manuscript is
[v1.1.0rc2](https://github.com/MannLabs/fasta-lake/releases/tag/v1.1.0rc2),
[archived separately on Zenodo](https://doi.org/10.5281/zenodo.23215189).

Code is under the [MIT licence](LICENSE). The proposed source-data licence is
[CC BY 4.0](LICENSE-DATA.md), subject to author confirmation of redistribution rights.
Report reproducible problems through [Issues](https://github.com/treitpeter/fastalake-paper/issues).
