# Conventions

Every unit agent follows this file exactly. If something here cannot be followed, write
the reason in the unit's `OPEN_ITEMS.md`; do not invent a variant.

## 1. Reference documents

| name | what | sha256 |
|---|---|---|
| v28d text | manuscript v28d, accepted text with line numbers (`v28d_accepted.txt`, 961 lines) | `7eec1945dc5348e94825d2563e999921ece1000dc913c014c4acbbdb6d365ff4` |
| v28d docx | manuscript v28d (`FastaLake_v28d.docx`) | `11168eee608e632845c88a8fb20e9105b76580c1b6364781b958fc1af6c6fd39` |
| v28f | v28d text plus the rebuilt figures queued for v28f | (set at release) |

`v28d_line` everywhere means the physical, 1-based line number in `v28d_accepted.txt`.
A range is written `123-125`. The manuscript is never copied into or edited from this
repository.

`LAB_ROOT` names the lab project directory that contains `FastaLake_v28d.docx`. Original
producer locations are written `LAB_ROOT/<relative path>`. Its absolute value never
appears in this repository.

## 2. Who may write where

| path | owner |
|---|---|
| `figures/<U>/`, `tables/<U>/`, `text_statistics/<U>/` | the agent assigned unit `<U>` only |
| `source_data/<U>/` | the agent assigned unit `<U>` only |
| `pipelines/<P>/` | the unit that creates it first (see below) |
| everything else (`README.md`, `CONVENTIONS.md`, `paper.py`, `tools/`, `style/`, `environment/`, `controlled/CONTROLLED_INPUTS.md`, `source_data/README.md`, `source_data/ZENODO_MANIFEST.tsv`, `docs/`, `.github/`) | the orchestrator only |

Pipelines are shared. A unit claims `pipelines/<P>/` with an atomic `mkdir` (without
`-p`); if that succeeds it writes `pipelines/<P>/OWNER` containing its unit name and
fills the directory. If the directory already exists, the unit does not write in it; it
references it from `UPSTREAM.json` and puts anything the pipeline lacks in its own
`OPEN_ITEMS.md` as `PIPELINE: <P>: <what is missing>`.

Requests for a shared file (a dependency, a controlled input, a Zenodo file) go in the
unit's `OPEN_ITEMS.md` with the prefixes in section 9. The orchestrator merges them.

Do not `git commit`. Do not modify any file outside this repository and its `_work`
sibling directory.

## 3. Names

- Units: `Fig1`..`Fig6`, `SFig01`..`SFig15`, `STab1`..`STab6`, `SData1`,
  `abstract_discussion`, `methods`, `supplementary_methods_1_5`,
  `supplementary_methods_6_7`, `supplementary_methods_8_11`.
- Panels: lower-case letters as printed (`a`, `b`, `d_inset` if the figure has an inset).
- Upstream analyses and pipelines: stable descriptive ids in lower-case `snake_case`, cohort or
  scope first and then the analysis, without internal codes or dates, for example
  `terminal_histidine_inference`. The internal lab code is kept only in the `Internal ID` line of
  the pipeline README and in `pipelines/INDEX.md`. The same id is used in `UPSTREAM.json`, `pipelines/<id>/` and
  every README. Before creating an id, look for an existing one in `pipelines/` and in
  other units' `UPSTREAM.json`; reuse it if it is the same analysis.
- Statistical units: exactly one of `acquisition`, `specimen`, `participant`, `patient`,
  `peptide`, `protein group`, `sequence`, `species`, `region`. Acquisitions are nested in
  specimens, specimens in participants.

## 4. Files of a figure unit, `figures/<U>/`

    README.md          purpose, panel table, Choices table, limitations
    <U>.ipynb          the notebook (one section per panel)
    panels.py          plotting and statistics functions
    NUMBERS.tsv        every printed number the unit covers
    STATS.tsv          every statistical test the unit covers
    ASSEMBLY.md        how the v28f composite was assembled from panels
    ARTWORK_NUMBERS.md only if the figure has schematics: every typed number traced
    PROVENANCE.md      original producer, job ids, hashes, Codex census reconciliation
    OPEN_ITEMS.md      everything not traced or not matching
    UPSTREAM.json      upstream analyses the unit relies on
    expected/          reference renders (inputs to pixel checks)
    artwork/           schematic panels as PDF, exactly as produced
    checks/            pixel-check receipts (*.json, *.diff.png) written by the notebook

Table units (`tables/<U>/`) have `README.md`, `<U>.ipynb` (regenerates the table from
`source_data/`), `NUMBERS.tsv`, `PROVENANCE.md`, `OPEN_ITEMS.md`, `UPSTREAM.json` and
`expected/` holding the table as published. Text-statistics units
(`text_statistics/<U>/`) have `README.md`, `<U>.ipynb`, `stats.py` (in place of
`panels.py`), `NUMBERS.tsv`, `STATS.tsv`, `OPEN_ITEMS.md`, `UPSTREAM.json`.

### 4.1 README.md

    # <Figure name as in the manuscript, e.g. Figure 2>

    <One paragraph: what the figure is for, in plain words.>

    ## Panels

    | panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
    |---|---|---|---|---|---|---|---|

    ## Choices

    | parameter | value | why | alternatives and where reported |
    |---|---|---|---|

    ## Known limitations

    - ...

"notebook cell" is the section heading in the notebook, for example `Panel b`. Schematic
panels write `artwork` in the test, cell and source_data columns and point to
`ARTWORK_NUMBERS.md`. Every value in the Choices table is a named constant in
`panels.py` (or `stats.py`) with the same value.

### 4.2 <U>.ipynb

1. First cell, markdown: what the figure shows (two to four sentences) and how to run it:
   `python paper.py build <U>` from the repository root.
2. Second cell, code: setup only. Paths are derived from the notebook location:

        from pathlib import Path
        import sys
        UNIT = Path.cwd()                 # the notebook runs with its own directory as cwd
        REPO = UNIT.parents[1]
        sys.path[:0] = [str(REPO), str(UNIT)]
        SD = REPO / "source_data" / UNIT.name
        OUT = REPO / "build" / UNIT.name   # rendered panels; gitignored
        OUT.mkdir(parents=True, exist_ok=True)
        from style import fl_style as fl
        import panels
        fl.init()

3. One section per panel, in figure order. Each section is a markdown cell headed
   `## Panel <x>` stating what is shown, which claim it supports (quote and v28d line),
   n with its unit, and the test; followed by code that only calls functions in
   `panels.py`, displays the panel and saves it to `OUT / "panel_<x>.png"` (and `.pdf`).
4. Final cell: runs the pixel checks for every quantitative panel (section 6) and prints
   every recomputed number in `NUMBERS.tsv` order.

Rules: reads only `source_data/` (and `FASTALAKE_CONTROLLED_DIR`, section 8); runs
headless in under 15 minutes on 4 cores; no random draw without a seed constant; stored
outputs contain no absolute paths (`tools/run_notebooks.py` strips them, but do not print
them in the first place); no `%`-magics that reach outside the repository.

### 4.3 panels.py

- Module docstring: which figure, which panels, which source_data files.
- Named constants at the top for every choice, each with a comment pointing to the
  Methods line, for example:

        SEED = 20260918            # Methods, v28d line 642
        N_PERMUTATIONS = 10_000    # Methods, v28d line 642
        Q_CUTOFF = 0.01            # Methods, v28d line 588
        LENGTH_WINDOW = (7, 30)    # Methods, v28d line 590
        ROSTER = "HSTOOL58"        # Methods, v28d line 571

- One function per panel (`panel_b(...)`) and one per statistic (`stat_b_...(...)`),
  NumPy-style docstrings with `Parameters` and `Returns`; the docstring of every
  statistics function names the statistical unit and n.
- Plot code is moved from the original producer as it is. Do not redraw or restyle.

### 4.4 NUMBERS.tsv

Tab-separated, this header exactly:

    v28d_line	quote	printed_value	recomputed_value	code	source_data	status

| column | content |
|---|---|
| `v28d_line` | line number or range in `v28d_accepted.txt` |
| `quote` | the exact text on that line around the number, enough to be unique |
| `printed_value` | the number exactly as printed, for example `0.318%`, `8,420`, `P = 3.2 × 10^-5` |
| `recomputed_value` | full-precision value printed by the notebook (no rounding) |
| `code` | `figures/<U>/panels.py:stat_b_fdp` |
| `source_data` | `source_data/<U>/<table>.tsv` |
| `status` | `MATCH`, `v28f_change:<new printed value>`, or `OPEN` |

One row per printed number, including numbers in figure legends. `tools/check_numbers.py`
checks each row; read its docstring for the agreement rule. `OPEN` means the number could
not be recomputed or traced; the reason goes in `OPEN_ITEMS.md`. Never fill
`recomputed_value` with anything the code did not compute.

### 4.5 STATS.tsv

One row per statistical test or interval. Tab-separated, this header exactly:

    claim	v28d_line	test	sided	unit	n	distinct_or_repeated	covariates	correction	family	family_size	estimate	uncertainty	uncertainty_construction	statistic	df	exact_P	effect_size	effect_size_how	draws	seed	software_version	code

| column | content |
|---|---|
| `claim` | the sentence the test supports, short |
| `test` | named test, for example `Wilcoxon signed-rank`, `PERMANOVA (Bray-Curtis)` |
| `sided` | `two-sided`, `one-sided (greater)`, `one-sided (less)`, or `NA` (interval only) |
| `unit` | statistical unit (section 3) at which the test is done |
| `n` | n per group, for example `58; 58` |
| `distinct_or_repeated` | `distinct` or `repeated: <what is repeated within what>` |
| `covariates` | covariates in the model, or `none` |
| `correction` | `none`, `Benjamini-Hochberg`, `Bonferroni`, ... |
| `family` | what the correction family is; `NA` if none |
| `family_size` | number of tests in the family |
| `estimate` | central value, with its name, for example `median difference 0.42` |
| `uncertainty` | for example `95% CI 0.31 to 0.55`, `IQR 3 to 7` |
| `uncertainty_construction` | for example `percentile bootstrap, 10,000 draws, resampling participants` |
| `statistic` | test statistic with its name, for example `W = 1211` |
| `df` | degrees of freedom or `NA` |
| `exact_P` | exact P as computed (not a bound) |
| `effect_size`, `effect_size_how` | value and how it was calculated |
| `draws`, `seed` | for permutation or bootstrap; `NA` otherwise |
| `software_version` | for example `scipy 1.16.3` |
| `code` | `figures/<U>/panels.py:stat_c_wilcoxon` |

Values come from the code. Where the Codex statistics register
(`STATISTICS_REGISTER.tsv`, `MISSING.md`) lists a gap the text does not state (seed, CI
construction, statistic, df, implementation, family size), fill it from the code and
note in `OPEN_ITEMS.md` as `TEXT-GAP:` that the manuscript does not state it. If the
code cannot answer it either, write `OPEN` in the cell and add an `OPEN_ITEMS.md` entry.

### 4.6 ASSEMBLY.md

How the v28f composite was made from panels, as an ordered list. Each step: input files
(repo-relative or `LAB_ROOT/...`), tool (script path and function, or `manual` with the
application), what changed, output file and its sha256. End with the sha256 of the
composite in `expected/<U>_v28f.png` and of the composite PDF it was rasterised from.

### 4.7 ARTWORK_NUMBERS.md (schematics only)

    | artwork file | typed text | value | table it comes from (repo path, row, column) | status |
    |---|---|---|---|---|

Every number typed in a schematic has a row. `status` is `MATCH` or `OPEN`.

### 4.8 PROVENANCE.md

Original producer locations are allowed only here, and only as `LAB_ROOT/...`.

    ## Panels
    | panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
    |---|---|---|---|---|---|

    ## Codex census reconciliation
    | panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
    |---|---|---|---|---|---|

    ## Source data
    | source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
    |---|---|---|---|

The Source Data package and `SOURCE_DATA_MAP_reanchored.tsv` are compared, not copied
blindly: each `source_data/<U>/` table states which package sheet it corresponds to and
whether the values agree.

### 4.9 OPEN_ITEMS.md

    # Open items: <U>

    - <PREFIX>: <one sentence> (<file, row or panel>)

Prefixes: `UNTRACED`, `MISMATCH`, `TEXT-GAP`, `PIXEL`, `ENV`, `CONTROLLED`, `ZENODO`,
`PIPELINE`, `AUTHOR`. If nothing is open, write `- none`.

### 4.10 UPSTREAM.json

    {
      "unit": "Fig2",
      "upstream": [
        {
          "id": "terminal_histidine_inference",
          "original_code_paths": ["LAB_ROOT/..."],
          "uses": "which outputs this unit uses, and for which panels",
          "pipeline_dir": "pipelines/terminal_histidine_inference"
        }
      ]
    }

`pipeline_dir` is `null` if the analysis is not (yet) in `pipelines/`.

## 5. Pipelines, `pipelines/<P>/`

    README.md     question; inputs with sha256; outputs; Slurm resources measured
                  (cores, peak RSS, wall time); run time; Choices table (section 4.1 format)
    OWNER         the unit that created the directory
    *.py, *.sh    the code, with module and function docstrings
    slurm/        launchers parametrised by environment variables
                  (FASTALAKE_DATA_DIR, FASTALAKE_OUT_DIR, ...); no absolute paths
    RECEIPTS.md   original Slurm job ids, output files (LAB_ROOT/...) and their sha256

Upstream code is moved as it ran. It is not rewritten, and a pipeline is not re-run unless
the task says so.

## 6. Pixel checks

- `expected/panel_<x>.png`: the panel as rendered by the ORIGINAL producer step (copied,
  not re-rendered), and `expected/<U>_v28f.png`: the v28f composite. Each at most 20 MB.
- The notebook's final cell calls

        from tools.check_pixels import compare
        compare(OUT / "panel_b.png", UNIT / "expected" / "panel_b.png",
                UNIT / "checks" / "panel_b", dpi=300)

  which writes `checks/panel_b.json` and `checks/panel_b.diff.png`.
- Compare at the producer's resolution. If the producer saved a PDF, compare PDFs at
  300 dpi.
- Pixel-exact (`IDENTICAL`) is the target. Text pixels can match only with the Söhne
  fonts (`FASTALAKE_FONT_DIR`). If fonts or a library version force a difference, keep
  the receipt and diff image, describe what differs in `OPEN_ITEMS.md` as `PIXEL:` and
  state that the check was visual. Never loosen the tolerance to hide a difference.
- Without `FASTALAKE_FONT_DIR`, the cell prints `pixel checks SKIPPED: fonts absent` and
  leaves committed receipts untouched (CI runs this way).

## 7. Numbers and values

- Never invent, round in advance, or copy a value from prose into a result. A value is
  either computed by repository code from `source_data/` or it is `OPEN`.
- Never subsample. Full data always.
- Anything beyond trivial runs on Slurm, never on the login node:

        srun -p p.hpcl92,p.hpcl8,p.hpcl93 --nice=10000 -c4 --mem=<measured>G -t <time> \
            <env>/bin/python ...

  The environment is `<work>/env` (see `environment/README.md`); `<work>` is the
  `fastalake-paper_v28_work` directory beside this repository. Logs and large
  intermediates go to `<work>/<U>/`, never into the repository.

## 8. Controlled inputs

PREDICT patient-level clinical covariates never enter the repository: no age, sex, MELD,
CRP, antibiotic exposure, visit dates, tube identifiers or raw patient identifiers, in any
table, notebook output or log. Public tables carry aggregates and de-identified
non-clinical per-patient values that a plot needs, keyed by pseudonymised ids.

A unit that needs a clinical input reads it as shown in `controlled/CONTROLLED_INPUTS.md`
and adds `CONTROLLED: <file> | <columns> | <sha256> | <notebook cell or function> | <purpose>`
to its `OPEN_ITEMS.md`; the orchestrator copies it into `CONTROLLED_INPUTS.md`.

Nothing from the old paper repository (`fasta-lake-paper_2026-10-06`) is copied if it
contains absolute lab paths, e-mail addresses, personal-genome references (the identifiers `tools/scan_release.py` rejects),
co-author review comments, vendored AlphaNovo code or patient clinical data.

## 9. Prose

Plain and short. No marketing, no bold, no em-dashes. Say what was done and what was
found. Numbers in prose carry their unit and n.

## 10. Before finishing a unit

    <env>/bin/python tools/scan_release.py          # 0 findings
    <env>/bin/python tools/check_numbers.py --unit <U>
    <env>/bin/python paper.py build <U>             # via srun if it takes more than a minute
    <env>/bin/python paper.py list                  # 0 missing files for <U>
