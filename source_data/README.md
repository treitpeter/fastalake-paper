# Source data

Every quantitative figure panel, table and text statistic is computed by repository code
from the files under this directory, and from nothing else (controlled clinical inputs
excepted, see `controlled/CONTROLLED_INPUTS.md`).

## Layout

    source_data/<Unit>/<table>.tsv
    source_data/<Unit>/README.md

`<Unit>` is the unit directory name (`Fig2`, `SFig04`, `STab1`, `methods`, ...). A table
used by several units lives under the unit that first needs it, and the others name that
path; it is never duplicated.

Tables are tab-separated UTF-8 with one header row, no index column and no comment lines.
Missing values are empty fields. Each file is at most 20 MB.

## Data dictionary

Each `source_data/<Unit>/README.md` describes every file in that directory:

    ## <table>.tsv
    sha256: <hex>
    rows: <n> (one row per <statistical or observational unit>)
    produced by: <pipeline or upstream analysis id from UPSTREAM.json>

    | column | unit | meaning |
    |---|---|---|

Units are physical units or counts (for example `count of peptides`, `%`, `log10
intensity`, `dimensionless`); identifiers have unit `id`.

## Larger files

Files above 20 MB, and the complete Source Data package, are deposited on Zenodo and
listed in `ZENODO_MANIFEST.tsv` (columns: `zenodo_path`, `bytes`, `sha256`, `unit`,
`description`, `used_by`). Code that needs such a file reads it from `build/zenodo/`
after download and checks its sha256.

Licence: CC BY 4.0 (`LICENSE-DATA.md`).
