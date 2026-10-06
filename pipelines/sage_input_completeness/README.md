# sage_input_completeness

Internal ID: A98 (previous folder name A98_sage_input_completeness)

## Question

This pipeline produces: Closing tag, declared spectra and highest scored scan of the 653 public-sequence and 58 healthy-stool searches (search_input_completeness.tsv, line 629); reconverted flag of the 306 PREDICT searches and end-of-file check of original mzMLs (lines 233, 303)

Used by: text_statistics/methods, text_statistics/supplementary_methods_1_5.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `sage_input_integrity.py` | Independent Sage input-coverage audit (read-only). |
| `check.sh` | 1) ending of every PREDICT mzML in 0_DATA |
| `run653.py` | Extended task list for the job-6048085 Sage input-coverage audit. |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `search_input_completeness.tsv` | not found on disk |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

No job id could be traced for this analysis; see OPEN_ITEMS.md.

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| `CAN` | `FL / 'canonical/publication_data'` | no comment in original code | not stated in code; see sage_input_integrity.py |
| `SCAN_RE` | `re.compile(r'scan=(\d+)')` | no comment in original code | not stated in code; see sage_input_integrity.py |
| `COUNT_RE` | `re.compile(rb'<spectrumList[^>]*count="(\d+)"')` | no comment in original code | not stated in code; see sage_input_integrity.py |
| `SOFT_RE` | `re.compile(rb'<software id="([^"]+)" version="([^"]+)"')` | no comment in original code | not stated in code; see sage_input_integrity.py |
| `A08_347` | `ROOT / 'final/public347/A08_347/results'` | no comment in original code | not stated in code; see run653.py |
| `ORIG` | `S.searches` | no comment in original code | not stated in code; see run653.py |
