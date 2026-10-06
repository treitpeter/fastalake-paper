# Environment

| item | value |
|---|---|
| operating system tested | Linux x86_64 (SUSE Linux Enterprise 15, kernel 6.4, glibc 2.38) |
| language | CPython 3.12.4 |
| direct dependencies | `requirements.txt` (pinned) |
| full lock | `requirements.lock.txt` (`pip freeze` of the tested environment, 66 packages) |
| non-standard hardware | none for the notebooks; upstream pipelines in `pipelines/` list their Slurm resources |
| typical install time | 66 s (4 CPU cores, 8 GB RAM, network install from PyPI, measured 2026-10-06) |

The versions are those of the interpreter that rendered the original figure panels, so
that pixel checks compare like with like. Some panels were produced by other interpreters;
where that matters, the unit's `PROVENANCE.md` names the version and the pixel receipt
records the difference.

## Install

    python3.12 -m venv .venv
    .venv/bin/python -m pip install --upgrade pip==25.2
    .venv/bin/python -m pip install -r environment/requirements.lock.txt
    .venv/bin/python -m ipykernel install --prefix .venv --name python3

Then, from the repository root:

    .venv/bin/python paper.py list
    .venv/bin/python paper.py build Fig2
    .venv/bin/python paper.py check

## Optional inputs

| variable | purpose |
|---|---|
| `FASTALAKE_FONT_DIR` | directory with the licensed Söhne fonts; needed only for pixel-exact text (see `style/README.md`) |
| `FASTALAKE_CONTROLLED_DIR` | directory with controlled clinical inputs (see `controlled/CONTROLLED_INPUTS.md`); cells that need them skip without it |
| `FASTALAKE_MANUSCRIPT_TXT` | accepted manuscript text; lets `tools/check_numbers.py` also verify quotes against line numbers |

## Adding a dependency

Unit agents do not edit these files. A unit that needs a package writes
`ENV: <package>==<version>, needed by <file:function>` in its `OPEN_ITEMS.md`; the
orchestrator adds it here and regenerates the lock.
