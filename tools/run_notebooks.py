#!/usr/bin/env python3
"""Execute repository notebooks headless, fail on any error, and clean stored outputs.

Each notebook runs with its own directory as the working directory and the repository
root on ``PYTHONPATH`` (so ``from style import fl_style`` and ``import panels`` work).
After a successful run, absolute paths in stored outputs are replaced: the repository
root by ``<repo>``, any other absolute path under a lab, user or temporary mount point
(the prefixes in ``_ABS``) by ``<path>``. The cleaned notebook is written back in place unless ``--no-write``.

Usage
-----
    python tools/run_notebooks.py                    # all notebooks
    python tools/run_notebooks.py figures/Fig2       # one unit (directory or .ipynb)
    python tools/run_notebooks.py --no-write --timeout 900
"""

import argparse
import os
import re
import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

REPO = Path(__file__).resolve().parents[1]

#: Notebook locations.
GLOBS = ["figures/*/*.ipynb", "tables/*/*.ipynb", "text_statistics/*/*.ipynb"]

#: Per-notebook wall-clock limit in seconds (CONVENTIONS.md: under 15 min).
DEFAULT_TIMEOUT = 900

_ABS = re.compile(r"(?<![A-Za-z0-9_.])/(?:fs|home|Users|Volumes|tmp|mpcdf|scratch)/[^\s'\"<>)\]]*")


def strip_paths(text):
    """Replace the repository root and any other absolute path in ``text``."""
    text = text.replace(str(REPO), "<repo>")
    return _ABS.sub("<path>", text)


def clean_outputs(nb):
    """Strip absolute paths from every stored output of ``nb`` in place.

    Returns
    -------
    int
        Number of output fields changed.
    """
    changed = 0
    for cell in nb.cells:
        for out in cell.get("outputs", []):
            for key in ("text", "evalue"):
                if key in out:
                    new = strip_paths(out[key])
                    changed += new != out[key]
                    out[key] = new
            if "traceback" in out:
                out["traceback"] = [strip_paths(t) for t in out["traceback"]]
            for mime, val in out.get("data", {}).items():
                if mime.startswith("text/") or mime == "application/json":
                    if isinstance(val, str):
                        new = strip_paths(val)
                        changed += new != val
                        out["data"][mime] = new
        cell.get("metadata", {}).pop("execution", None)
    return changed


def find(targets):
    """Resolve CLI targets (unit dirs or notebook files) to notebook paths."""
    if not targets:
        return sorted(p for g in GLOBS for p in REPO.glob(g) if ".ipynb_checkpoints" not in p.parts)
    out = []
    for t in targets:
        p = (REPO / t).resolve() if not Path(t).is_absolute() else Path(t)
        out += [p] if p.suffix == ".ipynb" else sorted(p.glob("*.ipynb"))
    return out


def run(path, timeout=DEFAULT_TIMEOUT, write=True, kernel="python3"):
    """Execute one notebook.

    Returns
    -------
    tuple
        ``(ok, seconds, message)``.
    """
    nb = nbformat.read(path, as_version=4)
    env_pp = os.pathsep.join([str(REPO), str(path.parent), os.environ.get("PYTHONPATH", "")])
    os.environ["PYTHONPATH"] = env_pp
    client = NotebookClient(nb, timeout=timeout, kernel_name=kernel,
                            resources={"metadata": {"path": str(path.parent)}})
    t0 = time.time()
    try:
        client.execute()
        ok, msg = True, "ok"
    except CellExecutionError as e:
        ok, msg = False, strip_paths(str(e))[-2000:]
    except Exception as e:  # timeouts, dead kernels
        ok, msg = False, f"{type(e).__name__}: {strip_paths(str(e))[-500:]}"
    secs = time.time() - t0
    if ok and secs > timeout:
        ok, msg = False, f"exceeded {timeout}s"
    if ok and write:
        clean_outputs(nb)
        nbformat.write(nb, path)
    return ok, secs, msg


def main(argv=None):
    """Command-line entry point."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("targets", nargs="*", help="unit directories or notebooks (default: all)")
    ap.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT)
    ap.add_argument("--no-write", action="store_true", help="do not write outputs back")
    ap.add_argument("--kernel", default="python3")
    a = ap.parse_args(argv)
    nbs = find(a.targets)
    if not nbs:
        print("no notebooks found")
        return 0
    failed = 0
    for p in nbs:
        ok, secs, msg = run(p, a.timeout, not a.no_write, a.kernel)
        failed += not ok
        print(f"{'PASS' if ok else 'FAIL'} {p.relative_to(REPO)} {secs:.0f}s" + ("" if ok else f"\n{msg}"))
    print(f"{len(nbs)} notebook(s), {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
