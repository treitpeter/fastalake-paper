"""Preflight: re-verify every input recorded in inputs/manifest.json (sha256 where recorded, bytes+mtime_ns
otherwise) plus the presence of the executables. Exit 1 on the first mismatch. MANIFEST env var overrides
the manifest path so the check can be tested against a tampered copy (negative control)."""
import os
import hashlib, json, os, sys
R = os.environ["FASTALAKE_OUT_DIR"]
MAN = os.environ.get("MANIFEST", R + "/inputs/manifest.json")
M = json.load(open(MAN))

def sha(p):
    """sha256 hex digest of a file."""
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()

fails = 0; n = 0
for section in ("C1", "C2"):
    for k, v in M[section].items():
        if not isinstance(v, dict) or "path" not in v:
            continue
        n += 1
        p = v["path"]
        if not os.path.exists(p):
            print(f"FAIL {section}.{k}: missing {p}"); fails += 1; continue
        st = os.stat(p)
        if st.st_size != v["bytes"]:
            print(f"FAIL {section}.{k}: bytes {st.st_size} != {v['bytes']}"); fails += 1; continue
        if "sha256" in v:
            d = sha(p)
            if d != v["sha256"]:
                print(f"FAIL {section}.{k}: sha256 {d} != {v['sha256']}"); fails += 1; continue
            print(f"ok   {section}.{k}: sha256 {d[:12]} bytes {st.st_size}")
        else:
            if st.st_mtime_ns != v["mtime_ns"]:
                print(f"FAIL {section}.{k}: mtime_ns {st.st_mtime_ns} != {v['mtime_ns']}"); fails += 1; continue
            print(f"ok   {section}.{k}: bytes {st.st_size} mtime_ns {st.st_mtime_ns} (no sha256 recorded; 61 GB)")
for exe in (R + "/bin/fasta_extractor", R + "/bin/parsimony_engine", os.environ.get("FASTALAKE_PYTHON", "python") + ""):
    n += 1
    if os.access(exe, os.X_OK):
        print(f"ok   executable {exe}")
    else:
        print(f"FAIL not executable {exe}"); fails += 1
print(f"PREFLIGHT {'FAIL' if fails else 'PASS'}: {n - fails}/{n} checks passed, manifest {MAN}")
sys.exit(1 if fails else 0)
