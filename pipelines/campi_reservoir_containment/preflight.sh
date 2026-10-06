#!/bin/bash
# R8 Task B preflight: every input path, size, mtime and (small-file) sha256 against inputs/manifest.json,
# the reference census, the reservoir header convention, the extractor binary, and a synthetic
# positive/negative control (overlap.py selftest: shuffled proteome must give 0 presence).
set -euo pipefail
D=${FASTALAKE_OUT_DIR}
PY=${FASTALAKE_PYTHON:-python3}
MAN=${R8B_MANIFEST:-$D/inputs/manifest.json}   # override only for the negative preflight test
fail() { echo "PREFLIGHT FAIL: $*" >&2; exit 1; }
[ -x "$PY" ] || fail "python missing: $PY"
"$PY" -c 'import sys; assert sys.version_info[:2]==(3,10), sys.version' || fail "python is not 3.10"
[ -f "$MAN" ] || fail "manifest missing"
[ -w "$D/results" ] && [ -w "$D/logs" ] || fail "results/ or logs/ not writable"

# 1. every manifest input: exists, bytes and mtime unchanged, sha256 recomputed for files < 100 MB
"$PY" - "$MAN" <<'EOF' || fail "manifest check"
import hashlib, json, os, sys
man = json.load(open(sys.argv[1]))
for key, e in man["inputs"].items():
    p = e["path"]
    if not os.path.exists(p): sys.exit(f"missing: {key} {p}")
    st = os.stat(p)
    if st.st_size != e["bytes"]: sys.exit(f"size changed: {key} {st.st_size} != {e['bytes']}")
    if st.st_mtime_ns != e["mtime_ns"]: sys.exit(f"mtime changed: {key}")
    if st.st_size < 100_000_000:
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()
        if e["sha256"] is None: sys.exit(f"small file without sha256 in manifest: {key}")
        if h != e["sha256"]: sys.exit(f"sha256 changed: {key}")
        print(f"  ok  {key}  {st.st_size} B  sha256 {h[:12]}")
    else:
        print(f"  ok  {key}  {st.st_size} B  (sha256 {'recorded' if e['sha256'] else 'pending pilot step 0'})")
print("manifest: all inputs present and unchanged")
EOF

# 2. reference census: 8 OS strings with the member counts, 116 nonmember entries, 29,673 total
REF=$("$PY" -c "import json;print(json.load(open('$MAN'))['inputs']['reference_fasta']['path'])")
MEM=$("$PY" -c "import json;print(json.load(open('$MAN'))['inputs']['members_tsv']['path'])")
n=$(grep -c '^>' "$REF"); [ "$n" -eq 29673 ] || fail "reference has $n entries, expected 29673"
nn=$(grep '^>' "$REF" | grep -vc 'OS=' || true); [ "$nn" -eq 116 ] || fail "nonmember entries $nn != 116"
while IFS=$'\t' read -r member os cnt alias; do
  [ "$member" = member ] && continue
  c=$(grep -c "OS=$os" "$REF" || true)
  [ "$c" -eq "$cnt" ] || fail "OS '$os': $c entries, members.tsv says $cnt"
done < "$MEM"
echo "reference: 29673 entries, 8 member OS strings match campi_members.tsv, 116 nonmembers"

# 3. reservoir header convention on the first bytes (cheap read) and the recorded stats
LAKE=$("$PY" -c "import json;print(json.load(open('$MAN'))['inputs']['reservoir_fasta']['path'])")
JH=$("$PY" -c "import json;print(json.load(open('$MAN'))['inputs']['joint_headers_tsv']['path'])")
head -c 4000 "$LAKE" | grep -q '^>MGYG' || fail "first reservoir headers are not UHGG (MGYG) as recorded"
head -c 4000 "$JH" | head -1 | grep -q $'^sha256_hash\tn_sources\tkept_header\tall_source_headers$' || fail "joint headers TSV header row differs"
grep -q 'Unique sequences (lake): 375,333,654' "$(dirname "$LAKE")/CAMPI_MOCK_lake_stats.txt" || fail "lake stats do not state 375,333,654"
echo "reservoir: header convention and recorded size confirmed"

# 4. the production extractor runs
EXT=$("$PY" -c "import json;print(json.load(open('$MAN'))['inputs']['extractor_binary']['path'])")
[ -x "$EXT" ] || fail "extractor not executable"
"$EXT" --help >/dev/null 2>&1 || fail "extractor --help failed"
echo "extractor: $EXT ok"

# 5. synthetic positive/negative control (shuffled proteome must give 0 exact and 0 tryptic presence)
T=${TMPDIR:-/tmp}/r8b_selftest_$$
mkdir -p "$T"
"$PY" -B "${PIPELINE_DIR}/overlap.py" selftest --out "$T" > "$T/selftest.log" 2>&1 || { cat "$T/selftest.log"; fail "selftest"; }
grep -q '^SELFTEST PASS' "$T/selftest.log" || { cat "$T/selftest.log"; fail "selftest did not pass"; }
grep 'selftest' "$T/selftest.log" | sed 's/^/  /'
rm -rf "$T"
echo "PREFLIGHT PASS $(date -u +%Y-%m-%dT%H:%M:%SZ) host=$(hostname) job=${SLURM_JOB_ID:-none}"
