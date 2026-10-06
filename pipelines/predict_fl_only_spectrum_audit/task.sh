#!/bin/bash
# One acquisition: two rank-5 Sage searches (public-only FL, specimen MG+human), then the R8A audit.
# usage: task.sh SAMPLE [RESULTS_DIR] [OUT_DIR]
set -uo pipefail
A=${FASTALAKE_OUT_DIR:?set FASTALAKE_OUT_DIR}
CODE=$(cd "$(dirname "$0")" && pwd)
SAGE=${FASTALAKE_SAGE:?set FASTALAKE_SAGE (Sage 0.14.6 binary)}
PY=${FASTALAKE_PYTHON:-python3}
S=$1; RES=${2:-$A/results}; OUT=${3:-$A/out}
mkdir -p "$OUT"
[ "$(sha256sum "$SAGE" | cut -d' ' -f1)" = "d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625" ] || { echo "FAIL sage binary hash"; exit 1; }
read -r MZ FL MG <<< "$($PY -c "
import json,sys
m=json.load(open('$A/inputs/manifest.json'))
t=[t for t in m['tasks'] if t['sample']=='$S']; assert len(t)==1, '$S not in manifest'
t=t[0]; print(t['mzml'], t['FL_fasta'], t['MG_fasta'])")"
for f in "$FL" "$MG" "$MZ"; do [ -e "$f" ] || { echo "FAIL $S: missing $f"; exit 1; }; done
[ -s "${OUT}/${S}.json" ] && { echo "already done $S"; exit 0; }
W=${RES}/${S}; mkdir -p "$W"
echo "sample=$S start=$(date -Is) host=$(hostname) job=${SLURM_JOB_ID:-none}"
for ARM in FL MG; do
  DBF=$([ "$ARM" = "FL" ] && echo "$FL" || echo "$MG")
  D=$W/$ARM; rm -rf "$D"; mkdir -p "$D"
  # identical to R8A_loss_psm5.sbatch: production search settings, lfq off, report_psms 5
  cat > "$D/config.json" <<CFG
{ "database": { "bucket_size": 32768,
    "enzyme": { "missed_cleavages": 2, "min_len": 5, "max_len": 50, "cleave_at": "KR", "restrict": "P" },
    "fragment_min_mz": 200, "fragment_max_mz": 2000, "peptide_min_mass": 500, "peptide_max_mass": 5000,
    "ion_kinds": ["b","y"], "min_ion_index": 2, "static_mods": { "C": 57.0215 },
    "variable_mods": { "M": [15.9949] }, "max_variable_mods": 2,
    "decoy_tag": "rev_", "generate_decoys": true, "fasta": "${DBF}" },
  "quant": { "tmt": null, "lfq": false },
  "precursor_tol": { "ppm": [-20, 20] }, "fragment_tol": { "ppm": [-20, 20] },
  "isotope_errors": [-1, 3], "deisotope": false, "chimera": false, "predict_rt": true,
  "wide_window": false, "min_peaks": 15, "max_peaks": 150, "min_matched_peaks": 4,
  "max_fragment_charge": 2, "report_psms": 5, "output_directory": "${D}",
  "mzml_paths": [ "${MZ}" ] }
CFG
  /usr/bin/time -v -o "$D/time.txt" "$SAGE" "$D/config.json" > "$D/sage.log" 2>&1 || { echo "FAIL SAGE $S $ARM"; exit 1; }
  $PY -c "
import json; j=json.load(open('$D/results.json'))
assert j['version']=='0.14.6' and j['report_psms']==5 and j['quant']['lfq'] is False, j
assert j['database']['fasta']=='$DBF' and j['mzml_paths']==['$MZ'], j" || { echo "FAIL results.json check $S $ARM"; exit 1; }
  grep -m1 'Maximum resident' "$D/time.txt" | sed "s/^/  $ARM /"
done
$PY "$CODE/audit.py" "$W/FL/results.sage.tsv" "$W/MG/results.sage.tsv" "$MG" "$S" "${OUT}/${S}.json.tmp" || { echo "FAIL audit $S"; exit 1; }
mv "${OUT}/${S}.json.tmp" "${OUT}/${S}.json"
gzip -f "$W/FL/results.sage.tsv" "$W/MG/results.sage.tsv"
rm -f "$W"/*/results.sage.pin
echo "done $S end=$(date -Is)"
