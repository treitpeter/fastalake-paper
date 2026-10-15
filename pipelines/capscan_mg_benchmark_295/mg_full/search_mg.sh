#!/bin/bash -l
# (b) same-spectrum MG Sage search. Copy of canonical/publication/arms/search_metag_hashed.sbatch
# (FL_MG_ARM=metag_hashed_human, production job 5674719). Changes, and only these:
#   - manifest = mg295/mg_full/manifest.tsv; row N+1 for task N (task 0 = control); 4th column
#     template_config replaces metag_solo/<sample>/config.json (participant 15 has no metag_solo
#     config; all 283 metag_solo CAPSCAN configs are identical apart from fasta/mzml/output, checked).
#   - OUTROOT = mg295/mg_full/out; code_commit line omitted (not run from canonical).
#   - after Sage: count accepted peptides with the frozen rule (count.py).
# Same Sage binary (${FASTALAKE_SAGE:-sage}, sha256 d3820543...), same CLI (--write-pin -o OUT).
set -uo pipefail
M=LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295
SAGE=${FASTALAKE_SAGE:-sage}
SAGE_SHA=d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625
PY=${FASTALAKE_PYTHON:-python3}
ARM=metag_hashed_human
MANIFEST=${M}/mg_full/manifest.tsv
OUTROOT=${M}/mg_full/out
[ "$(sha256sum ${SAGE} | cut -d' ' -f1)" = "${SAGE_SHA}" ] || { echo "FATAL: sage binary changed"; exit 1; }
IDX=${MG_TASK:-${SLURM_ARRAY_TASK_ID:-0}}
LINE=$(sed -n "$((IDX + 2))p" "${MANIFEST}")
[ -n "${LINE}" ] || { echo "FATAL: no manifest row for task ${IDX}"; exit 1; }
IFS=$'\t' read -r SAMPLE MZML FASTA OLDCFG <<< "${LINE}"
case "${SAMPLE}${MZML}${FASTA}${OLDCFG}" in *$'\r'*) echo "FATAL: CR in manifest"; exit 1 ;; esac
OUT=${OUTROOT}/${SAMPLE}
echo "task ${IDX} sample ${SAMPLE}"; echo "mzml ${MZML}"; echo "fasta ${FASTA}"; echo "template ${OLDCFG}"; echo "out ${OUT}"
[ -f "${MZML}" ]  || { echo "FATAL: mzML missing"; exit 1; }
[ -f "${FASTA}" ] || { echo "FATAL: FASTA missing"; exit 1; }
[ -f "${OLDCFG}" ] || { echo "FATAL: template config missing"; exit 1; }
rm -rf "${OUT}"; mkdir -p "${OUT}"
CFG=${OUT}/config.json
"${PY}" - "${OLDCFG}" "${FASTA}" "${MZML}" "${OUT}" "${CFG}" <<'PYEOF'
import json, sys
old, fasta, mzml, out, dst = sys.argv[1:6]
c = json.load(open(old))
c["database"]["fasta"] = fasta
c["mzml_paths"] = [mzml]
c["output_directory"] = out
json.dump(c, open(dst, "w"), indent=1)
print(f"  config copied from {old}")
print(f"  changed: database.fasta, mzml_paths, output_directory -- nothing else")
PYEOF
/usr/bin/time -v "${SAGE}" "${CFG}" --write-pin -o "${OUT}" 2> "${OUT}/sage.time"
RC=$?
N=$(awk -F'\t' 'NR>1 && $38<=0.01 && $8==1' "${OUT}/results.sage.tsv" 2>/dev/null | wc -l)
{
  printf 'key\tvalue\n'
  printf 'arm\t%s\n' "${ARM}"
  printf 'arm_description\tmetagenome with sha256 headers PLUS the full human SwissProt proteome (20,431 entries) appended, so the comparator can find host peptides that FastaLake can; nothing else changed\n'
  printf 'cohort\tCAPSCAN\n'
  printf 'sample\t%s\n' "${SAMPLE}"
  printf 'fastalake_used\tNO -- headers only, no lake, no evidence filter, no parsimony\n'
  printf 'database_fasta\t%s\n' "${FASTA}"
  printf 'database_sha256\t%s\n' "$(sha256sum "${FASTA}" | cut -d' ' -f1)"
  printf 'mzml\t%s\n' "${MZML}"
  printf 'search_config_source\tcopied verbatim from %s\n' "${OLDCFG}"
  printf 'slurm_job\t%s\n' "${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-none}}_${IDX}"
  printf 'psms_q0.01\t%s\n' "${N}"
  printf 'sage_exit\t%s\n' "${RC}"
} > "${OUT}/RUN_DESCRIPTION.tsv"
echo "psms at 1% peptide FDR: ${N}"
[ "${RC}" = 0 ] || exit "${RC}"
/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python -B ${M}/mg_full/count.py "${OUT}" "${SAMPLE}"
