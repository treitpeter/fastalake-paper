#!/bin/bash
set -u
D=${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/truncation_check
# 1) ending of every PREDICT mzML in 0_DATA
for f in ${FASTALAKE_DATA_DIR}/*.mzML; do
  if tail -c 200 "$f" | grep -q '</indexedmzML>\|</mzML>'; then s=COMPLETE; else s=TRUNCATED; fi
  printf "%s\t%s\t%s\n" "$(basename $f .mzML)" "$(stat -c %s "$f")" "$s"
done > $D/predict_mzml_endings.tsv
# 2) other copies of the five
for n in MicrobPred_13 MicrobPred_44 MicrobPred_90 MicrobPred_115 MicrobPred_181; do
  find ${FASTALAKE_DATA_DIR}/Peter_Microbpred ${FASTALAKE_DATA_DIR}/Peter -maxdepth 7 \( -name "*_SA_${n}.raw" -o -name "*_SA_${n}.mzML" -o -name "*_SA_${n}.mzML.gz" \) -printf "${n}\t%s\t%p\n" 2>/dev/null
done > $D/copies.tsv
echo DONE
