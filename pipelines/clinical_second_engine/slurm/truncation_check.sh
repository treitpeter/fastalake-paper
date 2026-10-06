#!/bin/bash
# Check whether the five PREDICT mzML files that failed in MSFragger are truncated, and list other copies.
# FASTALAKE_PREDICT_MZML_DIR: the PREDICT 0_DATA directory (one subfolder per acquisition).
# FASTALAKE_RAW_SEARCH_ROOTS: space-separated roots searched for other copies.
set -u
D=${FASTALAKE_OUT_DIR}/truncation_check; mkdir -p $D
# 1) ending of every PREDICT mzML in 0_DATA
for f in ${FASTALAKE_PREDICT_MZML_DIR}/*/*.mzML; do
  if tail -c 200 "$f" | grep -q '</indexedmzML>\|</mzML>'; then s=COMPLETE; else s=TRUNCATED; fi
  printf "%s\t%s\t%s\n" "$(basename $f .mzML)" "$(stat -c %s "$f")" "$s"
done > $D/predict_mzml_endings.tsv
# 2) other copies of the five
for n in MicrobPred_13 MicrobPred_44 MicrobPred_90 MicrobPred_115 MicrobPred_181; do
  find ${FASTALAKE_RAW_SEARCH_ROOTS} -maxdepth 7 \( -name "*_SA_${n}.raw" -o -name "*_SA_${n}.mzML" -o -name "*_SA_${n}.mzML.gz" \) -printf "${n}\t%s\t%p\n"
done > $D/copies.tsv
echo DONE
