#!/bin/bash -l
#SBATCH --job-name=s0_metag
# NO --partition DIRECTIVE: directives are parsed before the shell runs, so they
# cannot read ${FL_PARTITION}. It is passed at submit time -- see submit_cohort.sh
# and docs/SITE_CONFIG.md. A literal here would hard-code one cluster.
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=400G --time=12:00:00
#SBATCH --nice=10000
#SBATCH --output=${FASTALAKE_DATA_DIR}/%x_%A_%a.out
#SBATCH --error=${FASTALAKE_DATA_DIR}/%x_%A_%a.err
# ============================================================================
# STAGE 0 REBUILT: catalogues verbatim, metaG hashed, donor kept in the sidecar.
#
# Peter's rule, agreed 2026-08-01 and lost in the v26 rebuild:
# reference databases do not need hashing; metaG does.
#
# v26 stage 0 ran `--uniform-tag`, renaming EVERY sequence from EVERY source to
# `<COHORT>_<sha256hex>`, and v27 reuses that lake. Catalogue accessions are
# already globally unique, so hashing them destroys information for no gain.
# Only Prodigal contig ids needed it: `k141_1234_5` is unique within one assembly
# and collides across the 335 of them. And even that is now unnecessary, because
# prefix_assemblies.sbatch puts the sample into the header first, which removes
# the collision at source rather than papering over it with a hash.
#
# Everything that broke on 2026-08-02 traces to that flag: host/microbial
# classification (all 2,041,914 v26 CAPSCAN groups read "Microbial" against v25's
# Human 20,596 / Mixed 4,336), functional keying, 15.8% of members unannotatable,
# source attribution, and a 52 GB sidecar that exists to record what the
# accession would otherwise say for itself.
#
# WHY NOT JUST HASH THE ASSEMBLIES. Because all 335 share one source tag, so the
# sidecar records `CAPSCAN_MG  k141_0_1` with no way back to SAMN27256204 — the
# donor of a novel protein becomes unanswerable, which is exactly what the
# personal-barcode analysis asks. Prefixing keeps the donor in the accession, and
# as a bonus the result matches the 231,201 metaG keys the v25 eggNOG already
# holds, so 16.7% of sampled assembly accessions annotate with no re-keying.
#
# THE ONE ARGUMENT FOR --uniform-tag DOES NOT APPLY. Its help says per-source
# tags re-create a bias "because razor's default `alpha` tie-break picks the
# lexicographically smallest accession". v27 does not use alpha: chain/00_paths.sh
# line 119 sets TIEBREAK=hash-acc explicitly, ruled 2026-07-31, which hashes the
# accession before comparing. The bias is already prevented by the tie-break.
#
# WRITES TO A NEW DIRECTORY. projects/stage0_v26_20260731/ is untouched.
#
# Submit AFTER prefix_assemblies.sbatch:
#   sbatch --dependency=afterok:<prefix_jobid> --array=0-3 00_lake.sbatch
# ============================================================================
set -uo pipefail

REPO=${FASTALAKE_DATA_DIR}/FastaLake
OCEAN=${REPO}/projects/fasta_ocean_2026/scratch
UHGP=$OCEAN/uhgp-100.faa
GMGC=$OCEAN/GMGC_gut.faa
GMSC=$OCEAN/GMSC_gut.faa
SWISSPROT=${REPO}/databases/isoforms/human/uniprot/swissprot_human_reviewed_2026_03.fasta
CAPSCAN_DIR=${REPO}/canonical/publication_data/sources/assemblies_prefixed/CAPSCAN
PREDICT_ASM_DIR=${REPO}/canonical/publication_data/sources/assemblies_prefixed/PREDICT
# ONE code tree. This hardcoded ${REPO}/worktrees/v27 -- a second checkout outside the
# published repository -- which the 2026-08-03 WT fix in _config.sh did not reach because
# stage 0 does not source _config.sh.
LB=${FL_CODE_TREE:-${REPO}/canonical}/rust/fasta_extractor_v2/target/release/lake_builder

COHORTS=(CAPSCAN PREDICT HSTOOL SIHUMI)
# The array builds the four original cohorts. A cohort added later -- CAMPI_MOCK, whose
# lake is the catalogue alone -- is built ONE at a time by exporting COHORT, so it does not
# have to be wedged into an array index that other scripts count on. Without this branch an
# exported COHORT was silently ignored and `set -u` died on SLURM_ARRAY_TASK_ID instead,
# which is at least a loud failure rather than a lake built for the wrong cohort.
#   sbatch --array=0-3 00_lake.sbatch                        # the original four
#   sbatch --export=ALL,COHORT=CAMPI_MOCK 00_lake.sbatch     # one, by name
if [ -z "${COHORT:-}" ]; then
    [ -n "${SLURM_ARRAY_TASK_ID:-}" ] || { echo "FATAL: pass --array, or export COHORT=<name>"; exit 1; }
    COHORT=${COHORTS[${SLURM_ARRAY_TASK_ID}]}
fi

# --- OUTPUT ROOT ------------------------------------------------------------------
# 🔴 THIS LINE USED TO BE HARD-CODED to ${REPO}/canonical/publication_data/cohorts/
# ${COHORT}/0_lake, unconditionally. Stage 0 is the one stage that does NOT source
# _config.sh (see the LB comment above), so it never picked up FL_DATA_ROOT or
# FL_STAGE0_ROOT -- the two overrides every other stage honours.
#
# The consequence: building ANY ablation lake for a cohort that already has one --
# a metaG-free PREDICT lake, say -- would have OVERWRITTEN the published 118 GB lake
# in place, with a lake built from different sources, and said nothing. Every figure
# downstream would then be reading a lake nobody meant to build, and the only evidence
# would be an mtime. That is the destructive half of the "FL_STAGE0_ROOT alone
# overwrites cohorts/<C>" trap, sitting in the stage that writes the biggest file.
#
# Resolution order is now IDENTICAL to _config.sh, which resolves
#   S0=${FL_STAGE0_ROOT:-${DATA}/cohorts}/${COHORT}/0_lake   with DATA=${FL_DATA_ROOT:-...}
# so a lake built here lands exactly where the chain will look for it. Default is
# unchanged: with neither variable set this is the same path as before.
S0_ROOT=${FL_STAGE0_ROOT:-${FL_DATA_ROOT:-${REPO}/canonical/publication_data}/cohorts}
OUT=${S0_ROOT}/${COHORT}/0_lake
mkdir -p "$OUT" "${REPO}/canonical/publication_data/sources/stage0_build_logs"
LAKE=${OUT}/${COHORT}_lake.fasta
HEADERS=${OUT}/${COHORT}_joint_headers.tsv
STATS=${OUT}/${COHORT}_lake_stats.txt

# --- REFUSE TO CLOBBER AN EXISTING LAKE ---------------------------------------------
# The path fix above makes an accidental overwrite unlikely; this makes it impossible
# without saying so. A lake takes hours to build and is the input to everything, so
# overwriting one must be a decision, not a side effect of forgetting an export.
if [ -s "${LAKE}" ] && [ "${FL_LAKE_OVERWRITE:-0}" != "1" ]; then
    echo "FATAL: a lake already exists at ${LAKE}"
    echo "       $(ls -la "${LAKE}")"
    echo "       Building here would destroy it. Either point FL_STAGE0_ROOT / FL_DATA_ROOT"
    echo "       somewhere else, or set FL_LAKE_OVERWRITE=1 if you really mean to replace it."
    exit 1
fi

# TWO PROPERTIES ARE NEEDED AND THEY LIVE IN DIFFERENT PLACES.
#
#   1. SEQUENCE-STABLE ID. The same metaG sequence must carry the same id in
#      every sample and every study, or "is this the same protein in study X and
#      study Y" becomes a sequence comparison instead of a string match. This is
#      the whole purpose of hashing metaG, and it belongs in the ACCESSION:
#      CAPSCAN_MG_<sha256>.
#
#   2. DONOR PROVENANCE. Which subject contributed a novel protein. This belongs
#      in the SIDECAR, and it is what v26 lost: all 335 assemblies shared the tag
#      CAPSCAN_MG and Prodigal headers are bare k141_0_1, so the sidecar recorded
#      `CAPSCAN_MG  k141_0_1` and the donor was unrecoverable.
#
# An earlier version of this script tried to satisfy (2) by putting the sample in
# the accession and not hashing at all. That traded (1) away: an identical protein
# in two donors would get two different accessions. Wrong trade.
#
# Both are satisfied by doing the prefixing FIRST and then hashing. prefix_assemblies
# .sbatch rewrites the Prodigal headers to CAPSCAN_SAMN27256204_k141_0_1, so that is
# the string the sidecar preserves; --reassign-header-tags then derives the accession
# from the SEQUENCE. A protein carried by five donors ends up with ONE accession and
# FIVE sidecar entries naming all five.
# --- FL_NO_METAG=1 : THE METAGENOME-FREE ABLATION -----------------------------------
# Added 2026-08-06. PREDICT and CAPSCAN are the only cohorts whose lake carries their
# own metagenome (PREDICT_MG / CAPSCAN_MG at priority 4, 39,421,510 of PREDICT's
# 272,919,212 unique sequences -- 14.4%). Every published "FastaLake vs the matched
# metagenome" comparison therefore gives FastaLake the comparator's own assemblies PLUS
# three catalogues, so it is not a metaG-free result and cannot support a "without
# needing a metagenome" claim. This flag builds the arm that can.
#
# The metaG-free source list is EXACTLY the HSTOOL pattern -- UHGG, GMGC, GMSC,
# HUMAN_SP at priorities 1-4. HUMAN_SP is renumbered 5 -> 4 rather than leaving a gap
# at 4, so the priority sequence is contiguous and identical to the cohort that has
# always been built this way. Priority decides which source's header wins for a
# sequence seen more than once, so a gap is not obviously harmless and is not risked.
#
# METAG_TAG becomes empty, which switches off --reassign-header-tags and the stage-0
# gate's "metaG accessions are sequence-hashed" assertion -- both correct, because with
# no metaG in the lake there is nothing to hash and the assertion would fail a correct
# build. That branch is already conditional on METAG_TAG being non-empty.
_NO_MG=${FL_NO_METAG:-0}
case "$COHORT" in
  CAPSCAN) if [ "${_NO_MG}" = "1" ]; then
             SRC=(-s "UHGG:1:${UHGP}" -s "GMGC:2:${GMGC}" -s "GMSC:3:${GMSC}"
                  -s "HUMAN_SP:4:${SWISSPROT}")
             METAG_TAG=""
           else
             SRC=(-s "UHGG:1:${UHGP}" -s "GMGC:2:${GMGC}" -s "GMSC:3:${GMSC}"
                  -s "CAPSCAN_MG:4:${CAPSCAN_DIR}:*.faa"
                  -s "HUMAN_SP:5:${SWISSPROT}")
             METAG_TAG="CAPSCAN_MG"
           fi ;;
  PREDICT) if [ "${_NO_MG}" = "1" ]; then
             SRC=(-s "UHGG:1:${UHGP}" -s "GMGC:2:${GMGC}" -s "GMSC:3:${GMSC}"
                  -s "HUMAN_SP:4:${SWISSPROT}")
             METAG_TAG=""
           else
             SRC=(-s "UHGG:1:${UHGP}" -s "GMGC:2:${GMGC}" -s "GMSC:3:${GMSC}"
                  -s "PREDICT_MG:4:${PREDICT_ASM_DIR}:*.fasta"
                  -s "HUMAN_SP:5:${SWISSPROT}")
             METAG_TAG="PREDICT_MG"
           fi ;;
  HSTOOL)  SRC=(-s "UHGG:1:${UHGP}" -s "GMGC:2:${GMGC}" -s "GMSC:3:${GMSC}"
                -s "HUMAN_SP:4:${SWISSPROT}")
           METAG_TAG="" ;;
  # SIHUMIx carries NO human: human is this cohort's ENTRAPMENT database, so a
  # human sequence in the target would be target and entrapment at once.
  SIHUMI)  SRC=(-s "UHGG:1:${UHGP}" -s "GMGC:2:${GMGC}" -s "GMSC:3:${GMSC}")
           METAG_TAG="" ;;
  # CAMPI, both arms: THE BIG CATALOGUE AND NOTHING ELSE. Peter, 2026-08-04: "the whole
  # point of the CAMPI is that we try to figure out the species using our big catalog -
  # gmgc gmsc and uhgg", and "full lake - properly done - no predict or capscan metaG
  # added". So: no cohort metagenome, and no foreign metagenome either. The retired
  # campi_reference_2026 lake was 14.6% MuPr_ (MICROB-PREDICT) contigs under one uniform
  # CAMPI_<sha> tag -- another cohort's assemblies inside an external benchmark, with
  # every catalogue accession hashed away so no species could be read off a protein group.
  #
  # NO HUMAN in either arm. For the mock that is the SIHUMI rule: human is this cohort's
  # entrapment database, so a human sequence in the target would be target and entrapment
  # at once. For the faecal arm it is Peter's source list taken literally. ⚠ That differs
  # from HSTOOL, which does carry HUMAN_SP, and it means CAMPI_STOOL cannot report a host
  # share. Flagged rather than silently decided.
  #
  # Identical source list for both arms, so the two lakes are the same lake by
  # construction -- and identical to SIHUMI's 233,477,360.
  CAMPI_MOCK|CAMPI_STOOL)
           SRC=(-s "UHGG:1:${UHGP}" -s "GMGC:2:${GMGC}" -s "GMSC:3:${GMSC}")
           METAG_TAG="" ;;
  *) echo "FATAL: unknown cohort '$COHORT'"; exit 1 ;;
esac

# FL_NO_METAG on a cohort that never had a metagenome is a MISTAKE, not a no-op: the
# lake would be correct but the arm directory it lands in would claim an ablation that
# was not performed, and nothing downstream could tell the difference. Fail loudly.
if [ "${_NO_MG}" = "1" ]; then
    case "$COHORT" in
      PREDICT|CAPSCAN) ;;
      *) echo "FATAL: FL_NO_METAG=1 but ${COHORT}'s lake never contained a metagenome."
         echo "       There is nothing to ablate, and the output would misname itself."
         exit 1 ;;
    esac
fi

# Full, unfiltered UniProtKB augmentation (Peter, 2026-09-09). Both sections
# are mandatory: selecting only mock species or reviewed entries changes the question.
# Append priorities so all existing catalogue tie-breaking remains unchanged.
# FASTAs must be decompressed; lake_builder memory-maps plain FASTA inputs.
if [ -n "${FL_UNIPROT_DIR:-}" ]; then
    for section in sprot trembl; do
        [ -s "${FL_UNIPROT_DIR}/uniprot_${section}.fasta.gz.verified" ] || {
            echo "FATAL: UniProt ${section} download/checksum/decompression is not verified"; exit 1;
        }
    done
    _priority=0
    for _spec in "${SRC[@]}"; do
        [ "${_spec}" = "-s" ] && continue
        _p=${_spec#*:}; _p=${_p%%:*}
        [ "${_p}" -le "${_priority}" ] || _priority=${_p}
    done
    SRC+=(-s "UNIPROT_SP:$((_priority + 1)):${FL_UNIPROT_DIR}/uniprot_sprot.fasta"
          -s "UNIPROT_TR:$((_priority + 2)):${FL_UNIPROT_DIR}/uniprot_trembl.fasta")
    unset _priority _spec _p section
fi

echo "host=$(hostname) job=${SLURM_JOB_ID} cohort=${COHORT} start=$(date -Is)"
echo "FL_NO_METAG=${_NO_MG} (1 = metagenome-free ablation arm)"
echo "output root: ${OUT}"
echo "hashing: metaG only (${METAG_TAG:-none}); catalogues verbatim; donor preserved in the sidecar"

# --- preflight: fail before spending hours --------------------------------------
fail=0
[ -x "$LB" ] || { echo "FATAL: lake_builder missing: $LB"; exit 1; }
"$LB" --help 2>&1 | grep -q -- '--reassign-header-tags' \
  || { echo "FATAL: build lacks --reassign-header-tags"; exit 1; }
"$LB" --help 2>&1 | grep -q -- '--track-all-sources' \
  || { echo "FATAL: build lacks --track-all-sources"; exit 1; }
for spec in "${SRC[@]}"; do
  case "$spec" in -s) continue;; esac
  path=$(echo "$spec" | cut -d: -f3)
  if [ -e "$path" ]; then echo "  ok      $path"; else echo "  MISSING $path"; fail=1; fi
done
[ "$fail" -eq 0 ] || { echo "FATAL: inputs missing"; exit 1; }
echo "lake_builder sha256: $(sha256sum "$LB" | cut -d' ' -f1)"
printf '  source spec: %s\n' "${SRC[@]}"

# --- build ----------------------------------------------------------------------
REASSIGN=()
[ -n "$METAG_TAG" ] && REASSIGN=(--reassign-header-tags "$METAG_TAG")

# --- FL_DRY_RUN=1 : STOP HERE ------------------------------------------------------
# Everything above this line is the part that decides WHAT gets built and WHERE it
# lands: the cohort branch, the source list, the metaG ablation, the output root, the
# clobber guard, and the input-existence preflight. Everything below is the 1-3 hour
# build itself.
#
# So a dry run exercises the real submission path -- this script, these exports, this
# sbatch -- rather than a second copy of its logic. CLAUDE.md rule 8: a preflight that
# skips the submission path has not tested the submission. Duplicating the source-list
# logic into a separate checker would be the very defect this repository keeps hitting.
#
#   sbatch --export=ALL,COHORT=PREDICT,FL_NO_METAG=1,FL_DRY_RUN=1,... 00_lake.sbatch
if [ "${FL_DRY_RUN:-0}" = "1" ]; then
    echo "=== FL_DRY_RUN=1: stopping before lake_builder ==="
    echo "would run: ${LB} ${SRC[*]} ${REASSIGN[*]} --track-all-sources"
    echo "     -o ${LAKE} --output-headers ${HEADERS} --output-stats ${STATS} -t 64 --chunk-size 100000"
    echo "DRY RUN OK"
    exit 0
fi

echo "=== building $COHORT lake $(date '+%T') ==="
"$LB" "${SRC[@]}" \
    "${REASSIGN[@]}" \
    --track-all-sources \
    -o "$LAKE" --output-headers "$HEADERS" --output-stats "$STATS" \
    -t 64 --chunk-size 100000
rc=$?
# lake_builder can exit 141 (SIGPIPE) on otherwise-successful runs.
if [ "$rc" -ne 0 ] && [ "$rc" -ne 141 ]; then echo "FATAL: lake_builder exit $rc"; exit "$rc"; fi
echo "=== build finished $(date '+%T') (exit $rc) ==="

# --- STAGE 0 GATE, INVERTED FROM THE --uniform-tag VERSION -----------------------
# 🔴 The old gate asserted that EVERY accession is <COHORT>_<sha256> and that NO
# catalogue identifier survives. Under this design both are exactly backwards, so
# carrying it over unchanged would FAIL a correct build. Rewritten to assert the
# new invariant.
echo "=== STAGE 0 GATE: $COHORT ==="
G=0
chk() { if [ "$2" -eq 0 ]; then echo "  PASS  $1"; else echo "  FAIL  $1"; G=1; fi; }

N=$(grep -c '^>' "$LAKE"); echo "  sequences written: $N"
[ "$N" -gt 1000000 ]; chk "lake non-trivial ($N sequences)" $?

# 1. Catalogue accessions must SURVIVE VERBATIM. This is the point of the change,
#    and it is what makes references constant across cohorts and builds.
KEPT=$(grep '^>' "$LAKE" | grep -cE 'MGYG|GMGC10|GMSC10|sp\|')
echo "  catalogue accessions kept verbatim: $KEPT"
[ "$KEPT" -gt 100000 ]; chk "catalogue accessions survive ($KEPT)" $?

# 2. metaG accessions must be SEQUENCE-hashed, so the same protein carries the
#    same id in every sample and every study.
if [ -n "$METAG_TAG" ]; then
    HASHED=$(grep '^>' "$LAKE" | grep -cE "^>${METAG_TAG}_[0-9a-f]{64}$")
    echo "  metaG accessions sequence-hashed: $HASHED"
    [ "$HASHED" -gt 100000 ]; chk "metaG is ${METAG_TAG}_<sha256> ($HASHED)" $?
    LEAK=$(grep '^>' "$LAKE" | grep -cE '^>k[0-9]+_[0-9]+_[0-9]+')
    chk "no bare Prodigal contig id survives (offenders: $LEAK)" "$LEAK"

    # 2b. AND THE DONOR MUST BE RECOVERABLE. This is the check v26 never had, and
    #     its absence is why sample-level provenance was silently lost. The
    #     sidecar must carry the prefixed header, not the bare contig id.
    SAMPLED=$(head -2000000 "$HEADERS" | grep -c "${COHORT}_[^ 	|]*_k[0-9]" || true)
    echo "  sidecar rows carrying a sample-resolved metaG header: $SAMPLED (of 2M sampled)"
    [ "$SAMPLED" -gt 1000 ]; chk "donor recoverable from the sidecar ($SAMPLED)" $?
    BARE=$(head -2000000 "$HEADERS" | grep -cE "${METAG_TAG}	k[0-9]+_" || true)
    chk "sidecar holds no bare contig id (offenders: $BARE)" "$BARE"
else
    echo "  (no metaG source for this cohort)"
fi

# 3. Host must be readable from the accession, or classification breaks again.
#    Keyed on SRC, not on a cohort name. Keying on the name meant every cohort
#    except SIHUMI had to show host, so CAMPI_MOCK's correct catalogue-only lake
#    FAILED (job 5515060, 2026-08-04) and 26 downstream jobs died; and no
#    catalogue-only cohort had a zero-check at all, so host LEAKING into one
#    passed silently. Deriving the expectation from the cohort's own source list
#    fixes both, and a cohort added later cannot take the wrong branch.
#    Full UniProtKB also declares human coverage; it cannot serve as a human entrapment arm.
#    CAPSCAN/PREDICT/HSTOOL declare HUMAN_SP; SIHUMI and both CAMPI arms declare
#    none -- SIHUMI because human is its entrapment database, CAMPI because the
#    paper appends the human proteome to its GUT databases only (PMC8674281,
#    Methods). Controls: Fig2/analysis/gate_host_control.sh, all six cases.
HOST=$(grep '^>' "$LAKE" | grep -c '_HUMAN')
echo "  host accessions visible: $HOST"
if printf '%s\n' "${SRC[@]}" | grep -qE '^(HUMAN_SP|UNIPROT_SP):'; then
    [ "$HOST" -gt 1000 ]; chk "host declared in SRC and readable in the accession ($HOST)" $?
else
    [ "$HOST" -eq 0 ]; chk "$COHORT declares no host source and carries none ($HOST)" $?
fi

# 4. Provenance sidecar still written.
[ -s "$HEADERS" ]; chk "joint headers sidecar written" $?
[ -s "$HEADERS" ] && echo "  sidecar size: $(du -h "$HEADERS" | cut -f1)"

echo "=== GATE $( [ $G -eq 0 ] && echo PASS || echo FAIL ) ==="
[ "$G" -eq 0 ] || exit 1
echo "end=$(date -Is)"
