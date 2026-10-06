#!/bin/bash
# ============================================================================
# v27.1 chain — common environment. Sourced by every stage.
#
# v27.1 (2026-08-02): eggNOG annotation became a STAGE (11_annotate.sh) instead of a
# pre-existing file the chain hoped would cover its proteins. Peter: "run eggnog mapper
# for them - include this as a basic in the code as well v 27.1". Before this, four of
# five cohorts drew on funcdecoy_p05_eggnog.tsv, measured 2026-07-30 as covering 6.03%
# of selected UHGG and 0% of project metaG. CAPSCAN was the only cohort with a usable
# functional layer, solely because it had its own annotation.
#
# ONE script set for all SEVEN cohorts, parameterised by $COHORT. The canonical
# 2026-07-07 chain was four copied directories; copies drift, and a drifted copy
# is invisible. Here the cohort is the only variable.
#
#   COHORT in { PREDICT | CAPSCAN | SIHUMI | HSTOOL | SALIVA | CAMPI_MOCK | CAMPI_STOOL }
#
#   The five paper cohorts plus the two CAMPI benchmark cohorts. Every one of them now
#   reads its sample list from publication/cohorts/ and its lake from
#   publication_data/cohorts/<COHORT>/0_lake/, with no per-cohort exceptions.
#
# TWO THINGS THIS FILE EXISTS TO GET RIGHT
#
# 1. BINARIES COME FROM THE WORKTREE, NOT ${REPO}.
#    The canonical scripts point at ${REPO}/rust/*/target/release/*. ${REPO} is
#    branch v18.9-comment-resolution, where parsimony_engine still defaults to
#    --tiebreak alpha. v26 requires hash-acc (ruled 2026-07-31). Pointing at the
#    wrong binary silently runs the old tie-break and NOTHING downstream notices,
#    because both produce perfectly plausible output. So the path is the worktree
#    and the default is asserted below, not assumed.
#
# 2. THE HEADER ROW IS SKIPPED BY LINE NUMBER, NOT BY STRING.
#    The canonical scripts skip it with [ "${sample}" = "sample" ]. That works for
#    PREDICT, whose header is `sample`, and FAILS SILENTLY for CAPSCAN and SIHUMI,
#    whose header is `sample_id` -- the header row would be processed as a sample.
#    Every loop here starts at line 2.
# ============================================================================
: "${COHORT:?COHORT env required (PREDICT|CAPSCAN|SIHUMI|HSTOOL|SALIVA|CAMPI_MOCK|CAMPI_STOOL)}"

# REPO IS DERIVED FROM THIS FILE'S OWN LOCATION, NOT HARD-CODED.
#
# It used to be the absolute path below, which meant a CLONE COULD NEVER REFERENCE ITS
# OWN FILES. Measured 2026-08-04 (L7, job 5518506): cloning canonical and sourcing the
# clone's _config.sh resolved all 14 cohort inputs back to the ORIGINAL tree -- including
# the two CAMPI sample lists, which live inside the repository. A file being committed did
# not make it reachable from a clone, so "clone and run" was false even where the data
# shipped. That is invisible on this filesystem, because the hard-coded path exists and
# every check passes; it only appears for someone who is not us.
#
# This file is <repo>/publication/chain/_config.sh, so THE REPOSITORY ROOT IS TWO LEVELS
# UP. That root is `canonical/` here and is the clone directory in a clone -- which is the
# whole point: everything the repository owns is addressed relative to CANON and therefore
# follows the checkout.
#
# REPO (one further up) is kept ONLY for the inputs that genuinely live outside the
# repository -- the per-cohort sample lists under FastaLake/projects/ and saliva's lake
# under FastaLake/analyses/. In a clone those do not exist, and the chain should fail
# loudly on a missing file rather than silently reach back into someone else's tree, so
# REPO is derived from CANON and not from this file's depth.
_FL_CFG_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CANON=${FL_CANON_ROOT:-$(cd "${_FL_CFG_DIR}/../.." && pwd)}

# ---- SITE CONFIGURATION -------------------------------------------------------
# Values that differ between clusters are NOT hard-coded here. They live in
# fastalake.site.sh at the repository root, which is gitignored and created by copying
# fastalake.site.example.sh. See docs/SITE_CONFIG.md.
#
# This is what makes the repository portable: the chain says WHAT it needs -- a python,
# an mmseqs, a partition -- and the site says WHERE those are. Before this, 761 tracked
# lines named one project directory and 229 named partitions that exist on no other
# cluster, so cloning this gave you a pipeline that could not run anywhere.
#
# It is sourced BEFORE anything that uses it, and its absence is fatal rather than
# silently defaulted: a wrong path that happens to resolve produces plausible output from
# the wrong input, which is the failure this whole repository is built to avoid.
#
# 🔴 THE SITE FILE PROVIDES DEFAULTS. IT MUST NOT OVERRULE THE CALLER.
#
# Found 2026-08-06. fastalake.site.sh assigns unconditionally --
# `export FL_DATA_ROOT="..."` -- and it is sourced HERE, i.e. BEFORE line ~100 reads
# `DATA=${FL_DATA_ROOT:-...}`. So every documented submit-time override of a variable
# the site file happens to set was silently discarded:
#
#   sbatch --export=ALL,FL_DATA_ROOT=<arm root> ...   -> ignored, writes to the COHORT tree
#   sbatch --export=ALL,FL_SAGE=<0.15.0-beta.1> ...   -> ignored, runs 0.14.6
#
# The first is destructive: routing an ablation arm's output with FL_DATA_ROOT is the
# documented recipe, and instead of landing under arms/ it would have overwritten
# publication_data/cohorts/<COHORT>/ stage by stage. The second is the timsTOF recipe
# this very file documents at the FL_SAGE block below, which therefore could not work.
# Both fail SILENTLY and produce plausible output from the wrong input -- the exact
# failure mode the site-config comment above says this design exists to avoid.
#
# Fix: snapshot every FL_* variable the caller already set, source the site file for
# defaults, then restore the snapshot on top. The site file stays a plain list of
# assignments (it is gitignored and written per-cluster, so it cannot be relied on to
# use `:=`), and the precedence rule lives here, in tracked code, where it is reviewed.
_FL_SITE=${FL_SITE_CONFIG:-${CANON}/fastalake.site.sh}
_FL_CALLER_OVERRIDES=()
for _v in $(compgen -v | grep '^FL_' || true); do
    case "${_v}" in FL_SITE_CONFIG) continue ;; esac
    _FL_CALLER_OVERRIDES+=("${_v}=${!_v}")
done
unset _v
if [ -r "${_FL_SITE}" ]; then
    # shellcheck disable=SC1090
    . "${_FL_SITE}"
    for _kv in ${_FL_CALLER_OVERRIDES[@]+"${_FL_CALLER_OVERRIDES[@]}"}; do
        export "${_kv%%=*}=${_kv#*=}"
    done
    unset _kv
else
    echo "FATAL: no site configuration at ${_FL_SITE}" >&2
    echo "       cp fastalake.site.example.sh fastalake.site.sh  and fill it in." >&2
    echo "       See docs/SITE_CONFIG.md." >&2
    exit 1
fi
for _r in FL_PYTHON FL_MMSEQS FL_SAGE; do
    [ -n "${!_r:-}" ] || { echo "FATAL: ${_r} unset in ${_FL_SITE}" >&2; exit 1; }
done
unset _r
REPO=${FL_REPO_ROOT:-$(cd "${CANON}/.." && pwd)}
CHAIN_VERSION="27.1"   # bumped 2026-08-02: eggNOG is a stage, not a found file
# THE DATA ROOT. Everything the pipeline writes goes here, and it lives INSIDE the
# repository (gitignored) so there is exactly ONE directory to work in.
#
# This was `V27=${REPO}/v27`, and it survived three separate path sweeps because the
# variable hides the string: nothing matching "v27/cohorts" appears in any other file.
# Worse, line 165 does `mkdir -p "${ROOT}"`, so a check of the form `[ -d "$ROOT" ]`
# CREATES what it is checking for and then reports OK. That is how a broken root passed
# verification on 2026-08-03. A check that can manufacture its own success is not a check.
DATA=${FL_DATA_ROOT:-${CANON}/publication_data}
# THE CODE TREE. One copy, and it is canonical.
#
# This pointed at ${REPO}/worktrees/v27 -- a second checkout, outside the repository that
# gets published. So the committed chain executed binaries and Python that nobody reviewing
# canonical would ever see. Same defect as the duplicated _config.sh above.
#
# Verified 2026-08-03 before repointing:
#   SOURCE is byte-identical between the two trees for parsimony_engine, aggregation_engine
#   and fasta_export (diff -rq --exclude=target: no differences). rust-toolchain.toml pins
#   1.91.1 in both and the Cargo.lock files hash identically. The Python package differs
#   only by today's three additions (version 1.0.0, reproduce.py's _VERIF, the shipped
#   verification fixtures).
#
#   The BINARIES still differ, and one of them differs in code:
#     aggregation_engine  31/1,429,808 bytes, all .note.gnu.build-id + .rodata. .text IDENTICAL.
#     parsimony_engine    2,233/1,348,752 bytes, of which 84 are in .TEXT.
#   84 bytes of executable code from identical source and an identical pinned toolchain.
#   Neither binary embeds its build path or any src string and both are exactly the same
#   size, so inspection did not explain it and no explanation is asserted here.
#   publication/qc/verify_parsimony_binaries_agree.sbatch measures the thing that matters instead:
#   whether the two binaries produce byte-identical output over every sample of a cohort.
WT=${FL_CODE_TREE:-${CANON}}

# --- binaries: worktree build (job 5449159), NOT ${REPO} -------------------------
FASTA_EXTRACTOR=${WT}/rust/fasta_extractor_v2/target/release/fasta_extractor
FASTA_SORT=${WT}/rust/fasta_extractor_v2/target/release/fasta_sort_by_id
PARSIMONY=${WT}/rust/parsimony_engine/target/release/parsimony_engine
AGGREGATION=${WT}/rust/aggregation_engine/target/release/aggregation_engine
MMSEQS=${FL_MMSEQS}

# --- SAGE: PINNED, because it used to be whatever was on PATH ---------------------
# This was:
#     SAGE=$(command -v sage 2>/dev/null || ... || echo .../AlphaForge/sage/sage)
# so the search engine that produced a cohort's results was decided by the caller's
# PATH. That is defect (1) at the top of this file -- "pointing at the wrong binary
# silently runs the old behaviour and NOTHING downstream notices" -- applied to the
# one binary the entire depth readout comes from. It survived because all three
# fall-through targets happen to be the same 0.14.6 build (verified 2026-08-04:
# ${FASTALAKE_DATA_DIR}/sage and the AlphaForge copy are byte-identical). It
# worked by luck, and installing any other `sage` earlier on PATH would have
# silently changed every cohort's numbers with no record that anything moved.
#
# 0.14.6 IS THE PUBLISHED ENGINE. Every v27 number comes from it, so it is the
# default and it is asserted, not assumed.
SAGE=${FL_SAGE}
SAGE_EXPECTED_VERSION=${FL_SAGE_VERSION}

# --- 🔴 timsTOF / Bruker .d: 0.14.6 CANNOT QUANTIFY THEM --------------------------
# Established on ChemProt (Bruker DDA-PASEF .d), 2026-08-04, probes 1-4 under
# PhD/Peter_chemicalprot/rerun_latest_fastalake_2026-08/probe_sage_lfq/:
#
#   * SAGE gained Bruker .d reading in v0.14.4, documented in that release as
#     "ddaPASEF; NO MS1/LFQ support yet". 0.14.6 still has that limitation.
#   * Given a `quant` block, 0.14.6 on a .d writes an lfq.tsv containing ZERO ROWS
#     and exits 0. So `[ -s lfq.tsv ]` PASSES on a file with no quantification in
#     it. Any stage checking existence rather than rows will believe it worked.
#     -> 05_search.sh now counts rows. Do not weaken that check back to -s.
#   * LFQ on ion-mobility data arrives in v0.15.0-beta.1 (2025-04-07). On the same
#     .d that beta.1 reads fine, it reports "read 3284 ms1 spectra" and writes real
#     LFQ: 149-352 rows on sparse metaproteome samples, ~23,000 on HeLa, 6/6 samples.
#   * v0.15.0-beta.2 (2026-03-31) is NEWER and BROKEN for this input. It aborts in
#     under a second, before reading any spectrum:
#       panicked at timsrust-0.4.2/.../spectrum_reader/builder.rs:22:48
#       called `Result::unwrap()` on an `Err` value: IO(Os { code: 2, NotFound })
#     An strace of the whole process shows no failed open anywhere under the .d --
#     only /etc/ld.so.preload and /sys${FASTALAKE_DATA_DIR}/cpu.max, both benign -- so it is
#     not a missing file in the data. DO NOT "upgrade" a working beta.1 to beta.2.
#
# ION MOBILITY WAS NEVER THE OBSTACLE. SAGE fits a mobility model at rsq 0.93-0.94
# on this data, from the .d and from mzML alike. The obstacle is MS1: LFQ needs it,
# and a .d has it while a DDA-converted mzML frequently does not. ChemProt's two
# mzML sets are both MS2-ONLY -- what looked like 17,110 MS1 scans in them are
# <userParam name="ms level" value="1"/> tags inside <isolationWindow>, i.e.
# precursor metadata. Count MS1 SPECTRA before blaming the engine.
#
# For a timsTOF cohort, override explicitly and let the provenance record it:
#   sbatch --export=ALL,FL_SAGE=/path/to/sage-v0.15.0-beta.1/sage,FL_SAGE_VERSION=0.15.0-beta.1 05_search.sh
PY=${FL_PYTHON}
# fasta_lake is an editable/cwd-dependent install: `import fasta_lake` succeeds only
# when cwd happens to be ${REPO}. A SLURM job does not have that cwd, so stage 10
# would fail at the END of a 20-hour chain. Pin it explicitly.
export PYTHONPATH=${WT}${PYTHONPATH:+:${PYTHONPATH}}

# --- stage 0 output ---------------------------------------------------------------
# THE DEFAULT NOW POINTS AT THE LAKE THAT PRODUCED THE PUBLISHED v27 RESULTS.
#
# It used to default to projects/stage0_v26_20260731 (76 GB, built 13:27 on 2026-08-02)
# while every published v27 number came from the corrected lake (67 GB, 20:40) reached
# only by passing FL_STAGE0_ROOT explicitly. So `main` as committed did NOT reproduce the
# published data -- it would have built from the v26 lake, silently, and produced numbers
# nobody could match. That is the single worst kind of default.
#
# The old comment said to flip the default once no job was in flight against the old
# lake. Checked 2026-08-03: nothing was, so it is flipped.
#
# The lake also MOVED. It lived in v27/stage0_metag_only/<COHORT>/0_lake -- a name that was
# never true: that build pools UHGP-100, GMGC_gut, GMSC_gut, human SwissProt AND the
# prefixed metagenome assemblies. "metag" described what the rebuild FIXED (v26 gave all
# 335 assemblies the tag CAPSCAN_MG with bare Prodigal headers, losing donor provenance),
# not what the lake contains. A whole-file census measured metagenome at 14.44%.
# It now sits with the rest of its cohort's stages, so a cohort tree runs 0 -> 11 in
# one place instead of stage 0 living under a different parent with a false name.
#
# Override for an ablation arm:  sbatch --export=ALL,FL_STAGE0_ROOT=<root> ...
S0=${FL_STAGE0_ROOT:-${DATA}/cohorts}/${COHORT}/0_lake
LAKE=${S0}/${COHORT}_lake.fasta
LAKE_HEADERS=${S0}/${COHORT}_joint_headers.tsv
LAKE_STATS=${S0}/${COHORT}_lake_stats.txt

# --- per-cohort sample list -------------------------------------------------------
# Verified 2026-07-31: every col2 CSV and col3 mzML resolves, counts match
# DATA_MANIFEST.md. HSTOOL uses masters_circularity_F1 rather than
# healthy_eq1_an2026 because the latter has NO header row (its first line is data),
# which would break the +2 line arithmetic every array stage uses.
case "${COHORT}" in
  PREDICT) SAMPLE_LIST=${CANON}/publication/cohorts/predict_sample_list.tsv ;;
  CAPSCAN) SAMPLE_LIST=${CANON}/publication/cohorts/capscan_sample_list.tsv ;;
  # A5+B5 ONLY. E5 and F5 are failed injections: 769 / 799 precursors against A5's 17,986,
  # and SAGE logged "0 target proteins at 1% FDR" for E5, which survived earlier only because
  # 5 decoy PSMs won. OPEN_ITEMS B1 ruled 2026-07-28: report n = 2 and disclose the failures.
  # This must be fixed in the COHORT DEFINITION, not by selecting array indices -- the FDP
  # tool counts against this list and correctly refuses a subset ("A median FDP over a
  # subset is not the cohort's FDP").
  SIHUMI)  SAMPLE_LIST=${CANON}/publication/cohorts/sihumi_sample_list.tsv ;;
  HSTOOL)  SAMPLE_LIST=${CANON}/publication/cohorts/hstool_sample_list.tsv ;;
  # Added 2026-08-02. Saliva is cohort #2 in Peter's ordering and had NO v26 run --
  # it predates the chain's original four-cohort scope, not a defect in it.
  #
  # 🔴 ITS LAKE IS NOT THE SHARED CATALOGUE SET, AND THAT IS DELIBERATE.
  # The other four cohorts draw on gut catalogues (UHGG/GMGC/GMSC + human) built by
  # stage0_v26. Saliva is ORAL, so it uses its own environment-appropriate lake
  # (780,319 proteins, PXD003028 stage0_lake). That makes saliva a genuine instance
  # of generality axis 6 -- "reference collections from different microbial
  # environments" -- and it ALSO means saliva's absolute depth is NOT comparable to
  # the gut cohorts on one axis. Compare saliva to saliva; compare the METHOD across
  # cohorts, never the raw protein-group count.
  #
  # What this run DOES give on the v26 basis: hash-acc tiebreak (stage 6), cov-mode 0
  # clustering (stage 4), and the worktree binaries. What it does NOT give is a
  # source-blind stage-0 rebuild -- that is a separate, expensive decision, flagged
  # rather than silently assumed.
  SALIVA)
    SAMPLE_LIST=${CANON}/publication/cohorts/saliva_sample_list.tsv
    SALIVA_LAKE_OVERRIDE=${REPO}/projects/pride_e2e_test/PXD003028/stage0_lake/lake.fasta
    ;;
  # --- CAMPI (Van Den Bossche 2021, PXD023217) -------------------------------
  # 🔴 TWO COHORTS, ON PURPOSE, AND THEY MUST NOT BE MERGED.
  #
  # The S series is the SIHUMIx MOCK COMMUNITY (15 acquisitions: S01-S10, S11 in four
  # fractions, S12). The F series is HUMAN FAECAL (F05, F06). Stage 1 pools de novo
  # peptides across a cohort, so running them together would let faecal peptides select
  # sequences for the mock samples -- and the mock's whole purpose is a BLIND recovery of
  # a community whose composition is known. That is not a preference, it is the
  # difference between a test and a leak.
  #
  # 🔴 CORRECTED 2026-08-05. This block asserted the exact opposite -- "S is faecal, F is
  # the mock" -- and on that basis retracted a correct result. Three independent routes
  # settled it:
  #   1. CAMPI's own Figure 3 caption, PMC8674281: "the SIHUMIx sample S11 (A) ... the
  #      fecal sample F07 (B)".
  #   2. The data, no labels consulted -- SLURM 5523909, Fig2/analysis/which_series_is
  #      _the_mock.py: each acquisition's own de novo peptides tested for containment in
  #      CAMPI's eight-species proteome. Mean S 0.7562, mean F 0.1347; controls 1.0000
  #      for peptides cut from the proteome and 0.0000 reversed.
  #   3. The file counts -- CAMPI acquired 21 of each; we hold 15 S and 2 F.
  # So campi_reference_2026's S07/S09/S12 run was on the MOCK all along, and the
  # "species recall 8/8" retraction issued 2026-08-04 is WITHDRAWN. That number was
  # measured on the right samples. It has its own separate problem -- it rests on a
  # single peptide for L. plantarum -- which is why the recall panel prints per-member
  # evidence beside the tally. See Fig2/CAMPI_REPLICATION.md.
  #
  # THE LAKE IS THE BIG CATALOGUE AND NOTHING ELSE -- UHGG + GMGC + GMSC, accessions
  # verbatim. Peter, 2026-08-04: "the whole point of the CAMPI is that we try to figure out
  # the species using our big catalog". No cohort metagenome, and specifically no PREDICT
  # or CapScan assemblies: the retired campi_reference_2026 lake carried 14.6% MuPr_
  # (MICROB-PREDICT) contigs under a uniform CAMPI_<sha> tag, which is another cohort's
  # metagenome inside an external benchmark.
  CAMPI_MOCK)
    SAMPLE_LIST=${CANON}/publication/cohorts/campi_mock_sample_list.tsv ;;
  CAMPI_STOOL)
    SAMPLE_LIST=${CANON}/publication/cohorts/campi_stool_sample_list.tsv ;;
  *) echo "[FAIL] unknown COHORT '${COHORT}'" >&2; exit 1 ;;
esac

# --- FL_SAMPLE_LIST_OVERRIDE ------------------------------------------------------
# An ARM's cohort is not always the design cohort. The metaG baseline can only ever cover
# 296 of PREDICT's 306 acquisitions, because 10 have no assembly FASTA at all -- a data
# property, not a pipeline failure.
#
# 06_aggregate.sh refuses to run on an incomplete cohort, and that gate is correct: "a
# median FDP over a subset is not the cohort's FDP". The right way to aggregate an arm
# with a genuinely smaller n is therefore to hand it the list it actually covers, NOT to
# relax the gate. Relaxing it would let a silently-truncated run pass everywhere.
#
# The override is explicit and appears in the arm's RUN_DESCRIPTION, so a smaller n is
# always a stated decision rather than an accident.
if [ -n "${FL_SAMPLE_LIST_OVERRIDE:-}" ]; then
    [ -s "${FL_SAMPLE_LIST_OVERRIDE}" ] \
        || { echo "FATAL: FL_SAMPLE_LIST_OVERRIDE missing/empty: ${FL_SAMPLE_LIST_OVERRIDE}" >&2; exit 1; }
    SAMPLE_LIST=${FL_SAMPLE_LIST_OVERRIDE}
fi

# --- SALIVA lake override ---------------------------------------------------------
# Saliva is oral, so it does not draw on the gut catalogues that stage0_v26 built for
# the other four cohorts. Its lake is the v26 REBUILD of the four-arm study's v3 lake
# (human SwissProt + MGnify human-oral v1.0.1 protein_catalogue-100, 1,643,197 seqs),
# produced by analyses/saliva_v26_20260802/prep_lakes_v26.sbatch with the WORKTREE
# lake_builder. The 2026-08-01 original used the repo binary (f2855e3e, Jul 21)
# rather than the worktree one (59424a0e, Jul 31), which is exactly the trap this
# file's header documents -- that study never sourced this file, so it never got the
# guard. Hence the rebuild.
#
# Pointing SALIVA here rather than at stage0_v26 is deliberate and must be disclosed:
# saliva's ABSOLUTE depth is not comparable to the gut cohorts because the source
# catalogue differs by design. What IS comparable across all five is the METHOD and
# every measurement the bench harness and the entrapment stages make on it.
# SALIVA'S LAKE NOW LIVES WHERE EVERY OTHER COHORT'S DOES (2026-08-04).
#
# It was the one cohort reaching outside publication_data, to
# `analyses/saliva_v26_20260802/db/v3_*`, under filenames that named neither the cohort
# nor the stage. That made SALIVA the single exception to "one directory to work in" and
# it is why L7 reported saliva's lake as unreachable from a clone for a different reason
# than the other six.
#
# The three files were copied to the standard location and verified byte-identical
# (sha256 55eb3cae, f87825b9, faa9afd5), and the copy satisfies the chain's own gate:
# 1,643,197 FASTA records against 1,643,197 sidecar rows. The originals are left in
# place as history and are no longer read.
#
# What has NOT changed, and still must be disclosed: saliva's lake is human SwissProt +
# MGnify human-oral protein_catalogue-100, not the gut catalogues the other cohorts draw
# on. Its ABSOLUTE depth is therefore not comparable to theirs. What is comparable is the
# METHOD, and every measurement the bench harness and the entrapment stages make on it.
# Moving the file changes where it is read from, not what it contains.

ROOT=${DATA}/cohorts/${COHORT}
LOGS=${ROOT}/logs
mkdir -p "${ROOT}" "${LOGS}"

# FL_STRATEGY_OVERRIDE lets an arm change the inference RULE while holding everything
# else fixed. razor (default, MaxQuant analogue) | species_budget (groups by MGYG species
# prefix; only 75% of this evidence lake carries one) | uniform (keep all with >=1 peptide)
# | uniform_2pep (>=2 unique peptides -- the classic heuristic, and the rule Johannes
# argues removes the metagenome's gene-fragmentation artefact) | info_score | bayesian.
#
# The strategy is part of the output path (4_parsimony/${STRATEGY}_${TIEBREAK}) and of
# DATASET, so two strategies cannot overwrite each other and a file says which produced it.
# 🔴 THE DOCS AND THE BINARY SPELL THESE DIFFERENTLY. docs/INFERENCE_STRATEGIES.md and the
# --strategy help line both say `species_budget` / `uniform_2pep`, but clap's accepted
# values are `species-budget` / `uniform2pep`. Passing the documented spelling makes
# parsimony_engine exit 2 with "invalid value" -- which is how two sizing tests died on
# 2026-08-06. Accept either spelling from the caller and normalise to what the binary
# takes; keep the underscore form for paths and labels so output names stay readable.
STRATEGY=${FL_STRATEGY_OVERRIDE:-razor}
case "${STRATEGY}" in
  razor)                        STRATEGY_ARG=razor ;;
  species_budget|species-budget) STRATEGY=species_budget; STRATEGY_ARG=species-budget ;;
  uniform)                      STRATEGY_ARG=uniform ;;
  uniform_2pep|uniform2pep)     STRATEGY=uniform_2pep;  STRATEGY_ARG=uniform2pep ;;
  info_score|info-score)        STRATEGY=info_score;    STRATEGY_ARG=info-score ;;
  bayesian)                     STRATEGY_ARG=bayesian ;;
  *) echo "FATAL: FL_STRATEGY_OVERRIDE='${STRATEGY}' is not a parsimony_engine strategy" >&2; exit 1 ;;
esac
TIEBREAK=hash-acc          # ruled 2026-07-31; passed EXPLICITLY, never inherited
# --- CLUSTERING IDENTITY: now a variable, because 97% has never been justified --------
# Searched the whole tree 2026-08-06: there is NO identity sweep anywhere.
# publication/bench/SPEC.md:76 lists `dbsweep` as "not built", and the only evidence for
# any threshold is a single two-point comparison on SIHUMIx (STAGES.md:52, unclustered
# gives +4.8% peptides at no specificity cost). figures/cohort_art/cohorts/
# COHORT_FUNNELS.md:46 already warns that on SIHUMIx 97% collapses 1.70M -> 464k and
# "may be collapsing strain distinctions the cohort exists to resolve".
#
# So 97% is asserted, not measured. It is also the wrong instrument for the redundancy it
# is aimed at: Johannes_FL_problem_statement.docx section B0 measures 66.7% of
# FastaLake's ambiguous peptides as spanning MORE THAN ONE source catalogue -- the same
# gene entering the lake twice from two catalogues, landing below 97% identity, and
# surviving clustering as two entries. Clustering harder is the direct test of whether
# that ambiguity is bookkeeping rather than biology: peptides should stay flat while
# ambiguity drops.
#
# The output directory carries the threshold, so two sweeps cannot overwrite each other
# and a directory says which threshold produced it. CLUSTER_TAG is derived, never passed,
# so it cannot disagree with CLUSTER_ID.
CLUSTER_ID=${FL_CLUSTER_ID:-0.97}
CLUSTER_TAG=$(awk -v v="${CLUSTER_ID}" 'BEGIN{printf "%d", v*100 + 0.5}')
case "${CLUSTER_TAG}" in
  ''|*[!0-9]*) echo "FATAL: FL_CLUSTER_ID='${CLUSTER_ID}' is not a fraction like 0.97" >&2; exit 1 ;;
esac
[ "${CLUSTER_TAG}" -ge 1 ] && [ "${CLUSTER_TAG}" -le 100 ] \
  || { echo "FATAL: FL_CLUSTER_ID must be between 0.01 and 1.00, got '${CLUSTER_ID}'" >&2; exit 1; }
# Default 0.97 -> CLUSTER_TAG 97 -> 2_cluster97/cluster97_* : byte-identical to before.
CLUSTER_DIR=${ROOT}/2_cluster${CLUSTER_TAG}
CLUSTER_PREFIX=cluster${CLUSTER_TAG}
# 🔴 CLUSTER_REPS IS DERIVED AND NOT OVERRIDABLE, and it used to be both.
#
# `CLUSTER_REPS=${FL_CLUSTER_REPS:-...}` was added so stage 08 could point its
# conserved-peptide scan at a different target space. But THREE stages read this variable,
# not one:
#
#   02_cluster.sh:76   REP=${CLUSTER_REPS}      the file it checks and counts AFTER mmseqs
#   03_persample.sh:71 CLUSTER_LAKE=...         the lake every per-sample database is drawn from
#   08_entrapment_fdp.sh                        the intended consumer
#
# `--export=ALL` is how the chain submits, so exporting FL_CLUSTER_REPS once for stage 08
# left it set for the whole run. mmseqs would then write ${CLUSTER_DIR}/..._rep_seq.fasta
# while stage 02's `[ -s "${REP}" ]`, its count_headers log line and ${OUT}/.n_reps all
# described the OTHER file -- and the run passed. Stage 03 then built every per-sample
# database from foreign representatives: a different threshold, a different cohort, or a
# clustered file inside a run whose DATASET was never suffixed _noclust. Every stage-3 and
# stage-4 number, and the funnel rung the figures print, would describe a database this
# cohort never clustered, with no stage detecting it.
#
# The variable below is derived from CLUSTER_DIR and CLUSTER_PREFIX, so the comment above
# -- "CLUSTER_TAG is derived, never passed, so it cannot disagree with CLUSTER_ID" -- is
# true again. The override now has its own name and reaches only the stage that wanted it.
CLUSTER_REPS=${CLUSTER_DIR}/${CLUSTER_PREFIX}_rep_seq.fasta

# Stage 08's conserved-peptide target space, and ONLY stage 08's. Defaults to the
# representatives this build actually produced. When FL_CLUSTER_REPS names something else,
# stage 08 logs the substitution rather than absorbing it silently.
ENTRAP_REPS=${FL_CLUSTER_REPS:-${CLUSTER_REPS}}
CLUSTER_COV=0.8
CLUSTER_COVMODE=0          # what every published build actually used
MMSEQS_SHUFFLE=0           # respect input order; the sort in stage 03 is what makes it safe

# --- ABLATION SWITCHES: the two OPTIONAL stages, actually made optional ------------
# Added 2026-08-06 for the metagenome comparison arms.
#
# STAGES.md has called stages 2 (cluster) and 4 (parsimony) OPTIONAL since 2026-08-03,
# but the code could not skip either: 03_persample.sh hard-coded
# ${ROOT}/2_cluster97/cluster97_rep_seq.fasta and failed without it, and 05_search.sh
# hard-coded the parsimony FASTA. Documented-optional and actually-optional had drifted
# apart, which is worse than either, because the documentation is what gets believed.
#
# WHY THEY ARE NEEDED. Peter, 2026-08-06: "in order to be able to compare the metaG /
# fastalake we need to apply the exact same protein grouping to the two pipelines ...
# if we do not do the clustering step for the metaG - to make it comparable - we cannot
# cluster either". The metagenome arm searches raw Prodigal gene calls with no
# clustering, so any FastaLake-vs-metaG depth difference measured against a CLUSTERED
# FastaLake confounds the database source with the grouping rule. These switches remove
# that confound by matching the grouping.
#
#   FL_SKIP_CLUSTER=1     stage 3 reads 1_evidence/evidence_lake.fasta directly.
#                         Stage 1b (sort) is then unnecessary -- it exists only because
#                         MMseqs2 output depends on input order.
#   FL_SKIP_PARSIMONY=1   stage 5 searches the per-sample lake, NOT the parsimonious set.
#                         ⚠️ That is ~1.87M sequences per PREDICT sample against 8,147 --
#                         a 229x larger search space. Size the search job accordingly;
#                         05_search.sh's stock 16G/4h will not hold.
#
# Default 0 on both, so the published chain is bit-for-bit what it was.
SKIP_CLUSTER=${FL_SKIP_CLUSTER:-0}
SKIP_PARSIMONY=${FL_SKIP_PARSIMONY:-0}
case "${SKIP_CLUSTER}" in 0|1) ;; *) echo "FATAL: FL_SKIP_CLUSTER must be 0 or 1, got '${SKIP_CLUSTER}'" >&2; exit 1 ;; esac
case "${SKIP_PARSIMONY}" in 0|1) ;; *) echo "FATAL: FL_SKIP_PARSIMONY must be 0 or 1, got '${SKIP_PARSIMONY}'" >&2; exit 1 ;; esac

# The dataset label carries the ablation, so an output FILE says which arm produced it
# even after it has been copied somewhere else. The DIRECTORY layout is deliberately
# left alone -- arms are separated by FL_DATA_ROOT, and renaming
# ${STRATEGY}_${TIEBREAK} would ripple through stages 4, 5, 6 and every reader of them.
#
# ⚠️ if/fi, NOT `[ test ] && VAR=x`. That idiom evaluates to exit status 1 when the test
# is false, and it would be the last command of a top-level AND list, so `set -e` acts on
# it -- and submit_cohort.sh sources this file under `set -euo pipefail`. The driver would
# have died on the DEFAULT path, with no clustering ablation requested at all.
# --- MINIMUM DE NOVO PEPTIDE LENGTH -------------------------------------------------
# Added 2026-09-03. The length was written as a literal `--min-length 9` in THREE stages
# (01_evidence, 03_persample, 04_parsimony). The 2026-09-02 HSTOOL min-length sweep changed
# it in stage 1 only, so stages 3 and 4 of every "min-length L" arm would still have run
# at 9 -- a sweep that did not sweep. One variable, read by all three, so the three cannot
# disagree. Default 9 keeps the published chain bit-for-bit.
#
#   FL_MIN_LENGTH=10      Peter, 2026-09-03: in the I/L-merged reservoir the expected
#                         random count is ~2.74 at length 9 and ~0.195 at length 10, so
#                         the "below one expected occurrence" boundary is 10, not 9.
#                         Arm: publication_data/arms/10peptide_length_runs/
MIN_LENGTH=${FL_MIN_LENGTH:-9}
case "${MIN_LENGTH}" in
  ''|*[!0-9]*) echo "FATAL: FL_MIN_LENGTH must be an integer, got '${MIN_LENGTH}'" >&2; exit 1 ;;
esac
if [ "${MIN_LENGTH}" -lt 5 ] || [ "${MIN_LENGTH}" -gt 50 ]; then
    echo "FATAL: FL_MIN_LENGTH out of range [5,50]: ${MIN_LENGTH}" >&2; exit 1
fi

# --- DE NOVO SCORE CUT (top fraction per sample) ------------------------------------
# Added 2026-09-06. `--top-percent 0.30` was a literal in TWO stages (01_evidence,
# 03_persample); stage 04's parsimony_engine takes no cut (see its header). One variable,
# read by both, so the two cannot disagree. Default 0.30 keeps the published chain
# bit-for-bit. Written with two decimals so the dataset suffix is fixed-width (_tp005 ..
# _tp100) and sorts.
#
#   FL_TOP_PERCENT=1.00   Peter, 2026-09-06: "now we take top 30% - and i want to increase
#                         this ... make nice curves". Score sweep 0.05..1.00, HSTOOL ruled 58.
#                         Arm: publication_data/arms/score_sweep_runs/top<PCT>/
TOP_PERCENT=${FL_TOP_PERCENT:-0.30}
case "${TOP_PERCENT}" in
  0.[0-9][0-9]|1.00) ;;
  *) echo "FATAL: FL_TOP_PERCENT must be 0.01..1.00 with two decimals, got '${TOP_PERCENT}'" >&2; exit 1 ;;
esac
if [ "${TOP_PERCENT}" = "0.00" ]; then
    echo "FATAL: FL_TOP_PERCENT=0.00 keeps no peptides" >&2; exit 1
fi

DATASET=${COHORT}_V27_${STRATEGY}_${TIEBREAK}
if [ "${SKIP_CLUSTER}"   = "1" ]; then DATASET=${DATASET}_noclust; fi
if [ "${SKIP_PARSIMONY}" = "1" ]; then DATASET=${DATASET}_nopars; fi
if [ "${MIN_LENGTH}" != "9" ]; then DATASET=${DATASET}_ml${MIN_LENGTH}; fi
if [ "${TOP_PERCENT}" != "0.30" ]; then DATASET=${DATASET}_tp${TOP_PERCENT/./}; fi

# --- entrapment class sources (identical files to every prior entrapment run) -----
PFUR_REF=${FL_ENTRAP_PFUR}
PLANT_REF=${FL_ENTRAP_PLANT}
SHUFFLE_SEED=42

# --- Rust thread stack: set ONCE, here, so every stage inherits it ----------------
# SAGE dies with "thread '<unknown>' has overflowed its stack / fatal runtime error:
# stack overflow" (exit 134, SIGABRT) on large unfiltered databases. It is NOT an
# out-of-memory: the saliva v2 arm crashed at 93 GB RSS against a 500 G reservation,
# on 2026-08-01 (job 5459693_2) and again on 2026-08-02 (5478580_2) because the fix
# lived only in a one-off retry script instead of the code path everyone runs.
#
# Two settings, both required and neither sufficient alone:
#   ulimit -s        raises the MAIN thread's stack (the kernel-grown one)
#   RUST_MIN_STACK   raises every thread Rust SPAWNS -- ulimit does not reach those
#
# Harmless when unnecessary: it raises a ceiling, it does not reserve memory.
ulimit -s unlimited 2>/dev/null || ulimit -s 1048576 2>/dev/null || true
export RUST_MIN_STACK=${RUST_MIN_STACK:-1073741824}   # 1 GiB

log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${COHORT}] $*"; }
fail() { echo "[FAIL] [${COHORT}] $*" >&2; exit 1; }

n_samples() { echo $(( $(wc -l < "${SAMPLE_LIST}") - 1 )); }

# Row for array index $1 (0-based). Header is line 1, so data starts at line 2.
sample_row() { sed -n "$(( $1 + 2 ))p" "${SAMPLE_LIST}"; }

# --- CODE PROVENANCE GUARD --------------------------------------------------------
# Peter, 2026-08-03: "we have had A TON OF ISSUES V25 V26 due to non-committed branches
# and stupid shit like this. should NEVER HAPPEN AGAIN."
#
# What went wrong, concretely, on 2026-08-03: aggregation_engine was rebuilt and RUN before
# its source was committed. The binary that had aggregated all five cohorts (sha256
# 255f77f9) was overwritten by 2f6c8808, and because Rust builds here are not
# byte-reproducible -- rebuilding the same source yields a third hash, d9885fa1 -- the
# provenance of the shipped data could not be re-established by hashing afterwards. The
# evidence had to be reconstructed by re-running and comparing output.
#
# This makes that impossible to repeat. Every stage sources this file, so this is the one
# choke point that sees every run. For each binary it is given:
#   * the crate source must be COMMITTED (no dirty tracked files, no untracked .rs)
#   * the binary must be NEWER than the newest source file, or it is stale
#   * commit, binary sha256, and source state are appended to a manifest in the output dir
#
# FL_ALLOW_DIRTY=1 bypasses it and says so in the loudest terms available. A guard with no
# escape hatch gets commented out the first time it is inconvenient, and then it is worse
# than no guard because everyone believes it is still on.
assert_code_committed() {
    local out_dir="$1"; shift
    local manifest="${out_dir}/RUN_CODE_PROVENANCE.tsv"
    local bad=0

    mkdir -p "${out_dir}"
    [ -s "${manifest}" ] || printf 'utc\tstage\tbinary\tbinary_sha256\tcrate\tbranch\tcommit\tsource_state\n' > "${manifest}"

    for bin in "$@"; do
        [ -x "${bin}" ] || { echo "[FAIL] binary missing: ${bin}" >&2; bad=1; continue; }
        # crate root = two levels up from target/release/<bin>
        local crate; crate=$(cd "$(dirname "${bin}")/../.." 2>/dev/null && pwd)
        local src="${crate}/src"
        local state="clean"

        # THE COMMIT MUST COME FROM THE BINARY'S OWN WORKTREE, NOT FROM ${REPO}.
        # First version asked `git -C ${REPO}` and recorded 1ef0d37 -- the main checkout's
        # branch, v18.9-comment-resolution -- while the binary it was describing was built
        # in worktrees/v27. That is a commit which does not contain the code that ran, and
        # a wrong commit stamped on data is worse than no commit, because it will be
        # believed. Resolve from the crate directory so the answer is about this binary.
        local repo_git branch commit
        repo_git=$(git -C "${crate}" rev-parse --show-toplevel 2>/dev/null || echo "")
        if [ -z "${repo_git}" ]; then
            branch="NO-GIT"; commit="NO-GIT"
            state="NOT-IN-GIT"
            echo "[FAIL] ${crate} is not inside a git worktree; its provenance is unrecordable" >&2
            bad=1
        else
            branch=$(git -C "${crate}" rev-parse --abbrev-ref HEAD 2>/dev/null || echo DETACHED)
            commit=$(git -C "${crate}" rev-parse --short HEAD 2>/dev/null || echo NO-GIT)
        fi

        if [ -d "${src}" ] && [ -n "${repo_git}" ]; then
            # dirty tracked files, or untracked sources, in this crate
            local st; st=$(git -C "${crate}" status --porcelain -- "${crate}" 2>/dev/null | head -20)
            if [ -n "${st}" ]; then
                state="DIRTY"
                echo "[FAIL] ${crate#${repo_git}/} has UNCOMMITTED changes:" >&2
                echo "${st}" | sed 's/^/         /' >&2
                bad=1
            fi
            # stale binary: any source newer than the binary means it was not rebuilt
            local newer; newer=$(find "${src}" -name '*.rs' -newer "${bin}" 2>/dev/null | head -3)
            if [ -n "${newer}" ]; then
                state="${state},STALE-BINARY"
                echo "[FAIL] ${bin##*/} is OLDER than its source -- rebuild it:" >&2
                echo "${newer}" | sed "s|${repo_git}/|         |" >&2
                bad=1
            fi
        fi

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "${STAGE_NAME:-unknown}" \
            "${bin##*/}" "$(sha256sum "${bin}" | cut -d' ' -f1)" \
            "${crate#${repo_git:-/}/}" "${branch}" "${commit}" "${state}" >> "${manifest}"
    done

    if [ "${bad}" -ne 0 ]; then
        if [ "${FL_ALLOW_DIRTY:-0}" = "1" ]; then
            echo "[⚠ FL_ALLOW_DIRTY=1] RUNNING ANYWAY WITH UNCOMMITTED OR STALE CODE." >&2
            echo "[⚠] Whatever this run produces CANNOT be traced to a commit." >&2
            echo "[⚠] Do not ship, cite, or put in a figure anything it writes." >&2
        else
            echo "[FAIL] refusing to run: the code that would produce this data is not committed." >&2
            echo "       Commit it, or set FL_ALLOW_DIRTY=1 to override and forfeit provenance." >&2
            exit 1
        fi
    else
        echo "[ok] code provenance: commit ${commit}, all binaries built from committed source"
    fi
    echo "[ok] recorded -> ${manifest}"
}

count_headers()        { grep -c '^>' "$1" || true; }
count_entrap_headers() { grep -c '^>ENTRAP' "$1" || true; }
count_target_headers() { echo $(( $(count_headers "$1") - $(count_entrap_headers "$1") )); }

# --- the assertion that makes point 1 above real ---------------------------------
# Asserts BEHAVIOUR, not help text. --tiebreak is required with no default, so the
# binary must REFUSE to run when it is omitted. That is a stronger check than grepping
# for a default string: it proves no invocation can silently inherit a tie-break, which
# is the actual failure mode (the canonical script omitted the flag and got a binary
# that hardcodes alpha).
assert_tiebreak_explicit() {
    [ -x "${PARSIMONY}" ] || fail "parsimony_engine missing: ${PARSIMONY}"
    # 1. hash-acc must be an accepted value.
    "${PARSIMONY}" --help 2>&1 | grep -q 'hash-acc' \
        || fail "parsimony_engine at ${PARSIMONY} does not know hash-acc -- wrong build"
    # 2. omitting --tiebreak must fail. If this SUCCEEDS, a default exists and some
    #    script somewhere is inheriting it silently.
    if "${PARSIMONY}" --strategy razor --fasta /dev/null --peptides /dev/null \
                      --output-dir /dev/null --sample _probe >/dev/null 2>&1; then
        fail "parsimony_engine ACCEPTED a run with no --tiebreak: a silent default exists"
    fi
    # 3. and it must fail for that reason specifically, not for some unrelated one.
    #    Capture, do NOT pipe: every caller runs `set -o pipefail`, under which
    #    `binary ... | grep -q X` yields the binary's nonzero exit rather than grep's
    #    verdict, so the check would be false even when X matched.
    local out
    out=$("${PARSIMONY}" --strategy razor --fasta /dev/null --peptides /dev/null \
                         --output-dir /dev/null --sample _probe 2>&1 || true)
    grep -qiE 'tiebreak' <<<"${out}" \
        || fail "parsimony_engine rejected the probe run but not because --tiebreak was missing"
}
# Kept as an alias so any stage still calling the old name gets the new, stronger check
# rather than silently skipping it.
assert_tiebreak_default() { assert_tiebreak_explicit; }

# --- the same treatment for SAGE --------------------------------------------------
# The version is asserted rather than trusted, because SAGE's version decides whether
# a whole class of input can be quantified at all (see the timsTOF block above), and
# because until 2026-08-04 the binary was resolved from PATH. A run that silently used
# a different engine would produce perfectly plausible output.
#
# It reads the version from the BINARY, not from a variable someone set, so a wrong
# FL_SAGE cannot announce itself as the right one.
assert_sage_version() {
    [ -x "${SAGE}" ] || fail "SAGE missing: ${SAGE}"
    local got
    got=$("${SAGE}" --version 2>&1 | awk '{print $2}')
    [ -n "${got}" ] || fail "could not read a version from ${SAGE}"
    [ "${got}" = "${SAGE_EXPECTED_VERSION}" ] \
        || fail "SAGE version mismatch: ${SAGE} reports '${got}', expected '${SAGE_EXPECTED_VERSION}'.
       Set FL_SAGE and FL_SAGE_VERSION together, deliberately.
       0.14.6        = the published v27 engine; CANNOT do MS1/LFQ on Bruker .d.
       0.15.0-beta.1 = required for timsTOF .d LFQ.
       0.15.0-beta.2 = NEWER but ABORTS on .d -- do not use for Bruker input."
    echo "[ok] SAGE ${got}  ${SAGE}  sha256=$(sha256sum "${SAGE}" | cut -d' ' -f1)"
}

# Records the engine beside the data. assert_code_committed() covers the Rust binaries
# but never saw SAGE, so the one binary that decides every identification was the one
# with no provenance row.
record_sage_provenance() {
    local out_dir="$1"
    mkdir -p "${out_dir}"
    local manifest="${out_dir}/RUN_CODE_PROVENANCE.tsv"
    [ -s "${manifest}" ] || printf 'utc\tstage\tbinary\tbinary_sha256\tcrate\tbranch\tcommit\tsource_state\n' > "${manifest}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "${STAGE_NAME:-unknown}" \
        "sage" "$(sha256sum "${SAGE}" | cut -d' ' -f1)" \
        "EXTERNAL-BINARY" "-" "$("${SAGE}" --version 2>&1 | awk '{print $2}')" "pinned" \
        >> "${manifest}"
}
