# Terminal-residue counts, route overlap and local source ambiguity (A20)

Internal ID: A20 (previous folder name A20_terminal_routes)

Pipeline id `terminal_route_counts`. Used by: Fig6 (panels a to c).

## Question

For each acquisition, which peptides does each search route (tryptic, semi-enzymatic, fully non-specific, de novo top 30 percent) contribute, and what are their N-terminal, N-cleavage (P1) and C-terminal residues and local parent source? Answer: a recount of the accepted peptides of the existing searches (no search is repeated) followed by a collector that averages acquisitions to specimens and participants. Fig6 reads the counts and peptide-origin labels; the analysis was not re-run for the paper.

## Inputs

- The completed A05 (CapScan) and A13 (PREDICT) non-specific results, the production tryptic and archived semi-enzymatic CapScan results, the production-retained de novo prediction tables, and each acquisition searched FASTA (`inputs/manifest.json`, manifest sha256 in `RECEIPTS.md`).
- Declared specimen and participant mapping of the CapScan roster; PREDICT uses the tube mapping held in the controlled design table (not in the repository).

## Outputs

- `results/<cohort>/<acquisition>/COMPLETE.json` and per-acquisition route counts, residue profiles and local-origin annotations.
- `results/summary_<cohort>/`: `acquisition_design.tsv`, `residue_summary.tsv`, `per_acquisition_residues.tsv.gz`, `per_specimen_residues.tsv.gz`, `per_participant_region_residues.tsv.gz`, `per_acquisition_route_identity.tsv.gz`, `per_acquisition_context_bounds.tsv.gz`, `histidine_*` tables, `participant_exclusions.tsv`.

## Call chain and code

`search.py`, `digest.py`, `parameters.py`, `origin.py` (plus `lib.py`, `peptide.py` from A05) are the call chain of the per-acquisition stage; `collector.py` builds the summaries. `prep_chain.py` and `preflight_chain.py` built and checked the participant-15 run directories for A20, A24, A26, A27 and A28. `digest.py` is byte-identical to the copies in A15, A18 and A21.

| file | origin (LAB_ROOT relative) | sha256 at origin | sha256 here |
|---|---|---|---|
| `search.py` | `LAB_ROOT/MM_review_analysis/A20_4_7_terminal_routes/frozen/search.py` | `52a5b29a135817e3c68fe159f28a7cce090123853e9e12fd664f33e070dfcdfb` | `52a5b29a135817e3c68fe159f28a7cce090123853e9e12fd664f33e070dfcdfb` |
| `collector.py` | `LAB_ROOT/MM_review_analysis/A20_4_7_terminal_routes/frozen/collector.py` | `1be7bbef21521f4527a52beccbe6c77f86abbbac82daf6550d5d5cccc95108b2` | `1be7bbef21521f4527a52beccbe6c77f86abbbac82daf6550d5d5cccc95108b2` |
| `digest.py` | `LAB_ROOT/MM_review_analysis/A20_4_7_terminal_routes/frozen/digest.py` | `2c9d232df04b795e0cea3e7b30e2cee508710b96cb8b78df3f94744a7f307245` | `2c9d232df04b795e0cea3e7b30e2cee508710b96cb8b78df3f94744a7f307245` |
| `parameters.py` | `LAB_ROOT/MM_review_analysis/A20_4_7_terminal_routes/frozen/parameters.py` | `ecef25b4ba6520a256a994607010c8664d7b315a434a9109bfd56f88dd94db06` | `ecef25b4ba6520a256a994607010c8664d7b315a434a9109bfd56f88dd94db06` |
| `origin.py` | `LAB_ROOT/MM_review_analysis/A20_4_7_terminal_routes/frozen/origin.py` | `2e8a6c148d40cd083fe1db6bb44cbb9ceb1f04abbe8e30779fc5c18132085198` | `2e02ad8408aa428041200a01fd734a98f7c230732e69fdc919634d721de7a82b` |
| `prep_chain.py` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/scripts/prep_chain.py` | `1310c40e45aecc5558c19d66bfa00b39696d4334e81b32f8cca7533998501f81` | `1310c40e45aecc5558c19d66bfa00b39696d4334e81b32f8cca7533998501f81` |
| `preflight_chain.py` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/scripts/preflight_chain.py` | `03995877515d798b98528fe203b3dfb24b0251010e27f7685a9010b5c2d04b05` | `03995877515d798b98528fe203b3dfb24b0251010e27f7685a9010b5c2d04b05` |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| populations | primary 5 to 50 residues and nested 9 to 50, always separate | 9 to 50 is the de novo comparison window | no pooling across lengths |
| de novo route | positive score, 9 to 50 residues, not low-complexity, top 30 percent by score including ties (`lib.load_queries`); no database peptide q-value exists for these rows | the exact production selection | score >= 0.9 route is A26 |
| compatible local parents | K/R not before P, at most two missed cleavages, no initiator-M clipping; one enzyme-compatible end for semi; unrestricted for non-specific and de novo | a peptide can be reported as cleavage only where a parent supports it | A05 context summary clipped initiator M; superseded for the new annotations only |
| origin class | host, microbial, mixed, unknown; any set containing an unreadable or hash accession (`PREFIX_<64 hex>`) is Unknown; the readable-only classification is kept as a named sensitivity | a hash accession carries no source | `origin.py` |
| residue profiles | all 20 residues for N-first, N-cleavage P1 (equal weight per distinct eligible context, protein ends separate) and C-last; canonical profiles collapse I/L, bare-sequence and PSM profiles keep assigned I/L | spectra do not discriminate I from L | no residue selected post hoc |
| averaging unit | acquisitions to specimens to participant and region; PREDICT averages specimens per patient | no peptide or PSM pseudo-replication | participants without mapping stay acquisition-only |
| exact reproduction | tryptic and non-specific accepted identities recounted by two independent readers and required to agree | guards the recounts | none |
| resources | 1 CPU, 16 GB, 2 h per acquisition task; collectors 16 GB, 3 h | measured peak 0.6 GB per task and 2.7 GB for the collector | none |

## Resources and run time (measured)

Per-acquisition tasks (jobs 5927974, 5927975, 5927976; 634 tasks): 1 CPU, 16 GB requested, wall time median about 0.3 to 0.4 min, maximum 0.7 min, peak 0.6 GB. CapScan collector (job 5928114): 4.4 min, 2.7 GB. The PREDICT collector submitted in the same wave (5928115) was cancelled and replaced by the repaired-PREDICT collector of the recovery wave (job 5937595, 32 GB, 3.6 min, 2.3 GB; code `terminal_histidine_inference/collect_repaired_predict_20260915.py`). Participant-15 stage (job 6058103, 16 tasks): 2 GB, maximum 0.6 min.

Per-job accounting is in `RECEIPTS.md`.

## How to run

Every script takes a run directory as first argument and, for array stages, a task index as second argument. The run directory holds `inputs/manifest*.json` (task list with the path, size, modification time and sha256 of every input), `frozen/` (this pipeline's code, assembled by `capscan_nonspecific_search/slurm/assemble_run_dir.sh`) and `results/`.

    export FASTALAKE_PIPELINES_DIR=<repo>/pipelines
    export FASTALAKE_OUT_DIR=<run directory>      # contains inputs/manifest*.json
    export FASTALAKE_PYTHON=<python with the packages named below>
    sbatch --array=<first>-<last>%<concurrency> <flags from slurm/submit_example.sh> slurm/run_task.sbatch

The manifests list original instrument-derived inputs by absolute path and are not part of the repository (they name lab locations; their sha256 is in `RECEIPTS.md`). A re-run therefore needs a manifest rebuilt for the new location; the frozen code verifies every input against the manifest (path, size, modification time, sha256) and stops on any difference. Per the repository conventions this pipeline is not re-run for the paper; the code is the code that produced the outputs.

## Receipts and open items

`RECEIPTS.md` lists the original job ids and the sha256 of every output. `OPEN_ITEMS.md` lists what could not be traced.
