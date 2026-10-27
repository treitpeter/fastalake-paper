# Assembly: Supplementary Figure 12 (v30)

The v30 figure is one matplotlib canvas with six panels written by one script. No manual editing
took place. `LAB_ROOT/.../final/` abbreviates
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/`.

1. Draw. Inputs: the v29 build's exported source data
   `LAB_ROOT/.../final/suppfigs_n15/SF12/alloc347/build_n15a/source_data/` (here
   `source_data/SFig12/KO_allocation_coverage.tsv` (a), `healthy_stool_leader_diversity.tsv` (b),
   `source_data/Fig5/contributor_cells.tsv` (c, d; same rows and values as the build's
   `CapScan_function_contributors.tsv`, written with full float precision),
   `healthy_stool_grouping_leader_sensitivity.tsv` (e), `CapScan_grouping_leader_sensitivity.tsv`
   (f); the other four tables byte-identical to the build's). Tool:
   `LAB_ROOT/.../final/v30/figures/SFig12/build_sfig12.py` (sha256
   5fc81f1c74dc6f772ae7d50fc0621de3c67c077e09a72bad2d04d19c2b622902), repository
   `figures/SFig12/panels.py:render`. What changed relative to the v28f figure: e keeps only the
   four direct-peptide-assignment points and f only the `A33_CAPSCAN__simplified` line (original
   protein-grouping series and legends removed); a-d unchanged in content; on-canvas subpanel titles
   and the suptitle removed; fl_style Söhne, Okabe-Ito colours, per-participant ECDFs in one colour,
   x labels "Healthy stool"/"CapScan", 6 pt at 180 x 172 mm placed size. Outputs:
   `SupplementaryFigure12.pdf`, `SupplementaryFigure12.png` (300 dpi) and `plotted_summary.tsv`.
2. Placement in the manuscript: the PNG of step 1, unchanged.

Composite `expected/SFig12_v30.png`: sha256
e644ae023b54c93e6926e4041b7dfdfcb2fbd566f6d3b02a0d02e48f3a6c211c.
Composite PDF `expected/SFig12_v30.pdf`: sha256
fc58c6f14d523a36deefc4a15b11207748fc0c17eba6cfada064ebf8ede13c35 (drawn from the same figure
object; the PNG was not rasterised from it). Plotted summary `expected/plotted_summary_v30.tsv`:
sha256 548d5790055d8754f15f15d050ef4520802abe2f32000968e0e84ac9f910681b.

The notebook repeats step 1 from `source_data/`; PNG and PDF (at 300 dpi) are pixel-identical,
0 differing pixels (`checks/SFig12_v30.json`, `checks/SFig12_v30_pdf.json`), and the plotted
summary is byte-identical.

## v35 label change

Panels b and e: internal series names replaced by the series names of main Fig. 3 (`SERIES_LABEL`; Fig3/panels.py DESIGNS/DESIGN_LABELS): Interindividual -> Inter-individual, Intra -> Intra-individual, OnePiece -> Spatial sampling, TechnicalReplicates -> Technical replicates, WorkflowReplicates -> Workflow replicates. Panel e ticks read "<series> (primary)" / "<series> (well resolved)". Drawn order and plotted summary unchanged. Expected render replaced: `expected/SFig12_v35.png` sha256 4888ee3ca2b5be196ce4ba145eed05a8e9c2f739ea56a7049b56138ec60b1564, `expected/SFig12_v35.pdf` sha256 f2aa7f7342c489a4f1d41b463981f0dc53eaf93fdf69a70e95de7da9b45691e0. A control build of the unchanged code reproduced the previous expected PNG byte for byte.
