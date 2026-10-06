# Assembly: Supplementary Figure 9 (v30)

Panels a and b are one matplotlib canvas; no panel file was combined and no manual step took
place. `FIN` is
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final`.

1. 9-50-residue table. Inputs `source_data/SFig09/capscan_assembly_dna_yield_gain_by_acquisition.tsv`
   and `a08_benchmark282_old.tsv`. Tool `FIN/len950/SFig09/build/build.py` (here
   `panels.py:acquisitions_9_50`): replaces each acquisition's FastaLake and MG peptide counts
   by the 9-50 counts. Output `FIN/len950/SFig09/capscan_acquisitions_9_50.tsv`, sha256
   5bb9ec898849e853e93ef0d01da9b7eaebab64db2b1a075b89f8d9d44fe5ad39
   (`source_data/SFig09/capscan_acquisitions_9_50.tsv`).
2. Producer drawing (v29). Tool `FIN/len950/SFig09/build/build_ylim.py` (here
   `panels.py:figure_v29`): the v28d producer `figure_ab` call for call, with the panel-a gain
   axis lower limit lowered to -1.07% so the -0.03% specimen shows. Output
   `FIN/v29/figures/SFig09/SupplementaryFigure9.png`, sha256
   8f9050ddd3f7ed63005c85f4cc630596f641a2791c714992d018ef74313c77d1
   (`expected/SFig09_v29_reference.png`).
3. House-style restyle. Tool `FIN/v30/figures_style/SFig09/build_sfig09_style.py` (here
   `panels.py:render_v30`): checks the step 2 render against the v29 PNG pixel for pixel,
   redraws the same values at 180 x 91.2 mm with `fl_style` (Söhne, Okabe-Ito region colours,
   composed 10^n tick labels, panel letters, titles removed) and compares every plotted array
   exactly (`values_check.json`, 24 checks). Outputs `SupplementaryFigure9.pdf`, sha256
   dc579ef30da325e72260ca105e620779feea34966ffb329ca0b62105021a40f5, and `SupplementaryFigure9.png`
   (300 dpi), sha256 360f1f3e69ad7299cdcca3b26ea2d0744d0ff2394ad1b07604a4d3a6c0550fdf.
4. Manuscript insertion: the PNG from step 3 at 180 x 91.2 mm.

`expected/SFig09_v30.png`: sha256 360f1f3e69ad7299cdcca3b26ea2d0744d0ff2394ad1b07604a4d3a6c0550fdf.
`expected/SFig09_v30.pdf` (same canvas, same step): sha256
dc579ef30da325e72260ca105e620779feea34966ffb329ca0b62105021a40f5.

The repository render is IDENTICAL to both at 300 dpi (`checks/SFig09_v30.json`,
`checks/SFig09_v30_pdf.json`). The v28d/v28f composite (5-50 residues) is superseded and no
longer kept in `expected/`. This figure will be rebuilt once more for the 295-acquisition
CapScan MG benchmark; the inputs of steps 1-3 are parameters in `panels.py`.
