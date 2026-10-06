# Assembly: Supplementary Table 1

Supplementary Table 1 has no image and no schematic, so there is no ARTWORK_NUMBERS.md.
The v28f table files are written by script; nothing was edited by hand.

1. Inputs `LAB_ROOT/review/v15_integration/analyses/clustering_90_95_2026-09-29/figures/per_acquisition_identifications.tsv`
   (sha256 0187448f...1c06) and `.../per_acquisition_KO_coverage.tsv` (sha256 4dd098b1...6b23);
   in this repository `source_data/STab1/`. Tool: `LAB_ROOT/.../final/supp_table1/build_supp_table1.py`
   (Slurm job 6065679), ported unchanged as `tables/STab1/panels.py:build_table` and
   `write_table`. What changed: the two producer tables were joined one-to-one on
   acquisition and arm, renamed to 21 headed columns with units, and sorted by acquisition
   then arm order (unclustered, 99, 97, 95, 90, 75, 50%). Floats are written as their
   shortest round-trip text. Outputs:
   - `Supplementary_Table_1.tsv`, sha256 `4de2e9be55a9ec95587edb24a349803651ccccbce921f3ea1e5b414b0761c085`
   - `Supplementary_Table_1.xlsx` (one sheet `Supplementary_Table_1`, panes frozen at C2),
     sha256 `1bd83133aeca8a1cbe0b56970ad8c8a3e57e3a5b9a9e3a616e9dff39f7bdfd5e`
2. Legend: `LAB_ROOT/.../final/supp_table1/LEGEND.txt` (sha256 c0fb7167...8c8). The builder's docstring lists it as an
   output but the script does not write it; how it was written is not recorded; copied to `expected/LEGEND.txt`. Its first sentence is the v28d
   legend title (line 475); the rest is new text for v28f.
3. Placement in the Source Data package (pending, not done by this unit):
   `LAB_ROOT/.../final/supp_table1/REPLACE.md` describes replacing the 45-row
   `Supplementary_Tables/Supplementary_Table_1.xlsx` (sha256 9731b751...ba61) with the
   105-row workbook, adding the TSV as `Table1_support/Supplementary_Table_1_105_rows.tsv`
   and keeping the eight existing `Table1_support/*.tsv` files.

Regeneration in this repository (`python paper.py build STab1`) gives a TSV that is
byte-identical to step 1 and an XLSX that is cell-identical (2,226 cells, same sheet name
and frozen panes); XLSX bytes differ only because the file carries a creation timestamp.
Receipt: `checks/table.json`.
