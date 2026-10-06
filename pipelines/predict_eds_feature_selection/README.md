# predict_eds_feature_selection

Internal ID: A19 (previous folder name A19_eds_original_screen)

Which 20 features form the fixed EDS overview (Fig4 h, SFig10 b, STab4 Fixed_feature_key and 20-feature sheets)? This is the original screen that chose them; the later pipelines (A70, A71, A78) do not reselect.

Status: code copied as it ran (path edits only), not re-run in this repository. See RECEIPTS.md for original job ids and hashes.

## Inputs

Not traced. The selection code and its inputs could not be located in the lab tree.

## Outputs

`legacy_features.tsv` (the 20 selected features, as carried into the repaired analyses); see RECEIPTS.md.

## Slurm resources and run time (original logs, `sacct`)

Not recorded; no original job id could be traced.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| selection rule | not traced | the repaired analyses (A70) take the 20 features as fixed and do not reselect | none |

## Notes

No code is shipped because none could be traced to the output. The list is an input to A70 (`heatmap_features.tsv` was derived from it), not a result of any code in this repository. See OPEN_ITEMS.md.
