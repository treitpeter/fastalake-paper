# Source data: Abstract and Discussion statistics

This directory holds no table. Every table the unit `text_statistics/abstract_discussion`
reads lives with the unit that first needed it and is read there (`source_data/README.md`:
a table is never duplicated). The data dictionaries are in those directories' README files.

| table | sha256 | used for (v28d line) |
|---|---|---|
| source_data/Fig4/crossfit_per_acquisition.tsv | 9672eafdba626af69782af5e4aa3a3f975e0d0a3ed8582b7013c5ffd1f2a5eb8 | cross-fitted gains, columns `FL_pooled_9`, `MG_pooled_9`, `cohort`, `sampling_site`, `participant_sort_rank` (35, 183) |
| source_data/Fig1/acquisition_sizes_726.tsv | fd9ade3bb4ebd2ee6c129f69819b9697cf20abd7c10c77a12644bef25396e386 | database entry counts, `entries`; reservoir sizes, `reservoir_entries` (35, 181, 193) |
| source_data/Fig1/human_reference.json | 64363c011e5b760663f7735ac0f4a7740d686d0c05af9e5e5e5098fd5af5648c | human reference entries, `entries` (35, 181) |
| source_data/SFig15/compute_stage_summary.tsv | 37c1d94d821d5ba9defb5ff5c8495905ce93209a136952dec8b600b74816498f | lookup and Sage wall times, `wall_median_min` (193) |
| source_data/SFig15/compute_task_records.tsv | fb248054124262d7deff9a9e7653b85d6c4f055c08081cc8b1c29c15ff458c2a | per-acquisition selection wall times, `wall_minutes` (193) |
| source_data/Fig2/member_summary.tsv | bdf6d5f0ded8bc8455f95bd1ac8c4ab52af38294ac7030959a269715873591b8 | CAMPI members, `member` (35) |

No clinical column is read. PREDICT patients appear only as pseudonymised cluster ranks
(`participant_sort_rank`) inside the Fig4 table.
