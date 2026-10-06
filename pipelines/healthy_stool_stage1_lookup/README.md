# healthy_stool_stage1_lookup

Internal ID: none; previous folder name UNASSIGNED_reference58_stage1_lookup

## Question

This pipeline produces: stage1_stats.json: 233 M reservoir, 18,007,876 study-supported proteins (Fig1 panel a, Results lines 51, 53; Supp Methods line 567)

Used by: figures/Fig1, text_statistics/supplementary_methods_1_5.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `01_evidence.sh` | launcher |
| `01b_sort.sh` | launcher |
| `_config.sh` | ============================================================================ |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `stage1_stats.json` | not found on disk |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

| job | name | cores | requested memory | peak RSS (GB) | elapsed | state |
|---|---|---|---|---|---|---|
| 5872582 | gate_HSTOOL | 2 | 4G | 0.0 | 00:00:41 | COMPLETED |
| 5872583 | ev_HSTOOL | 64 | 100G | 73.7 | 00:05:15 | COMPLETED |
| 5872584_0 | ps_HSTOOL | 16 | 32G | 4.3 | 00:01:10 | COMPLETED |
| 5872584_1 | ps_HSTOOL | 16 | 32G | 4.3 | 00:01:12 | COMPLETED |
| 5872584_2 | ps_HSTOOL | 16 | 32G | 4.2 | 00:01:12 | COMPLETED |
| 5872584_3 | ps_HSTOOL | 16 | 32G | 4.3 | 00:01:13 | COMPLETED |
| 5872584_4 | ps_HSTOOL | 16 | 32G | 4.9 | 00:01:13 | COMPLETED |
| 5872584_5 | ps_HSTOOL | 16 | 32G | 4.9 | 00:01:19 | COMPLETED |
| 5872584_6 | ps_HSTOOL | 16 | 32G | 4.9 | 00:01:17 | COMPLETED |
| 5872584_7 | ps_HSTOOL | 16 | 32G | 1.5 | 00:00:47 | COMPLETED |
| 5872584_8 | ps_HSTOOL | 16 | 32G | 2.9 | 00:00:47 | COMPLETED |
| 5872584_9 | ps_HSTOOL | 16 | 32G | 1.9 | 00:00:46 | COMPLETED |
| 5872584_10 | ps_HSTOOL | 16 | 32G | 4.3 | 00:00:43 | COMPLETED |
| 5872584_11 | ps_HSTOOL | 16 | 32G | 4.9 | 00:00:44 | COMPLETED |
| 5872584_12 | ps_HSTOOL | 16 | 32G | 5.2 | 00:00:46 | COMPLETED |
| 5872584_13 | ps_HSTOOL | 16 | 32G | 4.2 | 00:00:49 | COMPLETED |
| 5872584_14 | ps_HSTOOL | 16 | 32G | 1.2 | 00:00:56 | COMPLETED |
| 5872584_15 | ps_HSTOOL | 16 | 32G | 5.2 | 00:00:45 | COMPLETED |
| 5872584_16 | ps_HSTOOL | 16 | 32G | 4.4 | 00:00:45 | COMPLETED |
| 5872584_17 | ps_HSTOOL | 16 | 32G | 4.3 | 00:00:45 | COMPLETED |
| 5872584_18 | ps_HSTOOL | 16 | 32G | 2.8 | 00:00:57 | COMPLETED |
| 5872584_19 | ps_HSTOOL | 16 | 32G | 3.6 | 00:00:49 | COMPLETED |
| 5872584_20 | ps_HSTOOL | 16 | 32G | 1.2 | 00:00:43 | COMPLETED |
| 5872584_21 | ps_HSTOOL | 16 | 32G | 5.3 | 00:00:45 | COMPLETED |
| 5872584_22 | ps_HSTOOL | 16 | 32G | 4.5 | 00:00:48 | COMPLETED |
| 5872584_23 | ps_HSTOOL | 16 | 32G | 4.5 | 00:00:43 | COMPLETED |
| 5872584_24 | ps_HSTOOL | 16 | 32G | 4.9 | 00:00:48 | COMPLETED |
| 5872584_25 | ps_HSTOOL | 16 | 32G | 1.6 | 00:00:51 | COMPLETED |
| 5872584_26 | ps_HSTOOL | 16 | 32G | 5.1 | 00:00:47 | COMPLETED |
| 5872584_27 | ps_HSTOOL | 16 | 32G | 4.7 | 00:00:50 | COMPLETED |
| 5872584_28 | ps_HSTOOL | 16 | 32G | 3.4 | 00:00:47 | COMPLETED |
| 5872584_29 | ps_HSTOOL | 16 | 32G | 4.9 | 00:00:46 | COMPLETED |
| 5872584_30 | ps_HSTOOL | 16 | 32G | 1.2 | 00:00:50 | COMPLETED |
| 5872584_31 | ps_HSTOOL | 16 | 32G | 1.3 | 00:00:51 | COMPLETED |
| 5872584_32 | ps_HSTOOL | 16 | 32G | 1.5 | 00:00:52 | COMPLETED |
| 5872584_33 | ps_HSTOOL | 16 | 32G | 1.5 | 00:00:57 | COMPLETED |
| 5872584_34 | ps_HSTOOL | 16 | 32G | 3.1 | 00:00:47 | COMPLETED |
| 5872584_35 | ps_HSTOOL | 16 | 32G | 5.0 | 00:00:44 | COMPLETED |
| 5872584_36 | ps_HSTOOL | 16 | 32G | 3.5 | 00:00:46 | COMPLETED |
| 5872584_37 | ps_HSTOOL | 16 | 32G | 3.1 | 00:00:47 | COMPLETED |
| 5872584_38 | ps_HSTOOL | 16 | 32G | 1.1 | 00:00:48 | COMPLETED |
| 5872584_39 | ps_HSTOOL | 16 | 32G | 5.2 | 00:00:45 | COMPLETED |
| 5872584_40 | ps_HSTOOL | 16 | 32G | 4.9 | 00:00:47 | COMPLETED |
| 5872584_41 | ps_HSTOOL | 16 | 32G | 3.9 | 00:00:49 | COMPLETED |
| 5872584_42 | ps_HSTOOL | 16 | 32G | 5.4 | 00:00:46 | COMPLETED |
| 5872584_43 | ps_HSTOOL | 16 | 32G | 5.0 | 00:00:47 | COMPLETED |
| 5872584_44 | ps_HSTOOL | 16 | 32G | 4.5 | 00:00:45 | COMPLETED |
| 5872584_45 | ps_HSTOOL | 16 | 32G | 5.4 | 00:00:45 | COMPLETED |
| 5872584_46 | ps_HSTOOL | 16 | 32G | 1.6 | 00:00:49 | COMPLETED |
| 5872584_47 | ps_HSTOOL | 16 | 32G | 4.5 | 00:00:51 | COMPLETED |
| 5872584_48 | ps_HSTOOL | 16 | 32G | 1.7 | 00:00:51 | COMPLETED |
| 5872584_49 | ps_HSTOOL | 16 | 32G | 5.8 | 00:00:46 | COMPLETED |
| 5872584_50 | ps_HSTOOL | 16 | 32G | 4.6 | 00:00:47 | COMPLETED |
| 5872584_51 | ps_HSTOOL | 16 | 32G | 5.1 | 00:00:47 | COMPLETED |
| 5872584_52 | ps_HSTOOL | 16 | 32G | 1.6 | 00:00:48 | COMPLETED |
| 5872584_53 | ps_HSTOOL | 16 | 32G | 5.4 | 00:00:49 | COMPLETED |
| 5872584_54 | ps_HSTOOL | 16 | 32G | 5.6 | 00:00:50 | COMPLETED |
| 5872584_55 | ps_HSTOOL | 16 | 32G | 3.9 | 00:00:57 | COMPLETED |
| 5872584_56 | ps_HSTOOL | 16 | 32G | 10.2 | 00:00:47 | COMPLETED |
| 5872584_57 | ps_HSTOOL | 16 | 32G | 6.4 | 00:00:47 | COMPLETED |
| 5872587 | agg_HSTOOL | 48 | 8G | 2.8 | 00:00:36 | COMPLETED |
| 5872585_0 | par_HSTOOL | 16 | 24G | 8.2 | 00:00:53 | COMPLETED |
| 5872585_1 | par_HSTOOL | 16 | 24G | 8.2 | 00:00:55 | COMPLETED |
| 5872585_2 | par_HSTOOL | 16 | 24G | 8.2 | 00:00:54 | COMPLETED |
| 5872585_3 | par_HSTOOL | 16 | 24G | 8.3 | 00:00:56 | COMPLETED |
| 5872585_4 | par_HSTOOL | 16 | 24G | 8.2 | 00:00:58 | COMPLETED |
| 5872585_5 | par_HSTOOL | 16 | 24G | 8.3 | 00:00:56 | COMPLETED |
| 5872585_6 | par_HSTOOL | 16 | 24G | 8.2 | 00:00:55 | COMPLETED |
| 5872585_7 | par_HSTOOL | 16 | 24G | 8.5 | 00:01:07 | COMPLETED |
| 5872585_8 | par_HSTOOL | 16 | 24G | 8.2 | 00:00:54 | COMPLETED |
| 5872585_9 | par_HSTOOL | 16 | 24G | 8.3 | 00:00:54 | COMPLETED |
| 5872585_10 | par_HSTOOL | 16 | 24G | 8.3 | 00:00:55 | COMPLETED |
| 5872585_11 | par_HSTOOL | 16 | 24G | 8.3 | 00:00:54 | COMPLETED |
| 5872585_12 | par_HSTOOL | 16 | 24G | 8.9 | 00:01:04 | COMPLETED |
| 5872585_13 | par_HSTOOL | 16 | 24G | 9.7 | 00:01:27 | COMPLETED |
| 5872585_14 | par_HSTOOL | 16 | 24G | 8.9 | 00:01:05 | COMPLETED |
| 5872585_15 | par_HSTOOL | 16 | 24G | 8.8 | 00:01:14 | COMPLETED |
| 5872585_16 | par_HSTOOL | 16 | 24G | 8.8 | 00:00:59 | COMPLETED |
| 5872585_17 | par_HSTOOL | 16 | 24G | 8.4 | 00:00:57 | COMPLETED |
| 5872585_18 | par_HSTOOL | 16 | 24G | 8.5 | 00:01:11 | COMPLETED |
| 5872585_19 | par_HSTOOL | 16 | 24G | 9.8 | 00:01:19 | COMPLETED |
| 5872585_20 | par_HSTOOL | 16 | 24G | 7.1 | 00:00:57 | COMPLETED |
| 5872585_21 | par_HSTOOL | 16 | 24G | 9.4 | 00:01:10 | COMPLETED |
| 5872585_22 | par_HSTOOL | 16 | 24G | 8.8 | 00:01:04 | COMPLETED |
| 5872585_23 | par_HSTOOL | 16 | 24G | 9.0 | 00:01:04 | COMPLETED |
| 5872585_24 | par_HSTOOL | 16 | 24G | 8.6 | 00:01:11 | COMPLETED |
| 5872585_25 | par_HSTOOL | 16 | 24G | 9.8 | 00:01:28 | COMPLETED |
| 5872585_26 | par_HSTOOL | 16 | 24G | 9.4 | 00:01:34 | COMPLETED |
| 5872585_27 | par_HSTOOL | 16 | 24G | 7.9 | 00:01:27 | COMPLETED |
| 5872585_28 | par_HSTOOL | 16 | 24G | 9.2 | 00:00:58 | COMPLETED |
| 5872585_29 | par_HSTOOL | 16 | 24G | 8.9 | 00:00:58 | COMPLETED |
| 5872585_30 | par_HSTOOL | 16 | 24G | 8.9 | 00:01:03 | COMPLETED |
| 5872585_31 | par_HSTOOL | 16 | 24G | 7.8 | 00:00:51 | COMPLETED |
| 5872585_32 | par_HSTOOL | 16 | 24G | 8.9 | 00:01:00 | COMPLETED |
| 5872585_33 | par_HSTOOL | 16 | 24G | 9.6 | 00:01:12 | COMPLETED |
| 5872585_34 | par_HSTOOL | 16 | 24G | 8.1 | 00:00:52 | COMPLETED |
| 5872585_35 | par_HSTOOL | 16 | 24G | 9.0 | 00:00:59 | COMPLETED |
| 5872585_36 | par_HSTOOL | 16 | 24G | 8.8 | 00:00:58 | COMPLETED |
| 5872585_37 | par_HSTOOL | 16 | 24G | 8.8 | 00:01:02 | COMPLETED |
| 5872585_38 | par_HSTOOL | 16 | 24G | 9.1 | 00:01:12 | COMPLETED |
| 5872585_39 | par_HSTOOL | 16 | 24G | 7.9 | 00:01:02 | COMPLETED |
| 5872585_40 | par_HSTOOL | 16 | 24G | 9.0 | 00:00:59 | COMPLETED |
| 5872585_41 | par_HSTOOL | 16 | 24G | 8.1 | 00:01:00 | COMPLETED |
| 5872585_42 | par_HSTOOL | 16 | 24G | 8.8 | 00:01:05 | COMPLETED |
| 5872585_43 | par_HSTOOL | 16 | 24G | 8.8 | 00:01:26 | COMPLETED |
| 5872585_44 | par_HSTOOL | 16 | 24G | 8.4 | 00:00:55 | COMPLETED |
| 5872585_45 | par_HSTOOL | 16 | 24G | 8.3 | 00:01:13 | COMPLETED |
| 5872585_46 | par_HSTOOL | 16 | 24G | 9.1 | 00:01:14 | COMPLETED |
| 5872585_47 | par_HSTOOL | 16 | 24G | 9.5 | 00:01:09 | COMPLETED |
| 5872585_48 | par_HSTOOL | 16 | 24G | 10.1 | 00:01:34 | COMPLETED |
| 5872585_49 | par_HSTOOL | 16 | 24G | 9.8 | 00:01:12 | COMPLETED |
| 5872585_50 | par_HSTOOL | 16 | 24G | 9.5 | 00:01:09 | COMPLETED |
| 5872585_51 | par_HSTOOL | 16 | 24G | 9.5 | 00:01:10 | COMPLETED |
| 5872585_52 | par_HSTOOL | 16 | 24G | 9.0 | 00:01:10 | COMPLETED |
| 5872585_53 | par_HSTOOL | 16 | 24G | 12.1 | 00:01:23 | COMPLETED |
| 5872585_54 | par_HSTOOL | 16 | 24G | 10.2 | 00:01:15 | COMPLETED |
| 5872585_55 | par_HSTOOL | 16 | 24G | 10.0 | 00:01:13 | COMPLETED |
| 5872585_56 | par_HSTOOL | 16 | 24G | 10.7 | 00:01:19 | COMPLETED |
| 5872585_57 | par_HSTOOL | 16 | 24G | 9.4 | 00:01:07 | COMPLETED |
| 5872586_0 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:27 | COMPLETED |
| 5872586_1 | sage_HSTOOL | 32 | 6G | 2.5 | 00:00:34 | COMPLETED |
| 5872586_2 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:26 | COMPLETED |
| 5872586_3 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:29 | COMPLETED |
| 5872586_4 | sage_HSTOOL | 32 | 6G | 2.8 | 00:00:34 | COMPLETED |
| 5872586_5 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:29 | COMPLETED |
| 5872586_6 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:26 | COMPLETED |
| 5872586_7 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:27 | COMPLETED |
| 5872586_8 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:30 | COMPLETED |
| 5872586_9 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:27 | COMPLETED |
| 5872586_10 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:27 | COMPLETED |
| 5872586_11 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:26 | COMPLETED |
| 5872586_12 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:30 | COMPLETED |
| 5872586_13 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:29 | COMPLETED |
| 5872586_14 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:28 | COMPLETED |
| 5872586_15 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:27 | COMPLETED |
| 5872586_16 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:28 | COMPLETED |
| 5872586_17 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:28 | COMPLETED |
| 5872586_18 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:30 | COMPLETED |
| 5872586_19 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:27 | COMPLETED |
| 5872586_20 | sage_HSTOOL | 32 | 6G | 2.0 | 00:00:33 | COMPLETED |
| 5872586_21 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:28 | COMPLETED |
| 5872586_22 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:29 | COMPLETED |
| 5872586_23 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:26 | COMPLETED |
| 5872586_24 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:30 | COMPLETED |
| 5872586_25 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:28 | COMPLETED |
| 5872586_26 | sage_HSTOOL | 32 | 6G | 2.1 | 00:00:32 | COMPLETED |
| 5872586_27 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:29 | COMPLETED |
| 5872586_28 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:28 | COMPLETED |
| 5872586_29 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:31 | COMPLETED |
| 5872586_30 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:31 | COMPLETED |
| 5872586_31 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:29 | COMPLETED |
| 5872586_32 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:31 | COMPLETED |
| 5872586_33 | sage_HSTOOL | 32 | 6G | 2.3 | 00:00:31 | COMPLETED |
| 5872586_34 | sage_HSTOOL | 32 | 6G | 2.3 | 00:00:33 | COMPLETED |
| 5872586_35 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:28 | COMPLETED |
| 5872586_36 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:29 | COMPLETED |
| 5872586_37 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:27 | COMPLETED |
| 5872586_38 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:26 | COMPLETED |
| 5872586_39 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:25 | COMPLETED |
| 5872586_40 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:31 | COMPLETED |
| 5872586_41 | sage_HSTOOL | 32 | 6G | 1.8 | 00:00:32 | COMPLETED |
| 5872586_42 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:28 | COMPLETED |
| 5872586_43 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:27 | COMPLETED |
| 5872586_44 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:26 | COMPLETED |
| 5872586_45 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:26 | COMPLETED |
| 5872586_46 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:31 | COMPLETED |
| 5872586_47 | sage_HSTOOL | 32 | 6G | 3.8 | 00:00:37 | COMPLETED |
| 5872586_48 | sage_HSTOOL | 32 | 6G | 3.3 | 00:00:39 | COMPLETED |
| 5872586_49 | sage_HSTOOL | 32 | 6G | 3.1 | 00:00:33 | COMPLETED |
| 5872586_50 | sage_HSTOOL | 32 | 6G | 3.4 | 00:00:37 | COMPLETED |
| 5872586_51 | sage_HSTOOL | 32 | 6G | 3.3 | 00:00:39 | COMPLETED |
| 5872586_52 | sage_HSTOOL | 32 | 6G | 3.1 | 00:00:33 | COMPLETED |
| 5872586_53 | sage_HSTOOL | 32 | 6G | 3.8 | 00:00:34 | COMPLETED |
| 5872586_54 | sage_HSTOOL | 32 | 6G | 3.2 | 00:00:37 | COMPLETED |
| 5872586_55 | sage_HSTOOL | 32 | 6G | 0.0 | 00:00:31 | COMPLETED |
| 5872586_56 | sage_HSTOOL | 32 | 6G | 3.5 | 00:00:41 | COMPLETED |
| 5872586_57 | sage_HSTOOL | 32 | 6G | 3.0 | 00:00:35 | COMPLETED |

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| none found as top-level constants | | | |
