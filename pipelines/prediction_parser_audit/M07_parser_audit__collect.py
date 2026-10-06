"""Collect exact M07 parser comparisons; fail on incomplete or changed inputs."""
from pathlib import Path
from collections import Counter,defaultdict
import csv,hashlib,json,datetime,subprocess
O=Path(__file__).resolve().parent
m=json.loads((O/"INPUT_MANIFEST.json").read_text())
rows=[json.loads((O/"results"/("%03d.json"%t["index"])).read_text()) for t in m["tasks"]]
checks={}
checks["all_710_exact_roster_keys"]=[(r["index"],r["cohort"],r["sample"]) for r in rows]==[(t["index"],t["cohort"],t["sample"]) for t in m["tasks"]]
checks["all_input_sizes_match"]=all(r["input_size_matches_manifest"] for r in rows)
checks["all_recorded_prior_hashes_match"]=all(r["input_hash_matches_prior"] is not False for r in rows)
checks["historical_current_lookup_rows_scores_exact"]=all(r["historical_current_lookup_rows_exact"] for r in rows)
checks["historical_current_selection_keys_scores_exact"]=all(r["historical_current_selection_exact"] for r in rows)
checks["all_710_header_and_csv_scans_clean"]=all(r["counts"]["csv_errors"]==0 and r["counts"]["unequal_field_count"]==0 for r in rows)
checks["score_fields_present_parseable_finite"]=all(r["counts"]["missing_score"]==0 and r["counts"]["unparseable_score"]==0 and r["counts"]["nonfinite_score"]==0 for r in rows)
checks["rank_zero_only"]=all(r["counts"]["nonzero_beam_rank"]==0 for r in rows)
checks["selected_fasta_ids_unique"]=all(not r["selected_duplicate_ids"] for r in rows)
checks["selected_fasta_sequences_uppercase_only"]=all(r["selected_nonuppercase_sequence_bytes"]==0 for r in rows)
checks["all_top30_queries_are_selection_keys"]=all(not r["lookup_top30_absent_from_selection"] for r in rows)
summaries=[];cohorts=defaultdict(list)
for r in rows:cohorts[r["cohort"]].append(r)
for cohort,rr in sorted(cohorts.items()):
 c=Counter()
 for r in rr:c.update(r["counts"])
 summaries.append(dict(cohort=cohort,files=len(rr),prediction_bytes=sum(r["prediction_bytes"] for r in rr),**dict(c),prior_hashes=sum(r["expected_sha256"] is not None for r in rr),lookup_equal_files=sum(r["historical_current_lookup_rows_exact"] for r in rr),selection_equal_files=sum(r["historical_current_selection_exact"] for r in rr),lookup_selection_query_mismatch_files=sum(bool(r["lookup_top30_absent_from_selection"]) for r in rr),top30_sequence_parser_difference_files=sum(bool(r["top30_queries_from_different_parsed_strings"]) for r in rr)))
def tsv(name,rr):
 if not rr:return
 keys=list(dict.fromkeys(k for r in rr for k in r))
 with (O/name).open("w",newline="") as f:
  w=csv.DictWriter(f,keys,delimiter="\t");w.writeheader();w.writerows(rr)
tsv("COHORT_SUMMARY.tsv",summaries)
tsv("ACQUISITION_SUMMARY.tsv",[dict(index=r["index"],cohort=r["cohort"],sample=r["sample"],prediction_sha256=r["prediction_sha256"],prediction_bytes=r["prediction_bytes"],prior_hash_matches=r["input_hash_matches_prior"],**r["counts"],lookup_rows_exact=r["historical_current_lookup_rows_exact"],selection_exact=r["historical_current_selection_exact"],top30_unique=r["historical_lookup_unique_top30"],selection_unique=r["historical_selection_unique"],selected_records=r["selected_records"],lookup_only_queries=len(r["lookup_top30_absent_from_selection"]),seconds=r["seconds"]) for r in rows])
special=[dict(index=r["index"],cohort=r["cohort"],sample=r["sample"],**c) for r in rows for c in r["parser_difference_cases"]]
tsv("ALL_ELIGIBLE_PARSER_DIFFERENCES.tsv",special)
out=dict(checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status="PASS" if all(checks.values()) else "FINDINGS",checks=checks,cohorts=summaries,total_rows=sum(r["counts"]["rows"] for r in rows),files=710,scope="Exact source-function input equivalence and selected-FASTA containment checks. Not a current-main full workflow rerun, candidate-FASTA uniqueness audit, or reconstruction of unknown upstream AlphaNovo inference states.")
(O/"VALIDATION.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))

