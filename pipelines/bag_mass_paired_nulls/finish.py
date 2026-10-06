from pathlib import Path
import csv,json,sys
HERE=Path(__file__).resolve().parent
m=json.loads((HERE/'manifest.json').read_text())
mode=sys.argv[1]
indices=m['paired_pilot_indices'] if mode=='gate' else range(len(m['tasks']))
rows=[];errors=[];receipts=[]
for i in indices:
 t=m['tasks'][i];p=Path(m['output'])/t['cohort']/t['sample']
 try:
  r=json.loads((p/'NULL_COMPARISON.json').read_text())
  retrieval=json.loads((p/'NULL_RETRIEVAL.json').read_text())
  assert r['status']=='NULL_DIAGNOSTIC_COMPLETE_REVIEW_REQUIRED'
  assert r['workflow_fdr']=='UNVERIFIED'
  assert len(r['rows'])==80
  for seed in (1,2):
   selected=[x for x in retrieval['paired_comparisons'] if x['arm'].startswith(f'shuffle_{seed}')]
   assert len(selected)==4 and len({x['queried'] for x in selected})==1
  for row in r['rows']:
   assert row['net']==row['accepted']-row['exact_accepted']==row['gained']-row['lost']
   rows.append(dict(index=i,cohort=t['cohort'],sample=t['sample'],**row))
  receipts.append(str(p/'NULL_COMPARISON.json'))
 except (OSError,ValueError,KeyError,AssertionError) as e:errors.append(dict(index=i,error=repr(e)))
result=dict(status=('FAIL' if errors else 'PASS') if mode=='gate' else ('INCOMPLETE' if errors else 'DIAGNOSTICS_COMPLETE'),errors=errors,completed=len(receipts),expected=len(indices),receipts=receipts,workflow_fdr='UNVERIFIED',scientific_review='PENDING',positive_result_required=False)
(HERE/('PILOT_GATE.json' if mode=='gate' else 'COLLECTION.json')).write_text(json.dumps(result,indent=2)+'\n')
if mode=='collect' and rows:
 with (HERE/'NULL_COMPARISONS.tsv').open('w') as h:
  w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
if mode=='gate' and errors:raise SystemExit('Pilot mechanical validation failed; see PILOT_GATE.json')
