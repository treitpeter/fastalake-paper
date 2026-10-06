from pathlib import Path
import sys,subprocess,json,datetime,math,re
from pipeline import W,check_frozen,dump,now,sha
from slurm_accounting import pilot_accounting

def sb(name,script,args,extra):
    cmd=['sbatch','--parsable','--partition=p.hpcl8,p.hpcl92,p.hpcl93','--nice=10000','--kill-on-invalid-dep=yes','--job-name='+name,'--chdir='+str(W),'--output='+str(W/'logs/%x_%A_%a.out'),'--error='+str(W/'logs/%x_%A_%a.err')]+extra+[str(W/script)]+list(map(str,args))
    job=subprocess.check_output(cmd,text=True).strip().split(';')[0]
    assert job.isdigit(),job
    record=dict(job_id=job,name=name,command=cmd,utc=now())
    with (W/'jobs/SUBMISSIONS.jsonl').open('a') as f:f.write(json.dumps(record)+'\n')
    return job
def launch():
    m=check_frozen();assert not (W/'jobs/SUBMITTED.json').exists()
    assert json.loads((W/'SELFTEST.json').read_text())['pass']
    inp=sb('FL_shfix_inputs','run.sbatch',['split',W],['--array=0-3%2','--cpus-per-task=1','--mem=4G','--time=02:00:00'])
    cache=sb('FL_shfix_cache','run.sbatch',['cache',W],['--array=0-7%1','--cpus-per-task=16','--mem=192G','--time=08:00:00','--dependency=afterok:'+inp])
    indices=[r['index'] for r in m['tasks'] if r['pilot']]
    pilot=sb('FL_shfix_pilot','run.sbatch',['worker',W],['--array='+','.join(map(str,indices))+'%4','--cpus-per-task=8','--mem=64G','--time=06:00:00','--dependency=afterok:'+cache])
    gate=sb('FL_shfix_gate','gate.sbatch',[W,pilot],['--cpus-per-task=1','--mem=2G','--time=00:20:00','--dependency=afterany:'+pilot])
    dump(W/'jobs/SUBMITTED.json',dict(utc=now(),inputs=inp,cache=cache,pilot=pilot,gate=gate,cancelled_superseded_launcher='6010380',rest_submitted=False,global_worker_throttle=4,manifest_sha256=sha(W/'MANIFEST.json')))
    print(json.dumps(json.loads((W/'jobs/SUBMITTED.json').read_text()),indent=2))
def gate(job):
    m=check_frozen();assert not (W/'jobs/REST_SUBMITTED.json').exists()
    pilot=[r for r in m['tasks'] if r['pilot']]
    raw=subprocess.check_output(['sacct','-j',job,'--format=JobID%64,State%32,ExitCode,MaxRSS','-P','--units=K'],text=True)
    (W/'jobs/PILOT_SACCT.tsv').write_text(raw)
    expected={r['index'] for r in pilot}
    states,rss_by_task=pilot_accounting(raw,job,expected)
    rss=list(rss_by_task.values())
    assert len(m['tasks'])==710 and len(pilot)==12
    assert [r['index'] for r in m['tasks']]==list(range(710))
    for co,cohort in m['cohorts'].items():
        tasks=[r for r in m['tasks'] if r['cohort']==co]
        assert len(tasks)==cohort['n'] and sum(r['pilot'] for r in tasks)==3
        receipt=W/f'inputs/{co}/COMPLETE.json'
        inputs=json.loads(receipt.read_text())
        assert inputs['n']==len(tasks)
        assert {r['sample'] for r in inputs['samples']}=={r['sample'] for r in tasks}
        for parity in ['even','odd']:
            cache=W/f'cache/{co}/{parity}'
            cj=json.loads((cache/'COMPLETE.json').read_text())
            assert cj['cohort']==co and cj['parity']==parity
            assert sha(receipt)==cj['input_receipt_sha256']
            assert sha(cache/'selected.fasta')==cj['fasta_sha256']
            for r in inputs['samples']:
                assert sha(W/f"inputs/{co}/{parity}/{r['sample']}_predictions.csv")==r['half_sha256'][parity]
    for r in pilot:
        out=W/f"outputs/{r['cohort']}/{r['sample']}"
        j=json.loads((out/'COMPLETE.json').read_text());assert j['full_reproduced'] is True
        assert (j['index'],j['cohort'],j['sample'],j['pilot'])==(r['index'],r['cohort'],r['sample'],True)
        assert json.loads((out/'FULL_REPRODUCTION.json').read_text())['equal'] is True
        assert set(j['metrics'])=={'FULL','HALF_SELF','HALF_DISJOINT'}
        p=json.loads((out/'PARITY.json').read_text());assert p['predictions_matched']==p['predictions'] and p['retained_odd_ms2']>0
        for arm in j['metrics']:
            assert sha(out/arm/'ACCEPTED.txt')==j['metrics'][arm]['accepted_sha256']
    memory=max(16,math.ceil(1.2*max(rss)+2))
    dump(W/'jobs/PILOT_GATE_PASS.json',dict(utc=now(),n=len(pilot),pilot_job=job,peak_rss_gib=max(rss),rss_gib_by_task=rss_by_task,rest_memory_gib=memory,manifest_sha256=sha(W/'MANIFEST.json'),scope='All rosters, eight cache hashes, half-input hashes, parity checks, native Sage q, FULL reproduction and success receipts passed'))
    indices=[r['index'] for r in m['tasks'] if not r['pilot']]
    # A fixed calendar timestamp prevents accidentally moving the repaired work to the next evening.
    rest=sb('FL_shfix_rest','run.sbatch',['worker',W],['--array='+','.join(map(str,indices))+'%4','--cpus-per-task=8','--mem='+str(memory)+'G','--time=06:00:00','--begin=2026-09-25T22:00:00'])
    dump(W/'jobs/REST_SUBMITTED.json',dict(utc=now(),rest=rest,evaluation=None,n=len(indices),memory_gib=memory))
    ev=sb('FL_shfix_summary','run.sbatch',['summary',W],['--cpus-per-task=1','--mem=2G','--time=00:20:00','--dependency=afterok:'+rest])
    dump(W/'jobs/REST_SUBMITTED.json',dict(utc=now(),rest=rest,evaluation=ev,n=len(indices),memory_gib=memory))
if __name__=='__main__':
    if sys.argv[1]=='launch':launch()
    elif sys.argv[1]=='gate':gate(sys.argv[2])
    else:raise ValueError(sys.argv[1])
