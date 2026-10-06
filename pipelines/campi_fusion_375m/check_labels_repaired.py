"""Read-only reservoir-presence audit; never an estimate of sample FDP.

Repairs: normalize BOTH sides; explicit q threshold; positives sampled from the
original reservoir rather than a fusion-selected FASTA; enforce positive gate.
Run --selftest, then --pilot, before --full. Original outputs are preserved.
"""
from pathlib import Path
import argparse, collections, csv, hashlib, json, random, sys, time
import ahocorasick

HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/'00_START_HERE_V7.md').exists())
FUS=ROOT/'00_PACKAGE_2026-09-17/cleanup_2026-09-22/campi_fusion_rerun'

def norm(s):return s.upper().replace('I','L')
def fasta(path,limit=None):
    with Path(path).open() as f:
        acc=None;seq=[];n=0
        for line in f:
            if line.startswith('>'):
                if acc is not None:
                    yield acc,''.join(seq)
                    n+=1
                    if limit is not None and n>=limit:return
                acc=line[1:].split()[0];seq=[]
            else:seq.append(line.strip())
        if acc is not None:yield acc,''.join(seq)

def matcher(patterns):
    a=ahocorasick.Automaton()
    for p in sorted(set(patterns)):a.add_word(norm(p),norm(p))
    a.make_automaton();return a

def scan_sequences(sequences,patterns):
    a=matcher(patterns);found=set()
    for s in sequences:
        for _,p in a.iter(norm(s)):found.add(p)
    return found

def read_controls(paths,q):
    out=set();all_rows=accepted=excluded=0;manifest=[]
    for path in paths:
        data=path.read_bytes();manifest.append({'path':str(path),'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data)})
        with path.open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                all_rows+=1
                if row['assignment'].strip()!='entrapment':continue
                value=float(row['peptide_q'])
                if value>q:excluded+=1;continue
                p=norm(row['peptide'])
                if 5<=len(p)<=50:out.add(p);accepted+=1
    return out,{'input_rows':all_rows,'accepted_entrapment_rows':accepted,'excluded_above_q':excluded,'inputs':manifest}

def shuffle(p,rng):
    chars=list(p)
    for _ in range(100):
        rng.shuffle(chars);x=''.join(chars)
        if x!=p:return x
    return None

def selftest():
    import tempfile
    # The old code normalizes queries but misses these raw I-containing positives.
    assert scan_sequences(['MAKPEPTIDEKGGIILIKQ'],['PEPTLDEK','LLLLK','ZZZZZ'])=={'PEPTLDEK','LLLLK'}
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)/'q.tsv';p.write_text('peptide\tpeptide_q\tassignment\nPEPTIDEK\t0.009\tentrapment\nOTHERK\t0.049\tentrapment\nSHAREDK\t0.001\tentrapment;target\n')
        peptides,meta=read_controls([p],.01)
        assert peptides=={'PEPTLDEK'} and meta['excluded_above_q']==1
    print('SELFTEST PASS: raw I/L normalization, q threshold, shared-label exclusion',flush=True)

def build(series,q):
    if series=='current375':
        reservoir=Path(json.loads((FUS/'inputs/PROVENANCE.json').read_text())['reservoir']['path'])
        files=sorted((FUS/'runs_v2/CAMPI_MOCK').glob('*/samples/*/independent_audit/accepted_sequences.tsv'))
        expected=375333654
    else:
        reservoir=ROOT/'canonical/publication_data/cohorts/CAMPI_MOCK/0_lake/CAMPI_MOCK_lake.fasta'
        files=sorted((ROOT/'analyses/fusion_entrapment_2026-09-07/runs_v2/CAMPI_MOCK').glob('*/samples/*/accepted_entrapment_queries.tsv'))
        expected=233477360
    assert len(files)==30,(series,len(files))
    entrap,meta=read_controls(files,q);assert entrap
    rng=random.Random(20260925);shuffled=set();unshufflable=0
    for p in sorted(entrap):
        s=shuffle(p,rng)
        if s is None:unshufflable+=1
        else:shuffled.add(s)
    # Select directly from the ORIGINAL unfused reference. Every positive is known present.
    originals=list(fasta(reservoir,10000));lengths=sorted(map(len,entrap));positive=set();positive_sources=[]
    for _ in range(4000):
        length=rng.choice(lengths)
        eligible=[i for i,(_,s) in enumerate(originals) if len(s)>=length]
        i=rng.choice(eligible);acc,s=originals[i];start=rng.randrange(len(s)-length+1)
        p=norm(s[start:start+length]);positive.add(p);positive_sources.append({'accession':acc,'start':start,'raw':s[start:start+length],'normalized':p})
    arms={'ENTRAP':entrap,'SHUFCTL':shuffled,'REALPOS':positive}
    return reservoir,expected,arms,dict(meta,unshufflable=unshufflable,positive_sources=positive_sources),originals

def summarize(arms,found):
    return [{'arm':arm,'band':f'{lo}-{hi}','present':sum(p in found for p in peps if lo<=len(p)<=hi),'total':sum(lo<=len(p)<=hi for p in peps)} for arm,peps in arms.items() for lo,hi in [(5,8),(9,50),(5,50)]]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--series',choices=['current375','historical233'],default='current375');ap.add_argument('--q',type=float,default=.01)
    ap.add_argument('--selftest',action='store_true');ap.add_argument('--pilot',action='store_true');ap.add_argument('--full',action='store_true');a=ap.parse_args();selftest()
    if a.selftest:return
    start=time.time();dest=HERE/a.series;dest.mkdir(exist_ok=True)
    reservoir,expected,arms,meta,originals=build(a.series,a.q)
    patterns=set().union(*arms.values());found=scan_sequences([s for _,s in originals],patterns)
    assert arms['REALPOS']<=found,'Positive gate failed on original-reservoir pilot'
    pil={'status':'PASS','series':a.series,'q':a.q,'reservoir':str(reservoir),'positive_total':len(arms['REALPOS']),'positive_recovered':len(arms['REALPOS']&found),'prefix_sequences':len(originals),'positive_with_raw_I':sum('I' in p['raw'] for p in meta['positive_sources']),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'arm_counts':{k:len(v) for k,v in arms.items()}}
    (dest/'PILOT.json').write_text(json.dumps(pil,indent=2));(dest/'INPUTS.json').write_text(json.dumps(meta,indent=2))
    print(json.dumps(pil),flush=True)
    if not a.full:return
    before=reservoir.stat();found=set();automaton=matcher(patterns);count=0
    for _,s in fasta(reservoir):
        count+=1
        for _,p in automaton.iter(norm(s)):found.add(p)
        if count%10000000==0:print(f'{count:,} sequences; {time.time()-start:.0f}s',flush=True)
    after=reservoir.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
    report={'status':'PASS' if count==expected and arms['REALPOS']<=found else 'FAIL','series':a.series,'q':a.q,'reservoir':str(reservoir),'reservoir_bytes':before.st_size,'reservoir_mtime_ns':before.st_mtime_ns,'reservoir_sequences':count,'expected_sequences':expected,'positive_total':len(arms['REALPOS']),'positive_recovered':len(arms['REALPOS']&found),'by_band':summarize(arms,found),'by_length':[{'arm':k,'length':n,'present':sum(p in found for p in v if len(p)==n),'total':sum(len(p)==n for p in v)} for k,v in arms.items() for n in range(5,51)],'seconds':time.time()-start,'whole_procedure_fdr_control_established':False,'meaning':'Reservoir presence only; neither sample presence nor corrected FDP. Distinct sequence denominator; not acquisition-weighted diagnostic denominator.','script_sha256':pil['script_sha256']}
    (dest/'RESULT.json').write_text(json.dumps(report,indent=2))
    with (dest/'PRESENCE.tsv').open('w') as f:
        f.write('arm\tpeptide\tpresent\n')
        for arm,peps in arms.items():
            for p in sorted(peps):f.write(f'{arm}\t{p}\t{int(p in found)}\n')
    assert report['status']=='PASS',report
    print(json.dumps(report['by_band']),flush=True)

if __name__=='__main__':main()
