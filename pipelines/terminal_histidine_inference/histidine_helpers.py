"""Protein-coordinate evidence and explicitly conditional sequence backgrounds."""
import hashlib
import itertools
from collections import Counter
import numpy as np

AA = 'ACDEFGHIKLMNPQRSTVWY'
MASS = dict(zip(AA, [71.03711,103.00919,115.02694,129.04259,147.0684,57.02146,
                    137.05891,113.08406,128.09496,113.08406,131.0405,114.04293,
                    97.05276,128.05858,156.1011,87.03203,101.04768,99.06841,
                    186.07932,163.06332]))
REGIONS = ['duodenum','jejunum','ileum','colon','stool']

def protein_key(seq):
    return hashlib.sha256(seq.replace('I','L').encode()).hexdigest()

def cuts(seq):
    return np.array([True]+[seq[i-1] in 'KR' and seq[i]!='P' for i in range(1,len(seq))]+[True])

def compatible(seq,start,end,enzyme):
    c=cuts(seq)
    return enzyme=='unspecific' or (int(c[start+1:end].sum())<=2 and
        ((c[start] and c[end]) if enzyme=='tryptic' else (c[start] or c[end])))

def windows(seq,positions):
    return [seq[max(0,int(p)-4):int(p)].rjust(4,'-')+seq[int(p):int(p)+4].ljust(4,'-') for p in positions]

def opportunities(seq,enzyme,side,mass_gate=True):
    """Internal boundaries allowing at least one 9-50 peptide on this side.

    The mass gate is 500-5000 Da with fixed C+57.0215 and 0-2 M+15.9949.
    No ionization/fragmentation detectability is inferred from this sequence null.
    """
    n=len(seq);ok=np.zeros(n+1,bool);c=cuts(seq);prefix=np.r_[0,np.cumsum(c)]
    unknown=np.r_[0,np.cumsum([a not in MASS for a in seq])]
    masses=np.r_[0.,np.cumsum([MASS.get(a,0)+(57.0215 if a=='C' else 0) for a in seq])]
    meth=np.r_[0,np.cumsum([a=='M' for a in seq])]
    for length in range(9,min(50,n)+1):
        starts=np.arange(n-length+1);ends=starts+length
        keep=(unknown[ends]-unknown[starts])==0
        if mass_gate:
            base=masses[ends]-masses[starts]+18.010565
            mods=np.minimum(2,meth[ends]-meth[starts])
            keep &= (base<=5000) & (base+mods*15.9949>=500)
        if enzyme!='unspecific':
            keep &= prefix[ends]-prefix[starts+1]<=2
            keep &= (c[starts]&c[ends]) if enzyme=='tryptic' else (c[starts]|c[ends])
        ok[(ends if side=='C' else starts)[keep]]=True
    ok[0]=ok[-1]=False
    return np.flatnonzero(ok)

def background(seq,enzyme,side):
    pp=opportunities(seq,enzyme,side)
    pp=np.array([p for p in pp if seq[p-1] in AA and seq[p-1] not in 'KR'],int)
    all_internal=np.array([p for p in range(1,len(seq)) if seq[p-1] in AA and seq[p-1] not in 'KR'],int)
    def counts(pos):
        mat=np.zeros((8,20),float)
        for w in windows(seq,pos):
            for j,a in enumerate(w):
                if a in AA:mat[j,AA.index(a)]+=1
        return mat
    # Motif backgrounds condition on P1=H for H-specific motifs; do not rediscover H.
    return dict(eligible=counts(pp),H=counts([p for p in pp if seq[p-1]=='H']),
                internal=counts(all_internal),eligible_positions=len(pp))

def brute_validate(seq,enzyme,side):
    expected=set()
    for start in range(len(seq)):
        for end in range(start+9,min(len(seq),start+50)+1):
            peptide=seq[start:end]
            if not set(peptide)<=set(AA):continue
            base=sum(MASS[a]+(57.0215 if a=='C' else 0) for a in peptide)+18.010565
            if not any(500<=base+v*15.9949<=5000 for v in range(min(2,peptide.count('M'))+1)):continue
            if compatible(seq,start,end,enzyme):expected.add(end if side=='C' else start)
    expected-={0,len(seq)}
    assert expected==set(opportunities(seq,enzyme,side)),(enzyme,side,seq)

_null_cache = {}
_boot_cache = {}

def signflip(delta):
    d=np.asarray(delta,float);d=d[np.isfinite(d)];n=len(d)
    if n<8:return dict(n=n,mean=float(d.mean()) if n else np.nan,p=1.,estimable=False)
    rng=np.random.default_rng(20260914+n)
    if n not in _null_cache:
        if n<=16:
            _null_cache[n]=((np.arange(2**n,dtype=np.uint32)[:,None]>>np.arange(n))&1).astype(np.int8)*2-1
        else:_null_cache[n]=rng.choice(np.array([-1,1],np.int8),size=(200000,n))
    signs=_null_cache[n]
    hits=int((np.abs(signs@d/n)>=abs(d.mean())-1e-12).sum())
    p=hits/len(signs) if n<=16 else (hits+1)/(len(signs)+1)
    if n not in _boot_cache:_boot_cache[n]=np.random.default_rng(20261014+n).multinomial(n,np.full(n,1/n),size=20000)
    ci=np.quantile(_boot_cache[n]@d/n,[.025,.975])
    return dict(n=n,mean=float(d.mean()),p=p,estimable=True,CI95_low=ci[0],CI95_high=ci[1],
                n_positive=int((d>0).sum()),n_negative=int((d<0).sum()))

def bh(p):
    p=np.asarray(p,float);order=np.argsort(p);v=np.empty(len(p))
    v[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return v
