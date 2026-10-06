"""FASTA reader and parent-context helpers shared by A25, A26, A27 and A28: fasta() reads a FASTA into an ordered dict with I replaced by L; context() returns the residues around a peptide occurrence."""
import hashlib
def fasta(path,audit):
    digest=hashlib.sha256();header=None;parts=[];n=0
    with open(path,'rb') as f:
        for raw in f:
            digest.update(raw);line=raw.strip()
            if line.startswith(b'>'):
                if header is not None:
                    n+=1;yield header,b''.join(parts).decode('ascii').upper()
                header=line[1:].decode('ascii');parts=[]
            elif line:parts.append(line)
    if header is not None:
        n+=1;yield header,b''.join(parts).decode('ascii').upper()
    audit.update(sha256=digest.hexdigest(),entries=n)
def context(s,start,end):
    n=start==0 or (s[start-1] in 'KR' and s[start]!='P')
    c=end==len(s) or (s[end-1] in 'KR' and s[end]!='P')
    missed=sum(s[i] in 'KR' and s[i+1]!='P' for i in range(start,end-1))
    return (s[start-1] if start else 'PROTEIN_END',s[end] if end<len(s) else 'PROTEIN_END',n,c,missed)
