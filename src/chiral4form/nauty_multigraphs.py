"""Streaming nauty/genbg enumeration of loopless 5-regular contraction multigraphs.

A loopless 5-regular multigraph on n five-form vertices has 5n/2 edges.
Encode each multiedge occurrence as a degree-2 vertex in a bipartite incidence
graph. Thus genbg with left degrees exactly 5 and right degrees exactly 2
generates the desired multigraphs up to colour-preserving isomorphism.
"""
from __future__ import annotations
import shutil,subprocess
from pathlib import Path

def require_genbg():
    exe=shutil.which("genbg")
    if exe is None:
        raise RuntimeError("nauty genbg is required. On macOS: brew install nauty")
    return exe

def decode_graph6(line):
    data=line.strip().encode("ascii") if isinstance(line,str) else line.strip()
    if not data: raise ValueError("empty graph6 record")
    if data[0]==126: raise ValueError("extended graph6 orders are not needed here")
    n=data[0]-63
    if not 0<=n<=62: raise ValueError("invalid graph6 order")
    bits=[]
    for c in data[1:]:
        v=c-63
        if not 0<=v<64: raise ValueError("invalid graph6 character")
        bits.extend((v>>k)&1 for k in (5,4,3,2,1,0))
    A=[[0]*n for _ in range(n)];q=0
    for j in range(1,n):
        for i in range(j):
            if q>=len(bits): raise ValueError("truncated graph6 record")
            if bits[q]: A[i][j]=A[j][i]=1
            q+=1
    return A

def collapse_incidence_graph(A,left_n):
    n=len(A);right=range(left_n,n)
    M=[[0]*left_n for _ in range(left_n)]
    if any(sum(A[i])!=5 for i in range(left_n)):
        raise ValueError("left incidence degrees are not all five")
    for r in right:
        nbr=[i for i in range(left_n) if A[r][i]]
        if len(nbr)!=2: raise ValueError("edge-vertex degree is not two")
        if any(A[r][s] for s in right): raise ValueError("right-right incidence edge")
        a,b=sorted(nbr)
        if a==b: raise ValueError("loop encountered")
        M[a][b]+=1;M[b][a]+=1
    if any(sum(row)!=5 for row in M): raise ValueError("collapsed graph is not 5-regular")
    return tuple(tuple(x for x in row) for row in M)

def graph_label(M):
    parts=[]
    for i in range(len(M)):
        for j in range(i+1,len(M)):
            if M[i][j]: parts.append(f"{i}-{j}^{M[i][j]}")
    return f"n{len(M)}["+",".join(parts)+"]"

def genbg_command(degree):
    if degree%2: raise ValueError("five-regular graph needs even vertex count")
    edges=5*degree//2
    return [require_genbg(),"-q","-l","-d5:2","-D5:2",str(degree),str(edges),f"{2*edges}:{2*edges}"]

def stream_multigraphs(degree):
    cmd=genbg_command(degree)
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,bufsize=1)
    try:
        assert proc.stdout is not None
        for line in proc.stdout:
            if line.strip():
                yield collapse_incidence_graph(decode_graph6(line),degree)
    finally:
        if proc.poll() is None: proc.terminate()
        try: proc.wait(timeout=5)
        except subprocess.TimeoutExpired: proc.kill()
