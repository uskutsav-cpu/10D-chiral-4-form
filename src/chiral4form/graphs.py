"""Exact graph evaluation and reverse differentiation, with bounded plans.

This is an independent evaluator, not an invariant-enumeration algorithm.
Graph slot order matches the predecessor's lexicographic edge convention.
"""
from __future__ import annotations
from functools import lru_cache
import re
import numpy as np
from .finite_field import integer
from .forms import constrained_gradient
from .tensors import as_field_array,binary_contract,DEFAULT_BUDGET


def parse_graph(label,valence=5):
    match=re.fullmatch(r'n(\d+)\[(.*)\]',label)
    if not match:
        raise ValueError('malformed graph label')
    n=int(match[1]);a=[[0]*n for _ in range(n)]
    if n<2 or n>32:
        raise ValueError('supported graph size is 2..32')
    for part in match[2].split(','):
        pattern=r'(\d+)-(\d+)\^(\d+)' if n>10 else r'(\d)(\d)\^(\d+)'
        m=re.fullmatch(pattern,part)
        if not m:
            raise ValueError('ambiguous or malformed edge record')
        i,j,w=map(int,m.groups())
        if not 0<=i<j<n or not 1<=w<=valence or a[i][j]:
            raise ValueError('invalid or repeated graph edge')
        a[i][j]=a[j][i]=w
    return validate_graph(a,valence)


def validate_graph(matrix,valence=5):
    a=tuple(tuple(integer(x) for x in row) for row in matrix);n=len(a)
    if not 2<=n<=32 or any(len(row)!=n for row in a):
        raise ValueError('invalid square graph')
    for i in range(n):
        if a[i][i] or sum(a[i])!=valence:
            raise ValueError('graph must have zero diagonal and uniform valence')
        for j in range(n):
            if a[i][j]!=a[j][i] or not 0<=a[i][j]<=valence:
                raise ValueError('invalid symmetric multiplicities')
    return a


def slots(matrix):
    n=len(matrix);labels=[[] for _ in range(n)];raised=[[] for _ in range(n)]
    edge=0
    for i in range(n):
        for j in range(i+1,n):
            for _ in range(matrix[i][j]):
                labels[i].append(edge)
                raised[j].append(len(labels[j]))
                labels[j].append(edge);edge+=1
    return tuple(map(tuple,labels)),tuple(map(tuple,raised))

@lru_cache(maxsize=256)
def plan(matrix,d):
    """Dynamic-programming tree minimizing peak tensor order, then work.

    Exact subset optimization is capped at 12 leaves; larger graphs use a
    deterministic greedy tree. This is an execution plan, not a proof claim.
    """
    labels,_=slots(matrix);n=len(labels)
    bounds={1<<i:tuple(labels[i]) for i in range(n)}
    trees={1<<i:i for i in range(n)}
    score={1<<i:(d**len(labels[i]),0) for i in range(n)}
    if n<=12:
        for mask in range(1,1<<n):
            if mask in trees:
                continue
            first=mask&-mask;rest=mask^first
            bounds[mask]=tuple(sorted(set(bounds[first])^set(bounds[rest])))
            best=None;besttree=None
            sub=(mask-1)&mask
            while sub:
                other=mask^sub
                if other and sub&first:
                    la,lb=bounds[sub],bounds[other]
                    output=d**len(bounds[mask]);work=d**len(set(la)|set(lb))
                    s=(max(score[sub][0],score[other][0],output),score[sub][1]+score[other][1]+work)
                    if best is None or s<best:
                        best=s;besttree=(trees[sub],trees[other])
                sub=(sub-1)&mask
            score[mask]=best;trees[mask]=besttree
        return trees[(1<<n)-1],score[(1<<n)-1]
    nodes=[(tuple(l),i,d**len(l),0) for i,l in enumerate(labels)]
    while len(nodes)>1:
        candidates=[]
        for i in range(len(nodes)):
            for j in range(i+1,len(nodes)):
                la,_,pa,wa=nodes[i];lb,_,pb,wb=nodes[j]
                out=tuple(sorted(set(la)^set(lb)));work=d**len(set(la)|set(lb))
                candidates.append(((max(pa,pb,d**len(out)),wa+wb+work),i,j,out))
        s,i,j,out=min(candidates)
        node=(out,(nodes[i][1],nodes[j][1]),s[0],s[1])
        nodes=[x for q,x in enumerate(nodes) if q not in (i,j)]+[node]
    return nodes[0][1],(nodes[0][2],nodes[0][3])



@lru_cache(maxsize=256)
def allocation_profile(matrix, d, gradient=False):
    """Conservative array-allocation estimate following the retained tape.

    This budgets arrays, including transpose/BLAS/reduction temporaries, not
    operating-system RSS, interpreter overhead, or vendor BLAS thread buffers.
    Every primitive separately checks its own allocation/work budget.
    """
    labels, _ = slots(matrix)
    tree, (peak_tensor, work) = plan(matrix, d)
    nodes = [{'labels':set(x), 'size':d**len(x), 'children':None} for x in labels]
    live = sum(x['size'] for x in nodes)
    peak = live
    def visit(t):
        nonlocal live, peak
        if isinstance(t, int):
            return t
        i, j = visit(t[0]), visit(t[1])
        out = nodes[i]['labels'] ^ nodes[j]['labels']
        size = d**len(out)
        peak = max(peak, live + 4*(nodes[i]['size'] + nodes[j]['size'] + size))
        nodes.append({'labels':out, 'size':size, 'children':(i,j)})
        live += size
        return len(nodes)-1
    root = visit(tree)
    if gradient:
        adj = {root:nodes[root]['size']}
        live += adj[root]
        for parent in range(root, len(matrix)-1, -1):
            i, j = nodes[parent]['children']
            for child, sibling in ((i,j),(j,i)):
                peak = max(peak, live + 4*(adj[parent] + nodes[sibling]['size'] + nodes[child]['size']))
                adj[child] = nodes[child]['size']
                live += adj[child]
        # Sum leaf adjoints and project/antisymmetrize the form derivative.
        peak = max(peak, live + 8*d**len(labels[0]))
    return {'reserved_bytes':8*peak, 'multiply_adds':work,
            'largest_tensor_elements':peak_tensor, 'gradient':bool(gradient)}


def evaluate_graph(matrix,form,p,*,gradient=False,budget=DEFAULT_BUDGET):
    form=as_field_array(form,p);k=form.ndim;d=form.shape[0]
    if form.shape!=(d,)*k:
        raise ValueError('form axes must have equal dimension')
    matrix=validate_graph(matrix,k);labels,raised=slots(matrix)
    tree,(peak,work)=plan(matrix,d)
    # Follow the actual retained tape; a peak tensor alone is not a memory bound.
    reserve=allocation_profile(matrix,d,gradient)['reserved_bytes']
    if reserve>budget.max_bytes or work>budget.max_multiply_adds:
        raise MemoryError(f'graph plan exceeds budget: reserved_bytes={reserve}, multiply_adds={work}')
    sign=np.ones(d,dtype=np.int64);sign[0]=-1
    nodes=[]
    def signed(x,axes):
        x=x.copy()
        for axis in axes:
            shape=[1]*k;shape[axis]=d
            x=x*sign.reshape(shape)%p
        return x
    for i in range(len(matrix)):
        nodes.append({'labels':labels[i],'value':signed(form,raised[i]),'children':None})
    def forward(t):
        if isinstance(t,int):
            return t
        i,j=forward(t[0]),forward(t[1]);a,b=nodes[i],nodes[j]
        out=tuple(sorted(set(a['labels'])^set(b['labels'])))
        value=binary_contract(a['value'],a['labels'],b['value'],b['labels'],out,p,budget)
        nodes.append({'labels':out,'value':value,'children':(i,j)})
        return len(nodes)-1
    root=forward(tree)
    if nodes[root]['labels']:
        raise ValueError('graph is not fully contracted')
    scalar=int(nodes[root]['value'])%p
    if not gradient:
        return scalar
    adj={root:np.array(1,dtype=np.int64)}
    for q in range(root,len(matrix)-1,-1):
        i,j=nodes[q]['children'];a,b=nodes[i],nodes[j]
        adj[i]=binary_contract(adj[q],nodes[q]['labels'],b['value'],b['labels'],a['labels'],p,budget)
        adj[j]=binary_contract(adj[q],nodes[q]['labels'],a['value'],a['labels'],b['labels'],p,budget)
    ambient=np.zeros_like(form)
    for i in range(len(matrix)):
        ambient=(ambient+signed(adj[i],raised[i]))%p
    return scalar,constrained_gradient(ambient,p,d,k)
