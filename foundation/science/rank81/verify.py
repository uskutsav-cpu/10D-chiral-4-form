#!/usr/bin/env python3
"""Recompute a rank81 theorem from graphs, coordinates and exact points only.

No stored values/Jacobians, candidate search, Hilbert counts, literature maps,
repository imports, files of intermediate arrays, or resume logic are used.
NumPy/opt_einsum provide bounded exact contractions and path selection.
"""
import argparse
from collections import Counter, OrderedDict
import hashlib
import itertools
import json
from pathlib import Path
import time
import numpy as np
import opt_einsum as oe

TUPLES=list(itertools.combinations(range(10),5))
POSITIONS={t:k for k,t in enumerate(TUPLES)}
COORDINATES=[k for k,t in enumerate(TUPLES) if 0 in t]
DIM=10


def parity(sequence):
    return (-1) ** sum(a > b for k, a in enumerate(sequence) for b in sequence[k + 1:])

def metric_product(indices):
    return -1 if 0 in indices else 1

def star_image(indices):
    complement = tuple(i for i in range(10) if i not in indices)
    return complement, parity(complement + indices) * metric_product(indices)

def compact_form(coordinates, prime):
    assert len(coordinates) == 126
    result = [0] * 252
    for a, k in zip(coordinates, COORDINATES):
        result[k] = int(a) % prime
        other, sign = star_image(TUPLES[k])
        result[POSITIONS[other]] = sign * int(a) % prime
    return result

def slot_records(item):
    n = item["degree"]
    slots = [[] for _ in range(n)]
    raised = [[] for _ in range(n)]
    for i, j, multiplicity in item["graph"]["edges"]:
        for k in range(1, multiplicity + 1):
            edge = (i, j, k)
            slots[i].append(edge)
            slots[j].append(edge)
            raised[i].append(False)
            raised[j].append(True)
    return slots, raised

def dense_form(coordinates, prime):
    tensor = np.zeros((10,) * 5, dtype=np.float64)
    for indices, value in zip(TUPLES, compact_form(coordinates, prime)):
        for ordered in itertools.permutations(indices):
            tensor[ordered] = parity(ordered) * value % prime
    return tensor

def terms(item):
    slots, raised = slot_records(item)
    labels = sorted({edge for group in slots for edge in group})
    mapping = {edge: k for k, edge in enumerate(labels)}
    return [[mapping[edge] for edge in group] for group in slots], raised

def plan(item):
    labels, raised = terms(item)
    equation = ",".join("".join(oe.get_symbol(k) for k in group) for group in labels) + "->"
    path, info = oe.contract_path(equation, *[(10,)*5 for _ in labels], shapes=True, optimize="dp")
    working = [list(group) for group in labels]
    largest_output = 10**5
    largest_dot_terms = 1
    work = 0
    for pair in path:
        assert len(pair) == 2
        i, j = pair
        left, right = working[i], working[j]
        shared = set(left) & set(right)
        output = [k for k in left if k not in shared] + [k for k in right if k not in shared]
        work += 10**len(set(left) | set(right))
        largest_output = max(largest_output, 10**len(output))
        largest_dot_terms = max(largest_dot_terms, 10**len(shared))
        for k in sorted(pair, reverse=True):
            working.pop(k)
        working.append(output)
    assert working == [[]]
    return {"id": item["id"], "path": [list(p) for p in path],
            "largest_output_elements": largest_output,
            "largest_dot_terms": largest_dot_terms,
            "forward_dual_multiply_terms": 3*work,
            "float64_exact_bound_at_largest_prime": largest_dot_terms * (32749-1)**2,
            "three_output_arrays_bytes": 3 * largest_output * 8}

def modular_dot(left, right, prime, counters):
    assert left.shape[1] == right.shape[0]
    bound = left.shape[1] * (prime - 1)**2
    assert bound < 2**53, ("STOP: exact float64 bound unavailable", bound)
    counters["dot_products"] += 1
    counters["largest_dot_bound"] = max(counters["largest_dot_bound"], bound)
    product = left @ right
    assert np.all(np.isfinite(product))
    assert np.all(product >= 0) and np.all(product <= bound)
    # This does not claim approximate arithmetic: every exact partial sum is
    # an integer no larger than the guard above, and is representable exactly.
    return np.remainder(product, prime)

class PointEngine:
    def __init__(self, coordinates, prime, direction=None, cache_limit=256*1024**2):
        self.prime = prime
        self.dual = direction is not None
        self.form = dense_form(coordinates, prime)
        self.tangent = dense_form(direction, prime) if self.dual else None
        self.inputs = {}
        self.ids = {}
        self.cache = OrderedDict()
        self.cache_bytes = 0
        self.limit = cache_limit
        self.stats = Counter()

    def identify(self, descriptor):
        if descriptor not in self.ids:
            self.ids[descriptor] = len(self.ids)
        return self.ids[descriptor]

    def initial(self, raised, permutation=None):
        permutation = tuple(range(5)) if permutation is None else tuple(permutation)
        key = ("input", tuple(raised), permutation)
        identity = self.identify(key)
        if identity not in self.inputs:
            arrays = [self.form.transpose(permutation)]
            if self.dual:
                arrays.append(self.tangent.transpose(permutation))
            signs = np.array([-1] + [1]*9, dtype=np.float64)
            for axis, flag in enumerate(raised):
                if flag:
                    shape = [1]*5
                    shape[axis] = 10
                    arrays = [np.remainder(array * signs.reshape(shape), self.prime) for array in arrays]
            self.inputs[identity] = tuple(arrays)
        return identity, self.inputs[identity]

    def contract_pair(self, left, right):
        la, ida, aa = left
        lb, idb, bb = right
        common = [label for label in la if label in lb]
        free_a = [label for label in la if label not in common]
        free_b = [label for label in lb if label not in common]
        order_a = tuple(la.index(label) for label in free_a + common)
        order_b = tuple(lb.index(label) for label in common + free_b)
        result_labels = free_a + free_b
        result_id = self.identify(("pair", ida, idb, order_a, order_b, len(common)))
        if result_id in self.cache:
            self.stats["cache_hits"] += 1
            self.cache.move_to_end(result_id)
            return result_labels, result_id, self.cache[result_id]
        shape_a = (10**len(free_a), 10**len(common))
        shape_b = (10**len(common), 10**len(free_b))
        a = aa[0].transpose(order_a).reshape(shape_a)
        b = bb[0].transpose(order_b).reshape(shape_b)
        value = modular_dot(a, b, self.prime, self.stats)
        arrays = [value]
        if self.dual:
            da = aa[1].transpose(order_a).reshape(shape_a)
            db = bb[1].transpose(order_b).reshape(shape_b)
            arrays.append(np.remainder(modular_dot(da, b, self.prime, self.stats)
                                       + modular_dot(a, db, self.prime, self.stats), self.prime))
        arrays = tuple(array.reshape((10,)*len(result_labels)) for array in arrays)
        size = sum(array.nbytes for array in arrays)
        if size <= self.limit:
            while self.cache and self.cache_bytes + size > self.limit:
                _, evicted = self.cache.popitem(last=False)
                self.cache_bytes -= sum(array.nbytes for array in evicted)
            self.cache[result_id] = arrays
            self.cache_bytes += size
            self.stats["max_cached_bytes"] = max(self.stats["max_cached_bytes"], self.cache_bytes)
        return result_labels, result_id, arrays

    def evaluate(self, item, selected_plan=None, slot_permutations=None):
        labels, raised = terms(item)
        working = []
        slot_permutations = slot_permutations or {}
        for v, (group, variance) in enumerate(zip(labels, raised)):
            identity, arrays = self.initial(variance, slot_permutations.get(v))
            working.append((group, identity, arrays))
        selected_plan = selected_plan or plan(item)
        for i, j in selected_plan["path"]:
            combined = self.contract_pair(working[i], working[j])
            for k in sorted([i, j], reverse=True):
                working.pop(k)
            working.append(combined)
        assert len(working) == 1 and working[0][0] == []
        return tuple(int(array) for array in working[0][2])

def sparse_coordinate_pullback():
    indices, signs = [], []
    for k in COORDINATES:
        pair, pair_sign = star_image(TUPLES[k])
        column_indices, column_signs = [], []
        for sorted_indices, sign in [(TUPLES[k],1),(pair,pair_sign)]:
            for permuted in itertools.permutations(sorted_indices):
                flat=0
                for i in permuted:
                    flat=10*flat+i
                column_indices.append(flat)
                column_signs.append(sign*parity(permuted))
        indices.append(column_indices)
        signs.append(column_signs)
    return np.array(indices,dtype=np.int64),np.array(signs,dtype=np.int64)

def reverse_gradient(engine,item,selected_plan=None):
    assert not engine.dual
    selected_plan=selected_plan or plan(item)
    labels,raised=terms(item)
    nodes=[]
    active=[]
    for group,variance in zip(labels,raised):
        identity,arrays=engine.initial(variance)
        nodes.append({"operand":(group,identity,arrays),"variance":variance})
        active.append(len(nodes)-1)
    for i,j in selected_plan["path"]:
        left,right=active[i],active[j]
        combined=engine.contract_pair(nodes[left]["operand"],nodes[right]["operand"])
        nodes.append({"operand":combined,"left":left,"right":right})
        for k in sorted([i,j],reverse=True):active.pop(k)
        active.append(len(nodes)-1)
    assert len(active)==1
    p=engine.prime
    adjoints={active[0]:np.array(1.0)}
    for index in range(len(nodes)-1,len(raised)-1,-1):
        node=nodes[index]
        la,_,aa=nodes[node["left"]]["operand"]
        lb,_,bb=nodes[node["right"]]["operand"]
        common=[label for label in la if label in lb]
        free_a=[label for label in la if label not in common]
        free_b=[label for label in lb if label not in common]
        oa=[la.index(label) for label in free_a+common]
        ob=[lb.index(label) for label in common+free_b]
        a=aa[0].transpose(oa).reshape(10**len(free_a),10**len(common))
        b=bb[0].transpose(ob).reshape(10**len(common),10**len(free_b))
        c=adjoints.pop(index).reshape(10**len(free_a),10**len(free_b))
        da=modular_dot(c,b.T,p,engine.stats).reshape((10,)*len(la)).transpose(np.argsort(oa))
        db=modular_dot(a.T,c,p,engine.stats).reshape((10,)*len(lb)).transpose(np.argsort(ob))
        # A contraction tree has each occurrence exactly once, even when the
        # forward cache happened to reuse an equal array at multiple nodes.
        assert node["left"] not in adjoints and node["right"] not in adjoints
        adjoints[node["left"]]=da
        adjoints[node["right"]]=db
    gradient=np.zeros((10,)*5,dtype=np.float64)
    metric=np.array([-1]+[1]*9,dtype=np.float64)
    for leaf,variance in enumerate(raised):
        g=adjoints[leaf]
        for axis,flag in enumerate(variance):
            if flag:
                shape=[1]*5
                shape[axis]=10
                g=np.remainder(g*metric.reshape(shape),p)
        gradient=np.remainder(gradient+g,p)
    compact=(gradient.ravel()[PULLBACK_INDICES].astype(np.int64)*PULLBACK_SIGNS).sum(axis=1)%p
    scalar=int(nodes[-1]["operand"][2][0])
    return scalar,[int(x) for x in compact]

def bareiss(matrix):
    a = [list(map(int, row)) for row in matrix]
    n = len(a)
    assert n and all(len(row) == n for row in a)
    prior = 1
    sign = 1
    swaps = 0
    for k in range(n - 1):
        if not a[k][k]:
            pivot_row = next((i for i in range(k + 1, n) if a[i][k]), None)
            if pivot_row is None:
                return 0, swaps
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
            swaps += 1
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                quotient, remainder = divmod(a[i][j] * pivot - a[i][k] * a[k][j], prior)
                assert remainder == 0, (k, i, j, "inexact Bareiss division")
                a[i][j] = quotient
            a[i][k] = 0
        prior = pivot
    return sign * a[-1][-1], swaps

PULLBACK_INDICES,PULLBACK_SIGNS=sparse_coordinate_pullback()

def elimination(matrix,p):
    a=[list(map(lambda x:int(x)%p,row)) for row in matrix]
    columns=[];det=1
    for c in range(len(a[0])):
        r=len(columns)
        pivot=next((i for i in range(r,len(a)) if a[i][c]),None)
        if pivot is None:continue
        if pivot!=r:a[r],a[pivot]=a[pivot],a[r];det=-det
        v=a[r][c];det=det*v%p;inverse=pow(v,-1,p)
        for i in range(r+1,len(a)):
            scale=a[i][c]*inverse%p
            if scale:a[i]=[(x-scale*y)%p for x,y in zip(a[i],a[r])]
        columns.append(c)
        if len(columns)==len(a):break
    return columns,det if len(columns)==len(a) else 0


def orbit_rows(coordinates,p):
    lookup=dict(zip([t for t in TUPLES if 0 in t],coordinates))
    def component(indices):
        if len(set(indices))<5:return 0
        t=tuple(sorted(indices));s=parity(indices)
        if 0 in t:return s*lookup[t]%p
        complement=tuple(i for i in range(10) if i not in t)
        return -s*parity(t+complement)*lookup[complement]%p
    rows=[]
    for i,j in itertools.combinations(range(10),2):
        reverse=1 if i==0 else -1
        full=[]
        for t in TUPLES:
            total=0
            for k,index in enumerate(t):
                v=list(t)
                if index==i:v[k]=j;total+=component(v)
                elif index==j:v[k]=i;total+=reverse*component(v)
            full.append(total%p)
        row=full[:126]
        assert full==compact_form(row,p)
        rows.append(row)
    return rows


def validate_input(data):
    assert data['schema']==1
    assert data['convention']=={'metric':[-1]+[1]*9,'coordinate_tuples':[list(t) for t in TUPLES if 0 in t],'hodge':'output-first epsilon; epsilon_0123456789=+1; F=sum A_I(e_I+star e_I)'}
    assert len(data['graphs'])==len({g['id'] for g in data['graphs']})==81
    expected_degrees=[4,6,6]+[8]*6+[10]*12+[12]*60
    assert [g['degree'] for g in data['graphs']]==expected_degrees
    seen_graphs=set()
    for item in data['graphs']:
        n=item['degree'];edges=item['graph']['edges']
        assert n==item['graph']['n'] and n>=4
        assert edges==sorted(edges) and len({(i,j) for i,j,m in edges})==len(edges)
        valences=[0]*n
        adjacency=[set() for _ in range(n)]
        for i,j,m in edges:
            assert all(type(x)is int for x in (i,j,m)) and 0<=i<j<n and 0<m<=4
            valences[i]+=m;valences[j]+=m
            adjacency[i].add(j);adjacency[j].add(i)
        assert valences==[5]*n
        graph_key=(n,tuple(tuple(edge) for edge in edges))
        assert graph_key not in seen_graphs
        seen_graphs.add(graph_key)
        seen={0};stack=[0]
        while stack:
            v=stack.pop()
            for w in adjacency[v]:
                if w not in seen:
                    seen.add(w);stack.append(w)
        assert len(seen)==n
    assert data['points']
    for point in data['points']:
        p=point['prime'];assert type(p)is int and 2<p<65536
        assert all(p%d for d in range(2,int(p**.5)+1))
        assert len(point['coordinates'])==126 and all(type(a)is int and 0<=a<p for a in point['coordinates'])
    for t in TUPLES:
        u,s=star_image(t);v,r=star_image(u)
        assert t==v and s*r==1


def main():
    if not __debug__:raise RuntimeError('Assertions must be enabled; do not use python -O')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();data=json.loads(args.input.read_text());validate_input(data)
    report={'schema':1,'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),'implementation_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'mode':'full fresh values, Jacobians, Bareiss minors and orbit action from minimal inputs','points':[]}
    for point in data['points']:
        start=time.monotonic();p=point['prime'];a=point['coordinates'];engine=PointEngine(a,p)
        values=[];rows=[]
        for i,item in enumerate(data['graphs']):
            value,row=reverse_gradient(engine,item)
            assert sum(x*y for x,y in zip(a,row))%p==item['degree']*value%p
            values.append(value);rows.append(row)
            print(f'point prime={p} row={i+1}/81 id={item["id"]}',flush=True)
        columns,det=elimination(rows,p);assert len(columns)==81 and det!=0
        integer,_=bareiss([[row[c] for c in columns] for row in rows]);assert integer%p==det
        orbit=orbit_rows(a,p);oc,od=elimination(orbit,p);assert len(oc)==45 and od!=0
        orbit_integer,_=bareiss([[row[c] for c in oc] for row in orbit]);assert orbit_integer%p==od
        assert all(sum(x*y for x,y in zip(row,tangent))%p==0 for row in rows for tangent in orbit)
        report['points'].append({**point,'values':values,'jacobian':rows,'pivot_columns':columns,'integer_minor_determinant':str(integer),'determinant_mod_p':det,'rank':81,'orbit_rows':orbit,'orbit_pivot_columns':oc,'orbit_integer_minor_determinant':str(orbit_integer),'orbit_determinant_mod_p':od,'orbit_rank':45,'seconds':time.monotonic()-start,'arithmetic_checks':dict(engine.stats)})
    report['passed']=True
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':True,'points':[{k:r[k] for k in ('prime','rank','determinant_mod_p','orbit_rank','orbit_determinant_mod_p','seconds')} for r in report['points']]},indent=2))
if __name__=='__main__':main()

