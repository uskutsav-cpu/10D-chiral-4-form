"""Explicit graph/product registry; never substitute 81 coordinates for a ring."""
from __future__ import annotations
from dataclasses import dataclass
import json
from pathlib import Path
import numpy as np
from .exact import rational
from .graphs import parse_graph,evaluate_graph,validate_graph
from .provenance import semantic_hash
from .tensors import DEFAULT_BUDGET

@dataclass(frozen=True)
class Invariant:
    id: str
    degree: int
    graph: tuple | None = None
    factors: tuple[str,...] = ()

class Registry:
    def __init__(self,items,degree_bases,metadata=None):
        self.items={x.id:x for x in items}
        if len(self.items)!=len(items):
            raise ValueError('duplicate invariant IDs')
        self.degree_bases={int(d):tuple(ids) for d,ids in degree_bases.items()}
        self.metadata=metadata or {}
        for item in items:
            if item.degree<2 or item.degree%2:
                raise ValueError('invalid invariant field degree')
            if (item.graph is None)==(not item.factors):
                raise ValueError('each invariant must be exactly one graph or product')
            if item.graph is not None:
                validate_graph(item.graph)
                if len(item.graph)!=item.degree:
                    raise ValueError('graph degree mismatch')
            else:
                if any(x not in self.items for x in item.factors):
                    raise ValueError('unknown product factor')
                if any(self.items[x].degree>=item.degree for x in item.factors):
                    raise ValueError('cyclic or non-lower product factor')
                if sum(self.items[x].degree for x in item.factors)!=item.degree:
                    raise ValueError('product degree mismatch')
        for degree,ids in self.degree_bases.items():
            if len(set(ids))!=len(ids) or any(i not in self.items or self.items[i].degree!=degree for i in ids):
                raise ValueError('invalid ordered homogeneous basis')

    def to_json(self):
        return {'schema':2,'metadata':self.metadata,'items':[{'id':x.id,'degree':x.degree,
            **({'graph_matrix':[list(row) for row in x.graph]} if x.graph is not None else {'factors':list(x.factors)})}
            for x in self.items.values()], 'degree_bases':{str(d):list(v) for d,v in self.degree_bases.items()}}

    @property
    def fingerprint(self):
        return semantic_hash(self.to_json())

    @classmethod
    def from_json(cls,payload):
        if payload.get('schema')!=2:
            raise ValueError('registry schema must be 2')
        items=[]
        for r in payload['items']:
            graph=parse_graph(r['graph']) if 'graph' in r else (validate_graph(r['graph_matrix']) if 'graph_matrix' in r else None)
            items.append(Invariant(r['id'],r['degree'],graph,tuple(r.get('factors',()))))
        return cls(items,payload['degree_bases'],payload.get('metadata'))

    @classmethod
    def load(cls,path):
        return cls.from_json(json.loads(Path(path).read_text()))

    def evaluate(self,item_id,form,p,*,gradient=False,cache=None,budget=DEFAULT_BUDGET):
        cache={} if cache is None else cache
        key=(item_id,gradient)
        if key in cache:
            return cache[key]
        item=self.items[item_id]
        if item.graph is not None:
            answer=evaluate_graph(item.graph,form,p,gradient=gradient,budget=budget)
        else:
            factors=[self.evaluate(i,form,p,gradient=gradient,cache=cache,budget=budget) for i in item.factors]
            values=[x[0] for x in factors] if gradient else factors
            value=1
            for v in values:
                value=value*v%p
            if gradient:
                derivative=np.zeros_like(form)
                for j,(_,g) in enumerate(factors):
                    c=1
                    for h,v in enumerate(values):
                        if h!=j:
                            c=c*v%p
                    derivative=(derivative+c*g)%p
                answer=(value,derivative)
            else:
                answer=value
        cache[key]=answer
        return answer

    def basis_values(self,degree,form,p,*,cache=None,budget=DEFAULT_BUDGET):
        cache={} if cache is None else cache
        return [self.evaluate(i,form,p,cache=cache,budget=budget) for i in self.degree_bases[degree]]


def import_predecessor_registry(root,max_degree=12):
    """Translate actual committed JSON formulas without importing upstream code."""
    if max_degree not in (4,6,8,10,12):
        raise ValueError('supported atlas cutoffs: 4,6,8,10,12')
    root=Path(root)
    lower=json.loads((root/'results/10d_order8.json').read_text())
    items=[Invariant(x['id'],int(x['order']),parse_graph(x['graph'])) for x in lower['generators'] if int(x['order'])<=max_degree]
    bases={d:[i.id for i in items if i.degree==d] for d in (4,6,8) if d<=max_degree}
    def product(name,d,factors):
        items.append(Invariant(name,d,None,tuple(factors)));return name
    if max_degree>=8:
        bases[8].append(product('I4_1^2',8,['I4_1','I4_1']))
    if max_degree>=10:
        raw=json.loads((root/'results/10d_order10.json').read_text())
        ten=[Invariant(x['id'],10,parse_graph(x['graph'])) for x in raw['generators']]
        if len(ten)!=12:
            raise ValueError('expected twelve explicit degree-10 graph formulas')
        items.extend(ten);bases[10]=[i.id for i in ten]
        bases[10]+=[product('I4_1*I6_1',10,['I4_1','I6_1']),product('I4_1*I6_2',10,['I4_1','I6_2'])]
    if max_degree>=12:
        raw=json.loads((root/'results/10d_order12.json').read_text())
        twelve=[Invariant(x['id'],12,parse_graph(x['graph'])) for x in raw['generators']]
        if len(twelve)!=62:
            raise ValueError('expected sixty-two explicit degree-12 graph formulas')
        items.extend(twelve)
        product('I4_1^3',12,['I4_1']*3)
        product('I6_1^2',12,['I6_1']*2)
        product('I6_1*I6_2',12,['I6_1','I6_2'])
        product('I6_2^2',12,['I6_2']*2)
        for i in range(1,7):
            product(f'I4_1*I8_{i}',12,['I4_1',f'I8_{i}'])
        bases[12]=[x['id'] for x in raw['degree12_basis']]
    return Registry(items,bases,{'source':'predecessor JSON graph formulas',
        'identity_scope':'explicit candidate bases; completeness is an external mathematical input'})
