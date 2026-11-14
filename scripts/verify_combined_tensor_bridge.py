#!/usr/bin/env python3
"""Reevaluate foundation rank-certificate values using the flow tensor engine."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import numpy as np
from chiral4form.registry import Registry
from chiral4form.forms import dense, hodge_compact, layout
from chiral4form.tensors import TensorBudget


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--certificate', type=Path, default=Path('reports/fresh-rank81.json'))
    ap.add_argument('--output', type=Path, default=Path('reports/tensor-bridge.json'))
    args=ap.parse_args()
    inp=Path('foundation/science/rank81/input.json')
    verifier=Path('foundation/science/rank81/verify.py')
    cert=json.loads(args.certificate.read_text()); inputs=json.loads(inp.read_text())
    if cert['input_sha256']!=hashlib.sha256(inp.read_bytes()).hexdigest() or cert['implementation_sha256']!=hashlib.sha256(verifier.read_bytes()).hexdigest():
        raise ValueError('rank certificate inputs or implementation have changed')
    spec=importlib.util.spec_from_file_location('independent_foundation_tensor', verifier)
    standalone=importlib.util.module_from_spec(spec);spec.loader.exec_module(standalone)
    reg=Registry.load('data/combined/registry.json')
    convention=inputs['convention']
    if convention['coordinate_tuples'] != [list(layout()['basis'][i]) for i in layout()['electric']]:
        raise ValueError('electric coordinate conventions differ')
    records=[]
    budget=TensorBudget(2*1024**3, 100_000_000_000, True)
    for point in cert['points']:
        p=point['prime'];electric=point['coordinates'];raw=np.zeros(252,dtype=np.int64)
        raw[layout()['electric']]=electric
        form=dense((raw+hodge_compact(raw,p))%p,p)
        if not np.array_equal(form,standalone.dense_form(electric,p)):
            raise ValueError('Hodge or dense-tensor conventions disagree')
        cache={}
        for i,item in enumerate(inputs['graphs']):
            inv=reg.items[item['id']]
            graph=[[0]*item['degree'] for _ in range(item['degree'])]
            for a,b,n in item['graph']['edges']:
                graph[a][b]=graph[b][a]=n
            if inv.graph != tuple(map(tuple,graph)):
                raise ValueError(f'graph definition mismatch: {item["id"]}')
            value=reg.evaluate(item['id'],form,p,cache=cache,budget=budget)
            if value!=point['values'][i]:
                raise ValueError(f'independent engines disagree: {item["id"]}')
            print(f'bridge prime={p} invariant={i+1}/81 {item["id"]}',flush=True)
        records.append({'prime':p,'independent_invariants_compared':len(inputs['graphs']), 'dense_tensor_and_graph_definitions_equal':True})
    report={'schema':1,'passed':True,'scope':'fresh flow-engine values compared to separately recomputed foundation certificate; same graph, metric, Hodge and electric coordinates','rank_certificate_sha256':hashlib.sha256(args.certificate.read_bytes()).hexdigest(),'points':records}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
