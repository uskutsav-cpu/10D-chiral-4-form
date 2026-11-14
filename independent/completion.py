#!/usr/bin/env python3
"""Exact homogeneous leading-block completion for the selected 78 extras.

This proves an invertible, constant leading block at each field degree. The
field-degree triangularity then supplies a global stage-by-stage construction.
It does not prove minimum nonlinear-control cardinality.
"""
import argparse
import json
from pathlib import Path
import sympy as s
from verify import ROOT, read, require, polynomial


def verify_completion(root):
    reg=read(root,'data/combined/registry.json')
    items={r['id']:r for r in reg['items']}
    graphs=[r['id'] for r in reg['items'] if 'graph_matrix' in r]
    symbols=s.symbols(f'g0:{len(graphs)}');functions=dict(zip(graphs,symbols))
    def expression(name):
        if name not in functions:
            functions[name]=s.prod(expression(f) for f in items[name]['factors'])
        return functions[name]
    for name in items:expression(name)
    model=read(root,'verification/degree10/rational_model.json')
    x=s.symbols('x0:24');zero=dict.fromkeys(x,0)
    for name in ('tr2','tr3','tr4','tr5'):
        values=[polynomial(q,x).xreplace(zero) for q in model['fields'][name]]
        functions[name]=s.expand(sum(c*expression(i) for c,i in zip(values,model['coordinate_ids'])))
    targets=read(root,'verification/all_orders/degree12_exact/exact_source_targets.json')
    tr6=next(r for r in targets['rows'] if r['source_row']==35)
    functions['tr6']=s.expand(sum(s.Rational(c)*expression(i) for c,i in zip(tr6['coordinates'],reg['degree_bases']['12'])))
    theorem=read(root,'verification/degree12/generalized_completion_theorem.json')
    extras=theorem['selected_extras']
    require(len(extras)==78 and len(set(extras))==78,'selected extra catalogue differs')
    require(all(i in graphs for i in extras),'extra is not a primitive graph')
    blocks={}
    for degree in (4,6,8,10,12):
        names=theorem['common_global_translation_bases'][str(degree)]
        ids=reg['degree_bases'][str(degree)]
        require(len(names)==len(ids),'leading basis is not square')
        basis_terms=[s.Poly(expression(i),*symbols).monoms() for i in ids]
        require(all(len(ms)==1 for ms in basis_terms),'homogeneous atlas element is not a primitive/product monomial')
        monomials=[ms[0] for ms in basis_terms]
        columns=[]
        for name in names:
            factors=name.split('*');q=s.S.One;weight=0
            for factor in factors:
                if factor.startswith('S:'):
                    item=factor[2:];require(item in extras,'unselected extra used')
                    q*=expression(item);weight+=items[item]['degree']
                else:
                    require(factor in ('tr2','tr3','tr4','tr5','tr6'),'unsupported stress generator')
                    q*=functions[factor];weight+=2*int(factor[2:])
            require(weight==degree,'leading generator degree mismatch')
            poly=s.Poly(s.expand(q),*symbols)
            require(set(poly.monoms())<=set(monomials),'generator is outside complete homogeneous basis')
            columns.append([poly.coeff_monomial(m) for m in monomials])
        matrix=s.Matrix.hstack(*map(s.Matrix,columns));det=matrix.det()
        require(det!=0,f'degree-{degree} completion block singular')
        blocks[str(degree)]={'dimension':len(ids),'generators':names,'determinant_QQ':str(det),
                            'matrix_QQ':[[str(q) for q in row] for row in matrix.tolist()]}
    return {'schema':1,'passed':True,'selected_extra_count':78,'ambient_dimension':96,
            'exact_constant_leading_blocks':blocks,
            'construction':'set coefficients in increasing field degree using the inverse constant block; higher-leading-degree controls preserve all lower coefficients',
            'structural_hypothesis':'HLS analytic derivative-free field-degree filtration; selected controls exclude tr1 and strictly lower coupling weight',
            'physics_inputs':'free stress trace coordinates through degree12; source identities independently certified by the full map runs',
            'minimality':'not claimed for nonlinear reachability or arbitrary additional scalars',
            'scope':'exact rational leading blocks plus degree-triangular global construction in the declared generalized polynomial class'}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--output',type=Path,default=Path('reports/generalized-completion-exact.json'));args=ap.parse_args()
    result=verify_completion(args.root);args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':True,'extras':78,'dimensions':{d:r['dimension'] for d,r in result['exact_constant_leading_blocks'].items()}}))

if __name__=='__main__':main()
