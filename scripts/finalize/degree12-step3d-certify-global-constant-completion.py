#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
from chiral4form.finite_field import determinant,matrix_rank
from chiral4form.fitting import fields_from_json
from chiral4form.normalization_audit import load_records,audit_model_records
from chiral4form.registry import Registry

DEGREES=(4,6,8,10,12)

def constant_value(poly):
    zero=(0,)*poly.n
    if any(powers!=zero for powers in poly.terms): return None
    return int(poly.terms.get(zero,0))

def greedy_basis(names,rows,p,ncols):
    chosen_names=[];chosen_rows=[];rank=0
    for name,row in zip(names,rows):
        nr=matrix_rank([*chosen_rows,row],p,ncols=ncols)
        if nr>rank:
            chosen_names.append(name);chosen_rows.append(row);rank=nr
            if rank==ncols: break
    return chosen_names,chosen_rows

def constant_layer_rows(ids,fields,catalogue,registry,degree):
    index={name:i for i,name in enumerate(ids)}
    lower=[index[name] for d in DEGREES if d<degree for name in registry.degree_bases[d]]
    cols=[index[name] for name in registry.degree_bases[degree]]
    cat_by_name={row['id']:row for row in catalogue}
    names=[];rows=[]
    for name in sorted(fields):
        cat=cat_by_name[name]
        if int(cat['leading_degree'])!=degree: continue
        field=fields[name]
        if any(field[j] for j in lower): continue
        vec=[];ok=True
        for j in cols:
            v=constant_value(field[j])
            if v is None: ok=False;break
            vec.append(v)
        if ok:
            names.append(name);rows.append(vec)
    return names,rows

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--run',type=Path,default=Path('runs/degree12-generalized-selected-verify'))
    p.add_argument('--discovery',type=Path,default=Path('verification/degree12/generalized_discovery.json'))
    p.add_argument('--output',type=Path,default=Path('verification/degree12/generalized_completion_theorem.json'))
    args=p.parse_args()
    discovery=json.loads(args.discovery.read_text());selected=list(discovery['selected_extras'])
    registry=Registry.from_json(json.loads((args.run/'registry.json').read_text()))
    records=load_records(sorted(args.run.glob('fields_prime*.json')));audit=audit_model_records(records)
    primes=tuple(sorted(records))
    if len(primes)<2: raise SystemExit('need at least two fresh verification primes')
    target={d:len(registry.degree_bases[d]) for d in DEGREES}
    common=None;per_prime={}
    print('===== STEP 3D GLOBAL CONSTANT-TRANSLATION CERTIFICATE =====')
    for pmod in primes:
        ids,fields=fields_from_json(records[pmod]);catalogue=records[pmod].get('generator_catalogue')
        if not catalogue: raise SystemExit(f'{pmod}: generator catalogue missing')
        dd={}
        for d in DEGREES:
            names,rows=constant_layer_rows(ids,fields,catalogue,registry,d)
            rank=matrix_rank(rows,pmod,ncols=target[d]);dd[d]=(names,rows,rank)
            print(pmod,f'd={d}: constant generators={len(names)} rank={rank}/{target[d]}')
            if rank!=target[d]:
                raise SystemExit(f'{pmod}: strong global translation criterion fails at degree {d}: {rank}/{target[d]}')
        if common is None:
            common={}
            for d in DEGREES:
                names,rows,_=dd[d];chosen,_=greedy_basis(names,rows,pmod,target[d])
                if len(chosen)!=target[d]: raise SystemExit(f'failed to choose d={d} basis')
                common[d]=chosen
        minors={}
        for d in DEGREES:
            names,rows,_=dd[d];by=dict(zip(names,rows));square=[by[name] for name in common[d]]
            det=determinant(square,pmod)
            if det==0: raise SystemExit(f'{pmod}: common d={d} translation minor vanished')
            minors[str(d)]=int(det)
        per_prime[pmod]={'constant_translation_ranks':{str(d):target[d] for d in DEGREES},'minor_residues':minors}
    theorem={'schema':2,'status':'characteristic_zero_global_generalized_completion_through_degree12',
      'normalization_audit':audit,'selected_extras':selected,'selected_count':len(selected),
      'candidate_class':discovery['candidate_class'],'fresh_verification_primes':list(primes),
      'homogeneous_dimensions':{str(d):target[d] for d in DEGREES},
      'common_global_translation_bases':{str(d):common[d] for d in DEGREES},
      'per_prime_certificates':{str(p):v for p,v in per_prime.items()},
      'characteristic_zero_rank_argument':'Full-size constant-translation minors are nonzero modulo good primes; the exact rational layer maps therefore have full rank.',
      'global_constructive_algorithm':{'degree_order':[4,6,8,10,12],
        'procedure':'At layer d solve the full-rank constant translation system. Chosen controls are identically zero below d; higher spillover is corrected at later layers.'},
      'full_polynomial_interaction_space_reachable_through_degree12':True,'global_generalized_completion':True,
      'scope':'finite polynomial interactions through field degree 12; arbitrary signed piecewise controls; selected primitive graph invariants admitted as genuine independent S factors in f(tau,S)',
      'not_claimed':['global cardinality minimality over arbitrary scalar combinations','all-orders completion','nonanalytic completion']}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(theorem,indent=2)+'\n')
    print('PASS: STEP 3D GLOBAL GENERALIZED COMPLETION');print('selected extras:',len(selected));print('homogeneous dimensions:',target);print('WROTE:',args.output)

if __name__=='__main__':main()
