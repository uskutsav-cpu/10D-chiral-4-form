#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
from chiral4form.finite_field import matrix_rank
from chiral4form.fitting import fields_from_json
from chiral4form.normalization_audit import load_records
from chiral4form.registry import Registry
DEGREES=(4,6,8,10,12)

def req(cat): return frozenset(f[2:] for f in cat['factors'] if f.startswith('S:'))
def constant_value(poly):
    zero=(0,)*poly.n
    if any(p!=zero for p in poly.terms): return None
    return int(poly.terms.get(zero,0))
def layer_rows(ids,fields,catalogue,registry,degree,selected):
    idx={n:i for i,n in enumerate(ids)}
    lower=[idx[n] for d in DEGREES if d<degree for n in registry.degree_bases[d]]
    cols=[idx[n] for n in registry.degree_bases[degree]]
    names=[];rows=[]
    for cat in catalogue:
        if int(cat['leading_degree'])!=degree or not req(cat)<=selected: continue
        field=fields[cat['id']]
        if any(field[j] for j in lower): continue
        vec=[];ok=True
        for j in cols:
            v=constant_value(field[j])
            if v is None: ok=False;break
            vec.append(v)
        if ok: names.append(cat['id']);rows.append(vec)
    return names,rows

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--run',type=Path,default=Path('runs/degree12-generalized-selected-verify'))
    p.add_argument('--completion',type=Path,default=Path('verification/degree12/generalized_completion_theorem.json'))
    p.add_argument('--output',type=Path,default=Path('verification/degree12/generalized_removal_minimality.json'))
    args=p.parse_args()
    comp=json.loads(args.completion.read_text());selected=tuple(comp['selected_extras'])
    target={int(d):int(n) for d,n in comp['homogeneous_dimensions'].items()}
    registry=Registry.from_json(json.loads((args.run/'registry.json').read_text()))
    records=load_records(sorted(args.run.glob('fields_prime*.json')));primes=tuple(sorted(records))
    results={};print('===== STEP 3E REMOVAL AUDIT =====')
    for removed in selected:
        remaining=frozenset(x for x in selected if x!=removed);per={};all_drop=True;count_ob=False;count_deg=None
        for pmod in primes:
            ids,fields=fields_from_json(records[pmod]);catalogue=records[pmod]['generator_catalogue'];degrows={}
            for d in DEGREES:
                names,rows=layer_rows(ids,fields,catalogue,registry,d,remaining)
                rank=matrix_rank(rows,pmod,ncols=target[d]);degrows[d]={'admissible_constant_generator_count':len(rows),'rank':rank,'target':target[d]}
                if len(rows)<target[d]: count_ob=True;count_deg=count_deg or d
            drops=[d for d in DEGREES if degrows[d]['rank']<target[d]]
            if not drops: all_drop=False
            per[pmod]={'degrees':{str(d):degrows[d] for d in DEGREES},'drop_degrees':drops}
        results[removed]={'all_fresh_primes_lose_strong_completion':all_drop,'exact_count_obstruction':count_ob,'exact_count_first_degree':count_deg,'per_prime':{str(p):v for p,v in per.items()}}
        tag='EXACT-COUNT NONREDUNDANT' if count_ob else ('FRESH-PRIME NONREDUNDANT' if all_drop else 'REMOVABLE UNDER STRONG CRITERION')
        print(removed,'=>',tag)
    exact=[n for n,r in results.items() if r['exact_count_obstruction']]
    nonred=[n for n,r in results.items() if r['all_fresh_primes_lose_strong_completion']]
    removable=[n for n,r in results.items() if not r['all_fresh_primes_lose_strong_completion']]
    out={'schema':2,'status':'generalized_selected_set_removal_audit_complete','selected_extras':list(selected),'fresh_verification_primes':list(primes),
         'removal_results':results,'exact_count_nonredundant_extras':exact,'fresh_prime_nonredundant_extras':nonred,
         'removable_under_strong_constant_translation_criterion':removable,
         'all_selected_extras_fresh_prime_nonredundant':len(nonred)==len(selected),
         'claim_scope':'removal audit for the strong layerwise global constant-translation completion criterion within the declared primitive graph-extra catalogue',
         'not_a_claim_of':['global minimality over arbitrary scalar linear combinations','minimality for every nonlinear dynamic completion','all-orders minimality']}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n')
    print('PASS: STEP 3E REMOVAL AUDIT COMPLETE');print('selected:',len(selected));print('exact count nonredundant:',len(exact));print('fresh-prime nonredundant:',len(nonred));print('removable:',len(removable));print('WROTE:',args.output)
if __name__=='__main__':main()
