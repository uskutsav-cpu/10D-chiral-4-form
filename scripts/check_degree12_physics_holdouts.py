#!/usr/bin/env python3
"""Fresh, disjoint finite-field physical checks of all 36 degree-12 sources.

A passed holdout is empirical validation of the map, not an identity proof.
"""
from pathlib import Path
import argparse
import hashlib
import json
from fractions import Fraction
from chiral4form import degree12_exact_targets as targets
from chiral4form.registry import Registry


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=Path('reports/degree12-physical-holdouts.json'))
    ap.add_argument('--primes',default='50021,50023')
    ap.add_argument('--seeds',default='2026111301,2026111302')
    args=ap.parse_args();reg=Registry.load('data/combined/registry.json')
    path=Path('verification/all_orders/degree12_exact/exact_source_targets.json')
    cert=json.loads(path.read_text());rows={int(r['source_row']):r for r in cert['rows']}
    descriptors=targets._all_source_descriptors();reports=[]
    for p in map(int,args.primes.split(',')):
        if p in cert['direct_physical_primes'] or p in cert['theorem_primes_verification_only']:
            raise ValueError('holdout prime overlaps source reconstruction/verification primes')
        substitution=targets._degree10_orbit_substitution_mod_p(reg,p)
        for seed in map(int,args.seeds.split(',')):
            form=targets._small_selfdual(seed,p)
            actual=targets._all_target_values_for_sample(reg,form,p,substitution,descriptors)
            basis=reg.basis_values(12,form,p,budget=targets.EXACT_D12_BUDGET)
            residuals=[]
            for i in range(36):
                coefficients=list(map(Fraction,rows[i]['coordinates']))
                expected=sum(q.numerator*pow(q.denominator,-1,p)*b for q,b in zip(coefficients,basis))%p
                residuals.append((actual[i]-expected)%p)
            record={'prime':p,'seed':seed,'source_residuals':residuals,'all_36_zero':not any(residuals)}
            reports.append(record)
            args.output.parent.mkdir(parents=True,exist_ok=True)
            payload={'schema':1,'status':'fresh_physical_holdouts_passed' if all(r['all_36_zero'] for r in reports) else 'physical_source_map_refuted',
                     'scope':'finite-field physical evaluations at fresh primes and seeds; not a characteristic-zero identity proof',
                     'target_certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'reports':reports}
            args.output.write_text(json.dumps(payload,indent=2)+'\n')
            print(f'prime={p} seed={seed} zero sources={36-sum(bool(x) for x in residuals)}/36',flush=True)
            if any(residuals):
                raise ValueError(f'physical degree-12 map REFUTED at p={p}, seed={seed}; failed rows={[i for i,x in enumerate(residuals) if x]}')

if __name__=='__main__':main()
