"""Lift matching sampled polynomial models with an UNUSED holdout prime.

The result is a candidate rational model with explicit coefficient bounds,
not an automatically proved characteristic-zero physics identity.
"""
from .exact import reconstruct_with_holdout
from .polynomial import Poly
from .fitting import fields_from_json,fields_to_json


def lift_models(records,bound):
    if len(records)<3:
        raise ValueError('at least two fit-prime models plus one held-out prime required')
    parsed=[fields_from_json(r) for r in records];ids=parsed[0][0]
    names=set(parsed[0][1]);primes=[]
    for names_ids,fields in parsed:
        if names_ids!=ids or set(fields)!=names:
            raise ValueError('coordinate or generator registries differ across primes')
        p=next(iter(fields.values()))[0].p
        if p is None or p in primes:
            raise ValueError('distinct prime-field models required')
        primes.append(p)
    result={};n=len(ids);count=0
    for name in sorted(names):
        coords=[]
        for i in range(n):
            monomials=sorted(set().union(*(f[name][i].terms for _,f in parsed)))
            terms={}
            for m in monomials:
                residues=[int(fields[name][i].terms.get(m,0)) for _,fields in parsed]
                lifted=reconstruct_with_holdout(residues[:-1],primes[:-1],residues[-1],primes[-1],bound)
                terms[m]=lifted['value'];count+=1
            coords.append(Poly(n,terms))
        result[name]=tuple(coords)
    output=fields_to_json(ids,result,'bounded_rational_candidate_with_unused_prime_validation_not_symbolic_physics_proof')
    output.update({'fit_primes':primes[:-1],'holdout_prime':primes[-1],'coefficient_height_bound':bound,
                   'coefficients_reconstructed':count,'generator_catalogue':records[0].get('generator_catalogue',[])})
    return output
