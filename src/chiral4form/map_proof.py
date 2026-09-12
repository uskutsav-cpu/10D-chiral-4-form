"""Bounded-residual certificate for the degree-eight physics/basis map.

Unlike an ordinary CRT fit, this uses a proved bound on *integer residuals at
fixed small integer tensors*, plus an injective evaluation matrix. The latter
uses the EXTERNAL characteristic-zero Hilbert upper bounds 1,2,7. The result
is conditional on those bounds, the declared HLS formula, and conventions.
See docs/DEGREE8_MAP_CERTIFICATE.md for the complete argument.
"""
from __future__ import annotations
from math import prod
from pathlib import Path
import numpy as np
from .forms import layout,dense,hodge_compact
from .registry import Registry
from .stress import evaluate_generators,residue
from .degree8_orbit import analytic_fields
from .finite_field import matrix_rank,_check_prime_field_modulus
from .provenance import atomic_json,semantic_hash
from .checkpoints import CheckpointStore
from .pipeline import engine_fingerprint

PRIMES=(30011,30013,30029,30047,30059,30071,30089)
UPPER={4:1,6:2,8:7}
DENOMINATOR=240**2


def bound_report():
    # One graph contraction has 5*n/2 independently summed labels.
    G={n:10**(5*n//2) for n in (4,6,8)}
    D4=4*10**5  # n vertices, with the differentiated vertex's 5 slots fixed.
    M=10**4
    B44=10**4*D4**2
    targets={
        'tr1_d8':60*G[8],
        'tr2_d8_c4_squared':2*10**2*M*25*B44+10*(2*G[4])**2,
        'tr2_d8_c6':2*10**2*M*4*G[6],
        'tr3_d8_c4':3*10**3*M**2*2*G[4],
        'tr4_d8':10**4*M**4,
        'tr1_squared_d8':400*G[4]**2,
        'tr1_tr2_d8':20*G[4]*10**2*M**2,
        'tr2_squared_d8':(10**2*M**2)**2}
    # Includes every candidate row; all coefficients have denominator dividing 3.
    candidate=400*G[8]
    bound=DENOMINATOR*(max(targets.values())+candidate)
    return {'tensor_component_bound':1,'graph_value_bounds':{str(k):v for k,v in G.items()},
            'quartic_gradient_component_bound':D4,'target_bounds':targets,
            'candidate_row_abs_bound':candidate,'common_denominator':DENOMINATOR,
            'integer_residual_abs_bound':bound,
            'derivation':'gradient antisymmetrization and Hodge projection are averages; 240 clears one gradient, 240^2 clears its quadratic stress term'}


def fixed_form(electric,p):
    if len(electric)!=126 or any(type(x) is not int or abs(x)>1 for x in electric):
        raise ValueError('certificate tensors require 126 signed integer electric components bounded by one')
    raw=np.zeros(252,dtype=np.int64);raw[layout()['electric']]=electric
    return dense((raw+hodge_compact(raw,p))%p,p)


def expected_keys():
    keys=set()
    ids,fields=analytic_fields()
    degrees=[4,6,6]+[8]*7
    for g,field in fields.items():
        for i,f in enumerate(field):
            for exponent in f.terms:
                m=tuple(name for name,power in zip(ids,exponent) for _ in range(power))
                keys.add((g,degrees[i],tuple(sorted(m))))
    # Terms permitted by degree counting but identically zero by self-duality.
    keys.update({('tr2',6,('I4_1',)),('tr2',8,('I6_1',)),('tr2',8,('I6_2',))})
    return sorted(keys)


def candidate_value(g,d,monomial,basis_values,reg,p):
    ids,fields=analytic_fields();m=tuple(monomial);answer=0
    for name,value in zip(reg.degree_bases[d],basis_values[str(d)]):
        f=fields[g][ids.index(name)]
        exponent=tuple(m.count(i) for i in ids)
        answer=(answer+residue(f.terms.get(exponent,0),p)*value)%p
    return answer


def sample(reg,electric,p):
    form=fixed_form(electric,p);cache={}
    targets,_=evaluate_generators(reg,form,p,8,(),cache=cache)
    bases={str(d):reg.basis_values(d,form,p,cache=cache) for d in UPPER}
    keys=expected_keys()
    if set(targets)-set(keys):
        raise AssertionError('unexpected coefficient monomial outside degree-eight envelope')
    rows=[]
    for g,d,m in keys:
        actual=targets.get((g,d,m),0);expected=candidate_value(g,d,m,bases,reg,p)
        rows.append({'generator':g,'degree':d,'monomial':list(m),'actual':actual,'candidate':expected,
                     'cleared_residual':DENOMINATOR*(actual-expected)%p})
    return {'prime':p,'electric':electric,'basis_values':bases,'residuals':rows}


def verify(record, *, reevaluate_tensors=False):
    """Recheck the explicit arithmetic certificate, not re-evaluate tensors."""
    try:
        bounds=bound_report();primes=record['primes'];points=record['electric_points']
        if record['bounds']!=bounds or record['external_upper_bounds']!={str(k):v for k,v in UPPER.items()}:
            return False
        if record.get('modulus_product') != prod(primes):
            return False
        if len(set(primes))!=len(primes) or prod(primes)<=2*bounds['integer_residual_abs_bound']:
            return False
        for p in primes:
            _check_prime_field_modulus(p)
            if DENOMINATOR%p==0 or p>65521:
                return False
        if len({tuple(x) for x in points})!=len(points):
            return False
        for point in points:
            if len(point)!=126 or any(type(x) is not int or abs(x)>1 for x in point):
                return False
        reg=Registry.from_json(record['registry'])
        if any(len(reg.degree_bases[k])!=n for k,n in UPPER.items()):
            return False
        samples=record['samples']
        if len(samples)!=len(primes)*len(points):
            return False
        bykey={}
        keys=expected_keys()
        for s in samples:
            key=(s['prime'],tuple(s['electric']))
            if key in bykey or s['prime'] not in primes or s['electric'] not in points:
                return False
            p=s['prime'];bykey[key]=s
            if any(len(s['basis_values'][str(d)])!=n for d,n in UPPER.items()):
                return False
            if [(r['generator'],r['degree'],tuple(r['monomial'])) for r in s['residuals']]!=keys:
                return False
            for r in s['residuals']:
                expected=candidate_value(r['generator'],r['degree'],r['monomial'],s['basis_values'],reg,p)
                if r['candidate']!=expected or r['cleared_residual']!=0 or (DENOMINATOR*(r['actual']-expected))%p:
                    return False
        # One injective evaluation matrix over QQ follows from this nonzero modular minor.
        first=primes[0]
        for d,n in UPPER.items():
            A=[bykey[(first,tuple(x))]['basis_values'][str(d)] for x in points]
            if matrix_rank(A,first)!=n:
                return False
        if reevaluate_tensors:
            for stored in samples:
                fresh=sample(reg,stored['electric'],stored['prime'])
                if fresh != stored:
                    return False
        clean={k:v for k,v in record.items() if k!='certificate_hash'}
        return record.get('certificate_hash')==semantic_hash(clean)
    except (KeyError,ValueError,TypeError,ZeroDivisionError):
        return False


def run(output,root=None,progress=print,prime_index=None):
    root=Path.cwd() if root is None else Path(root)
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    reg=Registry.load(root/'data/fixtures/low_degree_registry.json')
    rng=np.random.default_rng(202609121911)
    points=[rng.integers(-1,2,126).tolist() for _ in range(9)]
    bounds=bound_report()
    if prod(PRIMES)<=2*bounds['integer_residual_abs_bound']:
        raise ValueError('insufficient modulus for an integer residual proof')
    cache=CheckpointStore(output/'checkpoints',{'engine':engine_fingerprint(),'registry':reg.fingerprint,
                                               'bounds':bounds,'points':points})
    if prime_index is not None and (type(prime_index) is not int or not 0<=prime_index<len(PRIMES)):
        raise ValueError('proof prime-index must be 0..6')
    selected=PRIMES if prime_index is None else (PRIMES[prime_index],)
    samples=[]
    for p in selected:
        for i,electric in enumerate(points):
            s,cached=cache.run('bounded_map_sample',{'prime':p,'point':i},lambda:sample(reg,electric,p))
            if any(x['cleared_residual'] for x in s['residuals']):
                raise AssertionError(f'degree-eight map refuted at prime {p}, point {i}')
            samples.append(s)
            progress(f'map proof prime {p}: {i+1}/{len(points)} {"cached" if cached else "computed"}',flush=True)
    if prime_index is not None:
        samples=[]
        for p in PRIMES:
            for i in range(len(points)):
                saved=cache.load('bounded_map_sample',{'prime':p,'point':i})
                if saved is None:
                    partial={'status':'partial_proof_samples_not_a_completed_certificate',
                        'engine_fingerprint':engine_fingerprint(),'prime_index_completed':prime_index,
                        'certificate_hash':None}
                    atomic_json(output/f'partial_prime_{prime_index}.json',partial)
                    return partial
                samples.append(saved)
    record={'schema':1,'status':'bounded_integer_residual_certificate_conditional_on_external_Hilbert_bounds_and_HLS_conventions',
            'engine_fingerprint':engine_fingerprint(),'registry':reg.to_json(),
            'external_upper_bounds':{str(k):v for k,v in UPPER.items()},
            'external_source':'arXiv:2509.14350v2; homogeneous invariant dimensions in characteristic zero',
            'bounds':bounds,'primes':list(PRIMES),'modulus_product':prod(PRIMES),
            'electric_points':points,'samples':samples,
            'claim':'All declared degree-eight polynomial stress-generator basis identities over QQ.',
            'limitations':['external invariant-space upper bounds are inputs, not derived here',
                'HLS stress normalization and constrained-derivative conventions require expert review',
                'not an all-orders, nonanalytic, causality, or priority theorem']}
    record['certificate_hash']=semantic_hash(record)
    if not verify(record):
        raise AssertionError('bounded-residual certificate failed its independent arithmetic verifier')
    atomic_json(output/'map_certificate.json',record)
    atomic_json(output/'summary.json',{'status':record['status'],'certificate_hash':record['certificate_hash'],
        'samples':len(samples),'prime_count':len(PRIMES),'identity_rows_per_sample':len(expected_keys()),
        'integer_residual_bound':bounds['integer_residual_abs_bound'],'modulus_product':prod(PRIMES),
        'arithmetic_verification':True})
    return record
