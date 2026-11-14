#!/usr/bin/env python3
"""Self-contained bounded-CRT proof of the 36 degree-12 source identities.

Generates fresh 72-point physical evaluations; needs no ignored predecessor
run files and assumes no integrality of the full stress-map coefficients.
Derivative denominators are cleared explicitly. Checkpoints are hash-bound.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
import hashlib
import json
from math import prod
from pathlib import Path
import time
from chiral4form import degree12_exact_targets as t
from chiral4form.finite_field import matrix_rank, _check_prime_field_modulus
from chiral4form.provenance import atomic_json

SEEDS = tuple(202611130000 + j for j in range(72))
GRADIENT_DENOMINATOR = 240
TAU_DENOMINATOR = GRADIENT_DENOMINATOR**2


def parameters():
    reg = t.Registry.load(t.REGISTRY_PATH)
    traces = t._trace_bound_series(reg)
    substitution = t._degree10_orbit_substitution_qq(reg)
    rows = []
    for meta in t._all_source_descriptors():
        d, b, proof = t._target_bound(reg, meta, traces=traces, substitution=substitution)
        power = sum(int(token[2:]) for token in meta['generator'].split('*'))
        # A stress entry has at most two projected derivatives, each with
        # denominator dividing 2*5!=240. Multiplication of k entries therefore
        # clears with (240^2)^k. Identity/free terms divide the same clearer.
        clearer = TAU_DENOMINATOR**power
        rows.append({**meta, 'denominator': d*clearer, 'integer_bound': b*clearer,
                     'orbit_denominator': d, 'stress_denominator': clearer,
                     'target_absolute_bound_before_stress_clearer': b, 'bound_proof':proof})
    threshold = 2*max(t.GRAPH12_SAMPLE_BOUND, max(r['integer_bound'] for r in rows))
    primes=[];M=1
    for p in t.FULL_PRIMES:
        if any(r['denominator']%p==0 for r in rows):continue
        primes.append(p);M*=p
    for p in t._extra_prime_stream(set(primes)):
        if M > threshold:break
        if any(r['denominator']%p==0 for r in rows):continue
        primes.append(p);M*=p
    return reg, rows, tuple(primes)


def namespace(reg, rows, primes):
    files = sorted(Path('src/chiral4form').glob('*.py')) + [Path(__file__)]
    h=hashlib.sha256()
    for path in files:
        h.update(path.name.encode());h.update(path.read_bytes())
    body={'schema':1,'algorithm':'fresh_d12_direct_physics_explicit_derivative_denominators_v1',
          'engine_sha256':h.hexdigest(),'registry_sha256':hashlib.sha256(Path(t.REGISTRY_PATH).read_bytes()).hexdigest(),
          'target_sha256':hashlib.sha256(Path(t.CERT_PATH).read_bytes()).hexdigest(),
          'seeds':SEEDS,'rows':rows,'primes':primes,'basis_bound':t.GRAPH12_SAMPLE_BOUND}
    return hashlib.sha256(json.dumps(body,sort_keys=True).encode()).hexdigest(),body


def worker(output, p, digest):
    reg, rows, _ = parameters()
    path=Path(output)/f'prime{p}.json'
    data={'schema':1,'namespace_sha256':digest,'prime':p,'seeds':list(SEEDS),'samples':[]}
    if path.exists():
        data=json.loads(path.read_text())
        if data.get('namespace_sha256')!=digest or data.get('prime')!=p or data.get('seeds')!=list(SEEDS):
            raise ValueError(f'incompatible checkpoint: {path}; preserve it and use another output directory')
    samples=data['samples']
    if len(samples)>72:raise ValueError('too many checkpoint samples')
    for i,sample in enumerate(samples):
        if sample['seed']!=SEEDS[i] or len(sample['basis'])!=72 or len(sample['targets'])!=36:
            raise ValueError('invalid checkpoint sample structure')
    sub=t._degree10_orbit_substitution_mod_p(reg,p)
    for i in range(len(samples),72):
        start=time.monotonic();form=t._small_selfdual(SEEDS[i],p)
        values=t._all_target_values_for_sample(reg,form,p,sub,rows)
        basis=reg.basis_values(12,form,p,budget=t.EXACT_D12_BUDGET)
        samples.append({'seed':SEEDS[i],'basis':list(map(int,basis)),
                        'targets':[int(values[j]) for j in range(36)],'seconds':time.monotonic()-start})
        data['status']='complete' if len(samples)==72 else 'running'
        atomic_json(path,data)
        print(f'degree12 p={p} point={i+1}/72 seconds={samples[-1]["seconds"]:.2f}',flush=True)
    return p


def centered(residues, primes, bound):
    a,M=t.crt(residues,list(primes))
    if M<=2*bound:raise ValueError('insufficient CRT modulus')
    a=a-M if a>M//2 else a
    if abs(a)>bound:raise ValueError('lift exceeds proved bound')
    return a


def verify(output, *, write=True):
    reg, rows, primes=parameters();digest,config=namespace(reg,rows,primes)
    data={p:json.loads((Path(output)/f'prime{p}.json').read_text()) for p in primes}
    for p,record in data.items():
        _check_prime_field_modulus(p)
        if record['namespace_sha256']!=digest or record['seeds']!=list(SEEDS) or len(record['samples'])!=72:
            raise ValueError('incomplete/incompatible physical checkpoint')
    basis=[[centered([data[p]['samples'][i]['basis'][j] for p in primes],primes,t.GRAPH12_SAMPLE_BOUND) for j in range(72)] for i in range(72)]
    rank=matrix_rank(basis,primes[0],ncols=72)
    if rank!=72:raise ValueError('fresh evaluation points do not separate the 72-dimensional invariant space')
    candidate=json.loads(Path(t.CERT_PATH).read_text())
    coordinates={r['source_row']:list(map(Fraction,r['coordinates'])) for r in candidate['rows']}
    failures=[];target_numerators=[]
    for j,meta in enumerate(rows):
        D=meta['denominator'];B=meta['integer_bound'];values=[]
        for i in range(72):
            actual=centered([data[p]['samples'][i]['targets'][j]*(D%p)%p for p in primes],primes,B)
            expected=sum(q*b for q,b in zip(coordinates[j],basis[i]))*D
            if expected.denominator!=1 or actual!=expected:failures.append((j,i))
            values.append(str(actual))
        target_numerators.append(values)
    if failures:raise ValueError(f'physical source identities refuted: {failures[:12]}')
    char=json.loads(Path('foundation/science/character/independent-d5.json').read_text())
    upper=next(int(r['invariants']) for r in char['degrees'] if r['degree']==12)
    if upper!=72:raise ValueError('invariant-space dimension input changed')
    hashes={f'prime{p}.json':hashlib.sha256((Path(output)/f'prime{p}.json').read_bytes()).hexdigest() for p in primes}
    report={'schema':1,'status':'bounded_CRT_degree12_source_map_verified','namespace_sha256':digest,
            'configuration':config,'checkpoints_sha256':hashes,'prime_count':len(primes),'point_count':72,
            'source_count':36,'modulus_product':str(prod(primes)),'injective_basis_rank':rank,
            'independent_character_upper_bound':upper,'all_2592_exact_source_residuals_zero':True,
            'exact_basis_matrix':basis,'cleared_target_numerators':target_numerators,
            'integrality_assumption_for_full_map':False,
            'proof':'explicit denominator clearing + unique bounded integer lifts at fixed tensors + injective evaluations on the homogeneous invariant space',
            'assumptions':['declared Lorentzian HLS stress convention','projected derivative denominator divides 240','independent D5 character derivation applies to this real representation','exactness of the bounded tensor arithmetic'],
            'not_claimed':['all-orders orbit equality','nonanalytic ModMax reachability','physical admissibility','novelty']}
    if write:atomic_json(Path(output)/'map_certificate.json',report)
    print(json.dumps({'status':report['status'],'primes':len(primes),'points':72,'sources':36}),flush=True)
    return report


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=Path('reports/degree12-map'))
    ap.add_argument('--workers',type=int,default=2)
    ap.add_argument('--verify-only',action='store_true')
    args=ap.parse_args()
    if not 1<=args.workers<=4:raise ValueError('workers must be between 1 and 4')
    reg,rows,primes=parameters();digest,config=namespace(reg,rows,primes)
    args.output.mkdir(parents=True,exist_ok=True)
    if not args.verify_only:
        config_path=args.output/'configuration.json'
        if config_path.exists() and json.loads(config_path.read_text())['namespace_sha256']!=digest:
            raise ValueError('output has a different namespace; preserve it and choose another directory')
        atomic_json(config_path,{'namespace_sha256':digest,'configuration':config})
        print(f'fresh degree12 proof: {len(primes)} primes x 72 points x 36 sources; {args.workers} workers',flush=True)
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures=[pool.submit(worker,str(args.output),p,digest) for p in primes]
            for future in as_completed(futures):future.result()
    verify(args.output)

if __name__=='__main__':main()
