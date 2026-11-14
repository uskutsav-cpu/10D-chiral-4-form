#!/usr/bin/env python3
"""Replay degree-12 proof arithmetic without importing either tensor engine.

This checks the saved finite-field evidence and the declared bound arithmetic.
It does not independently derive stress bounds or regenerate tensor samples.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from math import prod
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
TARGET = 'verification/all_orders/degree12_exact/exact_source_targets.json'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def lift(residues, primes, bound):
    """Incremental CRT with unique centered bounded reconstruction."""
    require(len(residues) == len(primes) and bool(primes), 'CRT shape mismatch')
    value, modulus = 0, 1
    for residue, prime in zip(residues, primes):
        require(type(residue) is int and 0 <= residue < prime, 'noncanonical residue')
        value += modulus * (((residue-value)*pow(modulus, -1, prime)) % prime)
        modulus *= prime
    require(modulus > 2*bound, 'insufficient modulus for unique bounded lift')
    if value > modulus//2:
        value -= modulus
    require(abs(value) <= bound, 'reconstructed integer exceeds declared bound')
    return value


def rank_mod(matrix, prime):
    a = [[int(v) % prime for v in row] for row in matrix]
    rank = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(rank, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inverse = pow(a[rank][j], -1, prime)
        a[rank] = [(v*inverse) % prime for v in a[rank]]
        for i in range(rank+1, len(a)):
            factor = a[i][j]
            if factor:
                a[i] = [(v-factor*w) % prime for v,w in zip(a[i],a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def verify(root, directory):
    configuration = json.loads((directory/'configuration.json').read_text())
    config = configuration['configuration']
    digest = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
    require(digest == configuration['namespace_sha256'], 'configuration digest mismatch')
    require(config['registry_sha256'] == sha(root/'data/combined/registry.json'), 'registry changed')
    require(config['target_sha256'] == sha(root/TARGET), 'candidate identities changed')
    engine = hashlib.sha256()
    for path in sorted((root/'src/chiral4form').glob('*.py')) + [root/'scripts/prove_combined_degree12_map.py']:
        engine.update(path.name.encode()); engine.update(path.read_bytes())
    require(engine.hexdigest() == config['engine_sha256'], 'physical engine changed')
    primes, seeds, rows = config['primes'], config['seeds'], config['rows']
    require(len(set(primes)) == len(primes) and all(type(p) is int and sp.isprime(p) for p in primes), 'invalid primes')
    require(seeds == list(range(202611130000,202611130072)), 'unexpected physical sample seeds')
    require(len(rows) == 36 and [r['source_row'] for r in rows] == list(range(36)), 'source ordering mismatch')
    records = []
    for prime in primes:
        data = json.loads((directory/f'prime{prime}.json').read_text())
        require(data['schema'] == 1 and type(data['prime']) is int and data['prime'] == prime, 'prime/schema mismatch')
        require(data['namespace_sha256'] == digest and data['seeds'] == seeds, 'checkpoint namespace mismatch')
        require(data['status'] == 'complete' and len(data['samples']) == 72, 'incomplete checkpoint')
        for seed,sample in zip(seeds,data['samples']):
            require(sample['seed'] == seed, 'sample ordering mismatch')
            for key,size in [('basis',72),('targets',36)]:
                require(len(sample[key]) == size and all(type(v) is int and 0 <= v < prime for v in sample[key]), 'invalid physical vector')
        records.append(data)
    bound = config['basis_bound']
    require(bound == 10**30, 'unexpected graph bound')
    basis = [[lift([data['samples'][i]['basis'][j] for data in records],primes,bound) for j in range(72)] for i in range(72)]
    require(rank_mod(basis,primes[0]) == 72, 'evaluation matrix not injective')
    candidate = json.loads((root/TARGET).read_text())
    coordinates = {r['source_row']: list(map(Fraction,r['coordinates'])) for r in candidate['rows']}
    require(set(coordinates) == set(range(36)) and all(len(v) == 72 for v in coordinates.values()), 'candidate shape mismatch')
    for j,row in enumerate(rows):
        power = sum(int(token[2:]) for token in row['generator'].split('*'))
        require(row['stress_denominator'] == (240**2)**power, 'stress denominator mismatch')
        denominator = row['denominator']
        require(denominator == row['orbit_denominator']*row['stress_denominator'], 'combined denominator mismatch')
        require(row['integer_bound'] == row['target_absolute_bound_before_stress_clearer']*row['stress_denominator'], 'bound scaling mismatch')
        require(all(denominator % p for p in primes), 'bad denominator prime')
        for i in range(72):
            actual = lift([data['samples'][i]['targets'][j]*(denominator % p) % p for p,data in zip(primes,records)],primes,row['integer_bound'])
            expected = denominator*sum(q*v for q,v in zip(coordinates[j],basis[i]))
            require(expected.denominator == 1 and actual == expected, f'source identity residual: row {j}, point {i}')
    counts = json.loads((root/'foundation/science/character/independent-d5.json').read_text())['degrees']
    require(next(int(r['invariants']) for r in counts if r['degree'] == 12) == 72, 'character dimension mismatch')
    return {'schema':1,'passed':True,'namespace_sha256':digest,'prime_count':len(primes),'point_count':72,'source_count':36,
            'modulus_product':str(prod(primes)),'independent_modular_rank':72,'exact_source_comparisons':2592,
            'checkpoint_sha256':{f'prime{p}.json':sha(directory/f'prime{p}.json') for p in primes},
            'checker_sha256':sha(Path(__file__)),
            'scope':'independent saved-evidence CRT, denominator scaling, rank and rational identity arithmetic; imports neither tensor engine',
            'assumptions':['declared analytic stress formula and derived absolute bounds','saved tensor values faithfully produced by hash-bound engine','independent character bound applies to this representation'],
            'not_established':['fresh tensor evaluation','independent derivation of absolute bounds','all-orders classification','physical admissibility','novelty']}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--directory',type=Path,default=Path('reports/degree12-map'))
    ap.add_argument('--output',type=Path,default=Path('reports/degree12-map/independent-arithmetic.json'))
    args=ap.parse_args();report=verify(args.root,args.directory)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['passed','prime_count','exact_source_comparisons','scope']}))

if __name__ == '__main__':main()
