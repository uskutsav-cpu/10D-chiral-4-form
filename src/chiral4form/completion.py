"""Finite-catalogue completion using all mixed f(T,S) generator monomials."""
from itertools import combinations
from .fitting import fields_from_json
from .reachability import linear_invariant_hull,lie_rank


def search(record,candidates,target_rank,metric='lie_rank',max_depth=3,max_subsets=1024):
    ids,fields=fields_from_json(record);catalog={g['id']:g['factors'] for g in record.get('generator_catalogue',[])}
    if set(fields)-set(catalog):
        raise ValueError('full factor provenance is required for generalized flow completion')
    if metric not in ('lie_rank','linear_hull') or not 1<=target_rank<=len(ids):
        raise ValueError('invalid completion objective')
    if len(set(candidates))!=len(candidates):
        raise ValueError('candidate generator names must be distinct')
    available={x[2:] for factors in catalog.values() for x in factors if x.startswith('S:')}
    if set(candidates)-available:
        raise ValueError('missing genuine generalized-generator evaluations; seed augmentation is not accepted')
    tested=[]
    for k in range(len(candidates)+1):
        for subset in combinations(sorted(candidates),k):
            if len(tested)>=max_subsets:
                return {'status':'resource_limit','tested_subsets':tested,'global_minimality':False}
            chosen=set(subset)
            active=[field for name,field in fields.items() if all(not x.startswith('S:') or x[2:] in chosen for x in catalog[name])]
            rep=lie_rank(active,max_depth=max_depth) if metric=='lie_rank' else linear_invariant_hull(active)
            tested.append({'subset':list(subset),'rank':rep['rank']})
            if rep['rank']>=target_rank and (metric!='linear_hull' or rep['stabilized']):
                return {'status':'smallest_subset_meeting_declared_finite_test','metric':metric,'subset':list(subset),
                    'cardinality':k,'target_rank':target_rank,'tested_subsets':tested,'global_minimality':False,
                    'warning':'minimum for this catalogue and this test only; finite Lie depth cannot exclude deeper reachability for smaller sets'}
    return {'status':'no_subset_meets_declared_finite_test','tested_subsets':tested,'global_minimality':False}
