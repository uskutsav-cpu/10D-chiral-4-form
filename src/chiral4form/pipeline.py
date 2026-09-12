"""A resumable, sequential tensor -> fitting -> polynomial-model pipeline.

Large degree-12 jobs are opt-in and are NOT declared done from imported ranks.
Every stage has a content-addressed checkpoint and a scoped output statement.
"""
from __future__ import annotations
import json
import platform
from pathlib import Path
import numpy as np
import sympy
from .checkpoints import CheckpointStore
from .provenance import atomic_json,semantic_hash,sha256_file
from .finite_field import _check_prime_field_modulus,matrix_rank
from .registry import Registry,import_predecessor_registry
from .source import acquire_source
from .forms import random_selfdual,compact,hodge,inner
from .stress import evaluate_generators,coefficient_ids
from .tensors import TensorBudget
from .graphs import plan,allocation_profile
from .fitting import fit_degree,vector_fields_from_fits,fields_to_json
from .reachability import linear_invariant_hull,lie_rank
from .sextic import sextic_report,verify_sampled_sextic
from .degree8_orbit import report as orbit8_report,verify_sampled as verify_sampled8
from .localization import modmax_symbolic_report


def engine_fingerprint():
    root=Path(__file__).parent
    return semantic_hash({f.name:sha256_file(f) for f in sorted(root.glob('*.py'))})


def validate_config(config):
    if not isinstance(config, dict):
        raise ValueError('configuration must be a JSON object')
    for key in ('schema','max_degree','seed_start','fit_margin','holdouts','max_bytes','max_multiply_adds','lie_depth'):
        if key in config and type(config[key]) is not int:
            raise TypeError(f'{key} must be an integer')
    for key in ('exact_blas','flow_analysis'):
        if key in config and type(config[key]) is not bool:
            raise TypeError(f'{key} must be a boolean')
    if not isinstance(config.get('primes'), list) or not isinstance(config.get('extras', []), list):
        raise TypeError('primes and extras must be lists')
    if any(not isinstance(x, str) for x in config.get('extras', [])):
        raise TypeError('extra generator IDs must be strings')
    allowed={'schema','max_degree','primes','seed_start','fit_margin','holdouts','registry',
             'extras','max_bytes','max_multiply_adds','exact_blas','lie_depth','flow_analysis'}
    if set(config)-allowed:
        raise ValueError(f'unknown configuration keys: {sorted(set(config)-allowed)}')
    if config.get('schema')!=1 or config.get('max_degree') not in (4,6,8,10,12):
        raise ValueError('schema 1 and a supported even cutoff are required')
    primes=config.get('primes',[])
    if not primes or len(set(primes))!=len(primes):
        raise ValueError('nonempty distinct prime list required')
    for p in primes:
        _check_prime_field_modulus(p)
        if p<=5 or p>65521:
            raise ValueError('tensor backend needs 5 < prime <= 65521 for constrained derivatives')
    if config.get('fit_margin',2)<0 or config.get('holdouts',2)<2 or config.get('seed_start',0)<0:
        raise ValueError('invalid sampling policy; at least two holdout seeds required')
    if not 1<=config.get('lie_depth',3)<=8:
        raise ValueError('Lie depth must be 1..8; use the explicit analysis API for larger work')
    if len(set(config.get('extras',[])))!=len(config.get('extras',[])):
        raise ValueError('duplicate generalized invariant generator')
    return config


def prepare(config,root,source=None):
    validate_config(config);root=Path(root)
    if source is not None:
        source,manifest=acquire_source(source)
        registry=import_predecessor_registry(source,config['max_degree'])
    else:
        path=root/config.get('registry','data/fixtures/low_degree_registry.json')
        registry=Registry.load(path)
        manifest={'status':'bundled_curated_excerpt_not_full_upstream_snapshot','fixture_sha256':sha256_file(path),
                  'metadata':registry.metadata}
    needed=[d for d in (4,6,8,10,12) if d<=config['max_degree']]
    if any(d not in registry.degree_bases for d in needed):
        raise ValueError('the selected registry lacks the requested degree; freeze/import the actual predecessor source')
    budget=TensorBudget(config.get('max_bytes',512*1024**2),config.get('max_multiply_adds',5_000_000_000),config.get('exact_blas',True))
    return registry,manifest,budget


def execution_plan(config,root,source=None):
    reg,manifest,budget=prepare(config,root,source)
    degree=config['max_degree'];plans=[]
    # Through degree D, squared-derivative terms inside Tr(tau^k), k>=2,
    # need invariant gradients only through D-4.
    for i in reg.items.values():
        if i.degree<=degree and i.graph is not None:
            _,(peak,work)=plan(i.graph,10)
            gradient=i.degree<=degree-4
            reserve=allocation_profile(i.graph,10,gradient)['reserved_bytes']
            plans.append({'invariant':i.id,'degree':i.degree,'gradient_needed':gradient,
                          'reserved_bytes':reserve,'multiply_adds':work,
                          'within_budget':reserve<=budget.max_bytes and work<=budget.max_multiply_adds})
    largest=max(len(reg.degree_bases[d]) for d in reg.degree_bases if d<=degree)
    return {'status':'execution_plan_only','source':manifest,'registry_fingerprint':reg.fingerprint,
        'sample_count_per_prime':largest+config.get('fit_margin',2)+config.get('holdouts',2),
        'primes':config['primes'],'all_graphs_within_budget':all(x['within_budget'] for x in plans),'graphs':plans}


def compute_sample(reg,p,seed,degree,extras,budget):
    form=random_selfdual(seed,p);cache={}
    if not np.array_equal(hodge(form,p),form) or inner(form,form,p)!=0:
        raise AssertionError('self-duality sample gate failed')
    targets,catalog=evaluate_generators(reg,form,p,degree,extras,cache=cache,budget=budget)
    basis_values={str(d):reg.basis_values(d,form,p,cache=cache,budget=budget) for d in sorted(reg.degree_bases) if d<=degree}
    return {'seed':seed,'prime':p,'compact_selfdual_input':compact(form,p).tolist(),
        'basis_values':basis_values,
        'targets':[{'generator':g,'degree':d,'monomial':list(m),'value':v} for (g,d,m),v in sorted(targets.items())],
        'generator_catalogue':[{'id':g.id,'factors':list(g.factors),'leading_degree':g.leading_degree} for g in catalog]}


def run(config,root,output,source=None,progress=print):
    root,output=Path(root),Path(output);output.mkdir(parents=True,exist_ok=True)
    reg,manifest,budget=prepare(config,root,source)
    preflight=execution_plan(config,root,source)
    atomic_json(output/'execution_plan.json',preflight)
    if not preflight['all_graphs_within_budget']:
        failed=[x['invariant'] for x in preflight['graphs'] if not x['within_budget']]
        atomic_json(output/'status.json',{'state':'preflight_refused','over_budget_graphs':failed})
        raise MemoryError('graph resource preflight refused: '+', '.join(failed))
    environment={'python':platform.python_version(),'numpy':np.__version__,'sympy':sympy.__version__,
                 'platform':platform.platform()}
    namespace={'engine':engine_fingerprint(),'registry':reg.fingerprint,'config':config,
               'numpy_version':np.__version__,'source':manifest}
    key=semantic_hash(namespace)
    store=CheckpointStore(output/'checkpoints',namespace)
    atomic_json(output/'configuration.json',config)
    atomic_json(output/'source_manifest.json',manifest)
    atomic_json(output/'registry.json',reg.to_json())
    atomic_json(output/'environment.json',environment)
    degree=config['max_degree'];degrees=[d for d in sorted(reg.degree_bases) if d<=degree]
    fit_count=max(len(reg.degree_bases[d]) for d in degrees)+config.get('fit_margin',2)
    count=fit_count+config.get('holdouts',2);extras=tuple(config.get('extras',()))
    reports=[]
    run_status={'schema':1,'run_fingerprint':key,'state':'running','primes_finished':[],
                'degree':degree,'scientific_status':'computing_samples_not_verifying_imported_orbit_claims'}
    atomic_json(output/'status.json',run_status)
    try:
        for p in config['primes']:
            samples=[]
            for i in range(count):
                seed=config.get('seed_start',2026091201)+i
                sample,cached=store.run('tensor_sample',{'prime':p,'seed':seed},
                    lambda p=p,seed=seed:compute_sample(reg,p,seed,degree,extras,budget))
                samples.append(sample)
                progress(f'prime {p}: sample {i+1}/{count} {"cached" if cached else "computed"}',flush=True)
                run_status.update({'current_prime':p,'sample_finished':i+1,'samples_per_prime':count})
                atomic_json(output/'status.json',run_status)
            fits=[]
            for d in degrees:
                fitted,_=store.run('basis_fit',{'prime':p,'degree':d,'sample_hash':semantic_hash(samples)},
                    lambda d=d:fit_degree(samples,reg.degree_bases[d],d,p,fit_count))
                fits.append(fitted)
                atomic_json(output/f'fit_degree{d}_prime{p}.json',fitted)
            ids,fields=vector_fields_from_fits(fits)
            catalog=samples[0]['generator_catalogue']
            model=fields_to_json(ids,fields,'finite_field_samples_and_disjoint_holdouts; not a symbolic identity proof')
            model['generator_catalogue']=catalog
            atomic_json(output/f'fields_prime{p}.json',model)
            pure={name:field for name,field in fields.items() if not any(x.startswith('S:') for x in next(g['factors'] for g in catalog if g['id']==name))}
            checks={}
            if degree==6:
                checks['independent_analytic_sextic_map']=verify_sampled_sextic(pure,ids,p)
            if degree==8:
                checks['independent_analytic_degree8_map']=verify_sampled8(pure,ids,p)
            rep={'prime':p,'fit_count':fit_count,'holdouts':count-fit_count,
                 'basis_ranks':{str(d):len(reg.degree_bases[d]) for d in degrees},'cross_checks':checks}
            if config.get('flow_analysis',degree<=8):
                hull=linear_invariant_hull(list(fields.values()))
                lie=lie_rank(list(fields.values()),max_depth=config.get('lie_depth',3))
                rep.update({'linear_invariant_hull':hull,'finite_depth_lie_rank':lie})
                W=[[int(x)%p for x in row] for row in hull.get('basis',[])]
                rep['hull_projection_ranks']={str(d):matrix_rank([[row[ids.index(i)] for i in reg.degree_bases[d]] for row in W],p) for d in degrees}
            else:
                rep['flow_analysis']='not_run; coefficient fits are not a reachability classification'
            reports.append(rep);atomic_json(output/f'analysis_prime{p}.json',rep)
            run_status['primes_finished'].append(p);atomic_json(output/'status.json',run_status)
        atomic_json(output/'sextic_derivation.json',sextic_report())
        atomic_json(output/'degree8_reduced_model.json',orbit8_report())
        atomic_json(output/'modmax_known_identity.json',modmax_symbolic_report())
        summary={'schema':1,'state':'completed_declared_computations','run_fingerprint':key,'environment':environment,
            'degree':degree,'prime_reports':reports,'source_status':manifest['status'],
            'higher_degree_imported_dimensions_reproved':False,
            'warnings':['Linear invariant hull is not the nonlinear reachable set.',
                'Finite-field fits and held-out samples are not symbolic identities over QQ.',
                'The degree-eight analytic orbit report is conditional on its physics-to-basis map.',
                'ModMax output reproduces published algebra; alternative pure-stress reachability remains open.']}
        atomic_json(output/'summary.json',summary)
        run_status['state']='completed_declared_computations';atomic_json(output/'status.json',run_status)
        progress(f'COMPLETED declared degree-{degree} computations. Summary: {output/"summary.json"}',flush=True)
        return summary
    except BaseException as exc:
        run_status.update({'state':'failed_or_interrupted','error_type':type(exc).__name__,'error':str(exc)})
        atomic_json(output/'status.json',run_status)
        raise
