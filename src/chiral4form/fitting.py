"""Exact finite-field basis coordinates with disjoint fit and holdout samples.

Passing sampled tests is deliberately not called a symbolic polynomial
identity proof. The ambient basis's completeness is a separate proof input.
"""
from __future__ import annotations
from dataclasses import asdict
from .certificates import build_quotient_certificate
from .finite_field import independent_rows,solve_many,matmul
from .polynomial import Poly


def fit_degree(samples,basis_ids,degree,prime,fit_count):
    if any(s.get('prime') != prime for s in samples):
        raise ValueError('all samples must belong to the declared prime field')
    if not len(basis_ids)<=fit_count<len(samples):
        raise ValueError('need sufficient fit samples AND at least one disjoint holdout')
    if len({s['seed'] for s in samples})!=len(samples):
        raise ValueError('fit and holdout seed identities must be distinct')
    A=[s['basis_values'][str(degree)] for s in samples]
    if any(len(row)!=len(basis_ids) for row in A):
        raise ValueError('basis evaluation dimension mismatch')
    selected=independent_rows(A[:fit_count],prime)
    if len(selected)!=len(basis_ids):
        raise ValueError('fit sample basis is rank deficient; use additional fresh samples')
    selected=selected[:len(basis_ids)]
    keys=sorted({(g['generator'],tuple(g['monomial'])) for s in samples for g in s['targets'] if g['degree']==degree})
    dictionaries=[{(g['generator'],tuple(g['monomial'])):g['value'] for g in s['targets'] if g['degree']==degree} for s in samples]
    B=[[d.get(key,0) for key in keys] for d in dictionaries]
    if keys:
        X=solve_many([A[i] for i in selected],[B[i] for i in selected],prime)
        predicted=matmul(A,X,prime)
        for i,(actual,expected) in enumerate(zip(B,predicted)):
            if actual!=expected:
                raise ValueError(f'coordinate identity fails on {"holdout" if i>=fit_count else "fit"} seed {samples[i]["seed"]}')
        rows=[[X[j][i] for j in range(len(basis_ids))] for i in range(len(keys))]
    else:
        rows=[]
    certificate=build_quotient_certificate(rows,len(basis_ids),prime)
    return {'degree':degree,'prime':prime,'basis_ids':list(basis_ids),
        'fit_seed_indices':selected,'holdout_seed_indices':list(range(fit_count,len(samples))),
        'sample_seeds':[s['seed'] for s in samples], 'basis_evaluations':A,
        'target_evaluations':B,'all_holdouts_passed':True,
        'coordinates':[{'generator':key[0],'monomial':list(key[1]),'row':row} for key,row in zip(keys,rows)],
        'coefficient_row_span_certificate':asdict(certificate),
        'interpretation':'sample_validated_coefficients_not_a_nonlinear_reachable_set',
        'identity_proof_status':'finite_field_fit_and_disjoint_holdout_not_symbolic'}


def vector_fields_from_fits(fits):
    if not fits:
        raise ValueError('nonempty fitting result required')
    p=fits[0]['prime']
    if any(f['prime']!=p for f in fits):
        raise ValueError('cannot mix prime fields')
    ids=tuple(i for f in sorted(fits,key=lambda x:x['degree']) for i in f['basis_ids'])
    if len(set(ids))!=len(ids):
        raise ValueError('coordinate IDs must be unique across homogeneous pieces')
    index={name:i for i,name in enumerate(ids)};n=len(ids)
    generators=sorted({c['generator'] for f in fits for c in f['coordinates']})
    raw={g:[{} for _ in ids] for g in generators}
    for f in fits:
        for c in f['coordinates']:
            m=[0]*n
            for name in c['monomial']:
                if name not in index:
                    raise ValueError(f'coefficient variable {name} is absent from the registry')
                m[index[name]]+=1
            for name,value in zip(f['basis_ids'],c['row']):
                j=index[name];key=tuple(m)
                raw[c['generator']][j][key]=(raw[c['generator']][j].get(key,0)+value)%p
    return ids,{g:tuple(Poly(n,terms,p) for terms in coords) for g,coords in raw.items()}


def fields_to_json(ids,fields,status):
    return {'schema':1,'coordinate_ids':list(ids),'coefficient_status':status,
            'fields':{name:[f.to_json() for f in field] for name,field in fields.items()}}


def fields_from_json(record):
    ids=record['coordinate_ids'];fields={name:tuple(Poly.from_json(x) for x in field) for name,field in record['fields'].items()}
    if len(set(ids))!=len(ids) or not fields or any(len(f)!=len(ids) or any(x.n!=len(ids) for x in f) for f in fields.values()):
        raise ValueError('invalid polynomial vector-field model')
    primes={x.p for f in fields.values() for x in f}
    if len(primes)!=1:
        raise ValueError('all fields must use the same exact coefficient domain')
    return tuple(ids),fields
