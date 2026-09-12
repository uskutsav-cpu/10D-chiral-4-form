import json
from copy import deepcopy
from pathlib import Path
import pytest
from chiral4form.map_proof import verify,bound_report,expected_keys,fixed_form
from chiral4form.provenance import semantic_hash

ROOT=Path(__file__).resolve().parents[1]

@pytest.fixture(scope='module')
def certificate():
    return json.loads((ROOT/'verification/map8/map_certificate.json').read_text())

def test_packaged_bounded_residual_certificate(certificate):
    assert verify(certificate)
    assert len(certificate['samples'])==63
    assert certificate['modulus_product']>2*bound_report()['integer_residual_abs_bound']
    assert len(expected_keys())>=20

@pytest.mark.parametrize('kind',['residue','bound','point','prime','rank','hash'])
def test_forged_map_certificate_rejected(certificate,kind):
    bad=deepcopy(certificate)
    if kind=='residue':bad['samples'][0]['residuals'][0]['actual']+=1
    elif kind=='bound':bad['bounds']['integer_residual_abs_bound']=0
    elif kind=='point':bad['electric_points'][0][0]=2
    elif kind=='prime':bad['primes'][0]=bad['primes'][1]
    elif kind=='rank':
        for s in bad['samples']:s['basis_values']['8']=[0]*7
    else:bad['certificate_hash']='bad'
    if kind!='hash':
        bad['certificate_hash']=semantic_hash({k:v for k,v in bad.items() if k!='certificate_hash'})
    assert not verify(bad)

def test_fixed_tensors_cannot_change_between_primes(certificate):
    import numpy as np
    from chiral4form.forms import compact
    electric=certificate['electric_points'][0]
    a=compact(fixed_form(electric,30011),30011)
    b=compact(fixed_form(electric,30013),30013)
    a=np.where(a>15000,a-30011,a)
    b=np.where(b>15000,b-30013,b)
    assert np.array_equal(a,b)
    assert max(abs(int(x)) for x in a)<=1
