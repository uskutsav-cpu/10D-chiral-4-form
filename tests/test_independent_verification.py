"""Adversarial checks of the independent algebra verifier's trust boundaries."""
import importlib.util
import json
from pathlib import Path
from copy import deepcopy
import pytest
import sympy as sp
from chiral4form.registry_resolution import resolve_registry

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('independent_chiral_checker', ROOT/'independent/verify.py')
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)

@pytest.fixture(scope='module')
def model():
    return checker.degree10(ROOT)


def test_complete_registry_resolves_without_local_runs():
    reg, provenance = resolve_registry(12)
    assert {d:len(v) for d,v in reg.degree_bases.items()} == {4:1,6:2,8:7,10:14,12:72}
    assert provenance['path'] == 'data/combined/registry.json'


def test_invalid_explicit_registry_is_not_silently_replaced(tmp_path):
    p = tmp_path/'bad.json';p.write_text('{"schema": 0}')
    with pytest.raises(ValueError):
        resolve_registry(12, explicit_path=p)


def test_sign_error_in_octic_flow_is_detected(model):
    x, fields, *_ = model
    broken = deepcopy(fields)
    broken['tr2'][4] = -broken['tr2'][4]
    with pytest.raises(ValueError, match='Omega8 transport'):
        checker.check8(broken, x)


def test_omitted_generator_is_detected(model):
    x, fields, *_ = model
    broken = {k:v for k,v in fields.items() if k != 'tr4'}
    with pytest.raises(ValueError, match='catalogue'):
        checker.check8(broken, x)


def test_source_tail_corruption_is_detected(model):
    _, _, base, _, reduced = model
    source = checker.read(ROOT, 'verification/all_orders/degree12_exact/source_equation_space.json')
    source['all_source_rows'][0]['row'][72] = str(sp.Rational(source['all_source_rows'][0]['row'][72])+1)
    with pytest.raises(ValueError, match='tail derivation'):
        checker.degree12(ROOT, base, reduced, source_override=source)


def test_corrupt_kernel_is_rejected_even_with_passed_flags(model, monkeypatch):
    _, _, base, _, reduced = model
    original = checker.read
    def read(root, path):
        data = original(root, path)
        if path.endswith('/kernel_lift.json'):
            data['lift']['rational_rows'][0][0] = str(sp.Rational(data['lift']['rational_rows'][0][0])+1)
            data['exact'] = data['lift']['qq_verified'] = True
        return data
    monkeypatch.setattr(checker, 'read', read)
    with pytest.raises(ValueError, match='tangency residual'):
        checker.degree12(ROOT, base, reduced)


def test_constructive_endpoint_has_exact_signed_controls(model):
    _, _, base, _, reduced = model
    report, variables, fields, constraints, free, matrix, _ = checker.degree12(ROOT, base, reduced)
    target = ['-1','2/3','0','-4/5','1','2','-1/2','0','7/3','-1']
    result = checker.endpoint(variables, fields, constraints, free, matrix, target)
    assert result['all_68_constraints_zero']
    assert len(result['segments']) == 10
    assert report['vertical_determinant_QQ'] != '0'
    assert any(sp.Rational(segment['amount']) < 0 for segment in result['segments'])


def test_generalized_completion_cli_recomputes_all_exact_blocks(tmp_path):
    import subprocess
    import sys
    output=tmp_path/'completion.json'
    subprocess.run([sys.executable,'independent/completion.py','--output',str(output)],cwd=ROOT,check=True,capture_output=True,text=True)
    data=json.loads(output.read_text())
    assert data['selected_extra_count'] == 78
    assert {int(d): r['dimension'] for d,r in data['exact_constant_leading_blocks'].items()} == {4:1,6:2,8:7,10:14,12:72}
    assert all(sp.Rational(r['determinant_QQ']) != 0 for r in data['exact_constant_leading_blocks'].values())
    assert 'not claimed' in data['minimality']


def test_degree12_raw_checkpoint_rejects_prime_relabeling_and_noncanonical_values():
    import sys
    sys.path.insert(0,str(ROOT/'scripts'))
    try:
        from verify_degree12_checkpoint_inputs import validate_checkpoint
        from prove_combined_degree12_map import SEEDS
        record={'schema':1,'prime':101,'namespace_sha256':'unit-test','seeds':list(SEEDS),
                'status':'complete','samples':[{'seed':seed,'basis':[0]*72,'targets':[0]*36} for seed in SEEDS]}
        validate_checkpoint(record,101,'unit-test')
        with pytest.raises(ValueError,match='prime/schema'):
            validate_checkpoint(record,103,'unit-test')
        record['samples'][0]['targets'][0]=101
        with pytest.raises(ValueError,match='canonical'):
            validate_checkpoint(record,101,'unit-test')
    finally:
        sys.path.pop(0)


def test_classifier_rejects_explicit_octic_obstruction():
    result=checker.classify(ROOT,{'I8_2':'1'})
    assert result['classification']=='unreachable'
    assert result['Omega8']=='1'
    assert result['violated_coefficient_graph_relations_through_degree10']['I8_2']=='1'
    with pytest.raises(ValueError,match='rational strings'):
        checker.classify(ROOT,{'I4_1':0.1})
    with pytest.raises(ValueError,match='unknown coefficient'):
        checker.classify(ROOT,{'not_an_invariant':'1'})


def test_classifier_constructs_every_coefficient_of_reachable_interaction(model):
    _, _, base, graph10, reduced=model
    _, variables, fields, constraints, free, matrix, _=checker.degree12(ROOT,base,reduced)
    report=checker.attach_full_coefficients(checker.endpoint(variables,fields,constraints,free,matrix,['1','0','0','0','0','0','0','0','0','0']),ROOT,base,graph10)
    classified=checker.classify(ROOT,report['coefficients'])
    assert classified['classification']=='reachable'
    assert len(classified['coefficients'])==96
    assert classified['coefficients']==report['coefficients']
    assert classified['all_68_constraints_zero']
