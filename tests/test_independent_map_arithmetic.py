"""Adversarial tests for the engine-free saved-evidence arithmetic checker."""
import importlib.util
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('independent_map12',ROOT/'independent/map12_arithmetic.py')
checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker)


def test_centered_crt_recovers_signed_integers_and_rejects_bound_violation():
    primes=[101,103,107]
    for value in [-10000,-1,0,1,10000]:
        assert checker.lift([value%p for p in primes],primes,10000)==value
    with pytest.raises(ValueError,match='exceeds declared bound'):
        checker.lift([10001%p for p in primes],primes,10000)
    with pytest.raises(ValueError,match='insufficient modulus'):
        checker.lift([1],[101],51)
    with pytest.raises(ValueError,match='noncanonical'):
        checker.lift([101,1,1],primes,10000)


def test_independent_rank_distinguishes_integer_and_modular_dependence():
    assert checker.rank_mod([[0,1],[1,0]],101)==2
    assert checker.rank_mod([[1,2,3],[2,4,6]],101)==1
    assert checker.rank_mod([[101,0],[0,1]],101)==1
    assert checker.rank_mod([[101,0],[0,1]],103)==2


def test_configuration_mutation_is_rejected_before_success_flags(tmp_path):
    import json
    directory=tmp_path/'proof';directory.mkdir()
    (directory/'configuration.json').write_text(json.dumps({'configuration':{'schema':1,'passed':True},'namespace_sha256':'forged'}))
    with pytest.raises(ValueError,match='configuration digest mismatch'):
        checker.verify(ROOT,directory)


def test_prime_launcher_resumes_json_configuration_without_recomputing(tmp_path):
    import hashlib
    import shutil
    import subprocess
    import sys
    source=ROOT/'reports/degree12-map'
    for name in ['configuration.json','prime65521.json']:
        shutil.copyfile(source/name,tmp_path/name)
    checkpoint=tmp_path/'prime65521.json'
    original=hashlib.sha256(checkpoint.read_bytes()).hexdigest()
    subprocess.run([sys.executable,'scripts/run_degree12_prime.py','--prime','65521','--output',str(tmp_path)],cwd=ROOT,check=True,capture_output=True,text=True)
    assert hashlib.sha256(checkpoint.read_bytes()).hexdigest()==original
