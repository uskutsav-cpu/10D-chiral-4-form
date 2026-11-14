"""Small direct-contraction regressions for the independently read red brackets.

The patched operand is an arbitrary covariant six-slot tensor in dimension 3.
This isolates bracket placement and normalization; it does not pretend to be a
10D 1050 tensor. The oracle contracts all five operands as bounded int64 data,
without the production BracketProgram or modular contraction engine.
"""
import itertools
import math
import numpy as np
import pytest
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sdinv import published_degree10_invariants as published


def alternating(tensor, slots, prime):
    total=np.zeros_like(tensor)
    for perm in itertools.permutations(range(len(slots))):
        parity=(-1)**sum(perm[i]>perm[j] for i in range(len(perm)) for j in range(i+1,len(perm)))
        order=list(range(6))
        for k,slot in enumerate(slots): order[slot]=slots[perm[k]]
        total += parity*tensor.transpose(order)
    return total%prime*pow(math.factorial(len(slots)),-1,prime)%prime


def raised(tensor, slots, prime):
    result=tensor.copy()
    for slot in slots:
        axes=[slice(None)]*tensor.ndim;axes[slot]=0
        result[tuple(axes)]=-result[tuple(axes)]%prime
    return result


def scalar(word, tensors, prime):
    seen=set();ops=[]
    for term,tensor in zip(word.split(','),tensors):
        ops.append(raised(tensor,[i for i,c in enumerate(term) if c in seen],prime))
        seen.update(term)
    # At most 3**15 summands, each <=100**5: strictly below int64 maximum.
    assert 3**15*(prime-1)**5 < 2**63
    return int(np.einsum(word+'->',*ops,optimize='optimal')%prime)


@pytest.mark.parametrize('candidate', [10,11,12])
def test_source_red_program_matches_independent_direct_contraction(monkeypatch,candidate):
    prime=101
    Q=np.random.default_rng(90507+candidate).integers(0,prime,size=(3,)*6,dtype=np.int64)
    untouched=Q.copy()
    T=alternating(Q,(3,4,5),prime);B=alternating(Q,(4,5),prime)
    words={10:'abcdmn,abcxyz,dpqxzy,efghpq,efghmn',
           11:'rstabc,defabc,defijk,rslmni,tjmlnk',
           12:'rstabc,defabc,defijk,irlsmn,tjmkln'}
    factors=[B,Q,Q,Q,Q] if candidate==10 else [Q,Q,T,T,T]
    expected=scalar(words[candidate],factors,prime)
    monkeypatch.setattr(published,'composite_n1050',lambda *args:Q)
    monkeypatch.setattr(published,'_raise_axes',raised)
    actual=published.PUBLISHED_DEGREE10[f'P10_{candidate:02d}']['evaluator'](None,prime)
    assert actual==expected
    assert np.array_equal(Q,untouched)
    legacy=published.LEGACY_PUBLISHED_DEGREE10[f'P10_{candidate:02d}']['evaluator'](None,prime)
    expected_legacy=scalar(words[candidate],[Q]*5,prime)
    assert legacy==expected_legacy
    assert actual != legacy


def test_source_registry_promotes_resolved_red_readings():
    assert published.AMBIGUITY_VARIANTS=={}
    for number in [4,9,10,11,12]:
        name=f'P10_{number:02d}'
        assert published.PUBLISHED_DEGREE10[name]['source_verified']
        assert 'RED' in published.BRACKET_STAGES[name]
        assert 'RED' in published.PUBLISHED_DEGREE10[name]['brackets']
        assert not published.LEGACY_PUBLISHED_DEGREE10[name]['source_verified']


def test_legacy_registry_describes_its_actual_outer_brackets():
    for number in [10,11,12]:
        name=f"P10_{number:02d}"
        assert "legacy outer:" in published.LEGACY_PUBLISHED_DEGREE10[name]["brackets"]
        assert "no red trailing" in published.LEGACY_PUBLISHED_DEGREE10[name]["brackets"]
        assert "RED" in published.PUBLISHED_DEGREE10[name]["brackets"]
