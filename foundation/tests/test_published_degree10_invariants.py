"""Literal-source evaluators plus explicitly archived pre-audit projection evidence.

Source bracket extent is resolved from original arXiv v2 TeX and colored PDF.
Tests reading MAP below are historical artifact checks, not evidence for the
corrected-source registry. Fresh corrected-source map validation is separate.
"""

import json
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.modp import P, ALT_P, mod_einsum  # noqa: E402
from sdinv.forms import (  # noqa: E402
    selfdual_projector, to_dense, random_form, metric_signs)


def _transform_covariant_tensor(tensor, transformation, mod):
    """Apply one matrix per covariant index, reducing after every axis."""
    result = np.asarray(tensor, dtype=np.int64) % mod
    for axis in range(result.ndim):
        result = np.tensordot(transformation, result, axes=(1, axis))
        result = np.moveaxis(result, 0, axis) % mod
    return result


from sdinv.published_degree10_invariants import (  # noqa: E402
    NOT_IMPLEMENTED, PUBLISHED_DEGREE10, LEGACY_PUBLISHED_DEGREE10, evaluate_implemented)

ROOT = Path(__file__).resolve().parents[1]
MAP = ROOT / "results" / "intrinsic_candidates" / "published_degree10_map.json"


def _sample(prime, seed):
    projector = selfdual_projector(10, 5, True, prime)
    raw = random_form(10, 5, np.random.default_rng(seed), prime)
    return to_dense((projector @ raw) % prime, 10, 5, prime)


def test_implemented_and_unimplemented_partition_all_twelve():
    """Every one of the twelve must be either implemented or explicitly not.

    This guard exists so a candidate can never be silently dropped: the two
    sets must be disjoint and together cover P10_01..P10_12.
    """
    from sdinv.published_degree10_invariants import BRACKET_STAGES
    implemented = set(PUBLISHED_DEGREE10)
    blocked = set(NOT_IMPLEMENTED)
    assert implemented & blocked == set(), "a candidate is in both sets"
    assert implemented | blocked == {f"P10_{i:02d}" for i in range(1, 13)}
    # every implemented candidate above P10_02 must record its bracket stage,
    # because that is what decides whether the engine is needed
    for name in implemented:
        if name in ("P10_01", "P10_02"):
            continue
        assert name in BRACKET_STAGES, f"{name} has no recorded bracket stage"


def test_red_bracket_candidates_have_verified_source_stages():
    from sdinv.published_degree10_invariants import BRACKET_STAGES, AMBIGUITY_VARIANTS
    assert AMBIGUITY_VARIANTS == {}
    for name in ("P10_04", "P10_09", "P10_10", "P10_11", "P10_12"):
        assert "RED" in BRACKET_STAGES[name]
        assert "RED" in PUBLISHED_DEGREE10[name]["brackets"]
        assert PUBLISHED_DEGREE10[name]["source_verified"]
        assert not PUBLISHED_DEGREE10[name].get("ambiguity")


def test_p10_04_source_and_legacy_pair_readings_are_distinguished():
    """The RED symmetrisation must do something, and AMB-01 must be a real fork.

    If the two readings of the red bracket agreed, the ambiguity would be
    harmless and could be closed. This test records which of those worlds we
    are in, so the ambiguity is never quietly assumed away.
    """
    from sdinv.published_degree10_invariants import (
        p10_04_mm_m_red_n1050_n1050, p10_04_pairs)

    prime = P
    form = _sample(prime, 5)
    all4 = p10_04_mm_m_red_n1050_n1050(form, prime)
    pairs = p10_04_pairs(form, prime)
    assert isinstance(all4, int) and isinstance(pairs, int)
    # Record, do not assume: if these ever coincide, AMB-01 is closable and the
    # message below is what tells a future session so.
    if all4 == pairs:
        raise AssertionError(
            "The source and historical pair reading coincide at this sample")


def test_homogeneity_degree_ten():
    """Guards the int64 overflow that a bare np.einsum would reintroduce."""
    for prime in (P, ALT_P):
        form = _sample(prime, 20260729)
        base = evaluate_implemented(form, prime)
        for c in (2, 3, 5):
            scaled = evaluate_implemented((form * c) % prime, prime)
            for name, v in base.items():
                assert scaled[name] == (pow(c, 10, prime) * v) % prime, (
                    f"{name} not homogeneous of degree 10 at prime {prime}; a "
                    f"bare np.einsum on int64 silently wraps here")


def _boost(mod, t=7):
    """A hyperbolic rotation in the 0-1 plane over F_p: c^2 - s^2 = 1."""
    ti = pow(t, mod - 2, mod)
    half = pow(2, mod - 2, mod)
    c = ((t + ti) * half) % mod
    s = ((ti - t) * half) % mod
    assert (c * c - s * s) % mod == 1, "not a hyperbolic rotation"
    L = np.eye(10, dtype=np.int64)
    L[0, 0] = c; L[0, 1] = s; L[1, 0] = s; L[1, 1] = c
    eta = np.diag(metric_signs(10, True)).astype(np.int64) % mod
    assert np.array_equal((L.T @ eta @ L) % mod, eta % mod), "L not in SO(1,9)"
    return L


def test_every_published_candidate_is_boost_invariant():
    """THE test that was missing, and the one that catches metric misplacement.

    A published candidate is by construction a Lorentz scalar. If an index is
    raised on both ends of a contracted edge, that edge contracts with delta
    rather than eta. Under a pure ROTATION delta and eta agree on the spatial
    block, so the error is invisible; a BOOST mixes the timelike direction and
    exposes it.

    This is not hypothetical. P10_07 shipped with all six axes raised on both
    inner N factors, making the three alpha edges delta-contractions. It passed
    homogeneity, passed rotation invariance, and produced entirely plausible
    numbers. The only signals were this test -- absent at the time -- and the
    projection reporting `not_in_atlas_span`.
    """
    for prime in (P, ALT_P):
        form = _sample(prime, 5)
        boosted = _transform_covariant_tensor(form, _boost(prime), prime)
        base = evaluate_implemented(form, prime)
        moved = evaluate_implemented(boosted, prime)
        for name in base:
            assert base[name] == moved[name], (
                f"{name} is not boost invariant at prime {prime}: "
                f"{base[name]} != {moved[name]}. An index is almost certainly "
                f"raised on both ends of some contracted edge.")


def test_p10_07_alpha_edge_placement_is_free():
    """Either end of an alpha edge may carry the metric; the scalar is the same.

    This distinguishes a genuinely correct index placement from one that merely
    happens to be boost invariant. For a true tensor contraction, moving the
    metric from one end of an edge to the other cannot change the value.
    """
    from sdinv.stress import _raise_axes, composite_n1050, five_form_moment

    def variant(form, mod, raise_a, raise_b):
        _, mixed = five_form_moment(form, mod, "optimized")
        mm = (mixed @ mixed) % mod
        n = composite_n1050(form, mod, "optimized")
        return int(mod_einsum(
            "pqrstm,mk,pqrxyz,stkxzy->",
            [_raise_axes(n, (5,), mod), mm,
             _raise_axes(n, raise_a, mod), _raise_axes(n, raise_b, mod)],
            mod) % mod)

    for prime in (P, ALT_P):
        form = _sample(prime, 5)
        shipped = PUBLISHED_DEGREE10["P10_07"]["evaluator"](form, prime)
        flipped = variant(form, prime, (0, 1, 2, 3, 4, 5), (0, 1))
        assert shipped == flipped, (
            f"alpha-edge metric placement changed the value at {prime}: "
            f"{shipped} != {flipped}; the contraction is not a tensor")


def test_the_original_p10_07_placement_really_was_broken():
    """Guards the fix against silent reversion.

    Raising all six axes on both inner N factors -- the shipped bug -- must
    FAIL boost invariance. Without this, a regression to the old placement
    could pass the suite if the boost test above were ever weakened.
    """
    from sdinv.stress import _raise_axes, composite_n1050, five_form_moment

    prime = P
    form = _sample(prime, 5)
    boosted = _transform_covariant_tensor(form, _boost(prime), prime)

    def broken(f):
        _, mixed = five_form_moment(f, prime, "optimized")
        mm = (mixed @ mixed) % prime
        n = composite_n1050(f, prime, "optimized")
        up = _raise_axes(n, (0, 1, 2, 3, 4, 5), prime)
        return int(mod_einsum("pqrstm,mk,pqrxyz,stkxzy->",
                              [n, mm, up, up], prime) % prime)

    assert broken(form) != broken(boosted), (
        "the historical P10_07 placement is boost invariant after all, which "
        "would invalidate the recorded diagnosis")


def test_p10_05_and_p10_08_black_brackets_are_not_vacuous():
    """Dropping a BLACK antisymmetrisation must change the value.

    A bracket program that silently does nothing would let a wrong formula pass
    every other test in this file: homogeneity, boost invariance and atlas
    membership are all insensitive to whether the antisymmetriser fired. These
    two candidates are the ones whose brackets are applied explicitly rather
    than being supplied by `composite_n1050`, so they are exactly the ones
    where a no-op would go unnoticed.
    """
    from sdinv.index_symmetry_ops import BLACK, BracketOp, BracketProgram
    from sdinv.published_degree10_invariants import _outer
    from sdinv.stress import _raise_axes, composite_n1050, five_form_moment

    prime = P
    form = _sample(prime, 5)
    _, mixed = five_form_moment(form, prime, "optimized")
    mm = (mixed @ mixed) % prime
    n = composite_n1050(form, prime, "optimized")

    # --- P10_05 with and without its two black brackets -------------------
    n_up = _raise_axes(n, (0, 1, 2, 3, 4, 5), prime)
    raw5 = _outer(mm, mixed, prime)
    anti5 = BracketProgram(ops=[
        BracketOp("antisym", (0, 2), BLACK, True, "nu"),
        BracketOp("antisym", (1, 3), BLACK, True, "mu")]).apply(raw5, prime)
    spec = "nmvp,abcdmp,abcdnv->"
    with_b = int(mod_einsum(spec, [anti5, n, n_up], prime) % prime)
    without = int(mod_einsum(spec, [raw5, n, n_up], prime) % prime)
    assert with_b != without, (
        "P10_05 is unchanged by its black antisymmetrisations; the bracket "
        "program is a no-op here and the formula is not being tested")
    assert with_b == PUBLISHED_DEGREE10["P10_05"]["evaluator"](form, prime)

    # --- P10_08 with and without its black bracket ------------------------
    outer = _raise_axes(n, (4, 5), prime)
    n_a = _raise_axes(n, (0, 1, 2), prime)
    n_b = _raise_axes(n, (0, 3, 4, 5), prime)
    raw8 = _outer(mixed, mixed, prime)
    anti8 = BracketProgram(ops=[
        BracketOp("antisym", (0, 2), BLACK, True, "nu-mu")]).apply(raw8, prime)
    spec8 = "pqrsnm,ntmk,pqrxyz,stkxzy->"
    with_b8 = int(mod_einsum(spec8, [outer, anti8, n_a, n_b], prime) % prime)
    without8 = int(mod_einsum(spec8, [outer, raw8, n_a, n_b], prime) % prime)
    assert with_b8 == PUBLISHED_DEGREE10["P10_08"]["evaluator"](form, prime)

    # P10_08's black bracket is REDUNDANT, and that is a property of the
    # contraction rather than a defect. The check below states why, which is
    # stronger than asserting the bracket changes the value: the symmetric
    # part of the (nu, mu) pair contracts to exactly zero, so the rest of the
    # contraction already projects onto the antisymmetric part.
    sym8 = BracketProgram(ops=[
        BracketOp("sym", (0, 2), BLACK, True, "nu-mu")]).apply(raw8, prime)
    sym_value = int(mod_einsum(spec8, [outer, sym8, n_a, n_b], prime) % prime)
    assert sym_value == 0, (
        f"the (nu,mu)-symmetric part of P10_08 is {sym_value}, not 0. It was 0 "
        f"at both 32749 and 32719, which is the reason the black bracket is a "
        f"no-op there. If it is nonzero now, the bracket is doing real work "
        f"and the redundancy recorded in the index audit is wrong.")
    assert with_b8 == without8, (
        "P10_08's black bracket changed the value; given a vanishing symmetric "
        "part it must not")

    # ...and the engine is genuinely live on this tensor, so the redundancy is
    # a fact about the contraction and not a silently disabled bracket.
    other = BracketProgram(ops=[
        BracketOp("antisym", (0, 1), BLACK, True, "control")]).apply(raw8, prime)
    assert int(mod_einsum(spec8, [outer, other, n_a, n_b], prime)
               % prime) != without8, (
        "no antisymmetrisation of this tensor changes the result, so the "
        "bracket engine is not acting at all")


def test_p10_05_and_p10_08_index_mutation_is_detected():
    """Permuting one source index must change the invariant.

    Guards against a transcription in which two index labels are interchangeable
    -- if swapping them left the value alone, the einsum would not actually be
    encoding the published contraction pattern.
    """
    from sdinv.published_degree10_invariants import _outer
    from sdinv.stress import _raise_axes, composite_n1050, five_form_moment

    prime = P
    form = _sample(prime, 5)
    _, mixed = five_form_moment(form, prime, "optimized")
    n = composite_n1050(form, prime, "optimized")
    outer = _raise_axes(n, (4, 5), prime)
    n_a = _raise_axes(n, (0, 1, 2), prime)
    n_b = _raise_axes(n, (0, 3, 4, 5), prime)
    from sdinv.index_symmetry_ops import BLACK, BracketOp, BracketProgram
    block = BracketProgram(ops=[
        BracketOp("antisym", (0, 2), BLACK, True, "nu-mu")]).apply(
            _outer(mixed, mixed, prime), prime)

    base = int(mod_einsum("pqrsnm,ntmk,pqrxyz,stkxzy->",
                          [outer, block, n_a, n_b], prime) % prime)
    # the published alpha ordering is a1a2]a3 against a1a3]a2; flattening it to
    # a1a2]a3 / a1a2]a3 is the single most likely transcription slip
    mutated = int(mod_einsum("pqrsnm,ntmk,pqrxyz,stkxyz->",
                             [outer, block, n_a, n_b], prime) % prime)
    assert base != mutated, (
        "P10_08 does not distinguish the published alpha ordering a1a3]a2 from "
        "a1a2]a3, so the transcription of that detail is untested")


def test_archived_trM5_matches_the_M_only_result():
    """P10_01 = tr M^5 must land in the closure, as the M-only test found."""
    if not MAP.exists():
        return
    # MAP is archived pre-audit evidence and is bound to the legacy registry.
    assert all(spec["source_reading"] == "legacy-pre-audit"
               for spec in LEGACY_PUBLISHED_DEGREE10.values())
    with MAP.open() as stream:
        payload = json.load(stream)
    for record in payload["per_prime"].values():
        p1 = record["projections"]["P10_01"]
        assert p1["status"] == "solved"
        assert all(v == 0 for v in p1["quotient_vector"]), (
            "tr M^5 must project to zero in Q10, matching C-MONLY-01 via an "
            "independent code path")


def test_archived_published_projection_artifact_is_internally_sound():
    """Structural invariants that must hold whatever the rank turns out to be.

    Kept separate from the rank itself so that a change in the rank reads as a
    result rather than as a broken test.
    """
    if not MAP.exists():
        return
    # MAP is archived pre-audit evidence and is bound to the legacy registry.
    assert all(spec["source_reading"] == "legacy-pre-audit"
               for spec in LEGACY_PUBLISHED_DEGREE10.values())
    with MAP.open() as stream:
        payload = json.load(stream)
    assert payload["consistent"] is True, (
        "the Q10 rank disagrees between primes; a modular accident or a "
        "prime-dependent bug, and either way not a result")
    for record in payload["per_prime"].values():
        assert record["dim_Q10"] == 3
        for name, proj in record["projections"].items():
            assert proj["status"] == "solved", (
                f"{name} did not solve. A rank assembled from unsolved rows is "
                f"uninformative -- an unsolved candidate contributes 0 for a "
                f"reason that has nothing to do with the quotient. This is "
                f"exactly how the non-scalar P10_07 was caught.")


def test_archived_published_q10_rank_is_full():
    """The twelve published candidates span Q10: rank 3 of 3.

    This expectation was 0 until all twelve were implemented. The five simplest
    candidates -- P10_01 through P10_08 -- all project to the zero vector, and
    while only those were implemented the evidence pointed at the published
    list not reaching the quotient at all. It reaches it completely. The
    candidates that get there are exactly the four hardest to transcribe:
    P10_09 (a red bracket) and P10_10/11/12 (five N^(1050) blocks each).

    The lesson is worth keeping: a null result over a subset chosen for ease of
    implementation is not evidence about the whole set.
    """
    if not MAP.exists():
        return
    # MAP is archived pre-audit evidence and is bound to the legacy registry.
    assert all(spec["source_reading"] == "legacy-pre-audit"
               for spec in LEGACY_PUBLISHED_DEGREE10.values())
    with MAP.open() as stream:
        payload = json.load(stream)
    rank = payload["Q10_rank_from_implemented_published"]
    assert rank == 3, (
        f"published Q10 rank is {rank}, expected 3. If it dropped, a candidate "
        f"regressed; if it is None, the primes disagree.")


def test_archived_only_the_four_hardest_candidates_reach_the_quotient():
    """Pin which candidates carry the quotient, so a silent swap is caught.

    P10_01..P10_08 must project to zero and P10_09..P10_12 must not. If a
    transcription error moved a nonzero image onto a different candidate the
    rank would be unchanged and nothing else in this suite would notice.
    """
    if not MAP.exists():
        return
    # MAP is archived pre-audit evidence and is bound to the legacy registry.
    assert all(spec["source_reading"] == "legacy-pre-audit"
               for spec in LEGACY_PUBLISHED_DEGREE10.values())
    with MAP.open() as stream:
        payload = json.load(stream)
    for prime, record in payload["per_prime"].items():
        p = int(prime)
        for name, proj in record["projections"].items():
            if "[" in name or proj["status"] != "solved":
                continue
            index = int(name.split("_")[1])
            hit = any(v % p for v in proj["quotient_vector"])
            if index <= 8:
                assert not hit, (
                    f"{name} reaches Q10 at prime {prime}; it did not before, "
                    f"and that is a change in the result, not a test failure")
            else:
                assert hit, (
                    f"{name} no longer reaches Q10 at prime {prime}")
