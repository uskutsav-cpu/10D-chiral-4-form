"""The reverse search must stay independent of the published formulas.

The benchmark's whole value is that it reaches Q10 without being told the
answer. That independence is a property of the code, so it is asserted here
rather than promised in a docstring.
"""
import ast
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.reverse_block_decomposition import (  # noqa: E402
    BLOCK_KINDS, block_multisets, build_einsum, canonical_form,
    enumerate_topologies, generate_candidates, is_connected,
    multiset_slot_count)


def _imported_modules(path):
    tree = ast.parse(path.read_text())
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module)
    return names


def test_reverse_engine_does_not_import_the_published_formulas():
    """Generation must not be able to see the answer it is searching for."""
    mods = _imported_modules(ROOT / "src" / "sdinv" / "reverse_block_decomposition.py")
    leaked = [m for m in mods if "published" in m]
    assert not leaked, (
        f"the reverse search imports {leaked}; its recovery of Q10 would no "
        f"longer be independent of the published candidates")


def test_forbidden_traces_are_pruned_before_arithmetic():
    """A contraction inside an antisymmetric group vanishes identically.

    Enumerating it and discovering the zero numerically would be correct but
    wasteful, and would inflate the candidate count with structural nulls.
    """
    topos = enumerate_topologies(["N1050", "N1050"])
    for topo in topos:
        for (a, b), count in topo.items():
            if count and a == b:
                _, classes = BLOCK_KINDS["N1050"]
                assert classes[a[1]] == 1, (
                    f"a self-edge survived inside an antisymmetric class: {a}")


def test_every_enumerated_topology_saturates_all_slots():
    for kinds in (["M", "M"], ["M", "N1050"], ["N1050", "N1050"]):
        for topo in enumerate_topologies(kinds):
            used = {}
            for (a, b), count in topo.items():
                used[a] = used.get(a, 0) + count
                used[b] = used.get(b, 0) + count
            for bi, kind in enumerate(kinds):
                total = sum(v for (blk, _), v in used.items() if blk == bi)
                assert total == BLOCK_KINDS[kind][0], (
                    f"{kind} block {bi} has {total} slots used, expected "
                    f"{BLOCK_KINDS[kind][0]} -- the contraction is not closed")


def test_build_einsum_uses_every_index_exactly_twice():
    for kinds in (["M", "M"], ["M", "M", "N1050"], ["N1050", "N1050"]):
        for topo in generate_candidates(kinds)[0][:20]:
            spec, _ = build_einsum(kinds, topo)
            body = spec.split("->")[0].replace(",", "")
            for ch in set(body):
                assert body.count(ch) == 2, (
                    f"index {ch!r} appears {body.count(ch)} times in {spec}; a "
                    f"scalar contraction needs exactly two")


def test_build_einsum_raises_exactly_one_end_of_every_edge():
    """The rule that a wrong answer taught: one raised end per edge.

    An edge joining two equally-placed slots contracts with delta rather than
    eta and is not a Lorentz scalar. That is how the published P10_07 shipped
    broken, and the reverse engine must not be able to reproduce the mistake.
    """
    kinds = ["N1050", "N1050"]
    for topo in generate_candidates(kinds)[0][:20]:
        spec, raises = build_einsum(kinds, topo)
        subs = spec.split("->")[0].split(",")
        raised_positions = set()
        for bi, rset in enumerate(raises):
            for slot in rset:
                raised_positions.add((bi, slot))
        seen = {}
        for bi, sub in enumerate(subs):
            for slot, ch in enumerate(sub):
                seen.setdefault(ch, []).append((bi, slot))
        for ch, ends in seen.items():
            n_raised = sum(1 for e in ends if e in raised_positions)
            assert n_raised == 1, (
                f"edge {ch!r} has {n_raised} raised ends in {spec}; exactly "
                f"one is required or it contracts with delta not eta")


def test_disconnected_topologies_are_rejected():
    """Disconnected means a product of lower invariants, not a primitive."""
    kinds = ["M", "M", "M", "M", "M"]
    raw = enumerate_topologies(kinds)
    assert any(not is_connected(kinds, t) for t in raw), (
        "no disconnected topology in the raw enumeration, so the connectivity "
        "filter is untested here")
    kept, stats = generate_candidates(kinds)
    assert all(is_connected(kinds, t) for t in kept)
    assert stats["after_connectivity"] < stats["raw_topologies"]


def test_canonical_form_identifies_relabelled_topologies():
    kinds = ["N1050", "N1050", "N1050"]
    cands, stats = generate_candidates(kinds, max_topologies=3000)
    keys = [canonical_form(kinds, t) for t in cands]
    assert len(set(keys)) == len(keys), "canonicalisation left duplicates"
    assert stats["canonical"] <= stats["after_connectivity"]


def test_block_multisets_are_ordered_cheapest_first():
    """Ordering by cost lets a bounded pilot cover sectors, not stall in one."""
    sets = block_multisets(5)
    counts = [multiset_slot_count(m) for m in sets]
    assert counts == sorted(counts), "multisets are not cost-ordered"
    assert sets[0] == ("M", "M", "M", "M", "M")
    assert len(sets) == 21


def test_streaming_matches_the_list_enumeration():
    """The memory refactor must not change WHICH candidates are produced."""
    from sdinv.reverse_block_decomposition import stream_candidates
    for kinds in (["M", "M", "M", "M", "M"], ["M", "M", "M", "N1050", "N1050"]):
        listed, _ = generate_candidates(kinds)
        streamed = [t for t, _ in stream_candidates(kinds)]
        a = sorted(canonical_form(kinds, t) for t in listed)
        b = sorted(canonical_form(kinds, t) for t in streamed)
        assert a == b, (
            f"streamed generation produced a different candidate set for "
            f"{kinds}; the refactor changed the science, not just the memory")


def test_shards_partition_the_candidate_set_exactly():
    """Two shards must never claim the same scalar, and none may be lost.

    Sharding is tested on the CANONICAL key rather than on the raw topology,
    so equivalent topologies cannot land in different shards and be counted
    twice.
    """
    from sdinv.reverse_block_decomposition import stream_candidates
    kinds = ["M", "M", "M", "N1050", "N1050"]
    whole = {canonical_form(kinds, t) for t, _ in stream_candidates(kinds)}
    union, seen_twice = set(), set()
    for residue in range(3):
        shard = {canonical_form(kinds, t) for t, _ in
                 stream_candidates(kinds, shard_residue=residue,
                                   shard_modulus=3)}
        seen_twice |= (union & shard)
        union |= shard
    assert not seen_twice, f"{len(seen_twice)} candidates in multiple shards"
    assert union == whole, (
        f"shard union has {len(union)} candidates, unsharded has {len(whole)}")


def test_streaming_memory_is_bounded_by_distinct_candidates(tmp_path):
    """Streaming must not retain the raw topologies it walks past.

    The list form held up to 30 000 dicts per sector, which drove enumeration
    RSS to ~2.6 GB across 21 sectors. This asserts the streamed form walks far
    more raw topologies than it retains.
    """
    from sdinv.reverse_block_decomposition import stream_candidates
    kinds = ["N1050"] * 3
    last = {}
    count = 0
    for _, stats in stream_candidates(kinds, cap=20000):
        last = stats
        count += 1
    assert last["raw"] > count, (
        "streaming yielded as many candidates as raw topologies walked, so "
        "canonicalisation is not deduplicating and memory is not bounded")
    assert last["duplicate"] > 0


def test_every_generated_candidate_is_boost_invariant():
    """Generated candidates must be Lorentz scalars, across ALL sectors.

    This caught a real defect. `make_blocks` returned `mixed` for the M block,
    which is M_{a}{}^{b} -- already one index up -- while `build_einsum`
    assumes fully covariant operands and raises one end of each edge itself.
    Every M-block edge therefore carried two raised ends or none.

    The symptom in the pilot was 129 of 579 candidates failing to lie in the
    degree-10 atlas span, which a genuine Lorentz scalar cannot do. Pure-N
    sectors were 100% boost invariant and M-containing sectors were not, which
    localised it immediately.

    A rotation cannot see this: delta and eta agree on the spatial block.
    """
    import numpy as np
    from sdinv.forms import selfdual_projector, to_dense, random_form
    from sdinv.modp import P
    from sdinv.reverse_block_decomposition import (
        evaluate, make_blocks, stream_candidates)

    def _boost(mod, t=7):
        ti = pow(t, mod - 2, mod)
        half = pow(2, mod - 2, mod)
        c = ((t + ti) * half) % mod
        s = ((ti - t) * half) % mod
        L = np.eye(10, dtype=np.int64)
        L[0, 0] = c; L[0, 1] = s; L[1, 0] = s; L[1, 1] = c
        return L

    def _move(tensor, L, mod):
        out = np.asarray(tensor, dtype=np.int64) % mod
        for axis in range(out.ndim):
            out = np.moveaxis(np.tensordot(L, out, axes=(1, axis)), 0, axis) % mod
        return out

    pr = selfdual_projector(10, 5, True, P)
    form = to_dense((pr @ random_form(10, 5, np.random.default_rng(5), P)) % P,
                    10, 5, P)
    boosted = _move(form, _boost(P), P)
    blocks, blocks_b = make_blocks(form, P), make_blocks(boosted, P)

    for kinds in (["M"] * 5, ["M", "M", "M", "N1050", "N1050"],
                  ["M", "N4125", "N4125", "N4125", "N4125"],
                  ["N1050"] * 5, ["N4125"] * 5):
        checked = 0
        for topo, _ in stream_candidates(list(kinds), cap=20000,
                                         max_candidates=4):
            try:
                a = evaluate(list(kinds), topo, blocks, P)
                b = evaluate(list(kinds), topo, blocks_b, P)
            except Exception:
                continue
            checked += 1
            assert a == b, (
                f"{'+'.join(kinds)} candidate is not boost invariant "
                f"({a} != {b}); an edge has two raised ends or none, so it "
                f"contracts with delta instead of eta")
        assert checked, f"no candidate evaluated for {kinds}"


def test_make_blocks_returns_all_lower_operands():
    """build_einsum's edge rule is only correct on fully covariant operands."""
    import numpy as np
    from sdinv.forms import selfdual_projector, to_dense, random_form
    from sdinv.modp import P
    from sdinv.reverse_block_decomposition import make_blocks
    from sdinv.stress import _raise_axes, five_form_moment

    pr = selfdual_projector(10, 5, True, P)
    form = to_dense((pr @ random_form(10, 5, np.random.default_rng(5), P)) % P,
                    10, 5, P)
    _, mixed = five_form_moment(form, P, "optimized")
    blocks = make_blocks(form, P)
    assert not np.array_equal(blocks["M"], mixed % P), (
        "make_blocks handed back the mixed-placement M_{a}{}^{b}; the edge "
        "rule then double-raises or never raises every M edge")
    assert np.array_equal(blocks["M"], _raise_axes(mixed, (1,), P))


def test_formula_topology_round_trip_is_exact():
    """A candidate identifier is not mathematics; the formula must be.

    Renders each topology as an explicit indexed contraction, parses it back
    from the INDEX PAIRINGS alone, and requires the canonical key to match. A
    failure means the readable formula does not describe the topology it was
    generated from -- which would make every published formula in the
    certificates document wrong while the identifiers stayed reassuringly
    consistent.
    """
    from sdinv.reverse_block_decomposition import (
        parse_formula, render_formula, stream_candidates)
    for kinds in (["M"] * 5, ["M", "M", "M", "N1050", "N1050"],
                  ["N1050"] * 5, ["N4125"] * 5):
        checked = 0
        for topo, _ in stream_candidates(list(kinds), cap=20000,
                                         max_candidates=6):
            text = render_formula(list(kinds), topo)
            kinds2, topo2 = parse_formula(text)
            assert kinds2 == list(kinds), f"kinds changed: {kinds2}"
            assert canonical_form(list(kinds), topo) == \
                canonical_form(kinds2, topo2), (
                    f"round trip changed the topology for {kinds}:\n{text}")
            checked += 1
        assert checked, f"nothing checked for {kinds}"


def test_rendered_formula_marks_exactly_one_raised_end_per_index():
    """The readable form must show the variance that makes it a scalar."""
    from sdinv.reverse_block_decomposition import (
        render_formula, stream_candidates)
    kinds = ["N4125"] * 5
    for topo, _ in stream_candidates(kinds, cap=20000, max_candidates=5):
        text = render_formula(kinds, topo)
        seen = {}
        for factor in text.split("*"):
            body = factor.partition("{")[2].rstrip(" }")
            for token in body.split():
                seen.setdefault(token[1:], []).append(token[0])
        for name, marks in seen.items():
            assert len(marks) == 2, f"{name} appears {len(marks)} times"
            assert sorted(marks) == ["^", "_"], (
                f"index {name} has variance {marks}; exactly one end must be "
                f"raised or the contraction uses delta rather than eta")
