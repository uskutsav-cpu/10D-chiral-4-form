"""Correctness gates. Run these before trusting any 10D number."""
import hashlib
import importlib.util
import json
import sys, os, numpy as np, pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import sdinv.modp as modp
import sdinv.graphs as graphs_module
from sdinv.modp import P, ALT_P, mod_einsum, RankSieve
from sdinv.forms import (to_dense, random_form, metric_signs, check_star_squared,
                         selfdual_projector)
from sdinv.graphs import (canonical, enumerate_graphs, graph_from_label,
                          graph_from_record, graph_label, graph_to_record,
                          load_graph_catalog, validate_graph)
from sdinv.contract import (_slot_plan, _signed, value, jacobian_row,
                            jacobian_row_amputated, build_basis_flat,
                            build_compact_derivative_basis,
                            contraction_plan_cost,
                            contraction_plan_profile,
                            greedy_contraction_plan_cost,
                            planned_value,
                            value_and_jacobian_row)
from sdinv.catalog import (candidate_record, iter_graph_shard,
                           canonical_graph_id, write_graph_shard)
from sdinv.checkpoint import load_checkpoint, write_checkpoint
from sdinv.spinor_adapter import compare_column_spaces

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CATALOG = os.path.join(ROOT, "results", "10d_graph_catalog.json")
ORDER8_RESULT = os.path.join(ROOT, "results", "10d_order8.json")
ORDER10_RESULT = os.path.join(ROOT, "results", "10d_order10.json")
PIPELINE_SPEC = importlib.util.spec_from_file_location(
    "degree10_pipeline",
    os.path.join(ROOT, "scripts", "degree10_pipeline.py"),
)
degree10_pipeline = importlib.util.module_from_spec(PIPELINE_SPEC)
PIPELINE_SPEC.loader.exec_module(degree10_pipeline)


def _ops(M, Fd, PD, mod):
    slots, tails = _slot_plan(M, PD)
    s = metric_signs(Fd.shape[0], True) % mod
    sub = ",".join("".join(x) for x in slots) + "->"
    return sub, [_signed(Fd, tails[v], s, mod) for v in range(M.shape[0])]


def _transform_covariant_tensor(tensor, transformation, mod):
    """Apply one matrix per covariant index, reducing after every axis."""
    result = np.asarray(tensor, dtype=np.int64) % mod
    for axis in range(result.ndim):
        result = np.tensordot(
            transformation, result, axes=(1, axis))
        result = np.moveaxis(result, 0, axis) % mod
    return result


@pytest.mark.parametrize("n", [2, 4, 6])
def test_mod_einsum_matches_exact_bigint(n):
    """THE critical test. int64 overflow in einsum silently produces wrong
    numbers that look completely plausible. This caught three real bugs."""
    D, PD, mod = 6, 3, P
    Fd = to_dense(random_form(D, PD, np.random.default_rng(7), mod), D, PD, mod)
    for M in enumerate_graphs(n, PD):
        sub, ops = _ops(M, Fd, PD, mod)
        exact = int(np.einsum(sub, *[o.astype(object) for o in ops])) % mod
        assert int(mod_einsum(sub, ops, mod)) == exact


def test_float_blas_path_is_exact(monkeypatch):
    """The accelerated float64 branch must still be exact finite-field math."""
    rng = np.random.default_rng(17)
    a = rng.integers(0, P, size=(5, 7, 3), dtype=np.int64)
    b = rng.integers(0, P, size=(3, 7, 4), dtype=np.int64)
    subscripts = "abc,cbd->ad"
    exact = np.einsum(
        subscripts, a.astype(object), b.astype(object)) % P
    monkeypatch.setattr(modp, "FLOAT_BLAS_MIN_WORK", 0)
    got = mod_einsum(subscripts, [a, b], P)
    assert np.array_equal(got, exact.astype(np.int64))


def test_reverse_jacobian_matches_amputated_oracle():
    D, PD, mod = 6, 3, P
    Fd = to_dense(
        random_form(D, PD, np.random.default_rng(23), mod), D, PD, mod)
    basis = build_basis_flat(D, PD, None, mod)
    for n in (2, 4, 6):
        M = enumerate_graphs(n, PD)[0]
        fast = jacobian_row(M, Fd, basis, D, PD, True, mod)
        oracle = jacobian_row_amputated(M, Fd, basis, D, PD, True, mod)
        assert np.array_equal(fast, oracle)


def test_combined_value_jacobian_matches_independent_oracle():
    D, PD, mod = 6, 3, P
    Fd = to_dense(
        random_form(D, PD, np.random.default_rng(29), mod), D, PD, mod)
    basis = build_basis_flat(D, PD, None, mod)
    for n in (2, 4, 6):
        M = enumerate_graphs(n, PD)[0]
        fast_value, fast_row = value_and_jacobian_row(
            M, Fd, basis, D, PD, True, mod, backend="optimized")
        reference_value, reference_row = value_and_jacobian_row(
            M, Fd, basis, D, PD, True, mod, backend="reference")
        assert fast_value == reference_value == value(
            M, Fd, D, PD, True, mod)
        assert planned_value(M, Fd, D, PD, True, mod) == fast_value
        assert np.array_equal(fast_row, reference_row)


def test_compact_derivative_projection_matches_dense_basis():
    D, PD, mod = 6, 3, P
    Fd = to_dense(
        random_form(D, PD, np.random.default_rng(31), mod), D, PD, mod)
    dense_basis = build_basis_flat(D, PD, None, mod)
    compact_basis = build_compact_derivative_basis(
        D, PD, None, mod, independent=False)
    for M in enumerate_graphs(4, PD):
        dense_value, dense_row = value_and_jacobian_row(
            M, Fd, dense_basis, D, PD, True, mod)
        compact_value, compact_row = value_and_jacobian_row(
            M, Fd, compact_basis, D, PD, True, mod)
        assert dense_value == compact_value
        assert np.array_equal(dense_row, compact_row)


def test_global_contraction_plan_dominates_greedy_schedule_cost():
    for M in enumerate_graphs(6, 3):
        exact = contraction_plan_cost(M, 6, 3)
        greedy = greedy_contraction_plan_cost(M, 6, 3)
        assert exact <= greedy


def test_order12_graph_label_is_unambiguous_and_round_trips():
    M = np.zeros((12, 12), dtype=np.int64)
    for i, j, multiplicity in ((0, 4, 3), (0, 11, 2), (4, 10, 2)):
        M[i, j] = M[j, i] = multiplicity
    label = graph_label(M)
    assert label == "n12[0-4^3,0-11^2,4-10^2]"
    assert np.array_equal(graph_from_label(label), M)

    legacy = "n4[01^3,02^2,13^2,23^3]"
    assert graph_label(graph_from_label(legacy)) == legacy
    assert np.array_equal(
        graph_from_label("n4[0-1^3,0-2^2,1-3^2,2-3^3]"),
        graph_from_label(legacy),
    )


def test_contraction_memory_guard_rejects_before_execution():
    M = enumerate_graphs(4, 3)[0]
    profile = contraction_plan_profile(M, 6, 3)
    assert profile["largest_pair_work"] > 0
    assert profile["estimated_peak_bytes"] >= profile[
        "retained_forward_reverse_bytes"]
    Fd = to_dense(
        random_form(6, 3, np.random.default_rng(41), P), 6, 3, P)
    basis = build_basis_flat(6, 3, None, P)
    with pytest.raises(MemoryError, match="rejected before signed operands"):
        value_and_jacobian_row(
            M, Fd, basis, 6, 3, True, P, max_memory_bytes=1)


def test_rank_sieve_checkpoint_round_trip_is_deterministic(tmp_path):
    rows = [
        [0, 0, 0, 2, 0],
        [0, 3, 0, 1, 0],  # deliberately introduces a lower pivot later
        [5, 0, 0, 0, 0],
        [5, 3, 0, 3, 0],
        [0, 0, 7, 0, 0],
    ]
    uninterrupted = RankSieve(5, P)
    uninterrupted_results = [uninterrupted.add(row) for row in rows]

    partial = RankSieve(5, P)
    partial_results = [partial.add(row) for row in rows[:2]]
    identity = {"prime": P, "seed": 17, "catalog": "fixture"}
    checkpoint = tmp_path / "run.checkpoint.json"
    write_checkpoint(
        checkpoint, identity, {"sieve": partial.to_state(), "cursor": 2})
    restored_state = load_checkpoint(checkpoint, identity)
    restored = RankSieve.from_state(restored_state["sieve"])
    resumed_results = partial_results + [
        restored.add(row) for row in rows[restored_state["cursor"]:]
    ]

    assert resumed_results == uninterrupted_results
    assert restored.pivots == uninterrupted.pivots
    assert all(np.array_equal(a, b)
               for a, b in zip(restored.rows, uninterrupted.rows))
    with pytest.raises(ValueError, match="identity mismatch"):
        load_checkpoint(checkpoint, {**identity, "seed": 18})

    envelope = json.loads(checkpoint.read_text())
    envelope["payload"]["state"]["cursor"] = 999
    checkpoint.write_text(json.dumps(envelope))
    with pytest.raises(ValueError, match="checksum mismatch"):
        load_checkpoint(checkpoint, identity)


def test_discovery_candidate_staging_does_not_mutate_on_interruption():
    sieve = RankSieve(3, P)
    state = {
        "evaluated": 0,
        "cursor": {"id": ""},
        "generators": [],
    }

    def interrupted_generator(_rank):
        raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        degree10_pipeline._stage_discovery_candidate(
            state,
            sieve,
            np.array([1, 0, 0], dtype=np.int64),
            {"id": "candidate"},
            interrupted_generator,
        )
    assert sieve.rank == 0
    assert state["evaluated"] == 0
    assert state["cursor"] == {"id": ""}
    assert state["generators"] == []


def test_generation_executable_identity_follows_path(tmp_path, monkeypatch):
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()
    for directory, payload in ((first, "first"), (second, "second")):
        executable = directory / "geng"
        executable.write_text(f"#!/bin/sh\n# {payload}\n")
        executable.chmod(0o755)

    monkeypatch.setenv("PATH", str(first))
    first_identity = degree10_pipeline._executable_identity("geng")
    monkeypatch.setenv("PATH", str(second))
    second_identity = degree10_pipeline._executable_identity("geng")
    assert first_identity["path"] == str((first / "geng").resolve())
    assert second_identity["path"] == str((second / "geng").resolve())
    assert first_identity["sha256"] != second_identity["sha256"]


def test_streaming_catalog_round_trip_and_checksum(tmp_path):
    graphs = enumerate_graphs(4, 3)
    shard_a = tmp_path / "a.jsonl.gz"
    shard_b = tmp_path / "b.jsonl.gz"
    generation = {"order": 4, "residue": 0, "modulus": 1}
    manifest_a = write_graph_shard(shard_a, iter(graphs), generation)
    manifest_b = write_graph_shard(shard_b, iter(graphs), generation)

    loaded = list(iter_graph_shard(shard_a))
    assert manifest_a["count"] == manifest_b["count"] == len(graphs)
    assert manifest_a["logical_sha256"] == manifest_b["logical_sha256"]
    assert shard_a.read_bytes() == shard_b.read_bytes()
    assert [record["id"] for record, _ in loaded] == [
        candidate_record(M)["id"] for M in graphs]
    for (_, got), expected in zip(loaded, graphs):
        assert np.array_equal(got, expected)

    manifest_path = shard_a.with_suffix(shard_a.suffix + ".manifest.json")
    damaged = json.loads(manifest_path.read_text())
    damaged["logical_sha256"] = "0" * 64
    manifest_path.write_text(json.dumps(damaged))
    with pytest.raises(ValueError, match="logical SHA-256 mismatch"):
        list(iter_graph_shard(shard_a))


def test_catalog_id_supports_exact_small_graph_fallback(monkeypatch):
    M = enumerate_graphs(4, 3)[0]
    monkeypatch.setattr(graphs_module, "pynauty", None)
    graph_id = canonical_graph_id(M)
    assert len(graph_id) == 64
    assert graph_id == canonical_graph_id(M)


def test_nauty_stream_reaps_producer_when_multig_is_missing(tmp_path):
    producer = tmp_path / "producer.sh"
    producer.write_text("#!/bin/sh\nprintf 'header\\n'\n")
    producer.chmod(0o755)
    missing = tmp_path / "missing-multig"
    stream = graphs_module.iter_graphs_nauty(
        4, 3, 3, geng=str(producer), multig=str(missing))
    with pytest.raises(FileNotFoundError):
        next(stream)


def test_order10_graph_record_is_unambiguous():
    M = np.zeros((10, 10), dtype=np.int64)
    for i in range(5):
        M[i, i + 5] = M[i + 5, i] = 1
    record = graph_to_record(M)
    assert record["order"] == 10
    assert len(record["upper_triangle"]) == 45
    assert np.array_equal(graph_from_record(record), M)


def test_saved_graph_id_validation_is_catalog_independent():
    M = enumerate_graphs(4, 3)[0]
    item = {
        "id": "fixture",
        "graph": graph_label(M),
        "graph_id": canonical_graph_id(M),
    }
    degree10_pipeline._validate_saved_graph_ids([item])
    with pytest.raises(ValueError, match="canonical ID mismatch"):
        degree10_pipeline._validate_saved_graph_ids([
            {**item, "graph_id": "0" * 64},
        ])


def test_benchmark_hardware_metadata_is_detected(monkeypatch):
    monkeypatch.setattr(degree10_pipeline.platform, "machine", lambda: "arch")
    monkeypatch.setattr(
        degree10_pipeline.platform, "processor", lambda: "processor")
    monkeypatch.setattr(
        degree10_pipeline.platform, "platform", lambda: "TestOS")
    monkeypatch.setattr(degree10_pipeline.os, "cpu_count", lambda: 7)
    monkeypatch.setattr(
        degree10_pipeline, "_physical_memory_bytes", lambda: 12345)
    metadata = degree10_pipeline._hardware_metadata()
    assert metadata["architecture"] == "arch"
    assert metadata["processor"] == "processor"
    assert metadata["logical_cpu_count"] == 7
    assert metadata["memory_bytes"] == 12345
    assert metadata["operating_system"] == "TestOS"


def test_spinor_adapter_compares_column_spaces_not_column_names():
    trace = np.array([
        [1, 0],
        [0, 1],
        [1, 1],
        [2, 3],
    ], dtype=np.int64)
    change_of_basis = np.array([[2, 1], [1, 1]], dtype=np.int64)
    spinor = (trace @ change_of_basis) % P
    report = compare_column_spaces(trace, spinor, P)
    assert report["equal_column_spaces"]
    assert report["trace_rank"] == report["spinor_rank"] == 2

    outside = np.column_stack((spinor, [0, 0, 0, 1]))
    report = compare_column_spaces(trace, outside, P)
    assert not report["equal_column_spaces"]
    assert report["union_rank"] == 3


def test_contractions_are_lorentz_invariant():
    """Rotation by the 3-4-5 triple, exact over F_p. If a contraction is not
    invariant, the metric orientation is wrong and rank bounds are void."""
    D, PD, mod = 6, 3, P
    inv5 = pow(5, mod - 2, mod)
    c, s = (3 * inv5) % mod, (4 * inv5) % mod
    R = np.eye(D, dtype=np.int64)
    R[1, 1] = c; R[1, 2] = s; R[2, 1] = (-s) % mod; R[2, 2] = c
    Fd = to_dense(random_form(D, PD, np.random.default_rng(7), mod), D, PD, mod)
    Fr = _transform_covariant_tensor(Fd, R, mod)
    for n in [2, 4, 6]:
        for M in enumerate_graphs(n, PD):
            assert value(M, Fd, D, PD, True, mod) == value(M, Fr, D, PD, True, mod)


def test_6d_reproduces_five_invariants():
    """Elamaran-Ferko-Scarlett arXiv:2512.23750: 5 invariants, 1,2,1,1."""
    D, PD, mod = 6, 3, P
    Fd = to_dense(random_form(D, PD, np.random.default_rng(11), mod), D, PD, mod)
    B = build_basis_flat(D, PD, None, mod)
    sieve, pattern = RankSieve(B.shape[0], mod), []
    for n in [2, 4, 6, 8]:
        before = sieve.rank
        for M in enumerate_graphs(n, PD):
            sieve.add(jacobian_row(M, Fd, B, D, PD, True, mod))
        pattern.append(sieve.rank - before)
    assert sieve.rank == 5 and pattern == [1, 2, 1, 1]


def test_10d_selfduality_wellposed():
    assert check_star_squared(10, 5, True, P) == 1
    Pr = selfdual_projector(10, 5, True, P)
    assert np.array_equal((Pr @ Pr) % P, Pr % P)
    s = RankSieve(Pr.shape[1], P)
    for r in Pr:
        s.add(r)
    assert s.rank == 126


def test_10d_quadratic_invariant_vanishes():
    """F ^ F = 0 for odd-degree forms; self-duality gives F ^ *F = F ^ F,
    so F.F = 0 identically. Free correctness check on the projector."""
    Pr = selfdual_projector(10, 5, True, P)
    Fv = (Pr @ random_form(10, 5, np.random.default_rng(1), P)) % P
    Fd = to_dense(Fv, 10, 5, P)
    M = enumerate_graphs(2, 5, max_mult=5)[0]
    assert value(M, Fd, 10, 5, True, P) == 0


def test_exact_catalog_is_complete_and_collision_free():
    """The committed nauty catalog is exact at orders 4, 6, and 8."""
    with open(CATALOG) as stream:
        raw = json.load(stream)
    assert raw["generator"]["software"] == "nauty gtools 2.9.3"
    assert {n: raw["orders"][n]["count"] for n in ("4", "6", "8")} == {
        "4": 4,
        "6": 49,
        "8": 1689,
    }

    for order in (4, 6, 8):
        graphs = load_graph_catalog(CATALOG, order)
        certificates = set()
        for M in graphs:
            validate_graph(M, valence=5, max_mult=4)
            certificates.add(canonical(M))
        assert len(certificates) == len(graphs), (
            f"duplicate isomorphism class in order-{order} catalog")


def test_order8_basis_is_complete_under_two_primes():
    """Six new octic directions, matching the published Hilbert series.

    Cederwall et al. arXiv:2509.14350v2 give
      P(t) = 1 + t^4 + 2 t^6 + 7 t^8 + ...
    and factor it with (1-t^8)^-6. The seventh degree-8 scalar is I4^2,
    so six independent order-8 Jacobian directions are both necessary and
    sufficient for a complete octic generating set.
    """
    with open(ORDER8_RESULT) as stream:
        result = json.load(stream)
    generators = result["generators"]
    with open(CATALOG, "rb") as stream:
        assert hashlib.sha256(stream.read()).hexdigest() == (
            result["catalog_sha256"])
    assert [g["order"] for g in generators].count(8) == 6
    assert result["literature"]["new_generators"]["8"] == 6
    assert len(result["degree8_basis"]) == 7
    assert result["degree8_basis"][-1] == {
        "id": "I4_1^2",
        "kind": "composite",
        "expression": "I4_1^2",
    }
    for item in generators:
        catalog_graphs = load_graph_catalog(CATALOG, item["order"])
        assert graph_label(catalog_graphs[item["catalog_index"]]) == item["graph"]
    assert all(run["orders"]["8"]["rank"] == 9
               for run in result["runs"].values())

    D, PD = 10, 5
    for prime in (P, ALT_P):
        projector = selfdual_projector(D, PD, True, prime)
        Fd = to_dense(
            (projector @ random_form(
                D, PD, np.random.default_rng(20260727), prime)) % prime,
            D,
            PD,
            prime,
        )
        basis = build_compact_derivative_basis(
            D, PD, projector, prime, independent=True)
        sieve = RankSieve(basis.ncols, prime)
        increments = {4: 0, 6: 0, 8: 0}
        for item in generators:
            M = graph_from_label(item["graph"])
            assert sieve.add(jacobian_row(
                M, Fd, basis, D, PD, True, prime))
            increments[item["order"]] += 1
        assert increments == {4: 1, 6: 2, 8: 6}
        assert sieve.rank == 9


def test_saved_order10_result_has_complete_exact_evidence():
    with open(ORDER10_RESULT) as stream:
        result = json.load(stream)
    generators = result["generators"]
    assert len(generators) == 12
    assert len(result["degree10_basis"]) == 14
    assert result["catalog"]["candidate_count"] == 187392
    assert result["catalog"]["unique_canonical_ids"] == 187392
    assert result["discovery"]["evaluated"] == 132
    assert result["discovery"]["new_degree10"] == 12
    assert result["discovery"]["cumulative_rank"] == 21
    assert [item["id"] for item in result["degree10_basis"][-2:]] == [
        "I4_1*I6_1", "I4_1*I6_2"]

    ids = set()
    for item in generators:
        M = graph_from_label(item["graph"])
        validate_graph(M, valence=5, max_mult=4)
        assert canonical_graph_id(M) == item["graph_id"]
        ids.add(item["graph_id"])
    assert len(ids) == 12

    assert set(result["validation"]) == {str(P), str(ALT_P)}
    for run in result["validation"].values():
        assert len(run["jacobian_samples"]) == 3
        assert all(sample["cumulative_rank"] == 21
                   for sample in run["jacobian_samples"])
        assert run["value_space"]["rank"] == 14
        assert run["value_space"]["dimension"] == 14
        assert len(run["value_space"]["pivot_seeds"]) == 14
    assert result["validation"][str(P)]["optimized_reference_matches"] == [
        f"I10_{k}" for k in range(1, 13)]


def test_order10_basis_recomputes_rank21_at_alt_prime():
    with open(ORDER8_RESULT) as stream:
        lower = json.load(stream)["generators"]
    with open(ORDER10_RESULT) as stream:
        degree10 = json.load(stream)["generators"]

    D, PD, prime, seed = 10, 5, ALT_P, 20260731
    projector = selfdual_projector(D, PD, True, prime)
    Fd = to_dense(
        (projector @ random_form(
            D, PD, np.random.default_rng(seed), prime)) % prime,
        D,
        PD,
        prime,
    )
    basis = build_compact_derivative_basis(
        D, PD, projector, prime, independent=True)
    sieve = RankSieve(basis.ncols, prime)
    for item in lower:
        assert sieve.add(jacobian_row(
            graph_from_label(item["graph"]),
            Fd, basis, D, PD, True, prime))
    assert sieve.rank == 9
    for item in degree10:
        M = graph_from_label(item["graph"])
        scalar, row = value_and_jacobian_row(
            M, Fd, basis, D, PD, True, prime)
        assert scalar != 0
        assert sieve.add(row)
    assert sieve.rank == 21


def test_10d_contractions_survive_a_lorentz_boost():
    """A ROTATION cannot catch a wrong metric placement; a BOOST can.

    If raised/lowered indices are assigned per-tensor rather than per-edge,
    an edge joining two same-placement vertices contracts with delta instead
    of eta. Under a pure rotation delta and eta agree on the spatial block,
    so the error hides. A boost mixes the timelike direction and exposes it.

    Over F_p a boost in the 0-1 plane is any (c, s) with c^2 - s^2 = 1:
    take c = (t + 1/t)/2, s = (1/t - t)/2 for invertible t.
    """
    D, PD, mod = 10, 5, P
    t = 7
    ti = pow(t, mod - 2, mod)
    half = pow(2, mod - 2, mod)
    c = ((t + ti) * half) % mod
    s = ((ti - t) * half) % mod
    assert (c * c - s * s) % mod == 1, "not a hyperbolic rotation"

    L = np.eye(D, dtype=np.int64)
    L[0, 0] = c; L[0, 1] = s; L[1, 0] = s; L[1, 1] = c

    eta = np.diag(metric_signs(D, True)).astype(np.int64) % mod
    assert np.array_equal((L.T @ eta @ L) % mod, eta % mod), "L is not in SO(1,9)"

    Pr = selfdual_projector(D, PD, True, mod)
    Fd = to_dense((Pr @ random_form(D, PD, np.random.default_rng(5), mod)) % mod,
                  D, PD, mod)
    Fr = _transform_covariant_tensor(Fd, L, mod)

    for M in enumerate_graphs(4, PD, max_mult=PD - 1):
        a = value(M, Fd, D, PD, True, mod)
        b = value(M, Fr, D, PD, True, mod)
        assert a == b, f"{graph_label(M)} is not boost invariant: {a} != {b}"


def test_order4_is_exactly_one_invariant():
    """All four order-4 graphs are the SAME invariant, rescaled.

    Order 4 in 10D admits exactly 4 connected valence-5 multigraphs. Their
    values are pairwise proportional with ratios 1, 1/2, 1/4, 1/6, so the
    Jacobian rank is 1 -- there is exactly ONE independent invariant at
    order 4, not two or three. Checked under both primes: a spurious
    proportionality would not survive a change of modulus.
    """
    D, PD = 10, 5
    expected = [1, 2, 4, 6]          # value_k = value_0 / expected[k]

    for mod in (P, ALT_P):
        Pr = selfdual_projector(D, PD, True, mod)
        graphs = enumerate_graphs(4, PD, max_mult=PD - 1)
        assert len(graphs) == 4, f"expected 4 order-4 graphs, got {len(graphs)}"

        seen = None
        for seed in (1, 2, 3):
            Fd = to_dense((Pr @ random_form(D, PD, np.random.default_rng(seed), mod))
                          % mod, D, PD, mod)
            vals = [value(M, Fd, D, PD, True, mod) % mod for M in graphs]
            assert vals[0] != 0
            ratios = [(v * pow(vals[0], mod - 2, mod)) % mod for v in vals]
            if seen is None:
                seen = ratios
            assert ratios == seen, "ratios drifted between random points"

        for r, d in zip(seen, expected):
            assert (r * d) % mod == 1, f"ratio {r} is not 1/{d} mod {mod}"

        sieve = RankSieve(build_basis_flat(D, PD, Pr, mod).shape[0], mod)
        basis = build_basis_flat(D, PD, Pr, mod)
        for M in graphs:
            sieve.add(jacobian_row(M, Fd, basis, D, PD, True, mod))
        assert sieve.rank == 1, f"order-4 rank should be 1, got {sieve.rank}"
