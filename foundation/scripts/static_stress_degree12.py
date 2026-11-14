#!/usr/bin/env python3
"""Exact free-stress subalgebra map through five-form degree 12."""

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))

import degree12_pipeline as atlas  # noqa: E402
from sdinv.catalog import atomic_write_json  # noqa: E402
from sdinv.checkpoint import load_checkpoint, write_checkpoint  # noqa: E402
from sdinv.contract import value_and_jacobian_row  # noqa: E402
from sdinv.exactmap import (  # noqa: E402
    STRESS_TARGETS,
    fraction_record,
    rank_mod,
    reconstruct_vector,
    select_standard_complement,
    solve_full_column_rank,
)
from sdinv.graphs import graph_from_label, graph_label  # noqa: E402
from sdinv.invariant_registry import (  # noqa: E402
    load_verified_registry_through_degree12,
)
from sdinv.modp import inv  # noqa: E402
from sdinv.stress import five_form_moment, matrix_trace_power  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SEEDS = (20260729, 20260730, 20260731, 20260732)
DEFAULT_MAX_MEMORY = 2 * 1024 ** 3
FREE_STRESS_DENOMINATOR = 48
DEGREE12_TARGETS = (
    "tr_M6",
    "tr_M4*tr_M2",
    "tr_M3^2",
    "tr_M2^3",
)


def _file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _engine_sha256():
    digest = hashlib.sha256()
    for relative in (
        "scripts/static_stress_degree12.py",
        "src/sdinv/exactmap.py",
        "src/sdinv/invariant_registry.py",
        "src/sdinv/stress.py",
    ):
        digest.update(relative.encode())
        digest.update((ROOT / relative).read_bytes())
    digest.update(atlas._engine_sha256().encode())
    return digest.hexdigest()


def _trace_m6_graph():
    """The alternating six-M ring representing ``Tr(M^6)`` up to sign."""
    matrix = np.zeros((12, 12), dtype=np.int64)
    for index in range(6):
        left = 2 * index
        right = left + 1
        following = (left + 2) % 12
        matrix[left, right] = matrix[right, left] = 4
        matrix[right, following] = matrix[following, right] = 1
    return matrix


def _fraction_matrix(records):
    return [
        [Fraction(item["numerator"], item["denominator"]) for item in row]
        for row in records
    ]


def _standard_row(width, index, coefficient):
    return [
        Fraction(coefficient if column == index else 0)
        for column in range(width)
    ]


def _residue_vector(vector, prime):
    return [
        value.numerator * inv(value.denominator, prime) % prime
        for value in map(Fraction, vector)
    ]


def _evaluate_sample(prime, seed, basis, generators, lower,
                     max_memory_bytes):
    projector, derivative_basis, samples = atlas._contexts(prime, (seed,))
    del projector
    five_form = samples[0]
    lower_values, lower_rows, _ = atlas._evaluate_graphs(
        lower,
        five_form,
        derivative_basis,
        prime,
        max_memory_bytes,
    )
    product_rows = atlas._product_rows(lower_values, lower_rows, prime)
    basis_rows = [row for _, _, row in product_rows]
    basis_values = [value for _, value, _ in product_rows]

    for index, item in enumerate(generators, start=1):
        matrix = graph_from_label(item["graph"])
        scalar, row = value_and_jacobian_row(
            matrix,
            five_form,
            derivative_basis,
            10,
            5,
            True,
            prime,
            backend="optimized",
            max_memory_bytes=max_memory_bytes,
        )
        basis_values.append(int(scalar))
        basis_rows.append(row)
        if index % 10 == 0 or index == len(generators):
            print(
                f"prime {prime} seed {seed}: "
                f"degree-12 graph {index}/{len(generators)}",
                flush=True,
            )

    if len(basis_rows) != len(basis):
        raise AssertionError("degree-12 basis width changed")
    target_matrix = _trace_m6_graph()
    target_scalar, target_row = value_and_jacobian_row(
        target_matrix,
        five_form,
        derivative_basis,
        10,
        5,
        True,
        prime,
        backend="optimized",
        max_memory_bytes=max_memory_bytes,
    )
    _, m_mixed = five_form_moment(five_form, prime)
    direct_scalar = matrix_trace_power(m_mixed, 6, prime)
    if target_scalar == direct_scalar:
        target_sign = 1
    elif target_scalar == -direct_scalar % prime:
        target_sign = -1
    else:
        raise AssertionError("Tr(M^6) graph normalization is not a sign")
    target_row = target_sign * target_row % prime
    target_scalar = target_sign * target_scalar % prime
    if target_scalar != direct_scalar:
        raise AssertionError("Tr(M^6) scalar calibration failed")
    return {
        "seed": int(seed),
        "basis_values": [int(value) for value in basis_values],
        "basis_rows": [
            [int(value) for value in row] for row in basis_rows
        ],
        "tr_M6_value": int(direct_scalar),
        "tr_M6_row": [int(value) for value in target_row],
        "target_graph_sign": target_sign,
    }


def compute_certificate(args):
    prime = int(args.prime)
    seeds = tuple(int(seed) for seed in args.seeds)
    if len(seeds) < 4 or len(set(seeds)) != len(seeds):
        raise ValueError("at least four distinct samples are required")
    degree12_path = ROOT / "results" / "10d_order12.json"
    with degree12_path.open() as stream:
        degree12 = json.load(stream)
    generators = degree12["generators"]
    registry = load_verified_registry_through_degree12(ROOT)
    basis = [item.id for item in registry.basis(12)]
    lower, _ = atlas._lower_inventory(
        ROOT / "results" / "10d_order8.json",
        ROOT / "results" / "10d_order10.json",
    )
    identity = {
        "schema": 1,
        "prime": prime,
        "seeds": list(seeds),
        "basis": basis,
        "engine_sha256": _engine_sha256(),
        "degree12_sha256": _file_sha256(degree12_path),
    }
    checkpoint_path = (
        args.checkpoint
        if args.checkpoint is not None
        else ROOT / "work" / "static-degree12" / f"{prime}.checkpoint.json"
    )
    state = {"samples": []}
    if checkpoint_path.exists():
        state = load_checkpoint(checkpoint_path, identity)
        print(
            f"resuming prime {prime} from {len(state['samples'])} samples",
            flush=True,
        )
    completed_seeds = [sample["seed"] for sample in state["samples"]]
    if completed_seeds != list(seeds[:len(completed_seeds)]):
        raise ValueError("checkpoint sample prefix differs")

    started = time.perf_counter()
    for seed in seeds[len(state["samples"]):]:
        sample = _evaluate_sample(
            prime,
            seed,
            basis,
            generators,
            lower,
            args.max_memory_bytes,
        )
        state["samples"].append(sample)
        write_checkpoint(checkpoint_path, identity, state)
        print(
            f"prime {prime}: checkpointed seed {seed} "
            f"({len(state['samples'])}/{len(seeds)})",
            flush=True,
        )

    stacked_basis = np.asarray([
        np.concatenate([
            np.asarray(sample["basis_rows"][index], dtype=np.int64)
            for sample in state["samples"]
        ])
        for index in range(len(basis))
    ], dtype=np.int64) % prime
    stacked_target = np.concatenate([
        np.asarray(sample["tr_M6_row"], dtype=np.int64)
        for sample in state["samples"]
    ]) % prime
    basis_rank = rank_mod(stacked_basis, prime)
    if basis_rank != 72:
        raise AssertionError(
            f"stacked degree-12 basis rank {basis_rank}, expected 72")
    solution = solve_full_column_rank(
        stacked_basis.T, stacked_target, prime)
    if not np.array_equal(
        solution @ stacked_basis % prime, stacked_target
    ):
        raise AssertionError("Tr(M^6) stacked-gradient fit failed")

    width = len(basis)
    target_rows = [
        solution.tolist(),
        _residue_vector(
            _standard_row(width, basis.index("I4_1*I8_1"), 2),
            prime,
        ),
        _residue_vector(
            _standard_row(
                width, basis.index("I6_1^2"), Fraction(1024, 9)),
            prime,
        ),
        _residue_vector(
            _standard_row(width, basis.index("I4_1^3"), 8),
            prime,
        ),
    ]
    stress_rank = rank_mod(np.asarray(target_rows), prime)
    if stress_rank != 4:
        raise AssertionError(
            f"degree-12 stress rank {stress_rank}, expected 4")
    certificate = {
        "schema": 1,
        "claim": (
            "Exact degree-12 free-stress map on four stacked self-dual "
            "samples over one prime field."
        ),
        "prime": prime,
        "sample_seeds": list(seeds),
        "basis": basis,
        "basis_rank": basis_rank,
        "stacked_columns": int(stacked_basis.shape[1]),
        "target_graph": graph_label(_trace_m6_graph()),
        "target_graph_signs": [
            sample["target_graph_sign"] for sample in state["samples"]
        ],
        "targets": {
            name: [int(value) for value in row]
            for name, row in zip(DEGREE12_TARGETS, target_rows)
        },
        "stress_rank": stress_rank,
        "quotient_dimension": 72 - stress_rank,
        "engine_sha256": identity["engine_sha256"],
        "degree12_sha256": identity["degree12_sha256"],
        "seconds": round(time.perf_counter() - started, 6),
    }
    atomic_write_json(args.out, certificate)
    print(f"wrote {args.out}", flush=True)


def _rational_rank(rows):
    work = [list(map(Fraction, row)) for row in rows]
    if not work:
        return 0
    width = len(work[0])
    pivot = 0
    for column in range(width):
        selected = next(
            (row for row in range(pivot, len(work))
             if work[row][column]),
            None,
        )
        if selected is None:
            continue
        work[pivot], work[selected] = work[selected], work[pivot]
        scale = work[pivot][column]
        work[pivot] = [value / scale for value in work[pivot]]
        for row in range(len(work)):
            if row != pivot and work[row][column]:
                scale = work[row][column]
                work[row] = [
                    left - scale * right
                    for left, right in zip(work[row], work[pivot])
                ]
        pivot += 1
        if pivot == len(work):
            break
    return pivot


def _load_certificate(path):
    with Path(path).open() as stream:
        return json.load(stream)


def assemble(args):
    fit = [_load_certificate(path) for path in args.certificates]
    if len(fit) < 3:
        raise ValueError("at least three fit-prime certificates are required")
    primes = [int(certificate["prime"]) for certificate in fit]
    if len(set(primes)) != len(primes):
        raise ValueError("fit primes must be distinct")
    basis = fit[0]["basis"]
    engine = fit[0]["engine_sha256"]
    degree12_hash = fit[0]["degree12_sha256"]
    for certificate in fit:
        if (
            certificate["basis"] != basis
            or certificate["engine_sha256"] != engine
            or certificate["degree12_sha256"] != degree12_hash
            or certificate["basis_rank"] != 72
            or certificate["stress_rank"] != 4
        ):
            raise ValueError("incompatible degree-12 certificate set")

    modular_tr_m6 = [
        np.asarray(certificate["targets"]["tr_M6"], dtype=np.int64)
        for certificate in fit
    ]
    rational_tr_m6 = list(reconstruct_vector(modular_tr_m6, primes))
    width = len(basis)
    rational_rows = [
        rational_tr_m6,
        _standard_row(width, basis.index("I4_1*I8_1"), 2),
        _standard_row(
            width, basis.index("I6_1^2"), Fraction(1024, 9)),
        _standard_row(width, basis.index("I4_1^3"), 8),
    ]
    if _rational_rank(rational_rows) != 4:
        raise AssertionError("rational degree-12 stress rank is not four")

    holdout = None
    if args.validation_certificate is not None:
        holdout = _load_certificate(args.validation_certificate)
        holdout_prime = int(holdout["prime"])
        if holdout_prime in primes:
            raise ValueError("validation prime must be independent")
        if (
            holdout["basis"] != basis
            or holdout["engine_sha256"] != engine
            or holdout["degree12_sha256"] != degree12_hash
        ):
            raise ValueError("validation certificate is incompatible")
        expected = _residue_vector(rational_tr_m6, holdout_prime)
        if expected != holdout["targets"]["tr_M6"]:
            raise ValueError(
                "rational reconstruction failed the independent prime")

    complement_indices = select_standard_complement(rational_rows)
    if len(complement_indices) != 68:
        raise AssertionError("degree-12 complement dimension is not 68")
    degree12_entry = {
        "degree": 12,
        "full_dimension": 72,
        "stress_dimension": 4,
        "quotient_dimension": 68,
        "basis": basis,
        "stress_basis": [{
            "id": name,
            "coordinates": [
                fraction_record(value) for value in row
            ],
            "physical_free_stress_coordinates": [
                fraction_record(
                    value / FREE_STRESS_DENOMINATOR ** 6)
                for value in row
            ],
        } for name, row in zip(DEGREE12_TARGETS, rational_rows)],
        "complement_basis": [
            {
                "id": basis[index],
                "coordinate_index": index,
            }
            for index in complement_indices
        ],
        "exact_validation": {
            "fit_primes": primes,
            "sample_seeds": fit[0]["sample_seeds"],
            "independent_validation_prime": (
                int(holdout["prime"]) if holdout else None
            ),
            "rational_reconstruction_holdout_passed": holdout is not None,
            "engine_sha256": engine,
            "degree12_sha256": degree12_hash,
        },
    }

    with (ROOT / "results" / "stress_flow_exact_low_degree.json").open() as s:
        lower = json.load(s)
    degrees = {}
    for degree in (4, 6, 8, 10):
        degree_map = lower["degree_maps"][str(degree)]
        target_names = STRESS_TARGETS[degree]
        rows = [
            [
                Fraction(item["numerator"], item["denominator"])
                for item in degree_map["targets"][target]
            ]
            for target in target_names
        ]
        complement = select_standard_complement(rows)
        powers = degree // 2
        degrees[str(degree)] = {
            "degree": degree,
            "full_dimension": degree_map["basis_dimension"],
            "stress_dimension": len(rows),
            "quotient_dimension": (
                degree_map["basis_dimension"] - len(rows)
            ),
            "basis": degree_map["basis"],
            "stress_basis": [{
                "id": target,
                "coordinates": [
                    fraction_record(value) for value in row
                ],
                "physical_free_stress_coordinates": [
                    fraction_record(
                        value / FREE_STRESS_DENOMINATOR ** powers)
                    for value in row
                ],
            } for target, row in zip(target_names, rows)],
            "complement_basis": [
                {
                    "id": degree_map["basis"][index],
                    "coordinate_index": index,
                }
                for index in complement
            ],
            "exact_validation": {
                "primes": degree_map["primes"],
                "sample_count_per_prime": min(
                    run["sample_count"]
                    for run in degree_map["per_prime"].values()
                ),
            },
        }
    degrees["12"] = degree12_entry
    result = {
        "schema": 1,
        "claim": (
            "Exact homogeneous free-stress scalar subalgebra through "
            "five-form degree 12."
        ),
        "scope": (
            "Static scalar invariants of the free traceless INZ stress "
            "tensor T=M/48; interacting formal stress maps are separate."
        ),
        "columns": [
            "degree",
            "full_dimension",
            "stress_dimension",
            "quotient_dimension",
        ],
        "degrees": degrees,
        "dimension_rows": [{
            key: degrees[str(degree)][key]
            for key in (
                "degree",
                "full_dimension",
                "stress_dimension",
                "quotient_dimension",
            )
        } for degree in (4, 6, 8, 10, 12)],
        "reproduction": {
            "compute": (
                ".venv/bin/python scripts/static_stress_degree12.py "
                "compute --prime PRIME --out CERTIFICATE"
            ),
            "assemble": (
                ".venv/bin/python scripts/static_stress_degree12.py "
                "assemble --certificates CERTS --validation-certificate "
                "HOLDOUT --out results/stress_flow/dimension_table.json"
            ),
        },
    }
    atomic_write_json(args.out, result)
    print(f"wrote {args.out}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    compute = subparsers.add_parser("compute")
    compute.add_argument("--prime", type=int, required=True)
    compute.add_argument("--seeds", type=int, nargs="+", default=DEFAULT_SEEDS)
    compute.add_argument("--max-memory-bytes", type=int,
                         default=DEFAULT_MAX_MEMORY)
    compute.add_argument("--checkpoint", type=Path)
    compute.add_argument("--out", type=Path, required=True)
    compute.set_defaults(function=compute_certificate)

    assembly = subparsers.add_parser("assemble")
    assembly.add_argument("--certificates", type=Path, nargs="+",
                          required=True)
    assembly.add_argument("--validation-certificate", type=Path)
    assembly.add_argument("--out", type=Path, required=True)
    assembly.set_defaults(function=assemble)

    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()

