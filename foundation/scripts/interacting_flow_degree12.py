#!/usr/bin/env python3
"""Exact fully interacting stress-flow coefficient maps through degree 12."""

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

from sdinv.catalog import atomic_write_json  # noqa: E402
from sdinv.checkpoint import (  # noqa: E402
    load_checkpoint,
    load_checkpoint_payload,
    write_checkpoint,
)
from sdinv.contract import build_compact_derivative_basis  # noqa: E402
from sdinv.exactmap import (  # noqa: E402
    fraction_record,
    independent_row_indices,
    rank_mod,
    reconstruct_vector,
    sample_selfdual_five_form,
    solve_full_column_rank,
)
from sdinv.formal_flow import (  # noqa: E402
    TRACE_GENERATORS,
    flow_generator_directional_polynomials,
    flow_generator_polynomials,
    interaction_data,
)
from sdinv.forms import selfdual_projector, to_dense  # noqa: E402
from sdinv.invariant_registry import (  # noqa: E402
    load_verified_registry,
    load_verified_registry_through_degree12,
)
from sdinv.modp import inv  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
DEGREES = (4, 6, 8, 10, 12)
LOWER_FIT_SEEDS = tuple(range(20261101, 20261119))
LOWER_HOLDOUT_SEEDS = tuple(range(20261119, 20261123))
DEGREE12_HOLDOUTS_PER_SAMPLE = 2


def _file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _engine_sha256():
    digest = hashlib.sha256()
    for relative in (
        "scripts/interacting_flow_degree12.py",
        "src/sdinv/formal_flow.py",
        "src/sdinv/interaction.py",
        "src/sdinv/invariant_registry.py",
        "src/sdinv/stress.py",
    ):
        digest.update(relative.encode())
        digest.update((ROOT / relative).read_bytes())
    return digest.hexdigest()


def _target_id(generator, degree, coefficient_monomial):
    return (
        f"{generator}|d={int(degree)}|c="
        + ",".join(coefficient_monomial)
    )


def _target_schema(polynomials, allowed_degrees):
    generator_order = {
        generator.id: index
        for index, generator in enumerate(TRACE_GENERATORS)
    }
    records = []
    for generator, polynomial in polynomials.items():
        if generator == "tr_tau":
            continue
        for degree, coefficient_monomial in polynomial:
            if degree not in allowed_degrees:
                continue
            records.append({
                "id": _target_id(
                    generator, degree, coefficient_monomial),
                "generator": generator,
                "field_degree": int(degree),
                "coefficient_monomial": list(coefficient_monomial),
            })
    records.sort(key=lambda record: (
        record["field_degree"],
        generator_order[record["generator"]],
        record["coefficient_monomial"],
    ))
    if len({record["id"] for record in records}) != len(records):
        raise AssertionError("formal target IDs are not unique")
    return records


def _target_values(polynomials, schema):
    return {
        record["id"]: int(polynomials[record["generator"]].get(
            (
                record["field_degree"],
                tuple(record["coefficient_monomial"]),
            ),
            0,
        ))
        for record in schema
    }


def _static_context(prime, checkpoint_path, full_registry):
    payload = load_checkpoint_payload(checkpoint_path)
    identity = payload["identity"]
    state = payload["state"]
    if int(identity["prime"]) != int(prime):
        raise ValueError("static checkpoint prime differs")
    if identity["basis"] != [
        item.id for item in full_registry.basis(12)
    ]:
        raise ValueError("static checkpoint degree-12 basis differs")
    if len(state["samples"]) < 4:
        raise ValueError("static checkpoint has fewer than four samples")
    basis_rows = np.asarray([
        np.concatenate([
            np.asarray(sample["basis_rows"][index], dtype=np.int64)
            for sample in state["samples"]
        ])
        for index in range(72)
    ], dtype=np.int64) % prime
    if rank_mod(basis_rows, prime) != 72:
        raise ValueError("static checkpoint degree-12 basis rank changed")
    coordinate_count = len(state["samples"][0]["basis_rows"][0])
    if coordinate_count != 126:
        raise ValueError("static tangent coordinate count changed")

    round_robin = [
        sample_index * coordinate_count + coordinate_index
        for coordinate_index in range(coordinate_count)
        for sample_index in range(len(state["samples"]))
    ]
    selected_positions = independent_row_indices(
        basis_rows.T[round_robin], prime)
    fit_columns = [round_robin[index] for index in selected_positions]
    if len(fit_columns) != 72:
        raise ValueError("failed to select 72 degree-12 fit directions")
    fit_set = set(fit_columns)
    holdout_columns = []
    for sample_index in range(len(state["samples"])):
        available = [
            sample_index * coordinate_count + coordinate_index
            for coordinate_index in range(coordinate_count)
            if (
                sample_index * coordinate_count + coordinate_index
            ) not in fit_set
        ]
        holdout_columns.extend(
            available[:DEGREE12_HOLDOUTS_PER_SAMPLE])
    if len(holdout_columns) != (
        len(state["samples"]) * DEGREE12_HOLDOUTS_PER_SAMPLE
    ):
        raise ValueError("insufficient degree-12 holdout directions")
    return {
        "payload_sha256": _file_sha256(checkpoint_path),
        "static_engine_sha256": identity["engine_sha256"],
        "sample_seeds": [sample["seed"] for sample in state["samples"]],
        "basis_rows": basis_rows,
        "fit_columns": fit_columns,
        "holdout_columns": holdout_columns,
        "coordinate_count": coordinate_count,
    }


def _direction_context(prime):
    projector = selfdual_projector(10, 5, True, prime)
    basis = build_compact_derivative_basis(
        10, 5, projector, prime, independent=True)
    if basis.ncols != 126:
        raise ValueError("self-dual tangent dimension changed")
    return basis


def _compute_lower_sample(prime, seed, registry):
    five_form = sample_selfdual_five_form(seed, prime)
    values, derivatives = interaction_data(
        five_form,
        registry,
        prime,
        value_degrees=registry.degrees,
    )
    polynomials = flow_generator_polynomials(
        five_form, registry, prime, values, derivatives)
    schema = _target_schema(polynomials, (4, 6, 8, 10))
    return {
        "seed": int(seed),
        "basis_values": {
            str(degree): [
                int(values[item.id])
                for item in registry.basis(degree)
            ]
            for degree in (4, 6, 8, 10)
        },
        "target_schema": schema,
        "target_values": _target_values(polynomials, schema),
    }


def _compute_direction(prime, column, static, tangent_basis, registry):
    sample_index, direction_index = divmod(
        int(column), static["coordinate_count"])
    seed = static["sample_seeds"][sample_index]
    five_form = sample_selfdual_five_form(seed, prime)
    direction = to_dense(
        tangent_basis.directions[:, direction_index],
        10,
        5,
        prime,
    )
    polynomials, tangents, data = (
        flow_generator_directional_polynomials(
            five_form,
            direction,
            registry,
            prime,
            return_interaction_data=True,
        )
    )
    del polynomials
    schema = _target_schema(tangents, (12,))

    # The exact trace formula is independently enforced on every evaluated
    # tangent for the degrees whose interaction derivatives are materialized.
    value_tangents = data[1]
    trace_tangent = tangents["tr_tau"]
    for degree in (4, 6, 8):
        for item in registry.basis(degree):
            key = (degree, (item.id,))
            expected = (
                10 * (degree - 2) * value_tangents[item.id]
            ) % prime
            if trace_tangent.get(key, 0) != expected:
                raise AssertionError(
                    f"trace/homogeneity tangent failed for {item.id}")
    return {
        "column": int(column),
        "sample_index": sample_index,
        "sample_seed": int(seed),
        "direction_index": direction_index,
        "target_schema": schema,
        "target_tangents": _target_values(tangents, schema),
        "trace_homogeneity_checked": True,
    }


def _compatible_schema(records):
    if not records:
        raise ValueError("no computed records")
    schema = records[0]["target_schema"]
    for record in records:
        if record["target_schema"] != schema:
            raise ValueError("formal target schema changed across evaluations")
    return schema


def _solve_lower(prime, state, registry):
    fit_count = len(LOWER_FIT_SEEDS)
    records = state["lower_samples"]
    fit = records[:fit_count]
    holdout = records[fit_count:]
    if [record["seed"] for record in fit] != list(LOWER_FIT_SEEDS):
        raise ValueError("lower fit samples changed")
    if [record["seed"] for record in holdout] != list(
        LOWER_HOLDOUT_SEEDS
    ):
        raise ValueError("lower holdout samples changed")
    schema = _compatible_schema(records)
    solutions = []
    for degree in (4, 6, 8, 10):
        basis = [item.id for item in registry.basis(degree)]
        fit_basis = np.asarray([
            record["basis_values"][str(degree)] for record in fit
        ], dtype=np.int64) % prime
        if rank_mod(fit_basis, prime) != len(basis):
            raise ValueError(
                f"lower degree-{degree} fit basis is not full rank")
        holdout_basis = np.asarray([
            record["basis_values"][str(degree)] for record in holdout
        ], dtype=np.int64) % prime
        for target in (
            item for item in schema if item["field_degree"] == degree
        ):
            fit_target = np.asarray([
                record["target_values"][target["id"]]
                for record in fit
            ], dtype=np.int64)
            coordinates = solve_full_column_rank(
                fit_basis, fit_target, prime)
            holdout_target = np.asarray([
                record["target_values"][target["id"]]
                for record in holdout
            ], dtype=np.int64)
            if not np.array_equal(
                holdout_basis @ coordinates % prime,
                holdout_target % prime,
            ):
                raise ValueError(
                    f"lower holdout failed for {target['id']}")
            solutions.append({
                **target,
                "basis": basis,
                "coordinates": [
                    int(value) for value in coordinates],
                "fit_rank": len(basis),
                "holdout_passed": True,
            })
    return solutions


def _solve_degree12(prime, state, static, full_registry):
    records = state["directions"]
    by_column = {record["column"]: record for record in records}
    fit = [by_column[column] for column in static["fit_columns"]]
    holdout = [by_column[column] for column in static["holdout_columns"]]
    schema = _compatible_schema(records)
    basis = [item.id for item in full_registry.basis(12)]
    fit_basis = static["basis_rows"][:, static["fit_columns"]].T
    if rank_mod(fit_basis, prime) != 72:
        raise ValueError("selected degree-12 fit matrix is singular")
    holdout_basis = static["basis_rows"][
        :, static["holdout_columns"]].T
    solutions = []
    for target in schema:
        fit_target = np.asarray([
            record["target_tangents"][target["id"]]
            for record in fit
        ], dtype=np.int64)
        coordinates = solve_full_column_rank(
            fit_basis, fit_target, prime)
        holdout_target = np.asarray([
            record["target_tangents"][target["id"]]
            for record in holdout
        ], dtype=np.int64)
        if not np.array_equal(
            holdout_basis @ coordinates % prime,
            holdout_target % prime,
        ):
            raise ValueError(
                f"degree-12 holdout failed for {target['id']}")
        solutions.append({
            **target,
            "basis": basis,
            "coordinates": [int(value) for value in coordinates],
            "fit_rank": 72,
            "holdout_passed": True,
        })
    return solutions


def _analytic_trace_targets(prime, full_registry):
    targets = []
    for degree in DEGREES:
        basis = [item.id for item in full_registry.basis(degree)]
        for index, item_id in enumerate(basis):
            coordinates = [0] * len(basis)
            coordinates[index] = 10 * (degree - 2) % prime
            targets.append({
                "id": _target_id("tr_tau", degree, (item_id,)),
                "generator": "tr_tau",
                "field_degree": degree,
                "coefficient_monomial": [item_id],
                "basis": basis,
                "coordinates": coordinates,
                "fit_rank": len(basis),
                "holdout_passed": True,
                "derivation": "Tr(tau)[V_d]=10*(d-2)*V_d",
            })
    return targets


def compute(args):
    prime = int(args.prime)
    lower_registry = load_verified_registry(ROOT)
    full_registry = load_verified_registry_through_degree12(ROOT)
    static = _static_context(
        prime, args.static_checkpoint, full_registry)
    identity = {
        "schema": 1,
        "prime": prime,
        "engine_sha256": _engine_sha256(),
        "degree12_sha256": _file_sha256(
            ROOT / "results" / "10d_order12.json"),
        "static_checkpoint_sha256": static["payload_sha256"],
        "static_engine_sha256": static["static_engine_sha256"],
        "lower_fit_seeds": list(LOWER_FIT_SEEDS),
        "lower_holdout_seeds": list(LOWER_HOLDOUT_SEEDS),
        "degree12_sample_seeds": static["sample_seeds"],
        "degree12_fit_columns": static["fit_columns"],
        "degree12_holdout_columns": static["holdout_columns"],
    }
    state = {"lower_samples": [], "directions": []}
    if args.checkpoint.exists():
        state = load_checkpoint(args.checkpoint, identity)
        print(
            f"resuming prime {prime}: "
            f"{len(state['lower_samples'])} lower samples, "
            f"{len(state['directions'])} directions",
            flush=True,
        )
    started = time.perf_counter()
    lower_seeds = LOWER_FIT_SEEDS + LOWER_HOLDOUT_SEEDS
    completed_lower = [item["seed"] for item in state["lower_samples"]]
    if completed_lower != list(lower_seeds[:len(completed_lower)]):
        raise ValueError("lower sample checkpoint prefix changed")
    for seed in lower_seeds[len(completed_lower):]:
        state["lower_samples"].append(
            _compute_lower_sample(prime, seed, lower_registry))
        write_checkpoint(args.checkpoint, identity, state)
        print(
            f"prime {prime}: lower sample "
            f"{len(state['lower_samples'])}/{len(lower_seeds)}",
            flush=True,
        )

    columns = static["fit_columns"] + static["holdout_columns"]
    completed_columns = [item["column"] for item in state["directions"]]
    if completed_columns != columns[:len(completed_columns)]:
        raise ValueError("direction checkpoint prefix changed")
    tangent_basis = _direction_context(prime)
    for column in columns[len(completed_columns):]:
        state["directions"].append(_compute_direction(
            prime, column, static, tangent_basis, lower_registry))
        write_checkpoint(args.checkpoint, identity, state)
        print(
            f"prime {prime}: degree-12 direction "
            f"{len(state['directions'])}/{len(columns)}",
            flush=True,
        )

    targets = (
        _analytic_trace_targets(prime, full_registry)
        + _solve_lower(prime, state, lower_registry)
        + _solve_degree12(prime, state, static, full_registry)
    )
    generator_order = {
        generator.id: index
        for index, generator in enumerate(TRACE_GENERATORS)
    }
    targets.sort(key=lambda target: (
        target["field_degree"],
        generator_order[target["generator"]],
        target["coefficient_monomial"],
    ))
    certificate = {
        "schema": 1,
        "claim": (
            "Exact coefficient reduction of every analytic scalar "
            "stress-tensor generator through five-form degree 12."
        ),
        "normalization": (
            "tau=48*T; arbitrary generator coefficient functions absorb "
            "the corresponding constant powers of 48."
        ),
        "prime": prime,
        "engine_sha256": identity["engine_sha256"],
        "degree12_sha256": identity["degree12_sha256"],
        "static_checkpoint_sha256": identity[
            "static_checkpoint_sha256"],
        "static_engine_sha256": identity["static_engine_sha256"],
        "lower_fit_seeds": list(LOWER_FIT_SEEDS),
        "lower_holdout_seeds": list(LOWER_HOLDOUT_SEEDS),
        "degree12_sample_seeds": static["sample_seeds"],
        "degree12_fit_direction_count": len(static["fit_columns"]),
        "degree12_holdout_direction_count": len(
            static["holdout_columns"]),
        "degree12_fit_columns": static["fit_columns"],
        "degree12_holdout_columns": static["holdout_columns"],
        "trace_generators": [{
            "id": generator.id,
            "trace_powers": list(generator.trace_powers),
            "leading_field_degree": generator.leading_field_degree,
            "physical_T_rescaling": (
                f"product(Tr(tau^k))=48^"
                f"{sum(generator.trace_powers)}"
                "*product(Tr(T^k))"
            ),
        } for generator in TRACE_GENERATORS],
        "basis_dimensions": {
            str(degree): len(full_registry.basis(degree))
            for degree in DEGREES
        },
        "target_count": len(targets),
        "target_count_by_degree": {
            str(degree): sum(
                target["field_degree"] == degree for target in targets)
            for degree in DEGREES
        },
        "targets": targets,
        "all_holdouts_passed": all(
            target["holdout_passed"] for target in targets),
        "seconds": round(time.perf_counter() - started, 6),
    }
    atomic_write_json(args.out, certificate)
    print(f"wrote {args.out}", flush=True)


def _load(path):
    with Path(path).open() as stream:
        return json.load(stream)


def _target_signature(target):
    return (
        target["id"],
        target["generator"],
        target["field_degree"],
        target["coefficient_monomial"],
        target["basis"],
    )


def assemble(args):
    fit = [_load(path) for path in args.certificates]
    if len(fit) < 3:
        raise ValueError("at least three fit-prime certificates are required")
    primes = [int(certificate["prime"]) for certificate in fit]
    if len(set(primes)) != len(primes):
        raise ValueError("fit primes must be distinct")
    reference = fit[0]
    signatures = [
        _target_signature(target) for target in reference["targets"]]
    for certificate in fit:
        if (
            certificate["engine_sha256"] != reference["engine_sha256"]
            or certificate["degree12_sha256"]
            != reference["degree12_sha256"]
            or certificate["static_engine_sha256"]
            != reference["static_engine_sha256"]
            or certificate["trace_generators"]
            != reference["trace_generators"]
            or [
                _target_signature(target)
                for target in certificate["targets"]
            ] != signatures
            or not certificate["all_holdouts_passed"]
        ):
            raise ValueError("incompatible interacting-flow certificates")

    holdout = _load(args.validation_certificate)
    if int(holdout["prime"]) in primes:
        raise ValueError("validation prime must be independent")
    if (
        holdout["engine_sha256"] != reference["engine_sha256"]
        or [
            _target_signature(target)
            for target in holdout["targets"]
        ] != signatures
        or not holdout["all_holdouts_passed"]
    ):
        raise ValueError("incompatible interacting-flow holdout")

    rational_targets = []
    for index, signature in enumerate(signatures):
        modular = [
            np.asarray(
                certificate["targets"][index]["coordinates"],
                dtype=np.int64,
            )
            for certificate in fit
        ]
        rational = list(reconstruct_vector(modular, primes))
        holdout_prime = int(holdout["prime"])
        residue = [
            value.numerator * inv(value.denominator, holdout_prime)
            % holdout_prime
            for value in rational
        ]
        if residue != holdout["targets"][index]["coordinates"]:
            raise ValueError(
                f"independent-prime failure for {signature[0]}")
        source = reference["targets"][index]
        rational_targets.append({
            key: source[key]
            for key in (
                "id",
                "generator",
                "field_degree",
                "coefficient_monomial",
                "basis",
            )
        } | {
            "coordinates": [
                fraction_record(value) for value in rational],
            "derivation": source.get(
                "derivation", "exact finite-field coefficient reduction"),
        })

    result = {
        "schema": 1,
        "claim": (
            "Fully interacting formal coefficient equations for "
            "dV/dlambda=f(T,lambda), truncated at five-form degree 12."
        ),
        "normalization": reference["normalization"],
        "interaction_ansatz": (
            "V=sum_{d in {4,6,8,10,12}} sum_a c[d,a](lambda) "
            "I[d,a], including the verified product basis elements."
        ),
        "flow_equation": (
            "dot(c[d,a])=sum_G g_G(lambda) sum_C "
            "P[G,d,C,a]*C(c), with P recorded below."
        ),
        "trace_generators": reference["trace_generators"],
        "basis_dimensions": reference["basis_dimensions"],
        "target_count": len(rational_targets),
        "target_count_by_degree": reference["target_count_by_degree"],
        "targets": rational_targets,
        "exact_validation": {
            "fit_primes": primes,
            "independent_validation_prime": int(holdout["prime"]),
            "lower_fit_samples_per_prime": len(LOWER_FIT_SEEDS),
            "lower_holdout_samples_per_prime": len(
                LOWER_HOLDOUT_SEEDS),
            "degree12_generic_base_samples_per_prime": len(
                reference["degree12_sample_seeds"]),
            "degree12_fit_directions_per_prime": reference[
                "degree12_fit_direction_count"],
            "degree12_holdout_directions_per_prime": reference[
                "degree12_holdout_direction_count"],
            "all_modular_and_rational_holdouts_passed": True,
            "engine_sha256": reference["engine_sha256"],
            "degree12_sha256": reference["degree12_sha256"],
            "static_engine_sha256": reference[
                "static_engine_sha256"],
        },
    }
    atomic_write_json(args.out, result)
    print(f"wrote {args.out}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    compute_parser = subparsers.add_parser("compute")
    compute_parser.add_argument("--prime", type=int, required=True)
    compute_parser.add_argument(
        "--static-checkpoint", type=Path, required=True)
    compute_parser.add_argument("--checkpoint", type=Path, required=True)
    compute_parser.add_argument("--out", type=Path, required=True)
    compute_parser.set_defaults(function=compute)

    assemble_parser = subparsers.add_parser("assemble")
    assemble_parser.add_argument(
        "--certificates", type=Path, nargs="+", required=True)
    assemble_parser.add_argument(
        "--validation-certificate", type=Path, required=True)
    assemble_parser.add_argument("--out", type=Path, required=True)
    assemble_parser.set_defaults(function=assemble)

    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
