#!/usr/bin/env python3
"""Assemble interacting equations in the rational stress-adapted basis."""

import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.catalog import atomic_write_json  # noqa: E402
from sdinv.exactmap import (  # noqa: E402
    fraction_record,
    independent_row_indices,
    rank_mod,
    reconstruct_vector,
    solve_full_column_rank,
)
from sdinv.modp import inv  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
STATIC_TARGETS = (
    "tr_M6",
    "tr_M4*tr_M2",
    "tr_M3^2",
    "tr_M2^3",
)


def _load(path):
    with Path(path).open() as stream:
        return json.load(stream)


def _static_certificate(prime):
    return _load(
        ROOT / "results" / "stress_flow" / "certificates"
        / f"static_degree12_{prime}.json"
    )


def _adapted_basis_and_matrix(certificate, static):
    original = static["basis"]
    prime = int(static["prime"])
    candidates = [
        target for target in certificate["targets"]
        if (
            target["field_degree"] == 12
            and target["generator"] != "tr_tau"
        )
    ]
    # Put the four free/static rows first, then preserve formal target order.
    candidates = sorted(
        enumerate(candidates),
        key=lambda item: (
            0 if not item[1]["coefficient_monomial"] else 1,
            item[0],
        ),
    )
    candidate_targets = [target for _, target in candidates]
    candidate_rows = np.asarray([
        target["coordinates"] for target in candidate_targets
    ], dtype=np.int64) % prime
    selected_positions = independent_row_indices(
        candidate_rows, prime)
    forcing_targets = [
        candidate_targets[index] for index in selected_positions]
    current = candidate_rows[selected_positions]
    forcing_rank = rank_mod(current, prime)
    if forcing_rank != 21:
        raise ValueError(
            f"degree-12 new-forcing rank {forcing_rank}, expected 21")

    selected = []
    rank = forcing_rank
    for index in range(72):
        unit = np.zeros((1, 72), dtype=np.int64)
        unit[0, index] = 1
        next_rank = rank_mod(np.vstack((current, unit)), prime)
        if next_rank > rank:
            selected.append(index)
            current = np.vstack((current, unit))
            rank = next_rank
        if rank == 72:
            break
    if len(selected) != 51:
        raise ValueError("flow-adapted complement is not 51-dimensional")

    static_names = {
        "tr_tau6|d=12|c=": "tr_M6",
        "tr_tau2*tr_tau4|d=12|c=": "tr_M4*tr_M2",
        "tr_tau3^2|d=12|c=": "tr_M3^2",
        "tr_tau2^3|d=12|c=": "tr_M2^3",
    }
    definitions = []
    next_dynamic = 5
    for target in forcing_targets:
        item_id = static_names.get(target["id"])
        if item_id is None:
            item_id = f"F12_{next_dynamic:02d}"
            next_dynamic += 1
        definitions.append({
            "id": item_id,
            "kind": (
                "free_stress_coefficient"
                if target["id"] in static_names
                else "interacting_stress_coefficient"
            ),
            "source_target_id": target["id"],
        })
    if set(static_names.values()) - {
        definition["id"] for definition in definitions
    }:
        raise ValueError("the four static stress rows were not selected")
    definitions.extend({
        "id": original[index],
        "kind": "verified_atlas_complement",
        "original_atlas_index": index,
    } for index in selected)
    basis = [definition["id"] for definition in definitions]
    return (
        basis,
        current % prime,
        selected,
        definitions,
        [target["id"] for target in forcing_targets],
    )


def _target_id(generator, degree, coefficient_monomial):
    return (
        f"{generator}|d={degree}|c="
        + ",".join(coefficient_monomial)
    )


def _adapt_certificate(certificate):
    prime = int(certificate["prime"])
    static = _static_certificate(prime)
    if (
        static["engine_sha256"]
        != certificate["static_engine_sha256"]
    ):
        raise ValueError("dynamic/static engine mismatch")
    (
        adapted_basis,
        adapted_rows,
        complement,
        definitions,
        forcing_target_ids,
    ) = _adapted_basis_and_matrix(certificate, static)
    targets = []
    for target in certificate["targets"]:
        if target["field_degree"] != 12:
            targets.append(dict(target))
            continue
        if target["generator"] == "tr_tau":
            continue
        coordinates = solve_full_column_rank(
            adapted_rows.T,
            np.asarray(target["coordinates"], dtype=np.int64),
            prime,
        )
        adapted = dict(target)
        adapted["basis"] = adapted_basis
        adapted["coordinates"] = [
            int(value) for value in coordinates]
        targets.append(adapted)

    for index, item_id in enumerate(adapted_basis):
        coordinates = [0] * 72
        coordinates[index] = 100 % prime
        targets.append({
            "id": _target_id("tr_tau", 12, (item_id,)),
            "generator": "tr_tau",
            "field_degree": 12,
            "coefficient_monomial": [item_id],
            "basis": adapted_basis,
            "coordinates": coordinates,
            "fit_rank": 72,
            "holdout_passed": True,
            "derivation": "Tr(tau)[V_12]=100*V_12",
        })
    generator_order = {
        item["id"]: index
        for index, item in enumerate(certificate["trace_generators"])
    }
    targets.sort(key=lambda target: (
        target["field_degree"],
        generator_order[target["generator"]],
        target["coefficient_monomial"],
    ))
    result = dict(certificate)
    result["targets"] = targets
    result["target_count"] = len(targets)
    result["degree12_original_atlas_basis"] = static["basis"]
    result["degree12_flow_adapted_basis"] = adapted_basis
    result["degree12_flow_adapted_basis_definitions"] = definitions
    result["degree12_forcing_basis_target_ids"] = forcing_target_ids
    result["degree12_complement_indices"] = complement
    return result


def _signature(target):
    return (
        target["id"],
        target["generator"],
        target["field_degree"],
        target["coefficient_monomial"],
        target["basis"],
    )


def assemble(args):
    fit = [_adapt_certificate(_load(path))
           for path in args.certificates]
    if len(fit) < 3:
        raise ValueError("at least three fit primes are required")
    primes = [int(certificate["prime"]) for certificate in fit]
    if len(primes) != len(set(primes)):
        raise ValueError("fit primes must be distinct")
    holdout = _adapt_certificate(_load(args.validation_certificate))
    if int(holdout["prime"]) in primes:
        raise ValueError("validation prime must be independent")
    reference = fit[0]
    signatures = [_signature(target) for target in reference["targets"]]
    for certificate in fit + [holdout]:
        if (
            certificate["engine_sha256"]
            != reference["engine_sha256"]
            or certificate["degree12_sha256"]
            != reference["degree12_sha256"]
            or certificate["static_engine_sha256"]
            != reference["static_engine_sha256"]
            or certificate["trace_generators"]
            != reference["trace_generators"]
            or certificate["degree12_flow_adapted_basis_definitions"]
            != reference["degree12_flow_adapted_basis_definitions"]
            or [_signature(target) for target in certificate["targets"]]
            != signatures
            or not certificate["all_holdouts_passed"]
        ):
            raise ValueError("incompatible interacting certificates")

    holdout_prime = int(holdout["prime"])
    targets = []
    for index, signature in enumerate(signatures):
        rational = reconstruct_vector([
            np.asarray(
                certificate["targets"][index]["coordinates"],
                dtype=np.int64,
            )
            for certificate in fit
        ], primes)
        residue = [
            value.numerator * inv(value.denominator, holdout_prime)
            % holdout_prime
            for value in rational
        ]
        if residue != holdout["targets"][index]["coordinates"]:
            raise ValueError(
                f"independent-prime failure for {signature[0]}")
        source = reference["targets"][index]
        targets.append({
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
                "derivation",
                "exact finite-field coefficient reduction",
            ),
        })

    forcing_dimensions = {}
    forcing_quotients = {}
    for degree, full_dimension in (
        (4, 1), (6, 2), (8, 7), (10, 14), (12, 72)
    ):
        rows = np.asarray([
            target["coordinates"]
            for target in reference["targets"]
            if (
                target["field_degree"] == degree
                and target["generator"] != "tr_tau"
            )
        ], dtype=np.int64)
        dimension = rank_mod(rows, primes[0])
        for certificate in fit[1:] + [holdout]:
            check_rows = np.asarray([
                target["coordinates"]
                for target in certificate["targets"]
                if (
                    target["field_degree"] == degree
                    and target["generator"] != "tr_tau"
                )
            ], dtype=np.int64)
            if rank_mod(check_rows, int(certificate["prime"])) != dimension:
                raise ValueError(
                    f"new-forcing rank changes at degree {degree}")
        forcing_dimensions[str(degree)] = dimension
        forcing_quotients[str(degree)] = full_dimension - dimension

    result = {
        "schema": 2,
        "claim": (
            "Necessary and sufficient fully interacting coefficient "
            "equations for dV/dlambda=f(T,lambda), through five-form "
            "degree 12."
        ),
        "normalization": reference["normalization"],
        "interaction_ansatz": (
            "V=sum_d sum_a c[d,a](lambda) I[d,a], in the verified "
            "homogeneous bases, with the stress-adapted basis at degree 12."
        ),
        "flow_equation": (
            "dot(c[d,a])=sum_G g_G(lambda) sum_C "
            "P[G,d,C,a]*C(c); every nonzero P is recorded below."
        ),
        "coefficient_direction_count": 96,
        "generic_functional_rank": 81,
        "coordinate_count_explanation": (
            "The truncated polynomial vector space has "
            "1+2+7+14+72=96 homogeneous coefficients. The number 81 is "
            "the cumulative generic Jacobian/transcendence rank and is not "
            "a linear homogeneous basis dimension."
        ),
        "trace_generators": reference["trace_generators"],
        "flow_generator_count": len(reference["trace_generators"]),
        "basis_dimensions": reference["basis_dimensions"],
        "degree12_flow_adapted_basis": reference[
            "degree12_flow_adapted_basis"],
        "degree12_flow_adapted_basis_definitions": reference[
            "degree12_flow_adapted_basis_definitions"],
        "new_forcing_dimension_by_degree": forcing_dimensions,
        "new_forcing_quotient_dimension_by_degree": forcing_quotients,
        "new_forcing_definition": (
            "Span of all formal coefficient invariants from stress "
            "generators other than Tr(tau). Tr(tau) is excluded because it "
            "only propagates interaction coordinates already present in "
            "the seed."
        ),
        "target_count": len(targets),
        "target_count_by_degree": reference["target_count_by_degree"],
        "targets": targets,
        "exact_validation": {
            "fit_primes": primes,
            "independent_validation_prime": holdout_prime,
            "lower_fit_samples_per_prime": len(
                reference["lower_fit_seeds"]),
            "lower_holdout_samples_per_prime": len(
                reference["lower_holdout_seeds"]),
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
        "reproduction": (
            ".venv/bin/python "
            "scripts/assemble_interacting_stress_adapted.py "
            "--certificates FIT_CERTS --validation-certificate HOLDOUT "
            "--out results/stress_flow/closure_equations/"
            "interacting_through_degree12.json"
        ),
    }
    atomic_write_json(args.out, result)
    print(f"wrote {args.out}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--certificates", type=Path, nargs="+", required=True)
    parser.add_argument(
        "--validation-certificate", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assemble(args)


if __name__ == "__main__":
    main()
