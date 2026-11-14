#!/usr/bin/env python3
"""Assemble the static map in a rational stress-adapted degree-12 basis."""

import argparse
from fractions import Fraction
import json
from math import isqrt
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.catalog import atomic_write_json  # noqa: E402
from sdinv.exactmap import (  # noqa: E402
    STRESS_TARGETS,
    crt,
    fraction_record,
    rank_mod,
    rational_reconstruct,
    select_standard_complement,
)


ROOT = Path(__file__).resolve().parents[1]
FREE_STRESS_DENOMINATOR = 48
DEGREE12_TARGETS = (
    "tr_M6",
    "tr_M4*tr_M2",
    "tr_M3^2",
    "tr_M2^3",
)


def _load(path):
    with Path(path).open() as stream:
        return json.load(stream)


def _fraction_matrix(records):
    return [
        Fraction(record["numerator"], record["denominator"])
        for record in records
    ]


def _standard_row(width, index, coefficient):
    return [
        Fraction(coefficient if column == index else 0)
        for column in range(width)
    ]


def _modular_complement(rows, prime):
    rows = np.asarray(rows, dtype=np.int64) % prime
    selected = []
    current = rank_mod(rows, prime)
    for index in range(rows.shape[1]):
        unit = np.zeros((1, rows.shape[1]), dtype=np.int64)
        unit[0, index] = 1
        candidate = np.vstack(
            [rows]
            + [
                np.eye(rows.shape[1], dtype=np.int64)[
                    chosen:chosen + 1]
                for chosen in selected
            ]
            + [unit]
        )
        next_rank = rank_mod(candidate, prime)
        if next_rank > current:
            selected.append(index)
            current = next_rank
        if current == rows.shape[1]:
            break
    return selected


def assemble(args):
    certificates = [_load(path) for path in args.certificates]
    if len(certificates) < 3:
        raise ValueError("at least three prime certificates are required")
    primes = [int(certificate["prime"]) for certificate in certificates]
    if len(set(primes)) != len(primes):
        raise ValueError("certificate primes must be distinct")
    reference = certificates[0]
    basis = reference["basis"]
    complements = []
    for certificate in certificates:
        if (
            certificate["basis"] != basis
            or certificate["basis_rank"] != 72
            or certificate["stress_rank"] != 4
            or certificate["quotient_dimension"] != 68
            or certificate["engine_sha256"]
            != reference["engine_sha256"]
            or certificate["degree12_sha256"]
            != reference["degree12_sha256"]
        ):
            raise ValueError("incompatible static certificate set")
        rows = [
            certificate["targets"][target]
            for target in DEGREE12_TARGETS
        ]
        if rank_mod(np.asarray(rows), certificate["prime"]) != 4:
            raise ValueError("static stress rows lost rank")
        complements.append(
            _modular_complement(rows, certificate["prime"]))
    if not all(indices == complements[0] for indices in complements):
        raise ValueError("standard complement changes with characteristic")
    complement_indices = complements[0]
    if len(complement_indices) != 68:
        raise ValueError("degree-12 complement is not 68-dimensional")
    pivot_indices = sorted(set(range(72)) - set(complement_indices))

    rational_rows = {
        "tr_M4*tr_M2": _standard_row(
            72, basis.index("I4_1*I8_1"), 2),
        "tr_M3^2": _standard_row(
            72, basis.index("I6_1^2"), Fraction(1024, 9)),
        "tr_M2^3": _standard_row(
            72, basis.index("I4_1^3"), 8),
    }
    modular_tr_m6 = {
        str(certificate["prime"]):
        certificate["targets"]["tr_M6"]
        for certificate in certificates
    }

    failed_columns = []
    reconstructed_columns = {}
    modulus = 1
    for prime in primes:
        modulus *= prime
    for column in range(72):
        residue, _ = crt([
            certificate["targets"]["tr_M6"][column]
            for certificate in certificates
        ], primes)
        try:
            reconstructed_columns[str(column)] = fraction_record(
                rational_reconstruct(residue, modulus))
        except ValueError:
            failed_columns.append(column)
    if not failed_columns:
        raise ValueError(
            "stress-adapted fallback is unnecessary: all columns lifted")

    degree12_stress_basis = [{
        "id": "tr_M6",
        "definition": "Tr(M^6)",
        "coordinates_in_original_atlas_by_prime": modular_tr_m6,
        "rational_original_atlas_coordinates": None,
        "physical_free_stress_definition": "Tr(T_free^6)=Tr(M^6)/48^6",
    }]
    for target in DEGREE12_TARGETS[1:]:
        row = rational_rows[target]
        degree12_stress_basis.append({
            "id": target,
            "definition": target.replace("_", " "),
            "coordinates": [
                fraction_record(value) for value in row],
            "physical_free_stress_coordinates": [
                fraction_record(value / FREE_STRESS_DENOMINATOR ** 6)
                for value in row
            ],
        })

    degree12_entry = {
        "degree": 12,
        "full_dimension": 72,
        "stress_dimension": 4,
        "quotient_dimension": 68,
        "original_atlas_basis": basis,
        "stress_adapted_basis": (
            list(DEGREE12_TARGETS)
            + [basis[index] for index in complement_indices]
        ),
        "stress_basis": degree12_stress_basis,
        "complement_basis": [{
            "id": basis[index],
            "coordinate_index_in_original_atlas": index,
        } for index in complement_indices],
        "pivot_indices_in_original_atlas": pivot_indices,
        "pivot_ids_in_original_atlas": [
            basis[index] for index in pivot_indices],
        "exact_validation": {
            "primes": primes,
            "sample_seeds": reference["sample_seeds"],
            "prime_count": len(primes),
            "sample_count_per_prime": len(
                reference["sample_seeds"]),
            "basis_rank_every_prime": 72,
            "stress_rank_every_prime": 4,
            "same_standard_complement_every_prime": True,
            "engine_sha256": reference["engine_sha256"],
            "degree12_sha256": reference["degree12_sha256"],
        },
        "rational_lift_audit": {
            "crt_modulus": str(modulus),
            "uniqueness_bound": str(isqrt(modulus // 2)),
            "successfully_reconstructed_columns": reconstructed_columns,
            "failed_uniqueness_bound_columns": failed_columns,
            "policy": (
                "No rational coefficient is guessed outside the certified "
                "CRT uniqueness bound. Tr(M^6) is used as an intrinsic "
                "stress-adapted basis element instead."
            ),
        },
    }

    lower = _load(
        ROOT / "results" / "stress_flow_exact_low_degree.json")
    degrees = {}
    for degree in (4, 6, 8, 10):
        degree_map = lower["degree_maps"][str(degree)]
        target_names = STRESS_TARGETS[degree]
        rows = [
            _fraction_matrix(degree_map["targets"][target])
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
                    fraction_record(value) for value in row],
                "physical_free_stress_coordinates": [
                    fraction_record(
                        value / FREE_STRESS_DENOMINATOR ** powers)
                    for value in row
                ],
            } for target, row in zip(target_names, rows)],
            "complement_basis": [{
                "id": degree_map["basis"][index],
                "coordinate_index": index,
            } for index in complement],
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
        "schema": 2,
        "claim": (
            "Exact homogeneous free-stress scalar subalgebra through "
            "five-form degree 12, using a stress-adapted basis at degree 12."
        ),
        "scope": (
            "Static scalar invariants of the free traceless INZ stress "
            "tensor T=M/48; interacting formal stress maps are separate."
        ),
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
                ".venv/bin/python "
                "scripts/assemble_static_stress_adapted.py "
                "--certificates CERTIFICATES "
                "--out results/stress_flow/dimension_table.json"
            ),
        },
    }
    atomic_write_json(args.out, result)
    print(f"wrote {args.out}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--certificates", type=Path, nargs="+", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assemble(args)


if __name__ == "__main__":
    main()
