#!/usr/bin/env python3
"""Build the exact low-degree five-form/stress-tensor map.

The final artifact uses exact finite-field evaluations, CRT, and rational
reconstruction.  It never uses floating-point fitting.
"""

import argparse
from fractions import Fraction
import json
from math import isqrt, prod
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.exactmap import (  # noqa: E402
    DEFAULT_PRIMES,
    DEFAULT_SAMPLE_SEEDS,
    STRESS_TARGETS,
    compute_degree_map,
    fraction_record,
    rank_mod,
    rational_target_rows,
    sample_selfdual_five_form,
    select_standard_complement,
)
from sdinv.invariant_registry import (  # noqa: E402
    degree12_placeholder_items,
    degree12_product_items,
    load_verified_registry,
)
from sdinv.stress import (  # noqa: E402
    matrix_trace_power,
    modmax_stress,
    modmax_stress_square_formula,
    paper_i4_i8_i12,
    paper_i8_i12_expanded,
    stress_mixed,
    stress_v_i4,
    stress_v_i4_raw,
    symmetric_inner,
)


ROOT = Path(__file__).resolve().parents[1]
FREE_STRESS_DENOMINATOR = 48
PAPER_SOURCE_SHA256 = (
    "eb889f3dfd60280a95a907343e39ade60d62df6f56b4cd53c33cd445283ffb4e")
PAPER_PDF_SHA256 = (
    "c4f812035fb8d8d07aa8b0dba1a2e55de3e909a77f75397ee17f3fe21f4f5e90")
# The largest reconstructed numerator in the verified degree-10 map is below
# this conservative bound.  Rational reconstruction needs sqrt(modulus/2)
# to exceed both numerator and denominator bounds.
RECONSTRUCTION_BOUND = 300_000_000
MIN_SUPPORTED_PRIME = 11
MAX_SUPPORTED_PRIME = max(DEFAULT_PRIMES)
# These characteristics divide fixed rational coefficients in the verified
# maps.  The small ones are already excluded by MIN_SUPPORTED_PRIME.
SINGULAR_CHARACTERISTICS = frozenset((2, 3, 5, 7, 3871))


def _is_prime(value):
    value = int(value)
    if value < 2:
        return False
    if value == 2:
        return True
    if value % 2 == 0:
        return False
    for divisor in range(3, isqrt(value) + 1, 2):
        if value % divisor == 0:
            return False
    return True


def _validate_configuration(primes, sample_seeds):
    if len(set(primes)) != len(primes) or not all(
            _is_prime(prime) for prime in primes):
        raise ValueError("CRT moduli must be distinct primes")
    if any(
        prime < MIN_SUPPORTED_PRIME
        or prime > MAX_SUPPORTED_PRIME
        or prime in SINGULAR_CHARACTERISTICS
        for prime in primes
    ):
        raise ValueError(
            f"primes must be in the supported exact-arithmetic range "
            f"{MIN_SUPPORTED_PRIME}..{MAX_SUPPORTED_PRIME} and must not "
            "divide fixed map denominators"
        )
    reconstruction_modulus = prod(primes)
    if reconstruction_modulus <= 2 * RECONSTRUCTION_BOUND ** 2:
        raise ValueError(
            "the combined CRT modulus is too small for the verified "
            f"degree-10 reconstruction bound {RECONSTRUCTION_BOUND}; "
            "use at least four default-sized primes"
        )
    if len(sample_seeds) < 15 or len(set(sample_seeds)) != len(sample_seeds):
        raise ValueError(
            "at least 15 distinct samples are required: 14 degree-10 "
            "fit rows plus held-out validation"
        )


def _fraction(entry):
    return Fraction(entry["numerator"], entry["denominator"])


def _fraction_vector(entries):
    return [_fraction(entry) for entry in entries]


def _records(vector):
    return [fraction_record(value) for value in vector]


def _scale(vector, scalar):
    scalar = Fraction(scalar)
    return [scalar * value for value in vector]


def _add(*vectors):
    if not vectors:
        return []
    return [
        sum((vector[index] for vector in vectors), Fraction())
        for index in range(len(vectors[0]))
    ]


def _standard_vector(width, index, coefficient=1):
    return [
        Fraction(coefficient if position == index else 0)
        for position in range(width)
    ]


def _stress_normalization(target):
    powers = {
        "tr_M2": 2,
        "tr_M3": 3,
        "tr_M4": 4,
        "tr_M2^2": 4,
        "tr_M5": 5,
        "tr_M2*tr_M3": 5,
    }
    return Fraction(1, FREE_STRESS_DENOMINATOR ** powers[target])


def _decorate_map(degree_map):
    physical = {}
    for target in STRESS_TARGETS[degree_map["degree"]]:
        vector = _fraction_vector(degree_map["targets"][target])
        physical[target.replace("M", "T")] = _records(
            _scale(vector, _stress_normalization(target)))
    degree_map["physical_free_stress_monomials"] = physical
    return degree_map


def _stage1_reproduction(primes, seeds):
    records = {}
    for prime in primes:
        samples = []
        for seed in seeds:
            five_form = sample_selfdual_five_form(seed, prime)
            compact = paper_i4_i8_i12(five_form, prime)
            expanded = paper_i8_i12_expanded(five_form, prime)
            direct = stress_v_i4_raw(
                five_form, v=137, v_i=211, mod=prime)
            decomposed = stress_v_i4(
                five_form, v=137, v_i=211, mod=prime)
            stress = modmax_stress(five_form, b=173, mod=prime)
            square = modmax_stress_square_formula(
                five_form, b=173, mod=prime)
            samples.append({
                "seed": int(seed),
                "I4": compact["I4"],
                "tr_M3": compact["tr_M3"],
                "I8": compact["I8"],
                "I12": compact["I12"],
                "expanded_I8_matches": (
                    compact["I8"] == expanded["I8"]),
                "expanded_I12_matches": (
                    compact["I12"] == expanded["I12"]),
                "equation_3_3_equals_3_4": bool(
                    np.array_equal(direct, decomposed)),
                "stress_symmetric": bool(np.array_equal(stress, stress.T)),
                "stress_trace": matrix_trace_power(
                    stress_mixed(stress, prime), 1, prime),
                "stress_square_direct": symmetric_inner(
                    stress, stress, prime),
                "stress_square_formula": square["stress_square"],
                "root_term": square["root_term"],
                "delta_I8": square["delta_I8"],
                "delta_I12": square["delta_I12"],
            })
        records[str(prime)] = samples
    return records


def _subalgebra_report(degree_map):
    stress_targets = STRESS_TARGETS[degree_map["degree"]]
    stress_rows = rational_target_rows(degree_map)
    complement_indices = select_standard_complement(stress_rows)
    basis = degree_map["basis"]
    modular = {}
    for prime_text, run in degree_map["per_prime"].items():
        prime = int(prime_text)
        values = run["sample_values"]
        stress_values = np.column_stack([
            values["targets"][target] for target in stress_targets
        ])
        stress_rank = rank_mod(stress_values, prime)
        obstruction_ranks = {}
        for index in complement_indices:
            union = np.column_stack((
                stress_values,
                np.asarray(values["basis"], dtype=np.int64)[:, index],
            ))
            obstruction_ranks[basis[index]] = rank_mod(union, prime)
        modular[prime_text] = {
            "stress_value_rank": stress_rank,
            "stress_plus_obstruction_ranks": obstruction_ranks,
            "sample_seeds": run["sample_seeds"],
        }
    return {
        "degree": degree_map["degree"],
        "five_form_value_dimension": degree_map["basis_dimension"],
        "stress_generated_dimension": len(stress_rows),
        "stress_basis": [
            {
                "id": target,
                "coordinates": degree_map["targets"][target],
            }
            for target in stress_targets
        ],
        "complement_basis": [
            {
                "id": basis[index],
                "coordinates": _records(
                    _standard_vector(len(basis), index)),
            }
            for index in complement_indices
        ],
        "modular_nonmembership_certificates": modular,
    }


def _flow_pilot(maps):
    degree8 = maps[8]
    degree10 = maps[10]
    i8 = _fraction_vector(degree8["targets"]["paper_I8"])
    i4_squared = _fraction_vector(
        degree8["targets"]["tr_M2^2"])
    m2r = _fraction_vector(degree10["targets"]["tr_M2R"])
    m2m3 = _fraction_vector(degree10["targets"]["tr_M2*tr_M3"])

    tr_t2_degree8 = _add(
        _scale(i8, Fraction(1, 72)),
        _scale(i4_squared, Fraction(11, 2016)),
    )
    tr_t3_degree8 = _scale(
        i4_squared,
        Fraction(1, 8 * FREE_STRESS_DENOMINATOR ** 2),
    )
    tr_t3_degree10 = _add(
        m2r,
        _scale(m2m3, Fraction(-6, 7)),
    )
    tr_t3_degree10 = _scale(
        tr_t3_degree10,
        Fraction(1, FREE_STRESS_DENOMINATOR ** 2),
    )
    tr_t4_degree10 = _scale(
        m2m3,
        Fraction(1, 6 * FREE_STRESS_DENOMINATOR ** 3),
    )
    i4_i6_1 = _standard_vector(
        degree10["basis_dimension"],
        degree10["basis"].index("I4_1*I6_1"),
        Fraction(1, 2 * FREE_STRESS_DENOMINATOR ** 2),
    )
    i4_i6_2 = _standard_vector(
        degree10["basis_dimension"],
        degree10["basis"].index("I4_1*I6_2"),
        Fraction(1, 2 * FREE_STRESS_DENOMINATOR ** 2),
    )

    parameterizations = {}
    for degree, degree_map in maps.items():
        parameterizations[str(degree)] = {
            "basis": degree_map["basis"],
            "condition": (
                "The coefficient vector must be a rational linear "
                "combination of these rows."
            ),
            "rows": [
                {
                    "parameter": f"a_{degree}_{index}",
                    "stress_monomial": target,
                    "coordinates": degree_map["targets"][target],
                }
                for index, target in enumerate(
                    STRESS_TARGETS[degree], start=1)
            ],
        }

    return {
        "scope": (
            "Proven for the free-seed/linearized polynomial stress algebra "
            "through five-form degree 10. The displayed V=c4*I4 trace "
            "expansion includes all terms through degree 10 for that "
            "one-coupling ansatz. A full all-coupling nonlinear closure "
            "classification is not claimed."
        ),
        "linearized_general_ansatz": {
            "ansatz": (
                "dV/dlambda=sum_d sum_a c[d,a] I[d,a], d=4,6,8,10"
            ),
            "closed_parameterizations": parameterizations,
            "nontrivial_family_exists": True,
            "uniqueness": (
                "Closure fixes all quotient/complement components but "
                "leaves 1,1,2,2 stress parameters at degrees 4,6,8,10."
            ),
            "lowest_obstruction": {
                "degree": 6,
                "invariant": "I6_2",
                "reason": (
                    "The degree-6 stress space is span{tr_M3}="
                    "span{I6_1}; I6_2 raises the sampled rank from 1 to 2 "
                    "at every saved prime and sample set."
                ),
            },
        },
        "V_equals_c4_I4_trace_expansion": {
            "normalization": (
                "T=M/48 + c4*I4*identity/24 + "
                "c4^2*(-2*I4*M/7+R/3)"
            ),
            "TrT2": {
                "degree4": "tr_M2/48^2",
                "degree8_coefficient_of_c4^2": {
                    "basis": degree8["basis"],
                    "coordinates": _records(tr_t2_degree8),
                    "formula": "paper_I8/72+11*I4^2/2016",
                },
            },
            "TrT3": {
                "degree6": "tr_M3/48^3",
                "degree8_coefficient_of_c4": {
                    "basis": degree8["basis"],
                    "coordinates": _records(tr_t3_degree8),
                    "formula": "I4^2/(8*48^2)",
                },
                "degree10_coefficient_of_c4^2": {
                    "basis": degree10["basis"],
                    "coordinates": _records(tr_t3_degree10),
                    "formula": (
                        "(tr(M^2 R)-6 I4 tr(M^3)/7)/48^2"
                    ),
                },
            },
            "TrT4": {
                "degree8": "tr_M4/48^4",
                "degree10_coefficient_of_c4": {
                    "basis": degree10["basis"],
                    "coordinates": _records(tr_t4_degree10),
                    "formula": "I4 tr(M^3)/(6*48^3)",
                },
            },
            "TrT5": {"degree10": "tr_M5/48^5"},
            "TrT6_through_TrT10": (
                "Their leading five-form degrees are 12 through 20, so "
                "their truncations through degree 10 vanish."
            ),
        },
        "general_V6_linear_terms": {
            "TrT3_degree10_coefficient_of_c6_1": {
                "basis": degree10["basis"],
                "coordinates": _records(i4_i6_1),
            },
            "TrT3_degree10_coefficient_of_c6_2": {
                "basis": degree10["basis"],
                "coordinates": _records(i4_i6_2),
            },
            "formula": "I4*V6/(4*48^2)",
        },
    }


def build_result(primes, sample_seeds):
    registry = load_verified_registry(ROOT)
    maps = {}
    for degree in (4, 6, 8, 10):
        extras = {
            8: ("paper_I8",),
            10: ("tr_M2R",),
        }.get(degree, ())
        maps[degree] = _decorate_map(compute_degree_map(
            registry,
            degree,
            primes=primes,
            sample_seeds=sample_seeds,
            extra_targets=extras,
        ))
    subalgebra = {
        degree: _subalgebra_report(degree_map)
        for degree, degree_map in maps.items()
    }

    return {
        "schema": 1,
        "claim": (
            "Exact free-stress scalar map and linearized stress subalgebra "
            "through verified five-form degree 10, plus exact reproduction "
            "of the D=10 ModMax identities of arXiv:2509.14351v2."
        ),
        "scope_limits": [
            "No claim is made that the 81 functionally independent "
            "invariants generate the polynomial ring.",
            "The full nonlinear all-coupling flow-closure problem is not "
            "proved; the exact pilot scope is recorded separately.",
            "Degree-12 primitive formulas are not imported on this branch.",
        ],
        "paper": {
            "title": "On non-linear chiral 4-form theories in D=10",
            "arxiv": "https://arxiv.org/abs/2509.14351v2",
            "version": "v2",
            "source_sha256": PAPER_SOURCE_SHA256,
            "pdf_sha256": PAPER_PDF_SHA256,
            "hash_policy": (
                "Pinned arXiv v2 checksums; artifact regeneration does "
                "not depend on ignored local download files."
            ),
            "equations_reproduced": [
                "2.11-2.18",
                "2.33-2.36",
                "3.3-3.4",
                "3.10-3.16",
            ],
        },
        "conventions": {
            "metric": "diag(-1,+1,+1,+1,+1,+1,+1,+1,+1,+1)",
            "orientation": "epsilon^{01...9}=-1",
            "self_duality": "Lambda=*Lambda",
            "M": "M_mu^nu=Lambda_{mu rho(4)} Lambda^{nu rho(4)}",
            "free_stress": "T_{mu nu}=M_{mu nu}/(2*4!)=M/48",
            "flow_density": (
                "Equation (3.16) uses the 4!-rescaled INZ density; "
                "the action in (2.7) has an overall 1/4!."
            ),
        },
        "exact_arithmetic": {
            "primes": [int(prime) for prime in primes],
            "sample_seeds": [int(seed) for seed in sample_seeds],
            "method": (
                "Exact modular evaluation, deterministic full-rank fit "
                "samples, held-out sample validation, CRT, bounded rational "
                "reconstruction, and residue revalidation."
            ),
            "floating_point_used_for_claims": False,
        },
        "stage1_modmax_reproduction": _stage1_reproduction(
            primes[:3], sample_seeds[:3]),
        "stress_trace_registry": {
            "cayley_hamilton_limit": 10,
            "traceless": True,
            "traces": [
                {
                    "id": f"tr_T{power}",
                    "power": power,
                    "leading_five_form_degree": 2 * power,
                    "inside_verified_degree10_window": power <= 5,
                    "free_seed_target": (
                        f"tr_M{power}/{FREE_STRESS_DENOMINATOR}^{power}"
                    ),
                }
                for power in range(2, 11)
            ],
        },
        "degree_maps": {str(degree): value
                        for degree, value in maps.items()},
        "stress_subalgebra": {
            "lowest_nonzero_complement_degree": 6,
            "dimension_table": [
                {
                    "degree": degree,
                    "five_form_value_dimension": (
                        subalgebra[degree]["five_form_value_dimension"]),
                    "stress_generated_dimension": (
                        subalgebra[degree]["stress_generated_dimension"]),
                    "complement_dimension": len(
                        subalgebra[degree]["complement_basis"]),
                }
                for degree in (4, 6, 8, 10)
            ],
            "degrees": {str(degree): value
                        for degree, value in subalgebra.items()},
        },
        "flow_closure_pilot": _flow_pilot(maps),
        "degree12_interface": {
            "product_dimension": len(degree12_product_items()),
            "products": [item.id for item in degree12_product_items()],
            "primitive_slots": len(degree12_placeholder_items()),
            "primitive_slot_ids": [
                item.id for item in degree12_placeholder_items()],
            "modmax_I12_ready": True,
            "import_contract": (
                "Supply exactly 62 validated degree-12 InvariantItem "
                "objects using unambiguous graph_record formulas; the "
                "registry then exposes a 72-dimensional homogeneous basis."
            ),
        },
        "reproduction": {
            "command": (
                ".venv/bin/python scripts/stress_flow_pipeline.py "
                "--out results/stress_flow_exact_low_degree.json"
            ),
            "tests": ".venv/bin/python -m pytest tests -q",
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        default=str(ROOT / "results" / "stress_flow_exact_low_degree.json"),
    )
    parser.add_argument(
        "--primes", type=int, nargs="+", default=list(DEFAULT_PRIMES))
    parser.add_argument(
        "--sample-seeds",
        type=int,
        nargs="+",
        default=list(DEFAULT_SAMPLE_SEEDS),
    )
    args = parser.parse_args()
    try:
        _validate_configuration(args.primes, args.sample_seeds)
    except ValueError as exc:
        parser.error(str(exc))
    result = build_result(tuple(args.primes), tuple(args.sample_seeds))
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(f"wrote {output}")
    for row in result["stress_subalgebra"]["dimension_table"]:
        print(
            f"degree {row['degree']}: total "
            f"{row['five_form_value_dimension']}, stress "
            f"{row['stress_generated_dimension']}, complement "
            f"{row['complement_dimension']}")


if __name__ == "__main__":
    main()
