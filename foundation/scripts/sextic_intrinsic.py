#!/usr/bin/env python3
"""Build the intrinsic sextic change-of-basis certificate."""

import argparse
from fractions import Fraction
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.catalog import atomic_write_json  # noqa: E402
from sdinv.exactmap import (  # noqa: E402
    fraction_record,
    rank_mod,
    sample_selfdual_five_form,
)
from sdinv.invariant_registry import load_verified_registry  # noqa: E402
from sdinv.sextic import (  # noqa: E402
    INTRINSIC_TO_REGISTRY,
    REGISTRY_TO_INTRINSIC,
    paper_i6_1,
    paper_i6_2,
    paper_i6_2_direct,
    paper_i6_2_graph_expansion,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = (
    ROOT / "results" / "stress_flow" / "change_of_basis"
    / "sextic_intrinsic.json"
)
PRIMES = (32749, 32719, 32693)
SAMPLE_SEEDS = (20260921, 20260922, 20260923, 20260924)


def _record_matrix(matrix):
    return [
        [fraction_record(value) for value in row]
        for row in matrix
    ]


def _residue(value, prime):
    value = Fraction(value)
    return (
        value.numerator * pow(value.denominator, -1, prime) % prime
    )


def build_result():
    registry = load_verified_registry(ROOT)
    validations = {}
    for prime in PRIMES:
        samples = []
        intrinsic_values = []
        for seed in SAMPLE_SEEDS:
            five_form = sample_selfdual_five_form(seed, prime)
            registry_values = registry.evaluate_degree(
                6, five_form, prime)
            stress_value = paper_i6_1(five_form, prime)
            quotient_value = paper_i6_2(five_form, prime)
            direct_value = paper_i6_2_direct(five_form, prime)
            expected_stress = (
                _residue(Fraction(32, 3), prime)
                * registry_values[0]
            ) % prime
            expected_quotient = (
                _residue(Fraction(-1, 1125), prime)
                * registry_values[0]
                + _residue(Fraction(3, 125), prime)
                * registry_values[1]
            ) % prime
            if stress_value != expected_stress:
                raise AssertionError("Tr(M^3) change of basis failed")
            if quotient_value != expected_quotient:
                raise AssertionError("1050-cubic change of basis failed")
            if quotient_value != direct_value:
                raise AssertionError("1050 graph/direct formulas disagree")
            intrinsic_values.append((stress_value, quotient_value))
            samples.append({
                "seed": seed,
                "registry_values": {
                    "I6_1": int(registry_values[0]),
                    "I6_2": int(registry_values[1]),
                },
                "intrinsic_values": {
                    "Tr(M^3)": int(stress_value),
                    "K_1050": int(quotient_value),
                },
                "direct_1050_value": int(direct_value),
            })
        validations[str(prime)] = {
            "rank": rank_mod(np.asarray(intrinsic_values), prime),
            "samples": samples,
        }
        if validations[str(prime)]["rank"] != 2:
            raise AssertionError("intrinsic sextic basis lost rank")

    determinant = (
        INTRINSIC_TO_REGISTRY[0][0]
        * INTRINSIC_TO_REGISTRY[1][1]
        - INTRINSIC_TO_REGISTRY[0][1]
        * INTRINSIC_TO_REGISTRY[1][0]
    )
    return {
        "schema": 1,
        "claim": (
            "The intrinsic sextic basis (Tr(M^3), K_1050) is exactly "
            "related to the verified graph basis, and K_1050 represents "
            "the non-free-stress quotient class."
        ),
        "sources": [
            "https://arxiv.org/abs/2509.14351v2",
            "https://arxiv.org/abs/2509.14350v2",
            "results/10d_order8.json",
        ],
        "registry_basis": ["I6_1", "I6_2"],
        "intrinsic_basis": ["Tr(M^3)", "K_1050"],
        "definitions": {
            "Tr(M^3)": "M_mu^nu M_nu^rho M_rho^mu",
            "K_1050": (
                "N1050_[abc,de]f "
                "N1050^[abc,]_[ghi] N1050^[def,gi]h"
            ),
            "N1050": "Lambda^mn_[abc Lambda_de]fmn",
        },
        "intrinsic_to_registry": _record_matrix(
            INTRINSIC_TO_REGISTRY),
        "registry_to_intrinsic": _record_matrix(
            REGISTRY_TO_INTRINSIC),
        "determinant": fraction_record(determinant),
        "quotient": {
            "subspace": "span{Tr(M^3)}",
            "normalization": "[K_1050] maps to 1",
            "registry_coordinate_formula": (
                "q(c1*I6_1+c2*I6_2)=125*c2/3"
            ),
        },
        "K_1050_graph_expansion": [{
            "graph": item["graph"],
            "coefficient": fraction_record(item["coefficient"]),
        } for item in paper_i6_2_graph_expansion()],
        "spinor_status": {
            "available_basis": ["Sigma_1", "Sigma_2"],
            "source": (
                "arXiv:2509.14350v2, section F5 invariants in spinor "
                "formalism, equations labeled Sigma1 and S2"
            ),
            "Sigma_1": "Tr(H^3)",
            "Sigma_2": (
                "Omega^(a1a2,a3a4) I_(b1a1,a3b2) "
                "I_(c1a2,a4c2) F^(b1c1) F^(b2c2)"
            ),
            "caveat": (
                "The source proves this is a sextic basis but does not "
                "publish its exact change of basis to "
                "(Tr(M^3),K_1050)."
            ),
        },
        "validation": {
            "method": (
                "exact finite-field evaluation of the graph expansion, "
                "the independent dense 1050-tensor contraction, and the "
                "verified registry basis"
            ),
            "primes": list(PRIMES),
            "sample_seeds": list(SAMPLE_SEEDS),
            "runs": validations,
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    result = build_result()
    atomic_write_json(args.out, result)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
