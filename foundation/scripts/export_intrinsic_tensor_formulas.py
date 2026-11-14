#!/usr/bin/env python3
"""Stage 1: export Level-A intrinsic tensor formulas for the seven quotient
classes, with independent dense-vs-graph validation.

Level A is an explicit Einstein-index contraction of copies of the self-dual
five-form and the metric. It is coordinate-independent and fully specified,
which is what makes it an intrinsic representative. It is emphatically NOT
compact or canonical -- that is Level C and is not attempted here.

Validation is by a second evaluator that builds its einsum from the dummy
index assignment and applies the metric by raising one slot of each pair,
independent of the repository's slot-planner.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from sdinv.modp import P, ALT_P  # noqa: E402
from sdinv.forms import selfdual_projector, to_dense, random_form  # noqa: E402
from sdinv.contract import value  # noqa: E402
from sdinv.graph_to_tensor import (  # noqa: E402
    contraction_specification, dense_evaluate, parse_graph_label)

PRIMES = (32749, 32719, 32693)
SAMPLES = (20260729, 20260730, 20260731, 20260732)
FRESH = (20261201, 20261202)

# intrinsic provisional IDs, deliberately distinct from graph labels
CLASSES = {
    "Q10_A": ("I10_6", 10),
    "Q10_B": ("I10_7", 10),
    "Q10_C": ("I10_12", 10),
    "Q12_A": ("I12_59", 12),
    "Q12_B": ("I12_60", 12),
    "Q12_C": ("I12_61", 12),
    "Q12_D": ("I12_62", 12),
}


def graph_labels():
    labels = {}
    with (ROOT / "results" / "10d_order10.json").open() as stream:
        for g in json.load(stream)["generators"]:
            labels[g["id"]] = g["graph"]
    with (ROOT / "results" / "10d_order12.json").open() as stream:
        for g in json.load(stream)["generators"]:
            labels[g["id"]] = g["graph"]
    return labels


def sample(prime, seed):
    projector = selfdual_projector(10, 5, True, prime)
    raw = random_form(10, 5, np.random.default_rng(seed), prime)
    return to_dense((projector @ raw) % prime, 10, 5, prime)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path)
    ap.add_argument("--degree", type=int, choices=(10, 12))
    args = ap.parse_args()

    labels = graph_labels()
    selected = {k: v for k, v in CLASSES.items()
                if args.degree is None or v[1] == args.degree}

    records = {}
    all_ok = True
    for intrinsic_id, (graph_id, degree) in selected.items():
        label = labels[graph_id]
        spec = contraction_specification(label)
        matrix = parse_graph_label(label)

        checks = []
        for prime in PRIMES:
            for seed in SAMPLES + FRESH:
                dense_form = sample(prime, seed)
                g = value(matrix, dense_form, 10, 5, True, prime) % prime
                d = dense_evaluate(label, dense_form, True, prime)
                checks.append({
                    "prime": prime, "seed": seed,
                    "graph_evaluator": int(g), "dense_evaluator": int(d),
                    "agree": bool(d == g),
                    "fresh_sample": seed in FRESH,
                })
        agree = all(c["agree"] for c in checks)
        all_ok &= agree

        # homogeneity: F -> c F must scale by c^degree
        prime = PRIMES[0]
        base = sample(prime, SAMPLES[0])
        c = 7
        scaled = (base * c) % prime
        v1 = dense_evaluate(label, base, True, prime)
        v2 = dense_evaluate(label, scaled, True, prime)
        homog = (v2 % prime) == (pow(c, degree, prime) * v1) % prime

        records[intrinsic_id] = {
            "intrinsic_id": intrinsic_id,
            "graph_basis_label": graph_id,
            "graph": label,
            "field_degree": degree,
            "level": "A",
            "level_A_specification": spec,
            "validation": {
                "primes": list(PRIMES),
                "samples": list(SAMPLES),
                "fresh_samples": list(FRESH),
                "dense_vs_graph_checks": len(checks),
                "all_agree": agree,
                "homogeneity_degree_check": bool(homog),
                "detail": checks,
            },
            "level_B_status": "NOT DERIVED (M / N^(1050) / N^(4125) form)",
            "level_C_status": "NOT DERIVED (canonical compact form)",
            "caveat": (
                "Level A is explicit and coordinate-independent but is a "
                f"{degree}-fold contraction with "
                f"{spec['metric_factor_count']} metric factors. It must not "
                "be described as compact, canonical, or explanatory."),
        }
        print(f"{intrinsic_id} ({graph_id}, degree {degree}): "
              f"{len(checks)} checks, agree={agree}, homogeneous={homog}, "
              f"{spec['metric_factor_count']} metric factors")

    payload = {
        "schema": 1,
        "claim": ("Level-A explicit Einstein-index representatives for the "
                  "seven generalized-flow quotient classes, validated against "
                  "an independent dense evaluator."),
        "translation_rule": (
            "I = prod_v F_{slots[v]} x prod_pairs eta^{pair}; implemented with "
            "one einsum letter per contracted pair and the metric applied by "
            "raising one slot of each pair"),
        "classes": records,
        "all_validated": all_ok,
        "levels": {
            "A": "explicit F-index contraction — COMPLETE for all seven",
            "B": "M / N^(1050) / N^(4125) form — not derived",
            "C": "canonical compact form — not derived",
        },
    }
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
        print(f"\nwrote {args.out}")
    print(f"ALL VALIDATED: {all_ok}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
