#!/usr/bin/env python3
"""Fail-closed all-36 direct characteristic-zero degree-12 pipeline."""
from __future__ import annotations

import json
from pathlib import Path

from chiral4form.degree12_exact_targets import recover_exact_targets
from chiral4form.degree12_exact_source_rows import reconstruct_source_rows
from chiral4form.degree12_explicit_ideal import build_ideal
from chiral4form.degree12_exact_tangency import certify_tangency
from chiral4form.degree12_global_check import check_global
from chiral4form.degree12_all_orders_promotion import promote

ROOT = Path("verification/all_orders/degree12_exact")
INTERPOLATION = ROOT / "interpolation"
STALE = (
    ROOT / "exact_source_targets.json",
    ROOT / "exact_missing_targets.json",
    ROOT / "source_equation_space.json",
    ROOT / "kernel_lift.json",
    ROOT / "ideal.json",
    ROOT / "tangency.json",
    ROOT / "global_check.json",
    ROOT / "all_orders_theorem.json",
    ROOT / "target_recovery_status.json",
)


def _clear_stale_outputs():
    ROOT.mkdir(parents=True, exist_ok=True)
    if INTERPOLATION.exists() and not INTERPOLATION.is_dir():
        raise ValueError("degree-12 interpolation checkpoint path is not a directory")
    # Never delete anything under interpolation/: full basis, specialized B11,
    # and all-source per-sample checkpoints are expensive and resumable.
    for path in STALE:
        if path.exists():
            path.unlink()


def _gate(record, status):
    if record.get("status") != status:
        raise AssertionError(
            f"unexpected certificate status {record.get('status')!r}; expected {status!r}"
        )


def main():
    _clear_stale_outputs()
    print("=" * 78)
    print("DEGREE 12 HARDENED V2 — ALL 36 DIRECT QQ SOURCE ROWS")
    print("=" * 78)
    print("Eight full 72x72 checkpoints and all per-prime source checkpoints are preserved.")

    print("\n[1/6] Recover all 36 reduced physical source targets over QQ")
    targets = recover_exact_targets()
    _gate(targets, "exact_characteristic_zero_all_degree12_source_targets")
    print("PASS: 36/36 exact source targets")
    print("direct physical primes:", len(targets["direct_physical_primes"]))

    print("\n[2/6] Build all 36 exact 81D source equations and their QQ kernel")
    source, kernel = reconstruct_source_rows()
    _gate(source, "exact_characteristic_zero_degree12_source_equation_space")
    _gate(kernel, "exact_degree12_kernel_lift")
    print("PASS: equation count =", source["equation_count"])
    print("PASS: equation rank =", source["equation_rank"])
    print("PASS: kernel dimension =", kernel["kernel_dimension"])
    print("PASS: leading rank =", kernel["leading_rank"])
    print("free degree-12 columns =", kernel["free_degree12_columns"])

    print("\n[3/6] Build explicit 68-generator triangular ideal")
    ideal = build_ideal()
    _gate(ideal, "exact_characteristic_zero_degree12_triangular_ideal")
    print("PASS: constraints =", ideal["constraint_count"])
    print("PASS: fiber dimension =", ideal["fiber_dimension"])

    print("\n[4/6] Check exact tangency against all 36 physical source equations")
    tangency = certify_tangency()
    _gate(tangency, "exact_characteristic_zero_degree12_ideal_tangency")
    print("PASS: all 36 exact source residuals zero =", tangency["all_36_exact_source_residuals_zero"])

    print("\n[5/6] Match the finite global 6+4=10 orbit")
    global_check = check_global()
    _gate(global_check, "exact_degree12_ideal_global_orbit_equality_through_degree12")
    print("PASS:", global_check["base_dimension"], "+", global_check["fiber_dimension"], "=", global_check["orbit_dimension"])

    print("\n[6/6] Promote finite I12 to the all-orders invariant cylinder")
    theorem = promote()
    _gate(theorem, "unconditional_all_orders_degree12_obstruction_ideal_theorem")
    if theorem.get("all_orders_cylinder_invariant") is not True:
        raise AssertionError("all-orders cylinder invariant gate is false")
    print("PASS: all-orders cylinder invariant = True")

    final_path = ROOT / "all_orders_theorem.json"
    final = json.loads(final_path.read_text())
    _gate(final, "unconditional_all_orders_degree12_obstruction_ideal_theorem")
    if final.get("all_orders_cylinder_invariant") is not True:
        raise AssertionError("written final theorem certificate is not affirmative")

    print("\n" + "=" * 78)
    print("PASS: DEGREE-12 PROGRAM COMPLETE")
    print("=" * 78)
    print("Final certificate:", final_path)


if __name__ == "__main__":
    main()
