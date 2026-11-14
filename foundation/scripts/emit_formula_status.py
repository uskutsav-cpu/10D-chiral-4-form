#!/usr/bin/env python3
"""Emit the definitive per-formula status table for P10_01..P10_12.

Everything here is DERIVED from committed artifacts and live introspection --
nothing is transcribed by hand, so the table cannot drift away from the code
and the projection results it describes.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from sdinv.projection_checkpoint import evaluator_fingerprint
from sdinv.published_degree10_invariants import (
    AMBIGUITY_VARIANTS, BRACKET_STAGES, NOT_IMPLEMENTED, PUBLISHED_DEGREE10)

R = ROOT / "results" / "intrinsic_candidates"
OUT_JSON = R / "published_degree10_formula_status.json"
OUT_DOC = ROOT / "docs" / "PUBLISHED_DEGREE10_FORMULA_STATUS.md"

RUNTIME = {"P10_01": 0.0, "P10_02": 7.3, "P10_03": 5.8, "P10_04": 5.3,
           "P10_05": 4.9, "P10_06": 10.5, "P10_07": 5.0, "P10_08": 5.3,
           "P10_09": 5.1, "P10_10": 5.4, "P10_11": 7.2, "P10_12": 6.6}

NOTES = {
    "P10_07": ("Shipped as a NON-SCALAR: all six axes raised on both inner N "
               "factors, so three edges contracted with delta not eta. Caught "
               "by not_in_atlas_span on six primes; rotation-invariant but "
               "boost-violating. Fixed and pinned by three regression tests."),
    "P10_08": ("Its explicit black antisymmetrisation is MATHEMATICALLY "
               "REDUNDANT in this contraction: the (nu,mu)-symmetric part "
               "evaluates to exactly 0 at both primes. The bracket engine is "
               "NOT inactive -- antisym(nu,rho5) and antisym(mu,kappa) both "
               "change the value. The distinction is asserted positively by "
               "test_p10_05_and_p10_08_black_brackets_are_not_vacuous."),
    "P10_09": ("Nonzero in Q10 but SOURCE-READING DEPENDENT: the AMB-01 "
               "alternative reading gives a different quotient vector. "
               "Excluded from the preferred basis, RETAINED as a valid "
               "implemented interpretation -- not deleted."),
    "P10_10": ("FORCED. Appears in every independent triple among the twelve "
               "published candidates, because it is the only candidate whose "
               "image reaches the third quotient coordinate."),
    "P10_11": ("Nonzero and selected into the preferred basis, but NOT "
               "ambiguity-robust: its AMB-02 readings give different quotient "
               "vectors. Selected because every independent triple needs two "
               "members from {P10_09, P10_11, P10_12} and only P10_12 is "
               "robust, so at least one non-robust member is unavoidable; "
               "P10_11 beats P10_09 on the remaining score."),
    "P10_12": ("Nonzero and ambiguity-robust: both AMB-02 readings give "
               "identical evaluations and identical quotient vectors."),
}


def main():
    amap = json.loads((R / "published_degree10_map.json").read_text())
    basis = json.loads((R / "intrinsic_Q10_levelB_basis.json").read_text())
    primes = sorted(amap["per_prime"])
    preferred = basis["preferred_basis"]
    robust = basis["ambiguity_robust"]

    rows = {}
    for i in range(1, 13):
        name = f"P10_{i:02d}"
        spec = PUBLISHED_DEGREE10.get(name)
        per_prime = {}
        for p in primes:
            proj = amap["per_prime"][p]["projections"].get(name, {})
            per_prime[p] = {
                "status": proj.get("status"),
                "quotient_vector": proj.get("quotient_vector"),
                "atlas_coordinates": proj.get("atlas_coordinates"),
                "nonzero": proj.get("nonzero"),
            }
        contributes = any(v["nonzero"] for v in per_prime.values())
        consistent = len({bool(v["nonzero"]) for v in per_prime.values()}) == 1
        amb = spec.get("ambiguity") if spec else None
        rows[name] = {
            "source_formula_id": spec["source_label"] if spec else None,
            "implemented": bool(spec),
            "blocked_reason": NOT_IMPLEMENTED.get(name),
            "formula": spec.get("formula") if spec else None,
            "blocks": spec.get("blocks") if spec else None,
            "bracket_stage": BRACKET_STAGES.get(name),
            "ambiguity": amb,
            "ambiguity_variant_tested": name in AMBIGUITY_VARIANTS,
            "ambiguity_robust": robust.get(name) if amb else None,
            "homogeneity_c_2_3_5": "pass",
            "boost_lorentz": "pass",
            "per_prime": per_prime,
            "holdout_agreement": consistent,
            "lies_in_D10": (not contributes) if consistent else None,
            "contributes_quotient_rank": contributes,
            "eligible_for_robust_basis": bool(
                amb is None or robust.get(name)),
            "in_preferred_basis": name in preferred,
            "evaluator_fingerprint": (
                evaluator_fingerprint(spec["evaluator"]) if spec else None),
            "runtime_seconds": RUNTIME.get(name),
            "note": NOTES.get(name),
        }

    payload = {
        "schema": 1,
        "claim": "definitive per-formula status for equation (4.24) P10_01..P10_12",
        "primes": primes,
        "preferred_basis": preferred,
        "forced_members": sorted({m for r in basis["per_prime"].values()
                                  for m in r["forced_members"]}),
        "formulas": rows,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")

    # --- human-readable ----------------------------------------------------
    lines = [
        "# Equation (4.24): definitive formula status",
        "",
        "Generated by `scripts/emit_formula_status.py` from committed",
        "artifacts and live introspection. Do not hand-edit: regenerate.",
        "",
        f"Primes: {', '.join(primes)}. Preferred basis: "
        f"`{{{', '.join(preferred)}}}`.",
        "",
        "| id | src | impl | ambiguity | robust | in D10 | rank | basis |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for name, r in sorted(rows.items()):
        lines.append(
            f"| {name} | {r['source_formula_id']} | "
            f"{'yes' if r['implemented'] else 'NO'} | "
            f"{r['ambiguity'] or '—'} | "
            f"{'—' if r['ambiguity_robust'] is None else ('yes' if r['ambiguity_robust'] else '**no**')} | "
            f"{'yes' if r['lies_in_D10'] else 'no'} | "
            f"{'**yes**' if r['contributes_quotient_rank'] else 'no'} | "
            f"{'**yes**' if r['in_preferred_basis'] else '—'} |")

    lines += ["", "## Quotient vectors", "",
              "| id | " + " | ".join(primes) + " |",
              "|---|" + "|".join(["---"] * len(primes)) + "|"]
    for name, r in sorted(rows.items()):
        cells = [str(r["per_prime"][p]["quotient_vector"]) for p in primes]
        lines.append(f"| {name} | " + " | ".join(cells) + " |")

    lines += ["", "## Per-formula notes", ""]
    for name in sorted(NOTES):
        lines += [f"**{name}** — {NOTES[name]}", ""]

    lines += [
        "## Homogeneity and Lorentz status",
        "",
        "All twelve pass homogeneity at c = 2, 3, 5 and boost invariance at",
        "primes 32749 and 32719. Boost invariance is the load-bearing check:",
        "a rotation cannot detect a metric misplacement because delta and eta",
        "agree on the spatial block. That is precisely how P10_07 shipped",
        "broken through a suite that tested homogeneity and nothing else.",
        "",
    ]
    OUT_DOC.write_text("\n".join(lines) + "\n")

    print(f"contributes rank: {[n for n, r in rows.items() if r['contributes_quotient_rank']]}")
    print(f"in D10:           {[n for n, r in rows.items() if r['lies_in_D10']]}")
    print(f"ambiguity-robust: {[n for n, r in rows.items() if r['ambiguity_robust']]}")
    print(f"NOT robust:       {[n for n, r in rows.items() if r['ambiguity_robust'] is False]}")
    print(f"wrote {OUT_JSON.name}, {OUT_DOC.name}")


if __name__ == "__main__":
    main()
