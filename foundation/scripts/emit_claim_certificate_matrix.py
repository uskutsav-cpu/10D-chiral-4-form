#!/usr/bin/env python3
"""Map every generated manuscript number to the artifact that produced it.

Each entry records the macro, its current value, the artifact it is read from,
the command that regenerates it, and the strength of the underlying evidence.

Strength is the point. The manuscript mixes five kinds of statement and a reader
cannot tell them apart from the typeset number alone:

    analytic          a theorem or a literature result, no computation
    exact-char0       exact over Q, with a certificate
    exact-modular     exact over F_p; a bound in characteristic zero unless the
                      spanning-set argument applies
    multi-prime       agreement across primes; probabilistic, not a proof
    numerical         floating point, with a tolerance
    limitation        a stated scope restriction rather than a result

    python scripts/emit_claim_certificate_matrix.py
"""

from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GEN = ROOT / "manuscript" / "generated" / "numbers.tex"
OUT_MD = ROOT / "audit" / "FINAL_CLAIM_CERTIFICATE_MATRIX.md"
OUT_JSON = ROOT / "audit" / "FINAL_CLAIM_CERTIFICATE_MATRIX.json"

#: macro prefix/name -> (artifact, regeneration command, strength, note)
SOURCES = [
    (r"^(dimAten|dimBten|dimGten|dimPten|dimDten|dimQten|capBG|capPB|capPG|"
     r"capDB|sumPB|nIncidencePrimes)$",
     "results/intrinsic_candidates/degree10_space_incidence.json",
     "python scripts/emit_degree10_space_incidence.py --write",
     "exact-modular",
     "A10/B10/G10/P10 are additionally pinned over Q by the spanning-set "
     "argument; D10 is settled separately, see dimDtenQ"),
    (r"^(dimBcapP|dimBplusP|publishedPrimitiveContent)$",
     "results/intrinsic_candidates/degree10_published_product_intersection.json",
     "python scripts/emit_degree10_space_incidence.py --write",
     "exact-modular",
     "characteristic-zero status in docs/B10_P10_INTERSECTION_STATUS.md"),
    (r"^(dimDtenQ|dimQtenQ|dimAtenQ|czSettled|czMinorSize|czFitPrimes|"
     r"czHoldoutPrime|czIntegerRows|czLiftedRows|czSweeps|czQtenRelation)$",
     "results/stress_flow/D10_characteristic_zero.json",
     "python scripts/d10_characteristic_zero.py",
     "exact-char0",
     "CRT lift validated at a held-out prime; explicit non-vanishing minor"),
    (r"^(bridgePrimes|nBridgePrimes|starSquared|dimSelfdual|dimAntiSelfdual|"
     r"bridgeForwardRank|dimGammaTraceless|gammaTraceRank)$",
     "spinor_trace_bridge/results/bridge_validation.json",
     "cd spinor_trace_bridge && python -m pytest",
     "exact-modular",
     "kernel and image by two-way span equality, not dimension counting"),
    (r"^(comparisonPrime|nCommonSamples|nSample|traceRankDeg|spinorRankDeg|"
     r"spansEqualDeg|holdoutDeg)",
     "verification/spinor_trace_comparison.json",
     "python spinor_trace_bridge/scripts/run_comparison.py",
     "exact-modular",
     "spans compared by two-way containment on a hashed common registry"),
    (r"^bp(Settled|CapQ|SumQ|FitPrimes|HoldoutPrime|FreshPrime|FreshSamples|CapRelation)$",
     "results/degree10/B10_P10_intersection_exact.json",
     "python scripts/b10_p10_characteristic_zero.py",
     "exact-char0",
     "CRT lift of the published coordinates validated at a held-out prime; the "
     "generator is verified at a further prime and samples used nowhere else"),
    (r"^matrix(Cells|Primes|Seeds|Rank|Agree|StableRows|StableCols|UnstableRows)$",
     "results/rank81/full_rank_matrix_publication_final.json",
     "python scripts/aggregate_rank_matrix.py",
     "exact-modular",
     "15 cells, 5 primes x 3 seeds, identical rank and identical pivot "
     "structure; corroboration, not the basis of the rank claim"),
    (r"^gTen(LoadBearing|CounterfactualQten|DerivedQten)$",
     "results/stress_flow/G10_publication_certificate.json",
     "python -m pytest tests/test_G10_trace_activation.py",
     "analytic",
     "the vanishing of the quadratic stress trace is derived, not computed; "
     "the counterfactual is a computation showing it is load-bearing"),
    (r"^dEight",
     "verification/degree8_span_equality.json",
     "python spinor_trace_bridge/scripts/run_degree8_span.py",
     "exact-modular",
     "four primes, holdout-validated change of basis, family ablation"),
    (r"^(nExactJacPrimes|exactJac|charZeroLowerBound|cumRankDeg)",
     "results/rank81/certificate.json",
     "python spinor_trace_bridge/scripts/run_rank81_certificate.py --archive PATH",
     "exact-char0",
     "lower bound only; integral basis makes the modular rank unconditional"),
    (r"^(minorSize|minorPrimes|nMinorPrimes|minorDetNonzero)$",
     "results/rank81/minor81_certificate.json",
     "python spinor_trace_bridge/scripts/run_minor81_certificate.py --archive PATH",
     "exact-char0",
     "a single prime suffices: an integer determinant vanishing mod p would "
     "vanish mod every p"),
    (r"^(archiveJacobianRank|archiveJacobianRunsAgree|nArchive)",
     "verification/spinor_archive_jacobian_exact.json",
     "python spinor_trace_bridge/scripts/run_archive_jacobian_exact.py",
     "exact-modular",
     "superseded by the full-selection certificate; retained as history"),
    (r"^(nFloat|floatNondegenerateRanks)",
     "verification/SPINOR_JACOBIAN_RUNS.json",
     "python spinor_trace_bridge/scripts/analyse_archived_jacobians.py",
     "numerical",
     "float64 with an explicit noise floor; evidence about METHOD, not rank"),
    (r"^(nTraceTests|nBridgeTests)$",
     "pytest --collect-only",
     "python manuscript/scripts/make_numbers.py",
     "limitation",
     "reproducibility metadata, not a scientific claim"),
    (r"^(dimFiveFormComponents|dimSelfDualModule|dimLorentzAlgebra|"
     r"genericFunctionalDim|dimSymSpinor|spinorDim|dimQtenMinusOne)$",
     "literature / elementary",
     "not computed here",
     "analytic",
     "126 - 45 = 81 is the cited generic-orbit argument; NOT ours"),
]


def classify(name: str):
    for pattern, artifact, cmd, strength, note in SOURCES:
        if re.search(pattern, name):
            return artifact, cmd, strength, note
    return None, None, "UNCLASSIFIED", "no source mapping"


def main() -> int:
    if not GEN.exists():
        print("generated/numbers.tex missing; run make_numbers.py")
        return 1
    macros = re.findall(r"\\newcommand\{\\(\w+)\}\{([^}]*)\}", GEN.read_text())

    rows, unclassified = [], []
    for name, value in macros:
        artifact, cmd, strength, note = classify(name)
        if strength == "UNCLASSIFIED":
            unclassified.append(name)
        rows.append({"macro": name, "value": value, "artifact": artifact,
                     "regenerate": cmd, "strength": strength, "note": note})

    from collections import Counter
    counts = Counter(r["strength"] for r in rows)

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps({
        "schema": 1,
        "generated_by": "scripts/emit_claim_certificate_matrix.py",
        "n_macros": len(rows),
        "by_strength": dict(counts),
        "unclassified": unclassified,
        "entries": rows,
    }, indent=1) + "\n")

    L = ["# Final claim-to-certificate matrix", "",
         "Generated by `scripts/emit_claim_certificate_matrix.py`. Do not edit.",
         "",
         "Every scientific number in the manuscript is a generated macro. This "
         "table maps each to its artifact and, more importantly, to the "
         "**strength** of the evidence, which the typeset number does not show.",
         "", "## Evidence strengths", "",
         "| strength | meaning |", "|---|---|",
         "| `analytic` | theorem or literature result; no computation here |",
         "| `exact-char0` | exact over `Q`, with a certificate |",
         "| `exact-modular` | exact over `F_p`; a characteristic-zero bound unless a spanning-set argument applies |",
         "| `multi-prime` | agreement across primes; probabilistic |",
         "| `numerical` | floating point with a tolerance |",
         "| `limitation` | scope metadata, not a result |",
         "", "## Counts", "", "| strength | macros |", "|---|---:|"]
    for k, v in sorted(counts.items()):
        L.append(f"| `{k}` | {v} |")
    if unclassified:
        L += ["", f"**{len(unclassified)} UNCLASSIFIED macros**: "
                  + ", ".join(f"`{u}`" for u in unclassified),
              "", "An unclassified macro is a number in the paper with no "
                  "recorded provenance. This must be empty before submission."]
    else:
        L += ["", "No unclassified macros: every number in the manuscript has a "
                  "recorded artifact and evidence strength."]

    L += ["", "## Entries", "",
          "| macro | value | strength | artifact |", "|---|---|---|---|"]
    for r in sorted(rows, key=lambda r: (r["strength"], r["macro"])):
        L.append(f"| `{r['macro']}` | `{r['value']}` | `{r['strength']}` | "
                 f"`{r['artifact']}` |")
    L += ["", "## Regeneration", ""]
    seen = {}
    for r in rows:
        if r["regenerate"] and r["regenerate"] not in seen:
            seen[r["regenerate"]] = r["artifact"]
    for cmd, art in seen.items():
        L.append(f"- `{art}`  \n  `{cmd}`")
    OUT_MD.write_text("\n".join(L) + "\n")

    print(f"{len(rows)} macros; " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    if unclassified:
        print(f"UNCLASSIFIED: {unclassified}")
        return 1
    print(f"wrote {OUT_MD.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
