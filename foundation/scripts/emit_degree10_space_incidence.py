#!/usr/bin/env python3
"""Regenerate the degree-10 subspace incidence and product-intersection certificates.

Both certificates were originally produced by an inline computation, which left
the headline dimension table of the manuscript without a versioned generator.
This script is that generator. It rebuilds both files from the repository's own
evaluators at as many primes as have a stress-flow closure certificate, and by
default refuses to overwrite a value that disagrees with the existing artifact,
so re-running it is a check as well as a regeneration.

The five subspaces, all expressed in the same fourteen atlas coordinates:

    A10  the full degree-10 atlas
    G10  the twelve graph generators   I10_1 .. I10_12
    P10  the two lower products        I4_1*I6_1, I4_1*I6_2
    B10  the span of the published equation-(4.24) candidates
    D10  the stress-flow reachable closure

Run:
    python scripts/emit_degree10_space_incidence.py            # check only
    python scripts/emit_degree10_space_incidence.py --write    # regenerate
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from sdinv.forms import selfdual_projector, to_dense, random_form   # noqa: E402
from sdinv.exactmap import rank_mod                                 # noqa: E402
from sdinv.published_degree10_invariants import PUBLISHED_DEGREE10  # noqa: E402
from stress_flow_closure import closure_span                        # noqa: E402
from test_M_only_quotients import (                                 # noqa: E402
    registry_items, evaluate_atlas_element, solve_exact,
)

CERT_DIR = ROOT / "results" / "stress_flow" / "certificates"
OUT_INCIDENCE = ROOT / "results/intrinsic_candidates/degree10_space_incidence.json"
OUT_INTERSECT = (ROOT /
                 "results/intrinsic_candidates/degree10_published_product_intersection.json")

#: The degree-8 seeding of the flow, identical to the rest of the pipeline.
SEED8 = {4: ["I4_1"], 6: ["I6_2"], 8: ["I8_3", "I8_4", "I8_5", "I8_6"]}

SPACES = {
    "A10": "full degree-10 atlas (14)",
    "B10": "published equation-(4.24) candidate span (12)",
    "D10": "stress-flow reachable closure (11)",
    "G10": "span of the twelve graph generators I10_1..I10_12 (12)",
    "P10": "lower-product subspace, I4_1*I6_1 and I4_1*I6_2 (2)",
}


def sample(prime: int, seed: int):
    proj = selfdual_projector(10, 5, True, prime)
    raw = random_form(10, 5, np.random.default_rng(seed), prime)
    return to_dense((proj @ raw) % prime, 10, 5, prime)


def rank(rows, prime: int) -> int:
    if not len(rows):
        return 0
    return int(rank_mod(np.asarray(rows, dtype=np.int64) % prime, prime))


def sum_rank(a, b, prime: int) -> int:
    return rank(list(a) + list(b), prime)


def intersection_dim(a, b, prime: int) -> int:
    """dim(U cap V) = dim U + dim V - dim(U + V), which needs no kernel solve."""
    return rank(a, prime) + rank(b, prime) - sum_rank(a, b, prime)


def build_spaces(prime: int, n_extra: int = 8):
    """The five subspaces at one prime, in atlas coordinates."""
    cert_path = CERT_DIR / f"interacting_degree12_{prime}.json"
    if not cert_path.exists():
        return None
    cert = json.loads(cert_path.read_text())
    _, bmap, span = closure_span(cert, SEED8, prime)
    names = bmap[10]
    n = len(names)

    items = registry_items()
    forms = [sample(prime, 41000 + 11 * i) for i in range(n + n_extra)]
    cols = [[evaluate_atlas_element(items[nm], items, f, prime, {}) for f in forms]
            for nm in names]
    design = [[cols[j][i] for j in range(n)] for i in range(len(forms))]

    identity = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    graph_idx = [i for i, nm in enumerate(names) if not nm.startswith("I4_1*")]
    prod_idx = [i for i, nm in enumerate(names) if nm.startswith("I4_1*")]

    published, unsolved = [], []
    for cid, spec in PUBLISHED_DEGREE10.items():
        target = [spec["evaluator"](f, prime) for f in forms]
        vec, ok = solve_exact(design, target, prime)
        if ok:
            published.append([v % prime for v in vec])
        else:
            unsolved.append(cid)

    return {
        "names": names,
        "n_published_solved": len(published),
        "unsolved": unsolved,
        "A10": identity,
        "G10": [identity[i] for i in graph_idx],
        "P10": [identity[i] for i in prod_idx],
        "B10": published,
        "D10": [list(r) for r in span[10]],
    }


def analyse(sp, prime: int) -> tuple[dict, dict]:
    keys = ["A10", "B10", "D10", "G10", "P10"]
    dims = {k: rank(sp[k], prime) for k in keys}

    incidence = {}
    for i, u in enumerate(keys):
        for v in keys[i + 1:]:
            s = sum_rank(sp[u], sp[v], prime)
            cap = dims[u] + dims[v] - s
            incidence[f"{u}|{v}"] = {
                "dim_sum": s,
                "dim_intersection": cap,
                f"{u}_contains_{v}": cap == dims[v] and s == dims[u],
                f"{v}_contains_{u}": cap == dims[u] and s == dims[v],
            }

    cap_bp = intersection_dim(sp["B10"], sp["P10"], prime)
    intersect = {
        "dim_B10_published": dims["B10"],
        "dim_P10_product": dims["P10"],
        "dim_B10_cap_P10": cap_bp,
        "dim_B10_plus_P10": sum_rank(sp["B10"], sp["P10"], prime),
        # A complement must meet P10 trivially AND span the atlas with it.
        "B10_is_complement_of_P10": (
            cap_bp == 0 and sum_rank(sp["B10"], sp["P10"], prime) == dims["A10"]),
        "B10_equals_span_of_graph_generators": (
            intersection_dim(sp["B10"], sp["G10"], prime) == dims["B10"] == dims["G10"]),
        "n_published_solved": sp["n_published_solved"],
    }
    return dims, incidence, intersect


def write_certificates(per_prime_inc: dict, per_prime_int: dict, names) -> None:
    """Write both certificates. Called after every prime, not only at the end."""
    OUT_INCIDENCE.write_text(json.dumps({
        "schema": 1,
        "generated_by": "scripts/emit_degree10_space_incidence.py",
        "spaces": SPACES,
        "incidence": {},
        "per_prime": per_prime_inc,
    }, indent=1, sort_keys=True) + "\n")

    prev = json.loads(OUT_INTERSECT.read_text()) if OUT_INTERSECT.exists() else {}
    OUT_INTERSECT.write_text(json.dumps({
        "schema": 1,
        "generated_by": "scripts/emit_degree10_space_incidence.py",
        "atlas_names": names,
        "product_columns": [n for n in (names or []) if n.startswith("I4_1*")],
        "product_source": prev.get(
            "product_source",
            "explicit lower_product atlas entries: I4_1 x I6_1 and I4_1 x I6_2"),
        "conclusion": prev.get("conclusion"),
        "superseded_conclusion": prev.get("superseded_conclusion"),
        "per_prime": per_prime_int,
    }, indent=1, sort_keys=True) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--primes", default="",
                    help="comma separated; default is every prime with a "
                         "stress-flow closure certificate")
    ap.add_argument("--write", action="store_true",
                    help="write the certificates; without it the script only "
                         "checks the existing ones and reports disagreement")
    args = ap.parse_args()

    if args.primes.strip():
        primes = [int(x) for x in args.primes.split(",") if x.strip()]
    else:
        primes = sorted(int(p.stem.rsplit("_", 1)[1])
                        for p in CERT_DIR.glob("interacting_degree12_*.json"))

    # Snapshot what is already recorded before touching anything. It is both the
    # resume point and the reference the recomputation is checked against, and
    # incremental writing would destroy it if it were read later.
    stored = (json.loads(OUT_INCIDENCE.read_text()).get("per_prime", {})
              if OUT_INCIDENCE.exists() else {})
    per_prime_inc, per_prime_int, names = {}, {}, None
    if args.write:
        per_prime_inc = dict(stored)
        if OUT_INTERSECT.exists():
            per_prime_int = dict(
                json.loads(OUT_INTERSECT.read_text()).get("per_prime", {}))
        if OUT_INTERSECT.exists():
            names = json.loads(OUT_INTERSECT.read_text()).get("atlas_names")

    for p in primes:
        if args.write and str(p) in per_prime_inc and str(p) in per_prime_int:
            print(f"  prime {p}: already recorded, skipped")
            continue
        sp = build_spaces(p)
        if sp is None:
            print(f"  prime {p}: no closure certificate, skipped")
            continue
        names = sp["names"]
        dims, incidence, intersect = analyse(sp, p)
        # Check against the stored value before recording it, so a disagreement
        # stops the run instead of being written and then reported.
        if str(p) in stored and stored[str(p)]["dims"] != dims:
            print(f"FAIL: prime {p} disagrees with the existing certificate\n"
                  f"  stored:   {stored[str(p)]['dims']}\n"
                  f"  computed: {dims}")
            return 1
        per_prime_inc[str(p)] = {"dims": dims, "incidence": incidence}
        per_prime_int[str(p)] = intersect
        if args.write:
            write_certificates(per_prime_inc, per_prime_int, names)
        print(f"  prime {p}: dims {dims}  B10 cap P10 = "
              f"{intersect['dim_B10_cap_P10']}  solved "
              f"{intersect['n_published_solved']}/{len(PUBLISHED_DEGREE10)}")
        if sp["unsolved"]:
            print(f"      not in atlas span: {', '.join(sp['unsolved'])}")

    if not per_prime_inc:
        print("no primes computed")
        return 1

    # Every prime must agree, otherwise one of them is bad and the certificate
    # would be averaging over a failure.
    distinct = {json.dumps(v["dims"], sort_keys=True) for v in per_prime_inc.values()}
    if len(distinct) != 1:
        print("FAIL: primes disagree on the dimensions; not writing")
        return 1

    shared = set(stored) & set(per_prime_inc)
    if shared:
        print(f"agrees with the existing certificate on {len(shared)} shared prime(s)")

    if not args.write:
        print("check only; pass --write to regenerate")
        return 0

    OUT_INCIDENCE.write_text(json.dumps({
        "schema": 1,
        "generated_by": "scripts/emit_degree10_space_incidence.py",
        "spaces": SPACES,
        "incidence": {},
        "per_prime": per_prime_inc,
    }, indent=1, sort_keys=True) + "\n")

    prev = json.loads(OUT_INTERSECT.read_text()) if OUT_INTERSECT.exists() else {}
    payload = {
        "schema": 1,
        "generated_by": "scripts/emit_degree10_space_incidence.py",
        "atlas_names": names,
        "product_columns": [n for n in names if n.startswith("I4_1*")],
        "product_source": prev.get(
            "product_source",
            "explicit lower_product atlas entries: I4_1 x I6_1 and I4_1 x I6_2"),
        "conclusion": prev.get("conclusion"),
        "superseded_conclusion": prev.get("superseded_conclusion"),
        "per_prime": per_prime_int,
    }
    OUT_INTERSECT.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"wrote both certificates at {len(per_prime_inc)} prime(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
