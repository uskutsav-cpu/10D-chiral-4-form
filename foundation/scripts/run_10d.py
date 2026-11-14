"""10D self-dual (chiral) 5-form F^+_{mu1..mu5}.

Target: 81 functionally independent invariants (Hutomo-Lechner-Sorokin,
arXiv:2509.14351). The degree-by-degree partition function in
arXiv:2509.14350v2 fixes the complete low-order target at 1 quartic,
2 sextic, and 6 new octic generators. This script finds an explicit
contraction-graph basis for all nine.

Do NOT trust output from this script unless scripts/run_6d.py passes first.
"""

import sys, os, time, json, argparse, hashlib
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from sdinv.modp import P, ALT_P, RankSieve
from sdinv.forms import (selfdual_projector, to_dense, random_form,
                         check_star_squared)
from sdinv.graphs import enumerate_graphs, graph_label, load_graph_catalog
from sdinv.contract import (build_basis_flat, contraction_plan_cost,
                            jacobian_row, value)

D, PDEG = 10, 5
TARGET = 81
KNOWN_NEW_GENERATORS = {4: 1, 6: 2, 8: 6}
DEFAULT_CATALOG = os.path.join(
    os.path.dirname(__file__), "..", "results", "10d_graph_catalog.json")


def preflight(mod=P):
    c = check_star_squared(D, PDEG, True, mod)
    assert c == 1, f"star^2 = {c}, expected +1; self-duality is ill-posed"
    Pr = selfdual_projector(D, PDEG, True, mod)
    assert np.array_equal((Pr @ Pr) % mod, Pr % mod), "projector not idempotent"
    s = RankSieve(Pr.shape[1], mod)
    for r in Pr:
        s.add(r)
    assert s.rank == 126, f"self-dual subspace has rank {s.rank}, expected 126"
    return Pr


def _graphs_for_order(order, catalog_path):
    if catalog_path and os.path.exists(catalog_path):
        return load_graph_catalog(catalog_path, order)
    if order <= 6:
        return enumerate_graphs(order, PDEG, max_mult=PDEG - 1)
    raise RuntimeError(
        "order 8 requires the exact graph catalog; run "
        "scripts/generate_graph_catalog.py first")


def run(orders, mod=P, seed=20260727, catalog_path=DEFAULT_CATALOG):
    Pr = preflight(mod)
    rng = np.random.default_rng(seed)
    Fd = to_dense((Pr @ random_form(D, PDEG, rng, mod)) % mod, D, PDEG, mod)

    # quadratic invariant must vanish identically: F ^ F = 0 for odd-degree
    # forms, and self-duality gives F ^ *F = F ^ F, hence F.F = 0.
    g2 = enumerate_graphs(2, PDEG, max_mult=PDEG)
    assert value(g2[0], Fd, D, PDEG, True, mod) == 0, \
        "F.F did not vanish -- self-dual projection is wrong"

    basis = build_basis_flat(D, PDEG, Pr, mod)
    sieve = RankSieve(basis.shape[0], mod)
    log = {
        "prime": mod,
        "seed": seed,
        "catalog": os.path.relpath(catalog_path) if catalog_path else None,
        "orders": {},
        "generators": [],
    }

    for n in orders:
        t0, before = time.time(), sieve.rank
        graphs = _graphs_for_order(n, catalog_path)
        indexed_graphs = list(enumerate(graphs))
        if n >= 8:
            indexed_graphs.sort(
                key=lambda item: (contraction_plan_cost(item[1], D, PDEG),
                                  item[0]))

        expected_new = KNOWN_NEW_GENERATORS.get(n)
        evaluated = 0
        for catalog_index, M in indexed_graphs:
            evaluated += 1
            if sieve.add(jacobian_row(M, Fd, basis, D, PDEG, True, mod)):
                generator = {
                    "id": f"I{n}_{sieve.rank - before}",
                    "order": n,
                    "catalog_index": catalog_index,
                    "graph": graph_label(M),
                }
                log["generators"].append(generator)
                print(f"    NEW  {generator['id']} = {generator['graph']}")
            if expected_new is not None and sieve.rank - before == expected_new:
                break

        found = sieve.rank - before
        log["orders"][str(n)] = {
            "catalog_graphs": len(graphs),
            "evaluated": evaluated,
            "expected_new": expected_new,
            "new": found,
            "rank": sieve.rank,
            "complete": expected_new is not None and found == expected_new,
            "seconds": round(time.time() - t0, 1),
        }
        print(f"  order {n}: evaluated {evaluated}/{len(graphs)} graphs -> "
              f"+{found} new (rank {sieve.rank}/{TARGET}) "
              f"[{log['orders'][str(n)]['seconds']}s]")
        if expected_new is not None and found != expected_new:
            raise RuntimeError(
                f"order {n}: found {found}, expected {expected_new} from "
                "the published partition function")

    return sieve.rank, log


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--orders", type=int, nargs="+", default=[4, 6, 8])
    ap.add_argument("--catalog", default=DEFAULT_CATALOG)
    ap.add_argument("--single-prime", action="store_true")
    ap.add_argument("--out", default="results/10d_order8.json")
    a = ap.parse_args()
    if any(order >= 10 for order in a.orders):
        # Preserve the historical low-degree regression runner. The generic
        # runner includes lower generators and products and keeps the two
        # rank notions separate, which is essential at degree 12.
        from run_degree import main as degree_main
        output = a.out if a.out != "results/10d_order8.json" else (
            f"results/degree{max(a.orders)}_verification.json")
        arguments = ["--degree", str(max(a.orders)), "--out", output]
        if a.single_prime:
            arguments += ["--primes", str(P)]
        raise SystemExit(degree_main(arguments))
    print(f"10D self-dual 5-form: 126 components, target {TARGET}")
    primes = [P] if a.single_prime else [P, ALT_P]
    logs = []
    for mod in primes:
        print(f"\nprime {mod}")
        _, log = run(a.orders, mod, catalog_path=a.catalog)
        logs.append(log)

    generator_labels = [
        [item["graph"] for item in log["generators"]] for log in logs]
    if len(logs) == 2:
        assert generator_labels[0] == generator_labels[1], (
            "the two primes selected different bases")

    with open(a.catalog, "rb") as stream:
        catalog_sha256 = hashlib.sha256(stream.read()).hexdigest()
    result = {
        "schema": 1,
        "claim": (
            "Complete through order 8: 1 quartic, 2 sextic, and 6 new "
            "octic generators; cumulative functional rank 9."
        ),
        "literature": {
            "source": "https://arxiv.org/abs/2509.14350v2",
            "partition_function_through_order_8": "1 + t^4 + 2 t^6 + 7 t^8",
            "new_generators": {"4": 1, "6": 2, "8": 6},
            "note": "The seventh degree-8 scalar is the product I4^2.",
        },
        "catalog": os.path.relpath(a.catalog),
        "catalog_sha256": catalog_sha256,
        "generators": logs[0]["generators"],
        "degree8_basis": [
            {
                "id": item["id"],
                "kind": "connected_generator",
                "graph": item["graph"],
            }
            for item in logs[0]["generators"] if item["order"] == 8
        ] + [{
            "id": "I4_1^2",
            "kind": "composite",
            "expression": "I4_1^2",
        }],
        "runs": {str(log["prime"]): log for log in logs},
    }
    parent = os.path.dirname(a.out)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(a.out, "w") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    highest_order = max(a.orders)
    print(f"\nrank through order {highest_order}: "
          f"{logs[0]['orders'][str(highest_order)]['rank']}")
