#!/usr/bin/env python3
"""Project the implemented equation-(4.24) candidates into Q10."""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from sdinv.forms import selfdual_projector, to_dense, random_form
from sdinv.contract import value
from sdinv.graphs import graph_from_label
from sdinv.exactmap import rank_mod
from sdinv.published_degree10_invariants import PUBLISHED_DEGREE10
from solve_intrinsic_quotients import rref, project
from stress_flow_closure import closure_span
from test_M_only_quotients import registry_items, evaluate_atlas_element, solve_exact

CERT = ROOT / "results" / "stress_flow" / "certificates"
PRIMES = (32749, 32719, 32693, 32771, 32713, 32717)
SEED8 = {4: ["I4_1"], 6: ["I6_2"], 8: ["I8_3", "I8_4", "I8_5", "I8_6"]}

def sample(prime, seed):
    Pr = selfdual_projector(10, 5, True, prime)
    return to_dense((Pr @ random_form(10, 5, np.random.default_rng(seed), prime)) % prime, 10, 5, prime)

def main():
    items = registry_items()
    out = {}
    for prime in PRIMES:
        path = CERT / f"interacting_degree12_{prime}.json"
        if not path.exists(): continue
        with path.open() as s: cert = json.load(s)
        _, bmap, span = closure_span(cert, SEED8, prime)
        names = bmap[10]
        ech, piv = rref(span[10], prime)
        free = [j for j in range(len(names)) if j not in set(piv)]
        n = len(names) + 8
        forms = [sample(prime, 41000 + 11*i) for i in range(n)]
        cols = [[evaluate_atlas_element(items[nm], items, f, prime, {}) for f in forms] for nm in names]
        A = [[cols[j][i] for j in range(len(names))] for i in range(n)]
        proj = {}
        for cid, spec in PUBLISHED_DEGREE10.items():
            b = [spec["evaluator"](f, prime) for f in forms]
            x, ok = solve_exact(A, b, prime)
            if not ok:
                proj[cid] = {"status": "not_in_atlas_span"}; continue
            q = project(x, ech, piv, free, prime)
            proj[cid] = {"status": "solved", "quotient_vector": q,
                         "nonzero": any(v % prime for v in q)}
        vecs = [p["quotient_vector"] for p in proj.values() if p.get("status") == "solved"]
        qr = rank_mod(np.asarray(vecs, dtype=np.int64) % prime, prime) if vecs else 0
        out[prime] = {"projections": proj, "Q10_rank_from_published": qr,
                      "dim_Q10": len(free)}
        print(f"  prime {prime}: Q10 rank from implemented P10 = {qr} / {len(free)}")
        for cid, p in proj.items():
            print(f"      {cid}: {p.get('status')} q={p.get('quotient_vector')} nonzero={p.get('nonzero')}")
    ranks = {v["Q10_rank_from_published"] for v in out.values()}
    payload = {"schema": 1, "degree": 10, "implemented": list(PUBLISHED_DEGREE10),
               "not_implemented": 10, "per_prime": {str(k): v for k, v in out.items()},
               "consistent": len(ranks) == 1,
               "Q10_rank_from_implemented_published": ranks.pop() if len(ranks) == 1 else None}
    p = ROOT / "results" / "intrinsic_candidates" / "published_degree10_map.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"wrote {p}")

if __name__ == "__main__":
    main()
