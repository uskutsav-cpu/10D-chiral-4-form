#!/usr/bin/env python3
"""Project the equation-(4.25) structures into Q12. Bounded prime set."""
import argparse, json, sys, time
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from sdinv.forms import selfdual_projector, to_dense, random_form
from sdinv.exactmap import rank_mod
from sdinv.published_degree12_invariants import PUBLISHED_DEGREE12
from sdinv.invariant_registry import load_verified_registry_through_degree12
from solve_intrinsic_quotients import rref, project
from stress_flow_closure import closure_span
from test_M_only_quotients import evaluate_atlas_element, solve_exact

CERT = ROOT / "results" / "stress_flow" / "certificates"
SEED = {4: ["I4_1"], 6: ["I6_2"], 8: ["I8_3","I8_4","I8_5","I8_6"],
        10: ["I10_6","I10_7","I10_12"]}

def sample(prime, seed):
    Pr = selfdual_projector(10, 5, True, prime)
    return to_dense((Pr @ random_form(10,5,np.random.default_rng(seed),prime))%prime,10,5,prime)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--primes", type=int, nargs="+", default=[32749, 32717])
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()
    reg = load_verified_registry_through_degree12(str(ROOT))
    items = {}
    for d in reg.degrees:
        for it in reg.basis(d): items[it.id] = it
    print(f"registry items: {len(items)}", flush=True)
    out = {}
    for prime in a.primes:
        t0 = time.time()
        with (CERT / f"interacting_degree12_{prime}.json").open() as s: cert = json.load(s)
        _, bmap, span = closure_span(cert, SEED, prime)
        names = bmap[12]
        ech, piv = rref(span[12], prime)
        free = [j for j in range(len(names)) if j not in set(piv)]
        n = len(names) + 8
        forms = [sample(prime, 52000 + 13*i) for i in range(n)]
        print(f"  prime {prime}: evaluating {len(names)} x {n} atlas values", flush=True)
        cols = []
        for k, nm in enumerate(names):
            cols.append([evaluate_atlas_element(items[nm], items, f, prime, {}) for f in forms])
            if (k+1) % 12 == 0: print(f"    {k+1}/{len(names)} cols  {time.time()-t0:.0f}s", flush=True)
        A = [[cols[j][i] for j in range(len(names))] for i in range(n)]
        proj = {}
        for cid, spec in PUBLISHED_DEGREE12.items():
            b = [spec["evaluator"](f, prime) for f in forms]
            x, ok = solve_exact(A, b, prime)
            if not ok: proj[cid] = {"status": "not_in_atlas_span"}; continue
            q = project(x, ech, piv, free, prime)
            proj[cid] = {"status": "solved", "quotient_vector": q,
                         "nonzero": any(v % prime for v in q)}
        vecs = [p["quotient_vector"] for p in proj.values() if p.get("status")=="solved"]
        qr = rank_mod(np.asarray(vecs,dtype=np.int64)%prime, prime) if vecs else 0
        out[prime] = {"projections": proj, "Q12_rank_from_published": qr,
                      "dim_Q12": len(free), "seconds": round(time.time()-t0,1)}
        print(f"  prime {prime}: Q12 rank = {qr} / {len(free)}  ({time.time()-t0:.0f}s)", flush=True)
        for cid,p in proj.items():
            print(f"      {cid}: {p.get('status')} nonzero={p.get('nonzero')}", flush=True)
    ranks = {v["Q12_rank_from_published"] for v in out.values()}
    payload = {"schema":1,"degree":12,"structures":list(PUBLISHED_DEGREE12),
               "per_prime":{str(k):v for k,v in out.items()},
               "consistent":len(ranks)==1,
               "Q12_rank_from_published": ranks.pop() if len(ranks)==1 else None,
               "note":"tr(M^6) rational graph coordinates are NOT reconstructed; only modular quotient statements are certified."}
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(json.dumps(payload, indent=1, sort_keys=True)+"\n")
        print(f"wrote {a.out}")

if __name__ == "__main__":
    main()
