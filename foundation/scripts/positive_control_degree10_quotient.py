#!/usr/bin/env python3
"""Positive control: the degree-10 projector must recover Q10 rank 3.

Why this exists
---------------
Every published equation-(4.24) candidate projected to the zero vector in Q10.
A rank-0 result is only meaningful if the projection pipeline is CAPABLE of
returning a nonzero rank -- a projector that returned zero for every input
would produce exactly the same table and mean nothing at all.

So this runs the known-good case through the identical code path. The Level-A
quotient representatives Q10_A, Q10_B, Q10_C were established independently, as
explicit F-index contraction graphs, and they span Q10 by construction. Pushing
them through the same `closure_span` -> `rref` -> `project` chain must give
three vectors of rank exactly 3.

If this prints 3, the zero result for the published candidates is a statement
about those candidates. If it printed 0, the zero result would be a statement
about a broken projector, and nothing else.

No form evaluation is needed: the Level-A representatives are atlas basis
elements, so their atlas coordinate vectors are unit vectors and the control is
pure exact linear algebra over F_p.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from sdinv.exactmap import rank_mod
from solve_intrinsic_quotients import rref, project
from stress_flow_closure import closure_span

CERT = ROOT / "results" / "stress_flow" / "certificates"
PRIMES = (32749, 32719, 32693, 32771, 32713, 32717)
SEED8 = {4: ["I4_1"], 6: ["I6_2"], 8: ["I8_3", "I8_4", "I8_5", "I8_6"]}
LEVEL_A = ROOT / "results" / "intrinsic_candidates" / "explicit_F_contractions.json"


def main():
    with LEVEL_A.open() as stream:
        classes = json.load(stream)["classes"]
    labels = [classes[k]["graph_basis_label"] for k in ("Q10_A", "Q10_B", "Q10_C")]
    print(f"Level-A Q10 representatives: {labels}")

    out = {}
    for prime in PRIMES:
        path = CERT / f"interacting_degree12_{prime}.json"
        if not path.exists():
            continue
        with path.open() as stream:
            cert = json.load(stream)
        _, bmap, span = closure_span(cert, SEED8, prime)
        names = bmap[10]
        ech, piv = rref(span[10], prime)
        free = [j for j in range(len(names)) if j not in set(piv)]

        vecs, missing = [], []
        for label in labels:
            if label not in names:
                missing.append(label)
                continue
            unit = [0] * len(names)
            unit[names.index(label)] = 1
            vecs.append(project(unit, ech, piv, free, prime))

        rank = (rank_mod(np.asarray(vecs, dtype=np.int64) % prime, prime)
                if vecs else 0)
        out[str(prime)] = {"rank": rank, "dim_Q10": len(free),
                           "vectors": vecs, "missing": missing}
        status = "OK" if rank == len(free) else "CONTROL FAILED"
        print(f"  prime {prime}: recovered Q10 rank {rank} / {len(free)}  "
              f"{status}" + (f"  missing={missing}" if missing else ""))

    ranks = {v["rank"] for v in out.values()}
    dims = {v["dim_Q10"] for v in out.values()}
    ok = len(ranks) == 1 and ranks == dims
    payload = {"schema": 1, "control": "Level-A representatives must span Q10",
               "labels": labels, "per_prime": out,
               "recovered_rank": sorted(ranks)[0] if len(ranks) == 1 else None,
               "control_passes": ok}
    target = (ROOT / "results" / "intrinsic_candidates"
              / "degree10_positive_control.json")
    target.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"\ncontrol_passes={ok}; wrote {target}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
