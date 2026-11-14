"""6D validation: generic 3-form H_{mu nu rho} in d=6.

Target (Elamaran-Ferko-Scarlett, arXiv:2512.23750): 5 independent
invariants, appearing 1,2,1,1 at orders 2,4,6,8.

THIS IS THE GATE. If this does not print 5, nothing downstream is
trustworthy. Run it after every change to the core.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.modp import P, ALT_P, RankSieve
from sdinv.forms import basis_tuples, to_dense, random_form
from sdinv.graphs import enumerate_graphs, graph_label
from sdinv.contract import jacobian_row, build_basis_flat

D = 6
PDEG = 3
ORDERS = [2, 4, 6, 8]


def run(mod=P, seed=20260727, verbose=True):
    rng = np.random.default_rng(seed)
    ncomp = len(basis_tuples(D, PDEG))
    Fvec = random_form(D, PDEG, rng, mod)
    Fd = to_dense(Fvec, D, PDEG, mod)
    basis = build_basis_flat(D, PDEG, None, mod)

    sieve = RankSieve(ncomp, mod)
    per_order = {}
    for n in ORDERS:
        before = sieve.rank
        graphs = enumerate_graphs(n, PDEG)
        kept = []
        for M in graphs:
            row = jacobian_row(M, Fd, basis, D, PDEG, True, mod)
            if sieve.add(row):
                kept.append(graph_label(M))
        per_order[n] = sieve.rank - before
        if verbose:
            print(f"  order {n}: {len(graphs):4d} graphs -> "
                  f"+{per_order[n]} new (running rank {sieve.rank})")
            for lab in kept:
                print(f"        NEW  {lab}")
    return sieve.rank, per_order


if __name__ == "__main__":
    print(f"6D generic 3-form, components = {len(basis_tuples(D, PDEG))}")
    print(f"[prime {P}]")
    r1, pat1 = run(P)
    print(f"[prime {ALT_P}] confirmation")
    r2, pat2 = run(ALT_P, verbose=False)

    print()
    print(f"independent invariants : {r1}   (confirm: {r2})")
    print(f"pattern by order       : "
          f"{[pat1[n] for n in ORDERS]} at orders {ORDERS}")
    print(f"expected               : 5 and [1, 2, 1, 1]")
    ok = (r1 == 5 and r2 == 5 and [pat1[n] for n in ORDERS] == [1, 2, 1, 1])
    print("RESULT: " + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)
