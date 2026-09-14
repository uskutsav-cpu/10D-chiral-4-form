# Maximal all-orders program: gates 1–11

This bundle is intentionally fail-closed. It automates all computations that
can be justified from the current code/data and emits explicit blocker
certificates where a new theorem is genuinely required.

1. Unconditional degree-10 bounded-residual physics/basis map.
2. Promote the exact degree-10 orbit ideal to an unconditional all-orders ideal.
3. Attempt an explicit degree-12 orbit ideal from the stable equation-space lift.
4. Record degree-12 tangency status.
5. Promote degree 12 only if exact ideal/tangency is genuinely available.
6. Enumerate/rank-select a 247-dimensional degree-14 contraction basis using nauty `genbg`.
7. Fit the degree-14 pure-stress map on fresh primes with disjoint holdouts.
8. Compute the degree-14 obstruction space and, when a compatible full-ring
   lower ideal exists, the quotient modulo lower-ideal prolongations.
9. If degree 14 introduces no new generator, repeat the falsification probe at
   degree 16 (published homogeneous dimension 1364).
10. Never infer finite generation from finitely many zero quotients. A
    machine-checkable structural induction/module certificate is required.
11. Never infer orbit equality from obstruction containment. A compatible
    formal control-lifting/sufficiency certificate is required.

Published external homogeneous dimensions used as theorem inputs:
4:1, 6:2, 8:7, 10:14, 12:72, 14:247, 16:1364.

The higher-basis enumerator represents each loopless 5-regular contraction
multigraph as a bipartite incidence graph. `genbg` canonically enumerates the
bicoloured incidence graphs with left degree 5 and edge-vertex degree 2.

Important: stages 6–9 are potentially very compute-intensive and are resumable
only where the underlying repo checkpoint layer applies. The basis rank
selection itself writes progress records and fails closed on resource limits.
