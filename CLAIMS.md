# Claim Ledger

Statuses:
- **IMPORTED** — reported by the predecessor repository; not independently recomputed here yet.
- **REPRODUCED** — independently reproduced in this repository under the same declared assumptions.
- **ESTABLISHED** — analytic/exact proof or exact certificate plus scope audit.
- **CONJECTURE** — supported but unproved.
- **OPEN** — target only.

| ID | Claim | Status | Required evidence before publication |
|---|---|---|---|
| C01 | Full graded dimensions through degree 12 are `(1,2,7,14,72)` | IMPORTED | frozen source hashes + independent atlas verification |
| C02 | Static stress dimensions through degree 12 are `(1,1,2,2,4)` | IMPORTED | exact stress-row rank certificates |
| C03 | Free-seed dynamic reachable dimensions are `(1,1,3,11,67)` | IMPORTED | exact closure certificates on fresh primes/seeds; characteristic-zero scope stated |
| C04 | Quotient dimensions are `(0,1,4,3,5)` | IMPORTED | follows only after C01 and C03 are certified in the same spaces |
| C05 | Sextic `K6` is transported but not created from `q6=0` | IMPORTED | analytic generator-exhaustion proof + exact regression |
| C06 | `q6_dot = 40 a(lambda) q6` in pinned normalization | IMPORTED | convention derivation + symbolic projection + numerical cross-check |
| C07 | A five-direction **seed** augmentation closes through degree 8 | IMPORTED | preserve wording: seed result only |
| C08 | The same five directions form a minimal generalized **generator** set | OPEN | genuine `f(T,S_i,lambda)` closure + removal minimality |
| C09 | `Q8` has a compact intrinsic tensor/representation characterization | OPEN | explicit tensor basis or dual invariant functionals |
| C10 | `Q10,Q12` admit basis-independent annihilator certificates | OPEN | exact dual-space certificates |
| C11 | The 10D ModMax-like model is pure-stress reachable | OPEN | explicit generator reconstruction OR falsified by obstruction |
| C12 | The 10D ModMax-like model is pure-stress obstructed | OPEN | explicit nonzero intrinsic obstruction OR falsified by reconstruction |
| C13 | A finite all-orders obstruction module controls reachability | CONJECTURE | proof or sharply delimited theorem |
| C14 | Dynamic near-saturation persists at degree 14 | OPEN | independent degree-14 probe |

## Forbidden shortcuts

The manuscript must not:
- call C03 an all-orders result;
- infer C08 from C07;
- call an arbitrary graph complement “intrinsic”;
- call modular equality a characteristic-zero equality without an upper-bound/lift argument;
- describe `K6` as absent from every pure-stress trajectory;
- apply polynomial theorems directly to nonanalytic conformal models without a localization argument.
