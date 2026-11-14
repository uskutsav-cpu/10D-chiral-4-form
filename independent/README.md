# Independent orbit verification

Run from the combined project root:

```sh
python independent/verify.py --output reports/independent-algebra.json
python independent/verify.py --target '2,-3/7,5,-11/13,1/3,-2/5,7/11,-1,4/9,2/3' --output reports/endpoint-example.json
```

The checker imports neither `chiral4form` nor the bundled foundation. It reads
coefficient/exponent JSON and independently implements the orbit calculation
using SymPy exact rational arithmetic. It checks degree-eight transport and
catalogue completion, all degree-ten tangency relations, all 36 degree-twelve
source rows and 68 kernel constraints, and the degree-twelve vertical matrix.
The endpoint command integrates terminating Lie series to produce ten signed
control segments reaching arbitrary rational values of six base and four free
fiber coordinates. Every resulting constraint is checked exactly.

The four fiber coordinates are the columns reported in the verification JSON;
currently they are `I4_1^3`, `I6_1^2`, `I4_1*I8_1`, and `I12_58` (indices
0, 1, 4, 67 of the homogeneous degree-twelve basis). These are **couplings** to
invariant polynomials, not values of the invariant polynomials at a field.

This verifies reduced-model algebra, not its tensor-physics derivation. Fresh
bounded physics-map evaluation is a separate command. The regression suite
contains sign, omitted-generator, source-tail and kernel corruption cases. No
formal proof-assistant verification is claimed.

Classify an arbitrary coefficient dictionary (unlisted coefficients are zero):

```sh
python independent/verify.py --classify examples/combined/unreachable-octic.json --output reports/unreachable-example.json
python independent/verify.py --target '1,0,0,0,0,0,0,0,0,0' --output reports/reachable-example.json
python independent/verify.py --classify reports/reachable-example.json --output reports/reachable-classification.json
```

Inputs use exact invariant IDs and rational strings or integers. Floating JSON
numbers are rejected. An unreachable interaction reports nonzero obstruction
residuals; a reachable interaction returns ten signed control segments and its
complete 96-coordinate polynomial. The classification is for the declared
finite reduced model; its tensor-physics interpretation has the separate map
proof obligations described in the science status.
