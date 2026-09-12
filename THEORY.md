# Theory and proof contract

## Fields, normalization, and allowed flows

Use `eta=diag(-,+,...,+)`, lower `epsilon_0123456789=+1`, real self-dual
`Lambda_5=*Lambda_5`, and the constrained anti-self-dual derivative
`D=dV/dLambda^upper`, normalized by `Lambda^I D_I=d V_d` on a homogeneous term.
The computation uses HLS arXiv:2509.14351v2 equation (2.33):

```
tau = 48 T = M - 25 B(D,D) + identity * sum_d (d-2)V_d
M_mu^nu = Lambda_mu,rho4 Lambda^nu,rho4
B(D,D)_mu^nu = D_mu,rho4 D^nu,rho4
Tr(tau) = sum_d 10(d-2)V_d.
```

The trace sign is obtained from (2.33). The printed sign in (2.36) differs;
confirm the convention with the authors rather than silently switching it.

The flow law is `dV/dlambda = sum_g u_g(lambda) g(tau[V])`, with arbitrary
signed scalar controls and `g` polynomial trace monomials. There are no field
derivatives, negative/fractional powers, field-dependent external controls,
or constant vacuum interaction in this declared class. The interaction starts
at field degree four. Analytic functions at the zero-field point give the same
finite generator catalogue after truncation. We do not infer positivity,
causality, supersymmetry, or a UV completion from these algebraic assumptions.

## Three different geometric objects

1. `A_d`: a finite vector space of homogeneous Lorentz invariant **functions
   of Lambda**. Its basis dimension uses separate invariant-theory input.
2. `c`: coefficients of those basis functions, coordinates on a finite
   **space of interactions**. Coupling-space vector fields are polynomial in c.
3. `Orbit(0)`: the set reachable from zero couplings by the allowed vector
   fields. This need not be a vector subspace and can satisfy nonlinear
   equations between coefficients of different field degrees.

The minimal linear flow-invariant hull of the free seed contains the orbit,
but can strictly exceed it. The standard counterexample is `x'=1,y'=x`:
the orbit satisfies `y=x^2/2` while its linear hull is two-dimensional.
Coefficient-by-coefficient forcing spans can overestimate even the linear
hull by treating dependent coefficient monomials as independent controls.

Never define an orbit quotient `A_d/R_d` without proving `R_d` is the desired
linear object. Linear annihilators describe a supplied row span only. A
nonlinear obstruction such as `b2+16*a^3` belongs to the coupling coordinate
ring, not to the ring of tensor invariants `I_d(Lambda)`.

## Exactness levels

- A nonzero finite-field minor proves a characteristic-zero rank lower bound
  only when the matrix is the reduction of the stated rational/integer map.
- Agreement at extra seeds or primes is strong regression evidence, not an
  unconditional identity proof.
- A rational reconstruction bound guarantees uniqueness among bounded
  candidates. It does not prove the unknown physical coefficients have that
  height. `lifting.py` therefore labels its output a candidate.
- `map_proof.py` instead bounds integer residual values at **fixed bounded
  tensors**, eliminates each residual by a sufficiently large CRT modulus,
  and uses an injective evaluation matrix with external Hilbert upper bounds.
- Exact ideal tangency and constructive controls prove statements about the
  supplied finite model. Extending them to the full physical problem requires
  the map, conventions, admissible control class, and any higher-order argument.

## Basis and minimality

The 81 generic functional invariants are not automatically a polynomial
presentation of the full ring. The registry keeps every homogeneous graph and
product slot rather than replacing the ring with 81 free variables.

A linear change of invariant basis induces a contragredient change of
couplings. `basis_change.py` transforms the full polynomial vector fields and
obstruction ideal; it does not declare individual graph pivots intrinsic.

Adding a new homogeneous scalar to `f(tau,S)` gives an independent generator
and its permitted mixed products. Adding it only to `V(0)` is a different
operation. Catalogue-minimal completion, removal nonredundancy, and global
minimality over arbitrary scalars are different claims.
