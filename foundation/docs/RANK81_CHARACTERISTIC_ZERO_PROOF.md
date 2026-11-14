# Characteristic-zero rank: statement and proof

## What is proven

**Lower bound, unconditional.** `rank_Q J >= 81`, certified by an explicit
`81 x 81` minor.

**Upper bound, from representation theory with a stated hypothesis.**
`rank <= 81`.

**Together.** The generic functional dimension is exactly `81` — with the two
halves resting on completely different kinds of argument, which is the point of
separating them.

## The lower bound

### Why the minor is integral

The Jacobian is not merely *computed* mod `p`; it is the **reduction of a genuine
integer matrix**. Two facts make that true:

1. `integral.integral_gamma_traceless_basis()` returns a basis of the
   126-dimensional gamma-traceless space over `Z`, with entries in `{-1,0,+1}`,
   annihilated exactly over `Z` by the integral gamma-trace constraints. No
   division and no prime enters its construction.
2. The sample point is an integer combination of that basis, and every candidate
   evaluator is a polynomial with integer coefficients in the components of the
   point (integral `sigma` matrices, integral invariant tensor, integral
   derivation structure).

So there is an integer matrix `J_Z` with `J_Z mod p = J_p` for every prime.

### Why one prime suffices

We need the minor to be **nonzero**, not its value. For an integer matrix,

    det(M_Z) = 0   =>   det(M_Z) = 0 (mod p)   for every prime p.

Contrapositive: a single prime at which `det(M_p) != 0` proves `det(M_Z) != 0`,
hence `rank_Z M_Z = 81`, hence `rank_Q J >= 81`.

This is why no Hadamard bound, no height estimate and no CRT reconstruction is
needed. Those would be required to recover the *value* of the determinant; they
are not required to certify that it is nonzero.

### The certificate

`results/rank81/minor81_certificate.json` records:

- the 81 candidate identifiers whose rows form the minor;
- the 81 coordinate column indices;
- the sample seed and the prime at which the minor was selected;
- `det = 20755` modulo `32749`, nonzero;
- the same determinant from **two independent routines** — modular LU with
  pivoting and fraction-free Bareiss elimination carried out mod `p` — which
  agree. These use genuinely different recurrences, so a mistake in one is
  unlikely to be mirrored in the other;
- the SHA-256 of the minor itself.

Additional primes are checked not because the rank argument needs them but
because agreement across primes also catches an indexing error in selecting the
minor, which a single prime cannot.

## The upper bound

This half is **not** ours and is not computational.

Let `V` be the 126-dimensional module of self-dual five-forms and `G = SO(1,9)`,
`dim g = 45`. If the generic stabiliser is trivial, the generic orbit has
dimension 45, so the generic fibre of the quotient map has dimension 45 and the
number of functionally independent invariants is

    dim V - dim(generic orbit) = 126 - 45 = 81.

**Hypothesis, stated because it is doing real work:** the generic stabiliser of a
self-dual five-form in `D=10` is trivial. This is the literature input
(Ref. arXiv:2509.14351), not something established here. If the generic
stabiliser were positive-dimensional the count would rise, and the upper bound
would fail.

We do not re-derive this. We cite it, and the manuscripts attribute the number
81 to it explicitly.

## What is NOT proven

**Genericity.** The computations establish rank 81 at specific, explicitly
recorded sample points. That the rank is 81 at a *generic* point does not follow
from finitely many points alone; it follows from the upper bound above together
with the lower bound here, and the upper bound is the cited analytic statement.

Formally: the rank function is lower semicontinuous, so rank `>= 81` on a
Zariski-open set containing the tested points; combined with `rank <= 81`
everywhere, the generic rank is 81. The semicontinuity step is standard; the
`<= 81` step is the citation.

**Independence of the 83 candidates.** Rank 81 with 83 candidates means there are
**two functional dependencies** among the selected functions. The candidates are
not algebraically independent and no claim here says they are.

**Anything about a different point.** A non-generic point can have smaller rank.
None was found among the points tested, but the search was not exhaustive and no
exceptional locus is characterised.

## Permitted wording

> An explicit 81x81 minor of the integral Jacobian is nonzero, certifying
> `rank_Q >= 81`. With the cited generic-stabiliser result giving `rank <= 81`,
> the generic functional dimension is exactly 81.

## Prohibited wording

> We proved the generic rank is 81. *(Only the lower half is ours.)*
> There are 83 algebraically independent invariants. *(Rank 81 with 83
> candidates means the opposite.)*
> Modular agreement proves the characteristic-zero identity. *(For the rank
> lower bound the argument is integrality, not agreement across primes.)*
