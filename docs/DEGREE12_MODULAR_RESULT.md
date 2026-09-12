# Degree-12 modular stress-flow obstruction structure

## Scope

This checkpoint concerns the polynomial pure-stress system through
field degree twelve after restriction to the previously certified
six-parameter degree-ten orbit.

It is exact over eight independently tested finite fields.

It is not yet a characteristic-zero theorem.

## Degree-12 obstruction ansatz

After restriction to the degree-ten orbit, the cumulative model has

- 6 certified lower orbit coordinates,
- 72 degree-twelve coordinates,

for 78 reduced coordinates.

A degree-twelve obstruction has tr1 eigenvalue 100.

The complete lower-coordinate weight-100 monomial catalogue consists of

- A^5
- A^3 B
- A^2 C
- A^2 D
- A B^2
- A U
- A V
- B C
- B D

Together with the 72 degree-twelve coordinates, this gives an
81-dimensional obstruction ansatz.

## Eight-prime result

The computation was independently performed over

30203
30211
30223
30241
30253
30259
30269
30271.

Every prime produced exactly the same structural result:

- nonzero coefficient equations: 36
- equation rank: 13
- obstruction nullity: 68
- degree-twelve leading rank: 68
- direct free-seed control rank: 10
- obstruction-variety dimension upper bound: 10

Therefore, over every tested finite field,

    10 <= dim O_{<=12} <= 78 - 68 = 10.

Hence the modular dimension sandwich closes at 10.

The same 68 leading pivot columns occur at all eight primes.

## Interpretation

This is strong prime-independent evidence that the pure-stress orbit
dimension grows

    4 -> 6 -> 10

through degrees

    8 -> 10 -> 12.

However, the generic characteristic-zero rational reconstruction has
not yet passed independent-prime validation.

Accordingly, the repository does NOT yet promote

    dim O_{<=12} = 10

to an exact QQ theorem.

## Rational reconstruction issue

Both the generic 96-coordinate lift and the model restricted to the
degree-ten orbit failed independent-prime rational validation.

The invariant obstruction structure itself nevertheless stabilizes
exactly across all eight tested primes.

The next step is to lift a canonical basis of the 68-dimensional
obstruction space rather than reconstructing the full vector field.

## Evidence

See

    verification/degree12/modular/

for the full eight-prime obstruction scan, modular kernels, pivot data,
prime computation summaries, reconstruction diagnostic, provenance and
hash manifest.
