# Degree-8 span equality

## Result

At four primes --- two fitting (`32749`, `32719`) and two holdout (`32717`,
`32713`) --- on a common registry of 78 hashed self-dual five-forms:

    trace rank = 7,   spinor rank = 7,   union rank = 7

with two-way span containment in both directions, and a change-of-basis map
fitted on the fitting samples and **validated on holdout samples that took no
part in the fit**.

Union rank is the load-bearing number. Two different 7-dimensional subspaces
would give equal ranks and a union of 8; the union being 7 is what makes this
span equality rather than dimension agreement.

## Which family supplies the seventh direction

Dropping each family in turn and re-ranking the spinor side:

| family dropped | remaining spinor rank |
|---|---:|
| port graphs | 7 |
| products (`I4^2`) | 7 |
| structured `Omega`/`Theta` | 7 |
| **tensor words** | **6** |

Only the tensor words are indispensable. This locates precisely the direction
that a port-graph-only family misses: it is not reachable by any single-graph
contraction of the invariant tensor, and the `Omega`/`Theta` structured family
does not supply it either.

This is the independent confirmation of an observation made earlier from the
other side. An exact modular enumeration restricted to port graphs reaches rank
6 at degree 8; the archived enumeration, run to candidate exhaustion with
structured candidates included, reaches 7. The two numbers were consistent but
neither on its own said *which* family mattered. The drop test does.

## Why the spinor row count varies with the prime

The number of spinor rows is 23 at the fitting primes and 21 at the holdout
primes, while every rank, containment and drop-test result is identical. The
cause is benign and worth recording so it is not mistaken for prime-dependence
of the result: port graphs are drawn at random, and a graph whose evaluation
vector vanishes identically modulo `p` is dropped rather than kept as a zero row.
Two more vanish at `32713` and `32717` than at `32719` and `32749`.

A zero row contributes nothing to a span, so dropping it cannot change a rank,
a containment or which family is indispensable -- and the four primes agreeing on
all of those while disagreeing on the row count is the evidence for that.

## Scope

- Exact over `F_p`; no floating point on either side.
- Four primes. Modular agreement is not a characteristic-zero identity and no
  claim here depends on it being one.
- Span equality is a statement about evaluation on this sample registry. It is
  established by containment in both directions, not by matching dimensions.
- It is **not** a candidate-by-candidate identity: the two sides span the same
  space, and individual candidates are not claimed to correspond one to one.

## Artifact

`verification/degree8_span_equality.json` --- per prime: ranks, containments,
family drop test, registry hash, change-of-basis shape and hash, holdout
validation flag, wall time.
