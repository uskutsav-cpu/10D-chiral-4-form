# The frame orientation was not pinned — found via prime 32707

## What was observed

A 15-cell rank matrix run left three cells failing, all at `p = 32707`, at every
seed:

    RuntimeError: image has dimension 0, not 126; cannot build a left inverse

Twelve cells at the other four primes completed at rank 81 with a clean schedule.
A prime-specific failure at every seed is not a sampling accident.

## What it was not

Not an exceptional prime. 32707 is the only one of the six with both `-1` and `2`
a non-residue, which looked like a promising explanation and is a coincidence.

Not a degenerate construction either: at 32707 the self-dual and anti-self-dual
bases both still have dimension 126, and the forward map still has rank 126.

## What it was

The forward map's roles are **swapped** at 32707:

| prime | rank of image of self-dual | rank of image of anti-self-dual |
|---|---:|---:|
| 32749 | 126 | 0 |
| 32719 | 126 | 0 |
| 32713 | 126 | 0 |
| **32707** | **0** | **126** |

The null frame is built by a congruence `L^T eta_null L = eta_lorentzian`, which
is determined only up to the **square-root branch** chosen inside `congruence`.
That branch is not cosmetic: it reverses the frame's orientation, the orientation
reverses the sign of the Hodge star, and the sign of the star decides which
eigenspace the gamma map annihilates.

The decisive experiment: flipping the branch at 32707 restores the expected
behaviour, and flipping it at a *working* prime breaks that prime in exactly the
same way.

    p=32707  as-built            self-dual image   0   anti 126
    p=32707  other sqrt branch   self-dual image 126   anti   0
    p=32749  as-built            self-dual image 126   anti   0
    p=32749  other sqrt branch   self-dual image   0   anti 126

So the construction contained an unpinned convention. At five of six primes
`sqrt_mod` happened to return the branch matching the stated convention; at
32707 it returned the other one.

## The fix

`orientation_normalised_L` tries both branches in a fixed order and keeps the one
under which the self-dual space survives — the convention the package states
everywhere. Deterministic, and it raises rather than guesses if neither branch
works, which would be a genuinely exceptional prime.

All six primes now behave identically, 32707 included.

## What this does and does not change

**Does not change any published number.** Every committed certificate was
computed at 32749, 32719, 32717, 32713 or 32693, and at each of those the old
code already selected the correct branch. The recomputed values are unchanged.

**Does change what the construction is.** Before, "kernel is the anti-self-dual
126, image is the gamma-traceless 126" was true at the primes used and would have
been false at others, with no error until a later step failed. The statement is
now true by construction at any prime where the construction exists at all.

## Regression tests

`spinor_trace_bridge/tests/test_orientation.py`:

- the self-dual space survives, and the anti-self-dual space is annihilated, at
  all six primes including 32707;
- a left inverse exists at all six;
- **the opposite branch is asserted to break it**, at 32707 and at 32749. Without
  that third test a future refactor could drop the normalisation and the suite
  would still pass at five of six primes.

## Credit

The failing cells came from a rank-matrix run in a separate working tree
(`~/Downloads/sdinv-jhep`). The diagnosis and fix are recorded here; the failure
itself was surfaced by that run.
