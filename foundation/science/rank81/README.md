# Standalone rank and orbit certificate

Copy only `verify.py`, `input.json` and `requirements.txt` into an empty directory.
With Python 3.13 and dependencies installed, run:

    python verify.py input.json --output fresh-certificate.json

The verifier reconstructs the self-dual form from 126 coordinates, evaluates
all 81 graph contractions and every entry of their 81 by 126 Jacobian,
selects a nonsingular minor, and computes its integer Bareiss determinant.
It separately constructs all 45 infinitesimal Lorentz actions, checks their
self-duality and annihilation by the Jacobian, and certifies an orbit minor.
No existing certificate is an input. Assertions must remain enabled.

`input.json` freezes only the convention, ordered graph definitions and one
exact point. Metric signature is (-,+,...,+); the 126 independent components
are the lexicographic five-index tuples containing zero. Complement signs
are explicitly reconstructed by the Hodge convention in the implementation.
Graph slots follow increasing neighbor order. Each edge contracts with the
inverse metric exactly once. These polynomial definitions have integer
coefficients, so a nonzero determinant modulo a prime proves a nonzero
characteristic-zero determinant polynomial. The Jacobian criterion proves
algebraic independence; the 45-dimensional orbit bounds transcendence degree
by 126-45=81. See `../notes/PHYSICS_AND_QUOTIENT.md` for local versus global scope.

The fresh cold run used p=50021, seed 12959865808634581237. Its rank-81 minor
is 21539 modulo p and its orbit-45 minor is 10547. The complete output is
`certificate.json`; cold execution took about 56 seconds on the recorded host.
The seed is provenance; the explicit coordinates are the actual inputs.

The contraction functions were copied from the previously independent audit;
`implementation_provenance.json` records their source hashes. The new run
reuses no numerical result. NumPy contractions use checked exact integer-valued
floating operations below 2^53 and modular reduction; opt_einsum chooses paths.
An in-memory 256 MB contraction cache lasts only for the current point. There
are no checkpoints, persistent intermediate arrays, repository imports,
search tools, literature maps or Hilbert data in this verifier.
