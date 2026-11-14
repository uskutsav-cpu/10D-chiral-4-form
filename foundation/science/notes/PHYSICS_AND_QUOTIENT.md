# Physics interpretation and limits of the invariant classification

## Algebraic field space

Work at one spacetime point with η=diag(−1,1,...,1), a chosen orientation,
and ε_0123456789=+1. For a covariant five-form,

    (*F)_{μ1...μ5} = ε_{μ1...μ5ν1...ν5} F^{ν1...ν5}/5!.

There are C(10,5)=252 alternating components. The Lorentzian Hodge identity is
*²=(−1)^(5*5+1)=+1. In an increasing-index basis, Hodge exchanges a five-tuple
with its distinct complementary five-tuple; hence its trace is zero and its
±1 eigenspaces each have 126 dimensions.

For every increasing tuple I containing 0, let A_I=F_I. Its complement J has
no 0, and self-duality fixes F_J=−sgn(J,I)A_I. Thus

    F(A) = sum_{I contains 0} A_I(e_I+*e_I),  dim V_+=C(9,4)=126.

This is an integral coordinate parametrization, not a gauge choice for the
potential. No differential constraint was used in this pointwise count.
The contraction F·F vanishes because F∧*F=F∧F=0 for an odd form.

## Lorentz scalars and the orbit dimension

A graph vertex means one copy of F. Every edge contracts one lower index
from each endpoint with η^{-1}; parallel edges give distinct pairs. Under a
Lorentz transformation, the two matrices on each edge cancel against the
inverse metric by Λ^TηΛ=η. Therefore each complete contraction is an SO(1,9)
scalar. Fixed slot orders retain the signs required by the alternating form.
Hodge self-duality is preserved by orientation-preserving Lorentz matrices.

The Lie algebra consists of the 36 spatial rotations and9 boosts, dimension 45.
For X in that algebra its action is the sum of X acting in each of the five
covariant tensor slots. The standalone verifier constructs this action on all
126 coordinates, with no assumption that all 45 directions are independent.
At its fresh point the 45×45 orbit minor is 10547 modulo 50021, so that integral
minor polynomial is not identically zero. The orbit rank is therefore 45
on a nonempty Zariski-open set. The Lie stabilizer is zero there; this does
not assert that the entire stabilizer group is trivial or exclude a finite
stabilizer.

Every invariant differential annihilates these 45 tangent directions, giving
rank at most 126−45=81. The same standalone calculation obtains a nonzero
81×81 Jacobian minor,21539 modulo 50021, from the 81 fixed graph contractions.
Both nonzero modular minors lift to characteristic-zero polynomial statements.
Thus the invariant polynomial algebra has transcendence degree 81. This proof
does not use the Hilbert calculation, source-paper maps, or candidate search.

## What the quotient means

Where the orbit action has rank 45 and the invariant map G=(G1,...,G81) has
rank 81, its differential has kernel equal to the orbit tangent space. The
submersion theorem and a transverse local slice give 81 coordinates for a
local space of Lorentz orbits. A smooth invariant function restricted to a
sufficiently small such neighborhood can be written as f(G1,...,G81).

This is local classification of algebraic field configurations up to a
change of Lorentz frame. Lorentz transformations are not being declared an
extra gauge symmetry that removes 45 propagating degrees of freedom. The
number 81 is neither the number of dynamical polarizations nor a global
orbit-separation theorem. Disconnected fibers, singular strata and global
quotient topology are not resolved by the nonzero-minor argument. Nor does
it say every polynomial invariant is a polynomial in this 81-list.

## Connection with the chiral four-form

A four-form potential has gauge transformation A4→A4+dΛ3 and field strength
F5=dA4. The Bianchi identity dF5=0 is a differential identity. In the free
chiral theory one imposes F5=*F5; then its ordinary source-free field equation
d*F5=0 follows from the Bianchi identity. Pointwise self-duality still leaves
126 algebraic components; the differential equations and gauge equivalences
are additional structure. For the free massless field the transverse
little-group eight-space gives C(8,4)=70 four-form polarizations before
chirality and 35 after the middle-form duality split. Neither is 81 or 126.

For a general field strength define F5^+=(F5+*F5)/2. In a nonlinear chiral
formulation, one must not identify the physical full F5 with F5^+ off shell
without its constitutive/self-duality equation. Some formulations instead
use an auxiliary self-dual five-form Λ5. The algebraic classification applies
to either self-dual variable, when that is the variable of the chosen
formulation. Interactions depending on F5 have no derivatives of the field
strength, although F5 itself contains derivatives of A4.

Hutomo–Lechner–Sorokin distinguish the physical full field strength from an
auxiliary self-dual variable in their INZ-type construction, and relate its
interaction function to other chiral formulations. See [sections 1–2 and 4–5
of their primary article](https://arxiv.org/html/2509.14351v2). The role of our
81 scalars is to supply explicit local arguments for Lorentz-invariant
functions of that self-dual variable; it is not to replace those formulation-
specific field equations or equivalence proofs.

## What has and has not been established physically

On the regular local quotient the 81 graph invariants are local building
blocks for Lorentz-invariant nondifferential functions of a self-dual
five-form. This statement has the analytic argument and exact witnesses
above. Promoting a chosen function to a chiral action requires its complete
formulation, gauge symmetries and constitutive equations. The invariant
classification alone proves neither PST gauge invariance, nonlinear
self-duality consistency, conformal homogeneity, duality conditions,
hyperbolicity, positivity nor global existence of solutions.

In particular, inserting an arbitrary function directly into a guessed
physical-field action is not justified by Lorentz scalarity. An established
auxiliary-field construction may permit a class of arbitrary invariant
potentials under its own hypotheses; those hypotheses must be carried with
that construction. No new interacting action or equations of motion are
claimed here.

The stress-flow and spinor investigations are frozen as related follow-up
work at commit b732ae66034762cd2a9fc55068851687ee56fcb1. Closing mentor gates
G-1...G-10, the physical interpretation of Q10, and the Tr(tau) degree
assumption is not a premise or claimed achievement of this rank 81 project.
Likewise a full invariant-ring presentation and syzygies are outside scope.
