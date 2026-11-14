# Physics convention audit

Primary source: Hutomo--Lechner--Sorokin, arXiv:2509.14351v2,
https://arxiv.org/html/2509.14351v2, especially equations (2.1), (2.33),
(2.36), and the distinction between the auxiliary form and physical field strength.

The self-dual variable used by the invariant atlas is the auxiliary five-form
Lambda (or its on-shell identification with the self-dual part of the physical
field strength). The full nonlinear physical field strength is not assumed to
satisfy linear self-duality. In the unified manuscript, F denotes this auxiliary
self-dual variable.

The implemented convention is fixed by the full stress formula (2.33):

    tau = 48 T = M - 25 B(DV,DV) + 1 * sum_d (d-2) V_d.

Here D is the constrained anti-self-dual derivative, the Lorentz metric is
(-,+,...,+), and contraction of the homogeneous term gives F.DV_d=d V_d.
The traces of M and B(DV,DV) vanish because each is constructed from an odd
self-dual or anti-self-dual middle form. Direct contraction therefore yields

    tr(tau) = 10 * sum_d (d-2) V_d.

This calculation fixes the trace sign from (2.33). The trace printed in the
paper's (2.36) has the opposite overall sign relative to that contraction;
it is not silently substituted for the implemented expression. The code and
manuscript identify the chosen equation explicitly. This is an algebraic
convention audit, not an independent derivation of the full metric variation.

With independently signed controls, replacing an entire generator by its
negative changes its control parameterization, not its reachable orbit.
That observation does not permit changing individual components or the
relative signs in the full stress tensor.

The paper's quartic Tr(M^2) equals 2*I4_1 in the graph basis. The two tensor
engines agree on graph definitions, Hodge conventions, electric coordinates
and the selected 81 invariant values. The independent character count then
applies to the same self-dual representation. This connects the invariant
foundation to the flow map's homogeneous dimension premise.

Nonanalytic square-root interactions, higher field-strength derivatives,
physical causality/unitarity, and supersymmetric/Type-IIB completion require
separate arguments. The current finite polynomial classification asserts
none of those results.
