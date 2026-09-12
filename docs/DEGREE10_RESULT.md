# Degree-10 nonlinear stress-flow orbit

## Scope

This result applies to the finite polynomial truncation through field degree
ten, starting from the free seed, with piecewise signed independent scalar
stress controls in the declared analytic generator class.

It is not an all-orders theorem and does not cover nonanalytic ModMax-type
generators, positivity restrictions, or causality/unitarity conditions.

## Starting coordinate space

The cumulative invariant-coordinate dimensions are

- degree 4: 1
- degree 6: 2
- degree 8: 7
- degree 10: 14

for 24 cumulative coordinates.

The already-certified degree-eight orbit imposes six relations, leaving an
18-dimensional reduced space.

## Exact degree-10 obstruction system

An exact characteristic-zero search in the complete weight-80 ansatz found
12 independent obstruction relations.

Their degree-ten leading parts have rank 12.

Thus the invariant variety has dimension at most

    18 - 12 = 6.

The pure-stress control distribution at the free seed has rank 6, so

    dim O_{<=10} = 6

locally.

## Constructive sufficiency

The twelve obstruction relations solve globally for twelve dependent
degree-ten coordinates in terms of six free parameters

    A = c_I4_1
    B = c_I6_1
    C = c_I8_1
    D = c_I4_1^2
    U = c_I10_1
    V = c_I4_1*I6_1.

Every stress vector field is tangent to this six-dimensional graph.

The lower four coordinates are reached sequentially with

    tr2
    tr3
    tr4
    tr2*tr2.

The remaining two degree-ten directions are constant translations generated
by

    tr2*tr3
    tr5,

whose translation determinant is

    -6185484288 / 27097

and is nonzero.

Therefore every point on the free-seed component of this six-dimensional
degree-ten obstruction variety is reachable by an explicit six-stage
piecewise pure-stress control law.

## Important distinction

The linear invariant hull has dimension 8.

It is not the nonlinear orbit.

The nonlinear orbit dimension through degree ten is 6 in the stated model.

## Evidence

See `verification/degree10/` for:

- five-prime rational reconstruction,
- a completely unused validation prime,
- the exact 12-relation obstruction certificate,
- the accessibility analysis,
- and the constructive global-reachability certificate.
