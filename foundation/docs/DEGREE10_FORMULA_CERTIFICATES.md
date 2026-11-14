# Degree-10 formula certificates: the pure `N^(4125)` Q10 basis

Three explicit, indexed, reproducible contractions recovered by the
formula-independent reverse search. A canonical identifier is not mathematics;
these are the formulas.

## 1. The block, and exactly how it is built from F

Every factor below is the same tensor, `N^(4125)`, materialised by
`sdinv.stress.composite_n4125` from equation (2.17):

    N^(4125)  =  N  -  5 * antisym_{(3,4,5)}( N^(1050) )  -  T^(54)

where

    N              `composite_n` — the unprojected six-index tensor from F
    N^(1050)       `composite_n1050`, equation (2.15):
                       N^(1050)_[abc,de]f  =  Lambda^{mn}_[abc Lambda_{de]fmn}
                   with axes (0,1,2,3,4) normalised-antisymmetrised
    T^(54)         `_n_trace_54`, the 54 trace term of equation (2.17)

**Bracket ordering matters and is already applied inside the block.** The red
antisymmetrisation over axes (3,4,5) is performed *after* the five-index black
antisymmetrisation carried by `N^(1050)`, per the source convention that red
brackets act upon black ones. Nothing in the formulas below adds any further
symmetrisation: **there are no explicit BLACK or RED bracket operations, no
permutation expansion, and no free normalisation constant.** Each formula is a
plain full contraction of five copies of `N^(4125)`.

Field degree: `N^(4125)` is degree 2 in `F`, five factors, so degree **10**.

## 2. Index conventions

- `_mu` denotes a covariant (lower) slot, `^mu` a contravariant (raised) slot.
- Every index appears exactly **twice**, once raised and once lowered. That is
  what makes each expression a Lorentz scalar: an index pair with both ends
  raised, or both lowered, contracts with `delta` instead of `eta` and is not
  invariant. Only a boost detects the difference, because `delta` and `eta`
  agree on the spatial block.
- Metric `eta = diag(+,-,-,-,-,-,-,-,-,-)`; raising is `_raise_axes`, which is
  an involution for this diagonal metric, so which end of a pair carries the
  metric is arbitrary and provably cannot change the value.
- `N4125[k]` is the k-th tensor factor. All five factors are the same tensor;
  the bracketed number only labels the slot groups.

## 3. The three formulas

### R1 — first rank-gaining candidate

    einsum   abcdef,abcdgh,efijkl,ghimno,jklmno->

    N4125[0]{_mu _nu _rho _sigma _tau _alpha}
  * N4125[1]{^mu ^nu ^rho ^sigma _beta _gamma}
  * N4125[2]{^tau ^alpha _delta _epsilon _zeta _eta}
  * N4125[3]{^beta ^gamma ^delta _theta _iota _kappa}
  * N4125[4]{^epsilon ^zeta ^eta ^theta ^iota ^kappa}

### R2 — second rank-gaining candidate

    einsum   abcdef,abcdgh,efijkl,gijmno,hklmno->

    N4125[0]{_mu _nu _rho _sigma _tau _alpha}
  * N4125[1]{^mu ^nu ^rho ^sigma _beta _gamma}
  * N4125[2]{^tau ^alpha _delta _epsilon _zeta _eta}
  * N4125[3]{^beta ^delta ^epsilon _theta _iota _kappa}
  * N4125[4]{^gamma ^zeta ^eta ^theta ^iota ^kappa}

### R3 — third rank-gaining candidate

    einsum   abcdef,abcdgh,egijkl,fhimno,jklmno->

    N4125[0]{_mu _nu _rho _sigma _tau _alpha}
  * N4125[1]{^mu ^nu ^rho ^sigma _beta _gamma}
  * N4125[2]{^tau ^beta _delta _epsilon _zeta _eta}
  * N4125[3]{^alpha ^gamma ^delta _theta _iota _kappa}
  * N4125[4]{^epsilon ^zeta ^eta ^theta ^iota ^kappa}

All three share the first two factors and the `[0]-[1]` four-index bridge; they
differ in how the remaining six indices of `[0]` and `[1]` are routed into the
last three factors. R1 and R3 differ only in that routing, and R2 splits the
`[3]`/`[4]` pair differently again.

## 4. Round-trip certificate

Each formula is generated from its topology and parsed back from the printed
index pairings alone. The parse must reproduce the identical canonical key:

    topology -> render_formula -> text -> parse_formula -> topology'
    canonical_form(topology) == canonical_form(topology')

Verified for all three, and enforced for four block sectors by
`test_formula_topology_round_trip_is_exact`. A formula that did not describe
its own topology would otherwise sit here looking entirely plausible while the
identifiers stayed reassuringly consistent.

`test_rendered_formula_marks_exactly_one_raised_end_per_index` additionally
requires every index in the rendered text to appear once raised and once
lowered.

## 5. What these formulas are certified to be

**Permitted**: "A three-element compact basis for Q10 was independently
recovered in the pure `N^(4125)`x5 sector and shown to span the same quotient
space as the graph-derived and published-formula bases."

**Not claimed**:

- not universally canonical — no universal minimality theorem exists;
- not the simplest structures reaching Q10 — no larger class was swept;
- not a resolution of which published bracket reading is intended — these
  formulas *avoid* the ambiguity, which is not the same as settling it;
- not a characteristic-zero identity — certification is modular.
