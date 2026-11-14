# Publication claim allocation

Two manuscripts exist. This file states which claim belongs to which, so that
neither is a duplicate publication of the other.

## The Letter (`manuscript/prl/`) --- NOT YET SUBMISSION-READY

**Allocated claim.** The degree-ten obstruction: `dim A10 = 14`,
`dim D10 = 11`, `dim Q10 = 3`, `P10 subset D10`, three representatives, removal
minimality, and the resulting requirement on any deformation scheme.

Nothing else. Method is compressed to what is needed to believe the theorem.

**Status.** The PRL claim gate does not pass (`audit/PRL_CLAIM_GATE.md`):
Impact and Broad interest fail on present evidence, Innovation is conditional on
a literature sweep. The draft is complete and labelled accordingly.

## The long form (`manuscript/`) --- complete

**Allocated claims.** Everything else, and the full method for the above:

- the complete degree-ten atlas and its construction;
- the full incidence table of `A10`, `B10`, `G10`, `P10`, `D10`;
- the correction that the published span is *not* a product complement;
- the pure-`N` compact basis and its certificates;
- the formula-independent reverse recovery;
- the bridge: real forms, exact Clifford, equivariance, left inverse;
- the exact Jacobian and the rank-81 certificate;
- degree-8 and degree-10 common-sample comparisons;
- implementation, reproduction and the degree-12 outlook.

## If both were ever released

The Letter's theorem would also appear in the long form, which is normal for a
Letter plus a companion methods paper --- but only if the long form is a
**companion preprint**, not a second journal submission of the same result.
Publication ethics on that point are a human decision and are recorded as
outstanding.

Given that the PRL gate does not pass, the current position is simpler: the long
form is the vehicle, and the Letter is a draft held back.

## Synchronisation

Both manuscripts read their numbers from the same generator
(`manuscript/scripts/make_numbers.py`), which writes into both trees. A value
cannot drift between them, and a missing artifact shows as a loud marker in both.
