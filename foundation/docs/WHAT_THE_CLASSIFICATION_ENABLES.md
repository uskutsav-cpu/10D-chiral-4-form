# What the classification enables

A short, concrete list. Each item is a computation that is now finite and exact,
and was not previously well posed.

| question | before | now |
|---|---|---|
| Is a given degree-ten contraction a new invariant, or a combination of known ones? | comparable only against twelve expressions of unknown completeness | expand in a certified 14-dimensional basis; exact coordinates |
| Is a given degree-ten term reachable by the stress flow? | no way to ask | compute its class in `Q10`; reachable iff the class vanishes |
| How many parameters does a general degree-ten deformation need? | unknown | 14, of which 11 are reachable and 3 are not |
| Are the published structures independent? | unstated | they span 12 dimensions, and contain one product direction |
| How many *non-product* degree-ten structures does the published list contain? | assumed 12 | **11** |
| Do the tensor and spinor descriptions agree at degree ten? | untested at the level of values | evaluation spans coincide on a common sample registry, holdout-validated |
| Is the numerical rank-81 evidence trustworthy? | one float64 sample, tolerance-sensitive | replaced by an exact modular lower bound with no tolerance |

## The one that matters most for future work

`P10 ⊂ D10`. Because every degree-ten product is reachable, the three missing
directions cannot be reached by composing lower-degree results. Any extension of
the construction that hopes to cover degree ten must introduce genuinely new
structure, not merely iterate.
