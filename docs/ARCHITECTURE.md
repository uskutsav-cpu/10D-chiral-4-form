# Architecture

```
source JSON formulas -> registry / source hashes
                                |
fixed self-dual tensor -> graph values + constrained reverse derivatives
                                |
                    HLS interacting tau[V]
                                |
              weighted trace / extra-generator products
                                |
           coefficient monomials + homogeneous scalar values
                                |
       exact basis fit with separate fit and holdout samples
                                |
          polynomial vector fields over GF(p), then optional QQ lift
                 /                 |                    \
        linear hull        finite Lie brackets      time jets
                 \                 |                    /
                   scoped obstruction / completion analysis
```

The independent exact degree-eight model additionally has a bounded-residual
map certificate, ideal tangency, constructive controls, and a catalogue
completion proof. These are separate from the generic sampled-fit path.

## Modules

- `finite_field`, `exact`, `certificates`: exact arithmetic and explicit witnesses.
- `forms`, `graphs`, `tensors`, `registry`: field representation/evaluation.
- `stress`, `weighted`: interacting stress and all eligible scalar generators.
- `fitting`, `lifting`: preserve nonlinear coupling coefficients; never turn
  each monomial into an independent control.
- `polynomial`, `reachability`, `basis_change`, `ideals`: finite control-model tools.
- `sextic`, `degree8_orbit`, `degree8_completion`, `map_proof`: concrete scoped results.
- `localization`: exact algebraic square-root patch, not a completed ModMax classifier.
- `source`, `provenance`, `checkpoints`, `pipeline`, `worker`, `cli`: operation and audit.

A pipeline sample stores its seed, prime, compact field input, basis values,
generator catalogue, coefficient monomials, and values. Fit artifacts retain
the evaluation matrix, selected independent fit rows, every holdout row,
coordinate maps, and row-space rank/annihilator witnesses. Polynomial models
use coefficient/exponent JSON, never evaluated Python expressions.

Failures generate failure receipts and do not become successful checkpoints.
Changed engine code, registry, source hashes, or configuration changes the
checkpoint namespace. A stale lock is never silently stolen. Independent
prime workers use distinct directories; avoid launching two workers onto the
same shard output.
