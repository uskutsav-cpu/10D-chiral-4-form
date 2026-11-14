# Source-overlap audit

Distinguishes **allowed mathematical overlap** from **copied implementation text**.
The distinction matters because the mathematics is not anyone's property and the
expression of it is.

## Import audit, both directions

| direction | result |
|---|---|
| shipped bridge package imports the archive | **zero** — the only occurrences of `sd5_invariants` in `spinor_trace_bridge/src/` are three prose references inside docstrings |
| tensor implementation imports the archive | **zero** |
| archive imports project code | **zero** |
| bridge imports the tensor implementation | **yes, by design** — four `importlib.import_module` calls in `traceside.py` for `sdinv.forms`, `sdinv.modp`, `sdinv.contract`, `sdinv.invariant_registry` |

Only one script in the whole tree loads the archive at all
(`run_float_jacobian_matrix.py`), and it does so at runtime from a path the user
supplies with `--archive`. Nothing importable is affected, and the shipped
package runs without the archive present.

## Mathematical overlap, declared

The bridge re-implements mathematics the archive also implements. This is
overlap of *content*, and it is unavoidable — there is only one Clifford algebra
of `so(10)`.

| construction | how the bridge obtained it | copied? |
|---|---|---|
| oscillator realisation (wedge/contraction operators) | from the mathematical description; standard exterior algebra | no |
| Chevalley pairing with the reversal sign | standard; the *consequence* (symmetry of `sigma`) is asserted by test, not inherited | no |
| port-graph contraction rule | from the described rule: `m` invariant-tensor nodes, `2m` edges | no |
| tensor words `tr(word in A, B)` | from the described construction of `M` and `N` and the derivation action | no |
| the four `sd5_spinor_degree8_*` contractions | the index specifications are read from the archive as **data**, exactly as one reads an equation from a paper | index specification only |
| the selected 83-candidate list | read as a JSON data file | data only |

## The one place where specification is read verbatim

The four `sd5_spinor_degree8_*` candidates are defined by explicit index
contractions. There is no way to evaluate "the same invariant" without using the
same index specification — that specification *is* the mathematical object. The
surrounding implementation (modular arithmetic, stepwise reduction, exact
interpolation for derivatives) shares nothing with the archive's float64
`opt_einsum` calls.

This is recorded rather than glossed because it is the closest thing to
borrowing in the whole project.

## Independent-derivation evidence

The strongest evidence that the bridge is a genuine re-derivation rather than a
transcription is that it **disagreed with itself first and had to be fixed by
mathematics, not by copying**:

- the spinor index placement was determined by requiring exact equivariance;
- the structured degree-8 evaluator initially overflowed `int64` and was caught
  by a homogeneity test, not by comparison with the archive.

And the strongest evidence that it computes the *same* objects: the bridge's
spinor-to-five-form map, built through the Lorentzian frame and an exact modular
frame congruence, agrees with the archive's independent float `pinv`-based map
**exactly, entry by entry, ratio 1, on all 25 800 nonzero components**. Two
unrelated constructions landing on the same map is a real cross-check.

## Conclusion

Implementation independence holds and is mechanically checked. **Clean-room
independence is not claimed** — no protocol prevented shared assumptions, and
the same author had both in view. An automated wording gate fails the manuscript
build if the phrase "clean-room" appears as an assertion.
