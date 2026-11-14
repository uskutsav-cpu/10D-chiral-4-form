# External spinor-backend interface

The attached calculation note is an independent numerical reference, not
trusted source code. The repository does not reconstruct or claim ownership of
that implementation.

When the mentor's actual source arrives, keep it in `third_party/<name>/` with
its original license and attribution, or on a separate branch if redistribution
is not permitted. Write a thin adapter implementing
`sdinv.spinor_adapter.SpinorInvariantBackend`:

```python
class MentorSpinorBackend:
    name = "mentor-spinor"
    attribution = "Author and license supplied with the original source"

    def evaluate_degree(self, five_form_components, degree, prime):
        # Input shape: (samples, 252), in lexicographic sorted 5-index order.
        # Return shape: (samples, number_of_spinor_invariants_at_degree).
        ...
```

`five_form_components` contains lower-index components in the order returned by
`sdinv.forms.basis_tuples(10, 5)`. Samples must already satisfy the Lorentzian
self-duality projector, and every returned entry must be reduced modulo
`prime`.

Compare trace and spinor bases with:

```python
from sdinv.spinor_adapter import compare_column_spaces

report = compare_column_spaces(trace_values, spinor_values, prime=32749)
assert report["equal_column_spaces"]
```

This compares exact finite-field column spaces, not column-by-column formulas.
It therefore permits a change of basis, rescaling, and a different ordering.
Repeat under the second prime 32719 and at degrees 4, 6, 8, and 10. Save the
sample seeds, prime, adapter attribution, source commit, and report with every
comparison.
