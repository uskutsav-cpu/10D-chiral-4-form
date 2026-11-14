# Adversarial validation

17 mutation tests, each tied to a defect class that either occurred in this
project or would have been silent if it had. Every one **injects** a specific
error and asserts a specific check catches it. A mutation test that merely runs
proves nothing; if one of these ever stops failing on the mutated input, the
guard it protects has stopped working, and the test says so in its own assertion
message.

## Defects reproduced from this project's own history

| test | the defect |
|---|---|
| `test_old_defect_int64_overflow_is_caught_by_homogeneity` | The structured degree-8 evaluator once reduced mod `p` only at the end of each einsum. With `p ~ 2^15` a four-step contraction reaches `~2^72` and wraps `int64` silently; the wrapped values look entirely ordinary. The test performs the *unreduced* contraction deliberately and asserts homogeneity fails — which is how it was found — then asserts the fixed path passes. |
| `test_double_raising_M_changes_the_tensor_word_values` | `M_mu{}^nu` is a **mixed** tensor. Raising the already-upper index gives a different object. The test builds the doubly-raised version and requires at least one tensor word to change. |
| `test_checkpoint_key_includes_the_formula_hash` | A cache keyed only by candidate id would serve a stale row after the candidate's definition changed. |

## Defects that would have been silent

| test | what it guards |
|---|---|
| `test_row_normalisation_of_a_modular_matrix_is_meaningless` | Over `F_p` every nonzero row scales to leading entry 1, so "normalisation" carries no information. This is the step that inflated the float64 rank to 83. |
| `test_composite_modulus_is_rejected_or_detected` | Fermat-based inversion is simply wrong for a composite modulus. |
| `test_mutating_the_raising_permutation_breaks_the_metric_check` | Index raising is verified against the metric rather than assumed. |
| `test_transposing_N_changes_mixed_traces` | The `N` variance convention is load-bearing, not decoration. |
| `test_wrong_hodge_sign_destroys_the_kernel_identification` | Inverting the duality sign would send self-dual forms to zero. |
| `test_scaling_gamma_normalisation_changes_the_forward_map` | A different Clifford normalisation is a different map. |
| `test_split_signature_is_not_euclidean` | Asserts the `(5,5)` eigenvalue split **and** that the frame diagonal is zero, so a Euclidean misreading cannot return. |
| `test_forward_inverse_mismatch_is_detected` | `inverse . forward` is the self-dual projector, **not** the identity — the identity would contradict the kernel being the anti-self-dual part. |
| `test_edge_mutation_changes_a_port_graph_value` | Swapping two half-edges must change the contraction. |
| `test_word_order_matters_for_non_cyclic_rearrangement` | `AABB` and `ABAB` are different invariants and must not collapse. |
| `test_necklaces_do_not_double_count_cyclic_rotations` | `AAAB` and `AABA` are the same trace; only one may be counted. |
| `test_candidate_reordering_does_not_change_the_rank` | Rank must not depend on schedule order. |
| `test_degree8_span_equality_is_not_mere_dimension_agreement` | Constructs two genuinely different 7-spaces and confirms the union test separates them, so the check applied to real data is known to have teeth. |

## Test totals

- bridge suite: **72 passed** (55 structural + 17 adversarial), at two primes
- tensor suite: **199 passed**
- run from clean bytecode

## What adversarial testing does not cover

Mutation testing shows the checks are sensitive to the errors injected. It does
not show they are sensitive to errors nobody thought of, and it is not a proof of
correctness. The independent evidence that carries more weight is that the
bridge's spinor-to-five-form map agrees with a separately written third-party
implementation exactly, entry by entry, on all 25 800 nonzero components.
