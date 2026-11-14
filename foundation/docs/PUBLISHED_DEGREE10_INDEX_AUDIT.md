# Degree-10 source reading and exact coordinate audit

The current implementation follows the original colored TeX and PDF of Cederwall, Hutomo, Kuzenko, Lechner and Sorokin, *Some remarks on invariants*, [arXiv:2509.14350v2, section 4.1.4](https://arxiv.org/html/2509.14350v2#S4.SS1.SSS4). The PDF page 25 displays the twelve candidates near equation (4.24); the source TeX labels are `I101` through `I1012`. The candidates do not have twelve separate equation numbers.

The independent source transcription was frozen before the phase 01 audit opened the implementation. That audit found missing red operations in candidates 10–12 and returned FAIL. The user then explicitly authorized a separate repair and verification phase. The original failure report and frozen artifacts are retained; their results are not relabeled as corrected-source evidence.

## Authoritative bracket programs

Let `Q[abcde;f]=N_[abc,de]f` denote the tensor returned by `composite_n1050`. Its first-five black antisymmetrizer is already applied. All brackets are normalized. Red operations follow the black operation.

| Candidate | Further source operation | Current default |
|---|---|---|
|1,2,3,6,7|No further source bracket beyond each block's definition|Unchanged|
|4|Red full symmetrization of the four free slots mu,nu,rho,lambda|`all4`|
|5|Black antisymmetrization on each of the two free M-product pairs|Unchanged|
|8|Black antisymmetrization on the indicated free M-product pair|Unchanged|
|9|Red symmetrization of nu,rho,lambda; mu is excluded|`all3`|
|10|Red antisymmetrization of the last pair on the first Q factor|`source`|
|11,12|Red antisymmetrization of the last triple on each of factors C,D,E|`source`|

Define `B[abcd;ef]=(Q[abcde;f]-Q[abcdf;e])/2` and `T[abc;def]=Alt_def Q[abcde;f]`. Then the corrected factors are `B,Q,Q,Q,Q` for candidate 10 and `Q,Q,T,T,T` for candidates 11 and 12. The triple projectors are applied to all three affected covariant factors before their selected axes are raised. Their exact contractions are preserved in the [manual transcription](../results/audit/source-review/manual-transcription-preimplementation.md).

The original TeX lines 1641–1653 contain the decisive red delimiters. This resolves the former “AMB-01” and “AMB-02” reading questions. The historical no-red defaults and pair-on-C-only variants do not implement those displayed brackets. They remain available through `LEGACY_PUBLISHED_DEGREE10` and `LEGACY_READING_VARIANTS` solely for archived reproduction; `AMBIGUITY_VARIANTS` is empty.

The [frozen index audit](../results/audit/frozen-literature/PUBLISHED_DEGREE10_INDEX_AUDIT.md) is retained byte for byte. Its sections 4 and 6 describe obsolete reading assumptions, and its claim that “AMB-02 is harmless” does not assess the corrected source programs. Its historical numerical values are not fresh verification evidence.

## Metric and normalization conventions

The metric is `diag(-1,+1,...,+1)`, epsilon_0123456789 is +1, and F=*F. Each dummy edge carries exactly one inverse metric. Moving that metric between the two ends does not change the contraction; accidentally using two raised ends with an ordinary delta contraction does.

`M_ab=F_aijkl F_b^ijkl` has no 1/4! prefactor, and `N_abc,def=F_abcij F_def^ij` has no 1/2! prefactor. The trace subtraction is

`P=N-5T-(9/28) Alt_abc Alt_def(g_ad g_be M_cf)`.

Each triple projector divides by 3!, and the double projector divides by 36. The source degree 8 maps use the literal displayed formulas: nearby prose about trace-subtracted irreducible tensors does not supply another unprinted projector. The hatted H8_2 definition is equation (4.19); equation (4.20) is its following reduction.

## Verification and proof boundary

Small direct-contraction tests independently implement the normalized source projectors on arbitrary six-slot arrays, and distinguish the corrected defaults from legacy programs. A fresh 10D source oracle independently constructs the self-dual form, M,N,Q,P, and all 24 displayed degree 8/degree10 scalars. A full 120-term first-five projector agrees with its ten-shuffle reduction. Those blocks and scalar values agree with the production implementations at a fresh pilot point. The source oracle does not import the production bracket, block, or scalar evaluator.

The [exact source certificate](../results/audit/source-review/exact-source-evaluation-certificate.json) records 14 new binary coordinate vectors evaluated at 7 primes. Rigorous denominator and absolute-value bounds make centered integer CRT recovery unique. The reconstructed exact graph evaluation matrices have ranks 7 and 14. Exact rational solves produce the [octic map](../results/order8_change_of_basis.json) and [corrected decic atlas](../results/order10_change_of_basis.json); a separately written verifier checks all 98 point files, 630 reconstructed integer entries, exact minors, matrix residuals, and ranks.

The corrected degree-10 source span has rank 12, its quotient by lower products has rank 11, and its intersection with the two-dimensional product span has dimension 1. In repository normalization the witness is

`J6 + (5/6)J5 - J2/20 - (11/1008)J1 = -(13/63) I4_1 I6_1`.

The octic forward and inverse 7×7 maps, including their product corrections, agree with the frozen artifact. The corrected decic rows 1–9 agree with the frozen map, while rows 10–12 change. Their difference has rank 3 even after discarding the two product columns. These are validation conclusions for the corrected displayed formulas, under the published dimension assumption stated below.

All 83 independently computed dense modular points also agree with the exact maps. The literal primary octic list has rank 6 with or without the quartic square; the literal hatted list has rank 5, or rank 6 after adjoining the square. The selected T1–T5,H1 list plus the square has rank 7. The hatted list is therefore not promoted to a basis merely because of the source's terminology.

The exact evaluation argument uses the invariant-space upper bounds dim degree 8<= 7 and dim degree 10<= 14 stated by the source Hilbert series in equation (4.2), together with the fact that these complete tensor contractions are invariant homogeneous polynomials. Exact full-rank graph evaluation matrices therefore establish bases of the corresponding invariant spaces. Evaluation is injective on each such space, so an exact rational solve establishes polynomial coordinates. This argument does not infer a global identity from finite samples alone, and it does not independently rederive the Hilbert series.

The repaired literal-source mapping is distinct from the graph-only rank81 certificate. A change to the displayed source candidates does not by itself invalidate the independent graph minor certificate.
