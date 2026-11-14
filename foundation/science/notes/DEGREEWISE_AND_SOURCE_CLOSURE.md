# Degreewise completeness and source identities

The independent D5 calculation now supplies the upper bounds that were external
premises of the earlier source audit. The earlier certificate is preserved
byte for byte in `../evidence/source-certificate.json`; its historical conditional
wording describes the evidence available then. This note discharges that
Hilbert-count premise using `../character/DERIVATION.md` and the complete tables.
Analytic invariance, standard character theory and the checked arithmetic
implementations remain the stated mathematical and computational foundations.

| Degree | Invariant dimension | Products of lower degrees | New basis classes |
|---|---:|---:|---:|
| 4 | 1 | 0 | 1 |
| 6 | 2 | 0 | 2 |
| 8 | 7 | 1 | 6 |
| 10 | 14 | 2 | 12 |
| 12 | 72 | 10 | 62 |

The ten degree-12 products are I4^3, the six I4*I8_i, and the three
symmetric products I6_i*I6_j. The prior independent stacked-gradient
certificate has 72 independent rows for those ten products and the 62
explicit graph polynomials. Its two minor residues are 43334 modulo 50021
and 12699 modulo 32749. A nonzero modular minor proves characteristic-zero
linear independence of those polynomials. Together with dim I12=72 this
proves a complete homogeneous basis, including a 62-dimensional quotient
by products. It does not assert 72 algebraically independent degree-12
polynomials. The original computation is reused, not rerun.

For degrees 8 and 10, the exact graph evaluation matrices have ranks 7 and
14. Since these equal the independently derived dimensions, evaluation is
injective on each full invariant space. Every transcribed source contraction
is an invariant in the respective space. Therefore the exact rational matrix
solutions in the source certificate are polynomial identities, not just
identities at sample points. Integer reconstruction is unique because its
CRT modulus exceeds twice every proved numerator bound. All 12 displayed
octic and all 12 displayed decic formulas and their exact rational maps
are retained in that certificate and the original source audit. No source
formula or matrix has been replaced during this closure step.

The corrected source convention applies the red antisymmetrizers after the
black ones and before raising indices: J10 antisymmetrizes the trailing pair
of its first Q, while J11 and J12 antisymmetrize the trailing triples of their
last three Q factors. Normalized antisymmetrizers are essential.

Let L be the span of the twelve corrected decics and P the two-dimensional
product span. Exact rational elimination gives dim L=12 and dim(L+P)=13,
so dim(L intersect P)=12+2-13=1. An explicit generator is

    J6 + (5/6)J5 - J2/20 - (11/1008)J1 = -(13/63) I4_1 I6_1.

Thus the source decics supply eleven primitive classes; even after adding
both products they miss one dimension of I10. This is a statement about the
displayed source family, not a failure of the fourteen-dimensional graph basis.
The twelve displayed octics have the exact maps already retained; the
selected six source octics together with I4^2 form the seven-dimensional basis.
The six primary octics alone span six dimensions and already contain I4^2;
the six hatted octics span five and reach six after adjoining I4^2.

The standalone rank-81 theorem has no dependency on these degree counts or
source identities. Full invariant-ring generation, higher syzygies, a rational
Hilbert series and a globally separating invariant set are not established.
Stress-flow and mentor spinor identifications remain related follow-up work,
not prerequisites for the completed classification scope.
