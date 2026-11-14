# Independent D5 character calculation

The calculation reads no graph, fitted matrix, target Hilbert coefficient or
published character table. It computes the whole symmetric-power weight
character, rather than determining the invariant count from graph ranks.

## The representation and its weights

Complexifying the vector representation gives weights ±e_i, i=1,...,5.
Enumerate the C(10,5)=252 exterior basis vectors of its fifth exterior power.
A weight has either five, three or one nonzero coordinates, each ±1:

- Five nonzero coordinates: all 32 sign choices, multiplicity 1.
- Three nonzero coordinates: choose the three single indices and signs,
  then one cancelling pair from the other two indices; multiplicity 2.
- One nonzero coordinate: choose two cancelling pairs among the other four
  indices; multiplicity 6.

Hodge duality splits this middle exterior power into two 126-dimensional
modules. An orientation reversal in a zero-weight coordinate exchanges the
Hodge summands while fixing the weight, so the multiplicities with zeros
split equally. The highest-weight vector of weight (1,1,1,1,1) belongs to one
summand. Its D5 Weyl orbit consists of the 16 five-sign vectors with even
negative parity. This determines the chosen summand's complete weight list:
16 full-support weights of multiplicity 1,80 support-three weights of
multiplicity 1,10 support-one weights of multiplicity 3. Their total is 126.

Its highest weight is 2ω5, Dynkin labels (0,0,0,0,2). The other Hodge summand
has 2ω4 and is obtained by reversing one coordinate. This outer automorphism
preserves singlet multiplicities, so using the dual representation in the
polynomial coordinate ring gives the same answer. The Python checker derives
the list from all 252 exterior basis vectors independently of the C++ support
formula, and verifies the Weyl dimension formula gives 126.

The real Lorentz representation has this complexification, up to the choice
of Hodge summand. Taking kernels of the Lie-algebra action commutes with
extension of scalars R→C. These tensor representations factor through SO,
and the identity component of the real Lorentz group is Zariski dense in the
connected complex SO(10). Thus the polynomial invariant dimensions computed
here apply to SO(1,9); this is not an assumption about arbitrary smooth
functions on disconnected global orbit spaces.

## Symmetric powers without decomposition or graph data

Let c(w) be the input multiplicities and m_n(λ) the multiplicity of a weight
λ in Sym^n(V). Differentiate the formal product

    H(t,x) = product_w (1-t*x^w)^(-c(w)).

Its logarithmic derivative gives the exact Newton recurrence

    n*m_n(λ) = sum_{k=1}^n sum_w c(w)*m_{n-k}(λ-k*w),
    m_0(0)=1.

Every weight coordinate in Sym^n(V) has absolute value at most n; the sum of
coordinates has parity n. The code enumerates every dominant weight within
these necessary bounds, including weights whose answer is zero. There is no
truncation by a guessed orbit support or by a literature coefficient.

D5 acts by arbitrary permutations and even sign changes. A unique dominant
representative is a≥b≥c≥d≥|e|, with the sign of e recording the product of
nonzero signs when there are no zeros. The multiplicity is constant on each
orbit. Replacing a lookup by this representative therefore loses no terms.
The character is stored only after division by n is checked exact and its
multiplicities are checked nonnegative.

At every degree the sum of multiplicities times their exact Weyl-orbit sizes
must equal C(125+n,n). The C++ engine and Python checker both verify this
full-dimension identity. The complete nonzero orbit tables are retained in
`independent-d5-weights.tsv`, including odd degrees needed in the recurrence.

## Extracting the trivial representation

For D5, the positive roots are e_i−e_j and e_i+e_j for i<j, and
ρ=(4,3,2,1,0). The Weyl character formula implies that the constant term of

    ch(Sym^n V) * product_{α>0}(1-x^(-α))

is the trivial multiplicity: a nontrivial dominant highest weight cannot
have λ+ρ Weyl-conjugate to ρ. Equivalently,

    h_n = sum_{w in W(D5)} det(w)*m_n(ρ-wρ).

The engine enumerates all 120 permutations and16 even sign changes,1920 terms.
The checker independently multiplies the 20 positive-root factors, obtains
1920 nonzero monomials with coefficients ±1, and extracts the same constant
terms from the saved characters. This independently checks denominator signs,
not merely the final invariant counts.

The standard theorem is stated in [Etingof, Lie Groups and Lie Algebras,
Theorem 27.15 and Corollary 27.16](https://math.mit.edu/~etingof/lnlg.pdf).
No paper-specific Hilbert coefficient enters this argument or implementation.

## Results and checks

| Degree | 0 | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 | 22 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Dimension |1|0|1|2|7|14|72|247|1364|6851|40170|227979|

Odd degrees vanish. The output was compared only afterwards with the partition function in section 4 of
[Some remarks on invariants](https://arxiv.org/html/2509.14350v2); all listed
coefficients through 22 agree. The result is a truncated Hilbert series, not
a rational expression for the entire series and not a minimal ring presentation.

C++ computes every character coefficient through 22 exactly using signed 128
integers. Python arbitrary-precision arithmetic independently verifies every
Newton-recurrence entry, including missing zero entries, through degree 12;
for each degree 13–22 it checks 128 deterministic sampled entries. It checks
all degrees' complete dimensions and invariant extractions. The exhaustive
higher-degree computation is C++; higher-degree Python recurrence coverage
is explicitly sampled. The retained conservative intermediate bound
1920*126*22*C(147,22) is less than 2^127, covering Newton sums, binomial updates,
dimension sums and the alternating Weyl sum before cancellation.

Assertions must remain enabled. These are exact computational proofs using
standard character theory, not formal proof-assistant verification.

## Reproduce

From this directory:

    clang++ -std=c++17 -O3 -Wall -Wextra d5_symmetric_character.cpp -o /tmp/d5_character
    /tmp/d5_character 22 /tmp/d5
    python check_character.py /tmp/d5 --full-through 12 --output /tmp/d5-check.json

Python uses only its standard library. Passing `--full-through 22` additionally
checks every higher-degree recurrence entry in Python, if desired. Exact
commands, compiler version and hashes are recorded in `../evidence/`.
