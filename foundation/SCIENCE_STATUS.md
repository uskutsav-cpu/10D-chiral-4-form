# Completed science scope — September 2026

The independent character calculation closes the remaining Hilbert-count
premise. This completion adds mathematics, computation and a physics derivation;
no manuscript was changed. Earlier completed tensor evaluations were reused
with their exact hashes, rather than rerun.

| Requested gate | Result and evidence |
|---|---|
| Independent invariant dimensions | D5 symmetric characters give 1, 2, 7, 14, 72 at degrees 4, 6, 8, 10, 12; series through 22 retained in [character derivation](science/character/DERIVATION.md). |
| Degree-12 completeness | Independent upper count 72 meets the audited rank-72 certificate: ten products and 62 graph classes are a homogeneous basis. |
| Final source mathematics | All 24 corrected formulas and rational maps remain frozen; the external Hilbert premise is now discharged. Exact decic ranks are 12, 2, 13, intersection 1. See [closure proof](science/notes/DEGREEWISE_AND_SOURCE_CLOSURE.md). |
| Canonical independent rank theorem | [One-file verifier](science/rank81/README.md) computes F, all 81 values, the full Jacobian, and an integer minor from only graph definitions, convention and coordinates. Fresh residue 21539 modulo 50021. |
| Independent generic orbit | Same fresh run constructs 45 Lorentz tangent directions, checks self-duality and Jacobian annihilation, and obtains orbit minor 10547 modulo 50021. Thus transcendence degree is 126−45=81. |
| Physics interpretation | [Derivation note](science/notes/PHYSICS_AND_QUOTIENT.md) explains 126 components, metric invariance, the local quotient, A4/F5, and the limits of the interaction claim. |
| Stress-flow and spinors | Preserved as related follow-up at b732ae6; mentor gates and Q10 interpretation are not claimed closed or required here. |
| Full ring and syzygies | Explicitly outside this classification scope. |

The C++ engine exhaustively computes all character entries through degree 22.
Independent Python arbitrary-precision verification exhaustively checks the
recurrence through degree 12, samples 128 entries at each higher degree, and
checks every degree's full dimension and invariant extraction. This coverage
is explicit; no formal proof-assistant verification is claimed.

The standalone tensor run was executed cold in a directory initially containing
only its program and input, in about 56 seconds. It uses no stored Jacobian,
search, checkpoint, source map or Hilbert data. Seven new adversarial/provenance
tests pass. Commands, exit codes and output hashes are under
[science/evidence/commands](science/evidence/commands).

Run `python science/verify_science.py` for hash binding and saved-arithmetic
closure checks. Full tensor recomputation and character regeneration have
separate commands in their respective README/derivation files. The saved-data
checker does not claim to perform those full computations.

The historical failed frozen audit and its repaired follow-up remain intact
in [INDEPENDENT_AUDIT.md](INDEPENDENT_AUDIT.md). Its external-Hilbert caveat is
superseded by this new derivation, not silently rewritten. Canonical old source
maps and both copied evidence certificates are unchanged byte for byte.

A portable science bundle and a smaller standalone rank bundle are frozen under
`release/science/`, with SHA-256 sidecars. The science bundle contains all new
character tables, programs, notes, minimal rank inputs and certificates, plus
the two previously verified exact certificates needed for logical closure.
Detailed original source transcriptions and evaluator histories remain in the
repository's preserved independent audit.
