"""The cardinality half of PO-08.

`docs/PROOF_OBLIGATIONS.md` proves that any set closing a degree has at least
`dim Q` elements, where `Q` is the quotient of the atlas by the reachable
subspace. The proof is analytic and basis-free; these tests do not attempt to
establish it. What they do is guard the two things the proposition takes as
input, either of which could silently drift in an artifact:

  * that the recorded deficit really is `dim A - dim D` at that degree, so the
    quantity the bound speaks about is the one the certificates record;
  * that the exhibited set really does close the degree, which is the
    hypothesis of the proposition.

Given both, the proposition says `|S| >= deficit`, and the exhibited sets meet
it with equality -- so they are of minimum cardinality in *any* basis. That last
step is what these tests then check, and it is the only thing here that would
fail if someone shrank an exhibited set without redoing the mathematics.
"""

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CERTS = {
    10: ROOT / "results/generalized_flow/degree10_missing_directions.json",
    12: ROOT / "results/generalized_flow/degree12_missing_directions.json",
}


def load(degree: int) -> dict:
    path = CERTS[degree]
    if not path.exists():
        pytest.skip(f"{path.relative_to(ROOT)} absent")
    return json.loads(path.read_text())


@pytest.mark.parametrize("degree", sorted(CERTS))
def test_deficit_is_the_quotient_dimension(degree):
    """The recorded deficit is the codimension of the reachable subspace."""
    cert = load(degree)
    key = str(degree)
    base = cert["base_closure"][key]
    full = cert["closure_with_all"][key]
    assert cert["deficit"] == full - base, (
        f"degree {degree}: deficit {cert['deficit']} does not equal "
        f"{full} - {base}; the bound would then be about a different number"
    )


@pytest.mark.parametrize("degree", sorted(CERTS))
def test_exhibited_set_closes_the_degree(degree):
    """The hypothesis of the proposition: D + span(S) = A."""
    cert = load(degree)
    key = str(degree)
    assert cert["degree_closed"] is True, f"degree {degree} is not recorded as closed"
    atlas_dim = {10: 14, 12: 72}[degree]
    assert cert["closure_with_all"][key] == atlas_dim, (
        f"degree {degree}: adjoining the exhibited set reaches "
        f"{cert['closure_with_all'][key]}, not the full atlas dimension "
        f"{atlas_dim}, so the proposition's hypothesis does not hold"
    )


@pytest.mark.parametrize("degree", sorted(CERTS))
def test_exhibited_set_meets_the_lower_bound(degree):
    """|S| >= dim Q always; here |S| == dim Q, so S is of minimum cardinality.

    This is the statement that is basis-independent. Removal-minimality is a
    different and stronger property and is NOT tested here -- it remains open
    under general GL, as PO-08 records.
    """
    cert = load(degree)
    exhibited = cert["missing_directions"]
    deficit = cert["deficit"]
    assert len(exhibited) >= deficit, (
        f"degree {degree}: {len(exhibited)} directions cannot close a "
        f"codimension-{deficit} gap; this would contradict the proposition, so "
        f"one of the artifacts is wrong"
    )
    assert len(exhibited) == deficit, (
        f"degree {degree}: {len(exhibited)} exhibited directions against a "
        f"deficit of {deficit}. The set still closes the degree but is no "
        f"longer of minimum cardinality, so the minimality wording in "
        f"docs/CLAIM_LEDGER.md (C-MIN-02) must be narrowed"
    )


@pytest.mark.parametrize("degree", sorted(CERTS))
def test_bound_is_confirmed_at_several_primes(degree):
    """The dimensions feeding the bound are modular, so record the primes."""
    cert = load(degree)
    assert len(cert["primes_confirmed"]) >= 2, (
        f"degree {degree}: the closure dimensions are confirmed at only "
        f"{cert['primes_confirmed']}"
    )
