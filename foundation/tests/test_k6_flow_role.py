"""The dynamical role of the intrinsic sextic quotient class K6.

The static map says degree six splits as a one-dimensional free-stress span
<J6> plus a one-dimensional quotient spanned by [K6]. That is a statement
about a *static* span and says nothing on its own about what a *flow* can
reach.

This module pins the dynamical statement, read off the exact interacting
certificates. At field degree six only three trace generators can contribute
(leading_field_degree <= 6): tr_tau, tr_tau2, tr_tau3. Every degree-six
target is therefore enumerated, and exactly one of them carries a nonzero
I6_2 coordinate -- the one whose coefficient monomial is I6_2 itself.

Consequence: writing q6 for the intrinsic quotient coordinate,

    d q6 / d lambda = 10*(6-2) * a_tr_tau(lambda) * q6,

a linear HOMOGENEOUS ODE, so q6(lambda) = q6(0) * exp(40 * int a_tr_tau).
A seed with q6(0) = 0 keeps q6 identically zero for all lambda, and a seed
with q6(0) != 0 keeps it nonzero. K6 is neither universally forbidden nor
spontaneously generated: it is *transported*, never *created*.
"""

import json
import os
from pathlib import Path
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE_DIR = ROOT / "results" / "stress_flow" / "certificates"
PRIMES = (32749, 32719, 32717, 32693)

# Tr(tau) acts on a homogeneous degree-d piece as 10*(d-2); at d=6 that is 40.
EULER_FACTOR_AT_DEGREE_6 = 40


def _certificates():
    out = []
    for prime in PRIMES:
        path = CERTIFICATE_DIR / f"interacting_degree12_{prime}.json"
        if not path.exists():
            continue
        with path.open() as stream:
            out.append(json.load(stream))
    assert len(out) >= 3, "need at least three interacting certificates"
    return out


def test_degree6_generator_set_is_exhaustive():
    """Only generators with leading_field_degree <= 6 can reach degree six."""
    for certificate in _certificates():
        eligible = {
            generator["id"]
            for generator in certificate["trace_generators"]
            if generator["leading_field_degree"] <= 6
        }
        assert eligible == {"tr_tau", "tr_tau2", "tr_tau3"}

        seen = {
            target["generator"]
            for target in certificate["targets"]
            if target["field_degree"] == 6
        }
        assert seen == eligible, (
            "a degree-six generator produced no target row; the degree-six "
            "analysis would then be incomplete")


def test_k6_is_transported_never_created():
    """Exactly one degree-six target has a nonzero K6 coordinate, and it is
    the one already proportional to K6."""
    for certificate in _certificates():
        prime = certificate["prime"]
        nonzero = {}
        for target in certificate["targets"]:
            if target["field_degree"] != 6:
                continue
            index = target["basis"].index("I6_2")
            value = target["coordinates"][index] % prime
            assert target["holdout_passed"] is True
            if value:
                nonzero[target["id"]] = value

        assert list(nonzero) == [f"tr_tau|d=6|c=I6_2"], (
            f"prime {prime}: K6 is reachable from a direction that does not "
            f"already contain it -- {sorted(nonzero)}")
        assert nonzero["tr_tau|d=6|c=I6_2"] == EULER_FACTOR_AT_DEGREE_6


def test_free_seed_never_reaches_k6():
    """From a seed with no K6 component the degree-six flow has none."""
    for certificate in _certificates():
        prime = certificate["prime"]
        # A seed free of I6_2 activates only these coefficient monomials.
        seeded = {"", "I4_1", "I6_1"}
        total = 0
        for target in certificate["targets"]:
            if target["field_degree"] != 6:
                continue
            monomial = "".join(target["coefficient_monomial"])
            if monomial not in seeded:
                continue
            index = target["basis"].index("I6_2")
            total = (total + target["coordinates"][index]) % prime
        assert total == 0, (
            f"prime {prime}: a K6-free seed generated K6 at degree six")


def test_pure_trM3_generator_carries_no_k6():
    """Tr(tau^3) on the free seed is pure J6: the falsification gate."""
    for certificate in _certificates():
        prime = certificate["prime"]
        row = next(
            target for target in certificate["targets"]
            if target["id"] == "tr_tau3|d=6|c="
        )
        i2 = row["basis"].index("I6_2")
        i1 = row["basis"].index("I6_1")
        assert row["coordinates"][i2] % prime == 0, "Tr(tau^3) leaked into K6"
        assert row["coordinates"][i1] % prime != 0, (
            "Tr(tau^3) vanished at degree six, which would make the "
            "no-creation statement vacuous rather than substantive")
