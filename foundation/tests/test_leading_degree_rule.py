"""PO-02: leading field degrees must match an externally stated rule.

The concern PO-02 records is that `leading_field_degree` is computed by the same
code that validates it, so its internal consistency proves nothing. This test
states the rule independently and checks the certificate against it:

    Tr(tau)    has leading field degree 4, not 2, because the free stress tensor
               is traceless and the degree-2 term vanishes;
    Tr(tau^k)  has leading field degree 2k for k >= 2;
    a product   has the sum of its factors' leading degrees.

The rule is not derived here -- deriving it is a physics statement flagged for
coauthor confirmation. What is established is that the certificate's assignments
are not self-referential: they agree with a rule written down separately, and a
single misassigned target would break the pattern.
"""

import json
import re
from collections import defaultdict
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CERTS = sorted((ROOT / "results/stress_flow/certificates").glob(
    "interacting_degree12_*.json"))


def leading_degree(generator: str) -> int:
    """Predicted leading field degree, from the rule alone."""
    total = 0
    for factor in generator.split("*"):
        m = re.fullmatch(r"tr_tau(\d*)(?:\^(\d+))?", factor)
        assert m, f"unparsed generator factor {factor!r}"
        k = int(m.group(1)) if m.group(1) else 1
        power = int(m.group(2)) if m.group(2) else 1
        total += (4 if k == 1 else 2 * k) * power
    return total


@pytest.mark.parametrize("cert", CERTS, ids=lambda p: p.stem.rsplit("_", 1)[1])
def test_first_appearance_matches_the_rule(cert):
    d = json.loads(cert.read_text())
    first = defaultdict(lambda: 10**9)
    for t in d["targets"]:
        first[t["generator"]] = min(first[t["generator"]], t["field_degree"])
    wrong = {g: (obs, leading_degree(g)) for g, obs in first.items()
             if obs != leading_degree(g)}
    assert not wrong, (
        "generators whose first appearance disagrees with the stated rule "
        f"(observed, predicted): {wrong}")


@pytest.mark.parametrize("cert", CERTS, ids=lambda p: p.stem.rsplit("_", 1)[1])
def test_trace_of_tau_starts_at_four_not_two(cert):
    """The one assignment that is not a bare degree count."""
    d = json.loads(cert.read_text())
    degs = [t["field_degree"] for t in d["targets"] if t["generator"] == "tr_tau"]
    assert degs, "tr_tau absent from the certificate"
    assert min(degs) == 4, (
        f"Tr(tau) first appears at degree {min(degs)}; the traceless free stress "
        f"tensor requires 4. If this is 2 the leading-degree bookkeeping is wrong "
        f"and every closure dimension built on it is suspect")


def test_all_certificates_agree_on_assignments():
    if len(CERTS) < 2:
        pytest.skip("need at least two certificates")
    ref = None
    for c in CERTS:
        d = json.loads(c.read_text())
        got = {t["id"]: t["field_degree"] for t in d["targets"]}
        if ref is None:
            ref = got
        else:
            assert got == ref, f"{c.name} assigns degrees differently"
