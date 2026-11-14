"""The Q10 package must not drift into overclaiming.

Three phrases would each turn a careful result into a false one, and all three
are easy to write by accident when summarising. This makes them fail loudly.
"""
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results" / "intrinsic_candidates"

REQUIRED = ("Preferred ambiguity-minimal Level-B basis among the twelve "
            "published degree-10 candidates under the documented "
            "deterministic simplicity rule.")

ARTIFACTS = [
    "intrinsic_Q10_levelB_basis.json", "Q10_basis_search.json",
    "Q10_canonicality_manifest.json", "Q10_levelA_levelB_map.json",
    "published_degree10_formula_status.json",
    "degree10_reverse_benchmark.json",
]


@pytest.mark.parametrize("name", ARTIFACTS)
def test_artifact_carries_the_frozen_basis_wording(name):
    path = R / name
    if not path.exists():
        pytest.skip(f"{name} not generated")
    payload = json.loads(path.read_text())
    assert payload.get("basis_wording") == REQUIRED, (
        f"{name} lost the frozen wording; the basis is ambiguity-MINIMAL, not "
        f"robust, because no ambiguity-robust triple exists at all")


def test_the_basis_is_never_called_robust_or_canonical():
    """`ambiguity-robust` is the dangerous one: it is provably false here.

    P10_10 is forced and only P10_12 of the remainder is robust, so two robust
    members cannot be found. Every independent triple contains a member whose
    quotient image moves when the source reading is resolved.
    """
    banned = ("ambiguity-robust basis", "universally canonical",
              "the unique compact basis")
    offenders = []
    for path in list((ROOT / "docs").glob("*.md")) + list(R.glob("*.json")):
        text = path.read_text().lower()
        for phrase in banned:
            idx = text.find(phrase.lower())
            while idx != -1:
                # A negated or explicitly-disclaimed mention is the CORRECT
                # statement and must stay allowed: "**not** universally
                # canonical", '"not_claimed": ["universally canonical"]'.
                # Punctuation is normalised first, or markdown bold and JSON
                # key underscores hide the negation from a literal search.
                raw = text[max(0, idx - 90):idx]
                window = " " + re.sub(r"[^a-z0-9]+", " ", raw) + " "
                if not any(f" {neg} " in window for neg in
                           ("no", "not", "never", "cannot", "forbidden",
                            "banned", "wrong", "not_claimed".replace("_", " "),
                            "avoid")):
                    offenders.append(f"{path.name}: ...{text[idx-60:idx+50]}...")
                idx = text.find(phrase.lower(), idx + 1)
    assert not offenders, "overclaiming wording found:\n" + "\n".join(offenders)


def test_ambiguity_facts_are_recorded_consistently():
    path = R / "intrinsic_Q10_levelB_basis.json"
    if not path.exists():
        pytest.skip("basis artifact not generated")
    facts = json.loads(path.read_text())["ambiguity_facts"]
    assert facts["no_fully_robust_published_triple"] is True
    assert "forced" in facts["P10_10"]
    for name in ("P10_09", "P10_11"):
        assert "source-reading dependent" in facts[name], (
            f"{name} must be recorded as source-reading dependent")
    assert "identical" in facts["P10_12"]
