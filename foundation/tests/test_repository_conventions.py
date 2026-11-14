"""Guards on repository-level conventions that break the documented workflow.

These are cheap static checks. They exist because the failure they cover is
silent at authoring time and only appears when someone runs the documented
command in a clean checkout.
"""
import configparser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _pytest_config():
    parser = configparser.ConfigParser()
    parser.read(ROOT / "pytest.ini")
    return parser["pytest"]


def test_collection_is_scoped_so_bare_pytest_works():
    """`python -m pytest` must not abort on a module-name collision.

    `scripts/` holds helper modules whose names start with `test_`, most
    importantly `scripts/test_M_only_quotients.py`, which is imported by name
    from several runners. Without scoping, bare `pytest` collects it next to
    `tests/test_M_only_quotients.py` and aborts the ENTIRE run during
    collection -- not one failing test, zero tests.
    """
    config = _pytest_config()
    assert config.get("testpaths", "").split() == ["tests"], (
        "pytest.ini must set testpaths = tests, or bare pytest collides on "
        "the duplicated test_M_only_quotients module name")
    assert "scripts" in config.get("norecursedirs", ""), (
        "norecursedirs must exclude scripts/, which testpaths alone does not "
        "cover for invocations like `pytest .`")


def test_no_two_collected_test_modules_share_a_basename():
    """The collision guarded above must not reappear INSIDE tests/.

    pytest imports test modules by basename when there is no __init__.py, so
    two files with the same name anywhere under tests/ would collide exactly
    the same way.
    """
    names = {}
    for path in (ROOT / "tests").rglob("test_*.py"):
        names.setdefault(path.name, []).append(path)
    clashes = {n: p for n, p in names.items() if len(p) > 1}
    assert not clashes, f"duplicate test module basenames under tests/: {clashes}"


def test_helper_scripts_named_like_tests_are_documented():
    """A `test_*.py` under scripts/ is a trap; require it to say so.

    If someone adds another one without the explanatory docstring, the next
    person to hit the collision has no way to know it is a helper rather than
    a test that mysteriously fails to run.
    """
    for path in (ROOT / "scripts").glob("test_*.py"):
        head = path.read_text()[:1200]
        assert "not a pytest" in head.lower() or "helper" in head.lower(), (
            f"{path.name} looks like a test but lives in scripts/. Add a "
            f"docstring saying it is a helper module, or rename it.")


def test_published_span_is_not_a_primitive_complement():
    """Pins the computed product decomposition against a tempting false claim.

    The published span B10 and the product subspace P10 both sit in the
    14-dimensional degree-10 atlas. Two decompositions both read 12+2, which
    invites identifying them. They are NOT the same: B10 meets P10 in exactly
    one dimension and B10+P10 is only 13-dimensional, so B10 has primitive
    content 11, not 12.
    """
    import json
    path = (ROOT / "results" / "intrinsic_candidates"
            / "degree10_published_product_intersection.json")
    if not path.exists():
        return
    payload = json.loads(path.read_text())
    for prime, rec in payload["per_prime"].items():
        assert rec["dim_B10_published"] == 12, prime
        assert rec["dim_P10_product"] == 2, prime
        assert rec["dim_B10_cap_P10"] == 1, (
            f"at {prime} the published span meets the product subspace in "
            f"{rec['dim_B10_cap_P10']} dimensions, not 1; the primitive "
            f"content claim (11) depends on this")
        assert rec["dim_B10_plus_P10"] == 13, prime
        assert rec["B10_is_complement_of_P10"] is False, (
            "B10 must NOT be recorded as a primitive complement")
