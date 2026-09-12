from pathlib import Path
import pytest
from chiral4form.registry import Registry

ROOT = Path(__file__).resolve().parents[1]

@pytest.fixture(scope="session")
def registry():
    return Registry.load(ROOT / "data/fixtures/low_degree_registry.json")
