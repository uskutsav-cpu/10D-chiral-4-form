from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_claim_ledger_records_refuted_five_generator_completion():
    text = (ROOT / "CLAIMS.md").read_text()
    line = next(line for line in text.splitlines() if line.startswith("| C08 "))
    assert "REFUTED" in line


def test_readme_marks_baseline_imported():
    text = (ROOT / "README.md").read_text().lower()
    assert "imported claims" in text
    assert "not independently" in text
