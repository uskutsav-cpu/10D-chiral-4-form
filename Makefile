.PHONY: test baseline demo lint

test:
	PYTHONPATH=src python -m pytest -q

baseline:
	PYTHONPATH=src python scripts/validate_baseline.py

demo:
	PYTHONPATH=src python scripts/obstruction_module.py data/baseline/demo_reachable_rows.json

lint:
	python -m compileall -q src scripts tests
