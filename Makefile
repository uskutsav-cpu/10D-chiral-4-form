PYTHON ?= python3
export PYTHONPATH := src
export OPENBLAS_NUM_THREADS ?= 1
export OMP_NUM_THREADS ?= 1
export VECLIB_MAXIMUM_THREADS ?= 1

.PHONY: verify test lint baseline demo plan smoke degree8 completion8 map-proof symbolic

verify: lint test baseline demo
	$(PYTHON) -m chiral4form verify-map8 verification/map8/map_certificate.json

lint:
	$(PYTHON) -m compileall -q src scripts tests

test:
	$(PYTHON) -m pytest -q

baseline:
	$(PYTHON) scripts/validate_baseline.py

demo:
	$(PYTHON) scripts/obstruction_module.py data/baseline/demo_reachable_rows.json

plan:
	$(PYTHON) -m chiral4form plan --config configs/degree8.json

smoke:
	$(PYTHON) -m chiral4form run --config configs/smoke.json --output runs/smoke

degree8:
	$(PYTHON) -m chiral4form run --config configs/degree8.json --output runs/degree8

completion8:
	$(PYTHON) -m chiral4form run --config configs/degree8_completion.json --output runs/degree8-generalized
	$(PYTHON) -m chiral4form completion8 --output runs/completion8-exact.json

map-proof:
	$(PYTHON) -m chiral4form prove-map8 --output runs/map8-proof

symbolic:
	$(PYTHON) -m chiral4form sextic --output runs/sextic-exact.json
	$(PYTHON) -m chiral4form degree8 --output runs/degree8-orbit-exact.json
	$(PYTHON) -m chiral4form completion8 --output runs/completion8-exact.json
	$(PYTHON) -m chiral4form modmax --output runs/modmax-known-identity.json
