# Reproducing the degree-10 result

## 1. Quick path

```bash
scripts/reproduce_Q10_levelB.sh QUICK
```

Verifies committed artifacts and recomputes what is cheap and exact. Does not
re-evaluate the atlas; `FULL` does that.

## 2. Environment

| item | value |
|---|---|
| OS | macOS (Darwin 25.5.0), arm64 |
| Python | 3.13 (repo `.venv`) |
| hard dependencies | `numpy`, `pytest` |
| optional | `pypdf` (primary-source text extraction only) |
| memory | 8 GB; **one heavy job at a time** |

`ru_maxrss` is **bytes on Darwin, KiB on Linux**. `peak_rss_mb` decides by
`sys.platform`, never by magnitude — at ~1 GB the two units are
indistinguishable by value, and the magnitude heuristic reported 958576 MB.

## 3. Checkpoints must not live in iCloud

The canonical tree is under `~/Documents`, which is synced. Every runner
defaults its checkpoint root to local temp and honours `SDINV_CKPT_ROOT`:

```bash
export SDINV_CKPT_ROOT="${TMPDIR}/sdinv_ckpt"
```

Checkpoint validity is **semantic**, not commit-based: fingerprints are derived
from source, so a documentation commit preserves the atlas cache while a
reimplemented formula invalidates exactly its own units.

## 4. Full test suite

```bash
python -m pytest
```

Bare invocation works: `pytest.ini` scopes collection to `tests/`, past the
`scripts/test_M_only_quotients.py` ↔ `tests/test_M_only_quotients.py`
module-name collision that otherwise aborts the entire run during collection.

## 5. Regenerating the artifacts

```bash
# published candidates -> Q10   (checkpointed; cold ~15 min/prime, resumed ~1 min)
SDINV_PRIMES=32749,32717 python scripts/project_published_degree10_ckpt.py

# positive control: the projector must be able to return rank 3
python scripts/positive_control_degree10_quotient.py

# published Level-B basis selection and Level-A map
python scripts/build_Q10_levelB_basis.py
python scripts/emit_Q10_levelA_levelB_map.py
python scripts/emit_formula_status.py

# independent reverse search   (~17 min at --limit 15)
python scripts/reverse_engineer_degree10_benchmark.py --limit 15 \
       --topology-cap 30000 --max-rss-mb 2000

# pure-N4125 basis, maps, comparison
python scripts/validate_reverse_degree10_span.py
python scripts/emit_pure_N4125_basis.py
```

## 6. Measured costs

| operation | cost |
|---|---:|
| quadratic-block build, per sample | 1.6–4.4 s |
| one reverse contraction, blocks cached | ~119 ms |
| published evaluator (rebuilds blocks) | ~890 ms |
| reverse pilot, 259 candidates × 22 samples | 1025 s |
| published projection, 2 primes cold | 912 s |
| ...resumed after adding one evaluator | 55 s |
| full test suite | ~13 min |

## 7. Memory hazards, both real and measured

- **Evaluation, not enumeration, is now the memory ceiling.** Streaming fixed
  enumeration (2611 MB → 33 MB), but particular contraction paths still build
  large intermediates and the pilot peaks ~2.7 GB.
- A first corrected rerun **degraded from ~105 s to 1335 s per sample** as swap
  filled (2.0 GB of 3.0 GB). It was restarted with `--limit 15` and completed
  in 1025 s. Watch for per-sample time climbing; that is the signal.
- Never run two heavy jobs at once. Three pytest runs have been silently killed
  on this machine at ~60 MB free.

## 8. A trap worth knowing before you debug a "hung" job

A process whose stdout is redirected is **block buffered**. An empty log does
not mean an idle job. Check `ps -p <pid> -o time` for CPU accumulation and the
output artifact's mtime. A run killed *after* it wrote its artifact loses only
the buffered log — that happened here and briefly looked like lost work.

Also: `ps -Ao ... | grep <pattern>` can match the grep's own shell, and can
miss the target entirely. Check by **PID**.
