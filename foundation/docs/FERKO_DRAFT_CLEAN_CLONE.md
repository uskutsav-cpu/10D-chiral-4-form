# Clean-clone verification of the Ferko mentor-review draft

A fresh clone, outside iCloud, at the mentor-draft commit, with dependencies
installed from `requirements.txt` only. Every command below was run in that
clone and every exit code was recorded. Nothing here is carried over from the
working tree.

## Environment

| item | value |
|---|---|
| commit | `754d9c4bdb64187327d6a7cb7a5a07c17a1576aa` |
| branch | `publication/jhep-mentor-draft` |
| clone path | `/private/tmp/.../cleanclone/repo` (non-iCloud) |
| Python | 3.13.13 |
| NumPy | 2.5.1 |
| opt_einsum | 3.4.0 |
| pynauty | 2.8.8.1 |
| pytest | 9.1.1 |
| TeX engine | Tectonic (XeTeX), TeX Live 2026 bundle |
| host | macOS 25.5.0, arm64, 8 GiB RAM |

## Results

| step | command | exit | result |
|---|---|:--:|---|
| tensor suite | `pytest tests` (4 chunks) | 0,0,0,0 | **254 passed** (81+49+62+62), 1414 s |
| bridge suite | `pytest tests` in `spinor_trace_bridge/` | 0 | **86 passed**, 576 s |
| orientation | `test_orientation.py` (within bridge suite) | 0 | 14 passed |
| flow / G-10 | `test_G10_trace_activation.py`, `test_stress_flow*.py` | 0 | within the 254 |
| `D10`/`Q10` | `scripts/d10_characteristic_zero.py` | 0 | `A10=14`, `D10=11`, `Q10=3`; fixed point after 3 sweeps; 11×11 minor nonzero |
| `B10 ∩ P10` | `scripts/b10_p10_characteristic_zero.py` | 0 | `dim B10=12`, `dim P10=2`, `dim(B10+P10)=13`, **`dim(B10 ∩ P10)=1`**; STATUS exact |
| rank matrix | `scripts/aggregate_rank_matrix.py --expect-cells 15` | 0 | **15/15 cells, rank 81**, 81 stable pivot rows, 81 stable pivot columns, 0 unstable |
| rank-matrix guards | `--self-test` | 0 | 11 rejections fired; cells unmodified; order-independent |
| manuscript gates | `manuscript/jhep/scripts/check_draft.py` | 0 | **52 gates passed, 0 failed** |
| figures | `make_figures.py` | 0 | 10 figures, **byte-identical to committed** (empty `git diff`) |
| tables | `make_tables.py` | 0 | byte-identical except the environment table's commit field (see below) |
| compile | `tectonic -X compile main.tex` | 0 | **41 pages**, 0 errors, 0 overfull, 0 underfull, 0 undefined citations, 0 undefined references |
| package | `build_mentor_package.py` → `build_ferko_delivery.py` | 0 | 11 files |
| package, again | `build_ferko_delivery.py` a second time | 0 | `SHA256SUMS` **identical**; zip **byte-identical** |
| hashes | `shasum -a 256 -c SHA256SUMS` | 0 | **9/9 OK** |

The tensor suite was run as four sequential chunks rather than one invocation
purely to fit a per-command time limit; each chunk is a separate process with
its own recorded exit code, and the four cover all 26 files exactly once. The
same suite also completed as a single process in the working tree: 254 passed,
exit 0, 1411.78 s.

## Two artifacts that do not reproduce byte-for-byte, and why

**The compiled PDF does not.** Tectonic writes a fresh trailer `/ID` and
build-time metadata on every run. Three consecutive compilations of identical
source in the clean clone produced three distinct SHA-256 digests, differing in
about 7,000 bytes once the change propagates into the compressed streams;
`SOURCE_DATE_EPOCH=1735689600` does not suppress it. This is why the working
tree and the clean clone agree on the source archive and disagree on the PDF
hash. The PDF hash in the build manifest identifies a particular build; it does
not certify a reproducible one.

What *is* byte-identical across clean rebuilds: all ten figures, the generated
tables and macros, the LaTeX source archive, and the packaging step (verified by
building the delivery directory twice and diffing `SHA256SUMS`).

**`tab12_environment.tex` records the generating commit,** so regenerating it at
a later commit changes exactly one line. This is self-reference, not drift: the
table names the commit the artifacts were built from, which is necessarily the
parent of the commit that stores the rebuilt table. The clean clone's copy named
`754d9c4bdb64`; the committed copy named `038c95b9a3e6` until it was regenerated
and committed. No other field differs.

## One thing the clean clone found

`opt_einsum` was documented as optional, was absent from `requirements.txt`, and
was described as not installed for the certified runs. It is required. It
supplies both the contraction order and the intermediate-size and flop estimates
that `_modular_contract` enforces as a budget, and both checks sat inside
`if _oe is not None:`. In the documented environment neither budget was ever
evaluated and the code took an untested naive contraction order with no guard at
all. On this 8 GiB host the bridge suite died with `SIGKILL` and no traceback;
`StructuredDegree8.values` goes from 0.4 s at 0.27 GB peak RSS with opt_einsum
to an unbounded allocation without it. It is now a requirement, the fallback
raises instead of running unguarded, and the four documents that called it
optional are corrected. No computed value changed.

This is recorded because it is the second defect in this project found only by a
clean clone — the first being the undocumented Python 3.10 floor — and both were
invisible from a working tree that already had the right environment.
