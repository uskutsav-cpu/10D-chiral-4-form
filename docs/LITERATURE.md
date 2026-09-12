# Primary sources and implementation grounding

Consulted through web/connected source reads during this implementation.
This is a provenance list, not an exhaustive novelty review.

1. Hutomo, Lechner, Sorokin, **On non-linear chiral 4-form theories in D=10**,
   arXiv:2509.14351v2 (6 January 2026), JHEP 02 (2026) 147.
   https://arxiv.org/html/2509.14351v2
   Equation (2.33) supplies the interacting stress tensor. Equation (2.23)
   motivates the constrained derivative. Equations (3.12)-(3.16) supply the
   ModMax stress-square reproduction. The conclusion motivates the subclass
   question; our software does not re-claim the known qualitative failure of
   generic stress-only universality.
2. Cederwall, Hutomo, Kuzenko, Lechner, Sorokin, **Some remarks on invariants**,
   arXiv:2509.14350v2, J. Phys. A 59 (2026) 065203.
   https://arxiv.org/html/2509.14350v2
   Homogeneous invariant-count input through degree eight. Counts are external
   mathematical inputs, not results of this package's sparse polynomial code.
3. **Interacting Chiral Form Field Theories and TTbar-like Flows in Six and
   Higher Dimensions**, arXiv:2402.06947.
   https://arxiv.org/abs/2402.06947
   Lower-dimensional/chiral-flow context; no claim of novelty for the flow
   construction itself.
4. **Duality-Invariant Non-linear Electrodynamics and Stress Tensor Flows**,
   arXiv:2309.04253.
   https://arxiv.org/abs/2309.04253
   Related stress-flow program; a lower-dimensional result must not silently
   be promoted to a ten-dimensional theorem.
5. Predecessor source at commit
   `3ed32805b38ce34216b34888f6539e3538e90fb9`:
   https://github.com/uskutsav-cpu/selfdual-5form-invariants
   Actual graph formulas and declared conventions. The low-degree fixture
   records the original file's Git blob. A full source import additionally
   verifies source SHA-256 hashes and Git blobs automatically.

Software references:
- SymPy exact polynomial/Groebner documentation:
  https://docs.sympy.org/latest/modules/polys/reference.html
- Git bundle format and verification:
  https://git-scm.com/docs/git-bundle
- GitHub Actions Python build/test instructions:
  https://docs.github.com/en/actions/tutorials/build-and-test-code/python

Novelty questions remain: whether the exact degree-eight orbit graph,
nonlinear coupling obstruction, or finite-catalogue completion was already
stated in another formulation. Neither automated tests nor this source list
answers priority. Ask the mentor to review that before writing a novelty claim.
