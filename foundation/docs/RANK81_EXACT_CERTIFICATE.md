# Exact Jacobian rank certificate

## Result

At the sample point and prime recorded in `results/rank81/certificate.json`, the
exact analytic Jacobian of the selected spinor candidate family has

    rank = 81

computed over `F_p` with no finite differences, no step size, no rank tolerance
and no row normalisation.

Per-degree block ranks and cumulative ranks:

| degree | block rank | cumulative |
|---:|---:|---:|
| 4 | 1 | 1 |
| 6 | 2 | 3 |
| 8 | 6 | 9 |
| 10 | 12 | 21 |
| 12 | 61 | **81** |

Euler homogeneity `sum_r c_r dI/dc_r = deg * I` holds for **every** evaluated
candidate, exactly. No Jacobian row is zero.

## What this establishes, precisely

Four statements are kept apart. Merging them is the standard way this kind of
result gets overstated.

**1. Float64 evidence.** Not used. The finite-difference route is abandoned
entirely, not merely supplemented.

**2. Exact modular rank.** `rank_{F_p}(J) = 81` at the recorded point and prime.
This is a computation, not an estimate: over `F_p` an entry is zero or it is not.

**3. Characteristic-zero lower bound — unconditional.** The coordinate basis is
integral: `integral_gamma_traceless_basis()` returns entries in `{-1,0,+1}`
annihilated exactly over `Z`, and the sample point is an integer combination of
it. Every Jacobian is therefore the reduction of a genuine **integer** matrix,
and reduction can only drop rank:

    rank_{F_p}(J mod p)  <=  rank_Q(J)

So `rank_Q(J) >= 81`. This needs no probabilistic caveat. It is not a
Schwartz-Zippel argument and does not degrade with the number of primes tested.

**4. Generic rank.** *Not* established by this computation. The matching upper
bound `126 - dim so(10) = 126 - 45 = 81` is analytic and is a literature result.
Together the two pin the generic functional dimension at exactly 81 — but only
the lower half is ours, and it is a bound at specific points, not a statement
about a generic point.

## Permitted and prohibited wording

Permitted:

> The selected spinor candidate family has exact Jacobian rank 81 at the tested
> sample points, giving an unconditional characteristic-zero lower bound of 81.

Prohibited, and blocked by an automated gate:

> We discovered that there are 81 invariants.
> We proved the generic rank computationally.
> There are 83 algebraically independent invariants.

The last is worth stating explicitly: 83 candidates with rank 81 means there are
**functional dependencies** among the selected functions. The candidates are not
independent and the certificate does not claim they are.

## Scope: 83 of 83

The schedule is complete. Every one of the 83 candidates carries terminal status
`evaluated` with a value, a derivative row, an output hash, a runtime and a
peak-RSS figure; `evaluation_errors` is zero and no Jacobian row is zero.

An earlier revision of this document recorded `c046_portgraph_d12` as
`evaluation_error` after it exceeded the wall-clock budget. It was subsequently
evaluated under the dense-`I` contraction plan and the raised budget recorded in
`flop_limit`. The rank is unchanged, and could not have changed: the observed
rank already equalled the analytic upper bound, and adding a row to a matrix of
maximal rank cannot raise it. What changed is the *schedule* claim, which is now
complete rather than 82 of 83.

The family breakdown at the certified point is 70 port graphs, 9 tensor words
and 4 structured degree-8 candidates, all `evaluated`.

## Contraction plans

Two plans exist for the port-graph derivative and both are exact:

- **factorised** — the invariant tensor is split as
  `I = sum_mu sigma_mu (x) sigma_flip(mu)`, giving two small operands per node;
- **dense-I** — the node is kept whole as a 65536-entry tensor.

The factorised plan is usually far cheaper, but for a few degree-12 topologies
the doubled operand count forces a worse contraction order, and the dense plan
clears the budget where the factorised one does not. The two are required to
produce identical rows and a test asserts it.

Contraction order is chosen greedily, never by `opt_einsum`'s `auto` strategy:
`auto` switches to exact dynamic programming around twenty operands, and a
degree-10 amputated contraction has nineteen. The path *search* then runs away
before any array is allocated. That is what silently stalled two earlier runs.

## Reproduction

```
python spinor_trace_bridge/scripts/run_rank81_certificate.py \
    --archive PATH_TO_YOUR_COPY --fitting-primes 32749 --seeds 11
```

The run checkpoints every candidate row keyed by prime, seed, candidate id and
formula hash, so an interrupted run resumes rather than restarting. Including the
formula hash means a cached row is invalidated automatically if the candidate's
definition changes.
