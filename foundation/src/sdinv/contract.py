"""Evaluate a contraction graph, and build its Jacobian row.

METRIC. sum_mu F^mu G_mu = sum_mu s_mu F[mu] G[mu], so each edge needs
exactly ONE factor of the metric sign. Every edge is oriented and the sign
is applied at the tail. Applying it at both ends gives s^2 = 1 and silently
computes the Euclidean contraction -- a bug that will not announce itself.

JACOBIAN. I is a degree-n contraction of n copies of F, so it is
multilinear in the vertices and

    dI/dA^k = sum_v [ same graph, F at vertex v replaced by P(e_k) ].

The original implementation built each amputated tensor separately, costing
n tensor-network contractions per graph. The production path below contracts
the network once and reverse-differentiates that contraction tree. This gives
all n amputated tensors while reusing the forward intermediates. The older
amputated implementation remains as a small-case correctness oracle.

For 10D order 8, that reuse is the difference between an exhaustive run being
practical and taking many hours.
"""

import string
import itertools
from functools import lru_cache
import numpy as np

from .modp import P, RankSieve, mod_einsum
from .forms import basis_tuples, metric_signs, perm_sign, to_dense

LETTERS = string.ascii_letters


def _apply_sign(T, axis, s, mod):
    shape = [1] * T.ndim
    shape[axis] = len(s)
    return (T * s.reshape(shape)) % mod


def _slot_plan(M, valence):
    """Assign a letter to each edge-slot; record which axes are tails."""
    n = M.shape[0]
    slots = [[] for _ in range(n)]
    tails = [[] for _ in range(n)]
    li = 0
    for i in range(n):
        for j in range(i + 1, n):
            for _ in range(int(M[i, j])):
                L = LETTERS[li]
                li += 1
                tails[i].append(len(slots[i]))
                slots[i].append(L)
                slots[j].append(L)
    for v in range(n):
        assert len(slots[v]) == valence, "valence mismatch"
    return slots, tails


def _signed(T, axes, s, mod):
    for ax in axes:
        T = _apply_sign(T, ax, s, mod)
    return T


def evaluate(M, tensors, d, valence, lorentzian=True, mod=P):
    slots, tails = _slot_plan(M, valence)
    s = metric_signs(d, lorentzian) % mod
    ops = [_signed(T, tails[v], s, mod) for v, T in enumerate(tensors)]
    sub = ",".join("".join(sl) for sl in slots) + "->"
    return int(mod_einsum(sub, ops, mod))


def amputated(M, F_dense, v, d, valence, lorentzian=True, mod=P):
    """Contract every vertex except v, leaving v's slots free."""
    slots, tails = _slot_plan(M, valence)
    s = metric_signs(d, lorentzian) % mod
    ops, subs = [], []
    for u in range(M.shape[0]):
        if u == v:
            continue
        ops.append(_signed(F_dense, tails[u], s, mod))
        subs.append("".join(slots[u]))
    sub = ",".join(subs) + "->" + "".join(slots[v])
    A = mod_einsum(sub, ops, mod)
    return _signed(A, tails[v], s, mod)


def jacobian_row_amputated(M, F_dense, basis_flat, d, valence,
                           lorentzian=True, mod=P):
    """Independent reference Jacobian using separately amputated vertices."""
    n = M.shape[0]
    gradient = np.zeros_like(F_dense, dtype=np.int64)
    for v in range(n):
        A = amputated(M, F_dense, v, d, valence, lorentzian, mod).ravel()
        gradient = (gradient + A.reshape(F_dense.shape)) % mod
    return _project_gradient(basis_flat, gradient, mod)


@lru_cache(maxsize=256)
def _optimal_contraction_tree_cached(terms, operand_shapes):
    """Globally optimal binary tree for this small closed tensor network.

    A local greedy choice can create a disastrous later intermediate. There
    are at most twelve operands in the completed runs here. Boundary label sets
    are bitsets and proper subset scores are computed in numeric mask order.
    The score first minimises largest pairwise work, then total
    forward/reverse work.
    """
    n = len(terms)
    label_dims = {}
    label_order = list(dict.fromkeys("".join(terms)))
    label_bits = {label: 1 << k for k, label in enumerate(label_order)}
    label_sets = []
    for term, shape in zip(terms, operand_shapes):
        labels = sum(label_bits[label] for label in term)
        label_sets.append(labels)
        for label, size in zip(term, shape):
            old = label_dims.setdefault(label, size)
            assert old == size, f"dimension mismatch for index {label}"

    dimensions = [label_dims[label] for label in label_order]
    uniform_dimension = (
        dimensions[0]
        if dimensions and len(set(dimensions)) == 1
        else None
    )
    if uniform_dimension is not None:
        powers = [
            uniform_dimension ** k for k in range(len(dimensions) + 1)]
        work_for = lambda bits: powers[bits.bit_count()]
    else:
        work_cache = {0: 1}

        def work_for(bits):
            cached = work_cache.get(bits)
            if cached is not None:
                return cached
            work = 1
            remaining = bits
            while remaining:
                bit = remaining & -remaining
                work *= dimensions[bit.bit_length() - 1]
                remaining ^= bit
            work_cache[bits] = work
            return work

    boundary = [0 for _ in range(1 << n)]
    for mask in range(1, 1 << n):
        bit = mask & -mask
        v = bit.bit_length() - 1
        boundary[mask] = boundary[mask ^ bit] ^ label_sets[v]

    # score[mask] = (largest pair work, total fwd+reverse work).
    peak_score = [0 for _ in range(1 << n)]
    total_score = [0 for _ in range(1 << n)]
    split = [None for _ in range(1 << n)]
    for mask in range(1, 1 << n):
        if mask & (mask - 1) == 0:
            continue
        anchor = mask & -mask
        part = (mask - 1) & mask
        best_peak = best_total = None
        best_split = None
        while part:
            if part & anchor:
                other = mask ^ part
                if other:
                    pair_work = work_for(
                        boundary[part] | boundary[other])
                    candidate_peak = max(
                        peak_score[part], peak_score[other], pair_work)
                    candidate_total = (
                        total_score[part] + total_score[other]
                        + 3 * pair_work
                    )
                    if (
                        best_peak is None
                        or (candidate_peak, candidate_total)
                        < (best_peak, best_total)
                    ):
                        best_peak = candidate_peak
                        best_total = candidate_total
                        best_split = (part, other)
            part = (part - 1) & mask
        peak_score[mask] = best_peak
        total_score[mask] = best_total
        split[mask] = best_split

    full = (1 << n) - 1
    return split, boundary, (peak_score[full], total_score[full])


def _optimal_contraction_tree(terms, operand_shapes):
    return _optimal_contraction_tree_cached(
        tuple(terms), tuple(tuple(shape) for shape in operand_shapes))


def contraction_plan_cost(M, d, valence):
    """Return (largest_pair_work, total_forward_reverse_work) without running."""
    slots, _ = _slot_plan(M, valence)
    terms = ["".join(x) for x in slots]
    _, _, score = _optimal_contraction_tree(
        terms, [(d,) * valence for _ in terms])
    return score


def _contraction_plan_profile(terms, operand_shapes, itemsize=8):
    """Return an allocation-oriented profile of the exact contraction tree.

    ``largest_pair_work`` counts arithmetic terms and is intentionally
    separate from allocation size. The byte estimate includes every retained
    forward node, every reverse adjoint, and a conservative fourfold workspace
    for the largest pair. It is a pre-allocation guard, not an RSS claim:
    Python/NumPy allocator overhead and the input form live outside it.
    """
    split, boundary, score = _optimal_contraction_tree(terms, operand_shapes)
    n = len(terms)
    label_dims = {}
    label_order = list(dict.fromkeys("".join(terms)))
    for term, shape in zip(terms, operand_shapes):
        for label, size in zip(term, shape):
            old = label_dims.setdefault(label, int(size))
            if old != int(size):
                raise ValueError(f"dimension mismatch for index {label}")

    def elements(labels):
        count = 1
        remaining = labels
        while remaining:
            bit = remaining & -remaining
            count *= label_dims[label_order[bit.bit_length() - 1]]
            remaining ^= bit
        return count

    node_elements = []
    max_output_rank = 0
    max_pair_union_rank = 0
    max_pair_elements = 0

    def visit(mask):
        nonlocal max_output_rank, max_pair_union_rank, max_pair_elements
        if mask.bit_count() == 1:
            index = (mask & -mask).bit_length() - 1
            size = int(np.prod(operand_shapes[index], dtype=object))
            node_elements.append(size)
            max_output_rank = max(
                max_output_rank, len(operand_shapes[index]))
            return
        left, right = split[mask]
        visit(left)
        visit(right)
        union = boundary[left] | boundary[right]
        output = boundary[mask]
        left_elements = elements(boundary[left])
        right_elements = elements(boundary[right])
        output_elements = elements(output)
        node_elements.append(output_elements)
        max_output_rank = max(max_output_rank, output.bit_count())
        max_pair_union_rank = max(
            max_pair_union_rank, union.bit_count())
        max_pair_elements = max(
            max_pair_elements,
            left_elements + right_elements + output_elements,
        )

    visit((1 << n) - 1)
    retained_elements = sum(node_elements)
    retained_forward_reverse_bytes = (
        2 * retained_elements * int(itemsize))
    workspace_bytes = 4 * max_pair_elements * int(itemsize)
    return {
        "largest_pair_work": int(score[0]),
        "total_forward_reverse_work": int(score[1]),
        "max_output_rank": int(max_output_rank),
        "max_output_elements": int(max(node_elements, default=1)),
        "max_output_bytes": int(
            max(node_elements, default=1) * int(itemsize)),
        "max_pair_union_rank": int(max_pair_union_rank),
        "retained_node_elements": int(retained_elements),
        "retained_forward_reverse_bytes": int(
            retained_forward_reverse_bytes),
        "workspace_bytes": int(workspace_bytes),
        "estimated_peak_bytes": int(
            retained_forward_reverse_bytes + workspace_bytes),
    }


def contraction_plan_profile(M, d, valence, itemsize=8):
    """Profile the globally optimal plan without allocating its tensors."""
    slots, _ = _slot_plan(M, valence)
    terms = ["".join(x) for x in slots]
    return _contraction_plan_profile(
        terms, [(d,) * valence for _ in terms], itemsize)


def greedy_contraction_plan_cost(M, d, valence):
    """Cheap deterministic upper bound used only to schedule candidates.

    Production evaluation still uses the globally optimal dynamic program
    above. This O(n^3) heuristic is inexpensive enough to enrich an entire
    large catalog and puts narrow tensor networks first.
    """
    slots, _ = _slot_plan(M, valence)
    label_bits = {
        label: 1 << k
        for k, label in enumerate(dict.fromkeys("".join(
            "".join(term) for term in slots)))
    }
    boundaries = [
        sum(label_bits[label] for label in term)
        for term in slots
    ]
    peak = total = 0
    while len(boundaries) > 1:
        best = None
        for i in range(len(boundaries)):
            for j in range(i + 1, len(boundaries)):
                union = boundaries[i] | boundaries[j]
                result = boundaries[i] ^ boundaries[j]
                work = d ** union.bit_count()
                output = d ** result.bit_count()
                candidate = (work, output, i, j, result)
                if best is None or candidate[:4] < best[:4]:
                    best = candidate
        work, _, i, j, result = best
        peak = max(peak, work)
        total += 3 * work
        for k in sorted((i, j), reverse=True):
            boundaries.pop(k)
        boundaries.append(result)
    return peak, total


def _network_value_gradient(terms, operands, mod=P, plan_profile=None,
                            max_memory_bytes=None):
    """Reverse-differentiate a scalar pairwise tensor contraction.

    Each index label must occur exactly twice across the input terms, as it
    does for a fully contracted graph. Return the scalar value plus the
    derivative with respect to each operand in its original axis order.
    """
    if plan_profile is None:
        plan_profile = _contraction_plan_profile(
            terms, [np.shape(op) for op in operands])
    estimated = int(plan_profile["estimated_peak_bytes"])
    if max_memory_bytes is not None and estimated > int(max_memory_bytes):
        raise MemoryError(
            "exact contraction plan rejected before tensor allocation: "
            f"estimated {estimated} bytes exceeds limit "
            f"{int(max_memory_bytes)} bytes")

    nodes = [{"term": term, "value": np.asarray(op, dtype=np.int64) % mod,
              "left": None, "right": None}
             for term, op in zip(terms, operands)]
    split, boundary, _ = _optimal_contraction_tree(
        terms, [op.shape for op in operands])
    label_bits = {
        label: 1 << k
        for k, label in enumerate(dict.fromkeys("".join(terms)))
    }

    def forward(mask):
        if mask.bit_count() == 1:
            return (mask & -mask).bit_length() - 1
        left_mask, right_mask = split[mask]
        left, right = forward(left_mask), forward(right_mask)
        lt, rt = nodes[left]["term"], nodes[right]["term"]
        keep_labels = boundary[mask]
        keep = "".join(c for c in dict.fromkeys(lt + rt)
                       if keep_labels & label_bits[c])
        result = mod_einsum(
            f"{lt},{rt}->{keep}",
            [nodes[left]["value"], nodes[right]["value"]],
            mod,
        )
        parent = len(nodes)
        nodes.append({"term": keep, "value": result,
                      "left": left, "right": right})
        return parent

    root = forward((1 << len(operands)) - 1)
    assert nodes[root]["term"] == "", "network is not fully contracted"
    scalar = int(np.asarray(nodes[root]["value"]).item()) % mod
    adjoints = {root: np.array(1, dtype=np.int64)}

    for parent in range(root, len(operands) - 1, -1):
        node = nodes[parent]
        left, right = node["left"], node["right"]
        pt = node["term"]
        lt, rt = nodes[left]["term"], nodes[right]["term"]
        adj = adjoints[parent]
        adjoints[left] = mod_einsum(
            f"{pt},{rt}->{lt}", [adj, nodes[right]["value"]], mod)
        adjoints[right] = mod_einsum(
            f"{pt},{lt}->{rt}", [adj, nodes[left]["value"]], mod)

    return scalar, [adjoints[k] for k in range(len(operands))]


def _network_value(terms, operands, mod=P, max_memory_bytes=None):
    """Evaluate a closed network on the globally optimal binary tree."""
    profile = _contraction_plan_profile(
        terms, [np.shape(op) for op in operands])
    # The full forward/reverse estimate is conservative for this forward-only
    # path and therefore safe to reuse as its allocation guard.
    if (
        max_memory_bytes is not None
        and profile["estimated_peak_bytes"] > int(max_memory_bytes)
    ):
        raise MemoryError(
            "exact scalar contraction plan rejected before tensor allocation: "
            f"estimated {profile['estimated_peak_bytes']} bytes exceeds limit "
            f"{int(max_memory_bytes)} bytes")
    split, boundary, _ = _optimal_contraction_tree(
        terms, [np.shape(op) for op in operands])
    label_bits = {
        label: 1 << k
        for k, label in enumerate(dict.fromkeys("".join(terms)))
    }
    leaves = [
        np.asarray(op, dtype=np.int64) % mod for op in operands
    ]

    def forward(mask):
        if mask.bit_count() == 1:
            index = (mask & -mask).bit_length() - 1
            return terms[index], leaves[index]
        left_mask, right_mask = split[mask]
        left_term, left_value = forward(left_mask)
        right_term, right_value = forward(right_mask)
        keep_labels = boundary[mask]
        keep = "".join(
            label for label in dict.fromkeys(left_term + right_term)
            if keep_labels & label_bits[label]
        )
        result = mod_einsum(
            f"{left_term},{right_term}->{keep}",
            [left_value, right_value],
            mod,
        )
        return keep, result

    term, result = forward((1 << len(operands)) - 1)
    if term:
        raise ValueError("network is not fully contracted")
    return int(np.asarray(result).item()) % mod


def _network_gradient(terms, operands, mod=P):
    """Compatibility wrapper returning only reverse-mode gradients."""
    return _network_value_gradient(terms, operands, mod)[1]


def value_and_jacobian_row(M, F_dense, basis_flat, d, valence,
                           lorentzian=True, mod=P, backend="optimized",
                           max_memory_bytes=None):
    """Return an exact scalar value and Jacobian row from one graph.

    ``optimized`` shares one contraction tree between the forward scalar and
    reverse derivative. ``reference`` intentionally uses the independent
    amputated implementation and is a small-case correctness oracle.
    """
    if backend == "reference":
        return (
            value(M, F_dense, d, valence, lorentzian, mod),
            jacobian_row_amputated(
                M, F_dense, basis_flat, d, valence, lorentzian, mod),
        )
    if backend != "optimized":
        raise ValueError(f"unknown contraction backend: {backend!r}")

    slots, tails = _slot_plan(M, valence)
    s = metric_signs(d, lorentzian) % mod
    terms = ["".join(x) for x in slots]
    plan_profile = _contraction_plan_profile(
        terms, [(d,) * valence for _ in terms])
    if (
        max_memory_bytes is not None
        and plan_profile["estimated_peak_bytes"] > int(max_memory_bytes)
    ):
        raise MemoryError(
            "exact contraction plan rejected before signed operands: "
            f"estimated {plan_profile['estimated_peak_bytes']} bytes exceeds "
            f"limit {int(max_memory_bytes)} bytes")
    operands = [_signed(F_dense, tails[v], s, mod)
                for v in range(M.shape[0])]
    scalar, operand_gradients = _network_value_gradient(
        terms,
        operands,
        mod,
        plan_profile=plan_profile,
        max_memory_bytes=max_memory_bytes,
    )

    gradient = np.zeros_like(F_dense, dtype=np.int64)
    for v, grad in enumerate(operand_gradients):
        # The forward operand is signed(F), so the chain rule applies the
        # same diagonal metric signs once more to map dI/d(signed F) to dI/dF.
        gradient = (gradient + _signed(grad, tails[v], s, mod)) % mod
    return scalar, _project_gradient(basis_flat, gradient, mod)


def jacobian_row(M, F_dense, basis_flat, d, valence,
                 lorentzian=True, mod=P, backend="optimized"):
    """Jacobian row from the selected exact contraction backend."""
    return value_and_jacobian_row(
        M,
        F_dense,
        basis_flat,
        d,
        valence,
        lorentzian,
        mod,
        backend,
    )[1]


def value(M, F_dense, d, valence, lorentzian=True, mod=P):
    return evaluate(M, [F_dense] * M.shape[0], d, valence, lorentzian, mod)


def planned_value(M, F_dense, d, valence, lorentzian=True, mod=P,
                  max_memory_bytes=None):
    """Scalar value using the same globally optimal plan as the Jacobian."""
    slots, tails = _slot_plan(M, valence)
    s = metric_signs(d, lorentzian) % mod
    terms = ["".join(x) for x in slots]
    profile = _contraction_plan_profile(
        terms, [(d,) * valence for _ in terms])
    if (
        max_memory_bytes is not None
        and profile["estimated_peak_bytes"] > int(max_memory_bytes)
    ):
        raise MemoryError(
            "exact scalar contraction plan rejected before signed operands: "
            f"estimated {profile['estimated_peak_bytes']} bytes exceeds limit "
            f"{int(max_memory_bytes)} bytes")
    operands = [
        _signed(F_dense, tails[v], s, mod) for v in range(M.shape[0])
    ]
    return _network_value(
        terms, operands, mod, max_memory_bytes=max_memory_bytes)


def build_basis_flat(d, p_deg, projector=None, mod=P):
    """Flattened dense basis directions, stacked as rows."""
    tuples = basis_tuples(d, p_deg)
    N = len(tuples)
    rows = []
    for k in range(N):
        vec = np.zeros(N, dtype=np.int64)
        vec[k] = 1
        if projector is not None:
            vec = (projector @ vec) % mod
        rows.append(to_dense(vec, d, p_deg, mod).ravel())
    return np.array(rows, dtype=np.int64)


class CompactDerivativeBasis:
    """Project dense antisymmetric gradients without dense basis tensors.

    A dense p-form repeats every sorted component p! times with matching
    permutation signs. Pulling a dense gradient back to antisymmetric
    components therefore reduces exactly to signed orbit sums followed by a
    small compact matrix product.
    For the 10D self-dual 5-form this replaces a 252 x 100,000 array with a
    252 x 126 matrix.
    """

    def __init__(self, d, p_deg, directions, coordinate_indices, mod=P):
        self.d = int(d)
        self.p_deg = int(p_deg)
        self.mod = int(mod)
        self.tuples = basis_tuples(self.d, self.p_deg)
        self.directions = np.asarray(directions, dtype=np.int64) % self.mod
        self.coordinate_indices = [int(x) for x in coordinate_indices]
        expected = (len(self.tuples), len(self.coordinate_indices))
        if self.directions.shape != expected:
            raise ValueError(
                f"direction matrix has shape {self.directions.shape}, "
                f"expected {expected}")
        self.ncols = len(self.coordinate_indices)
        permutations = list(itertools.permutations(range(self.p_deg)))
        dense_indices = []
        signs = []
        for index in self.tuples:
            for permutation in permutations:
                dense_indices.append(tuple(index[k] for k in permutation))
                signs.append(perm_sign(permutation))
        self._flat_indices = np.ravel_multi_index(
            np.asarray(dense_indices, dtype=np.int64).T,
            (self.d,) * self.p_deg,
        )
        self._signs = np.asarray(signs, dtype=np.int64)

    def project(self, gradient):
        signed = gradient.ravel()[self._flat_indices] * self._signs
        compact = signed.reshape(len(self.tuples), -1).sum(axis=1) % self.mod
        return (self.directions.T @ compact) % self.mod


def _project_gradient(basis, gradient, mod):
    if isinstance(basis, CompactDerivativeBasis):
        if basis.mod != mod:
            raise ValueError("compact derivative basis uses a different prime")
        return basis.project(gradient)
    return (np.asarray(basis, dtype=np.int64) @ gradient.ravel()) % mod


def build_compact_derivative_basis(d, p_deg, projector=None, mod=P,
                                   independent=True):
    """Build a compact derivative basis, optionally removing redundancies."""
    component_count = len(basis_tuples(d, p_deg))
    if projector is None:
        directions = np.eye(component_count, dtype=np.int64)
    else:
        directions = np.asarray(projector, dtype=np.int64) % mod
        if directions.shape != (component_count, component_count):
            raise ValueError("projector has the wrong shape")

    if independent:
        sieve = RankSieve(component_count, mod)
        indices = [
            k for k in range(directions.shape[1])
            if sieve.add(directions[:, k])
        ]
        directions = directions[:, indices]
    else:
        indices = list(range(directions.shape[1]))
    return CompactDerivativeBasis(
        d, p_deg, directions, indices, mod)
