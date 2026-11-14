"""Independent reverse search: degree-10 scalars from quadratic-block topology.

Purpose
-------
The published equation-(4.24) candidates span Q10. This module asks the same
question from the other end, WITHOUT reference to those formulas: enumerate
compact contractions of quadratic blocks directly, evaluate them exactly, and
see whether the quotient is reached.

Independence is a hard constraint. Nothing here imports
`published_degree10_invariants`, and no published coordinate vector may be
consulted during generation. The published result is ground truth for
comparison AFTER a span has been recovered, never an input to the search.

Why the enumeration is tractable at all
---------------------------------------
A degree-10 scalar is five quadratic blocks, each degree 2 in F. Blocks are
`M` (2 slots), `N^(1050)` (6 slots) and `N^(4125)` (6 slots). A naive search
over perfect matchings of the raw index slots is hopeless: five `N` blocks give
30 slots and 29!! ~ 6.2e15 matchings.

But the slots are not distinguishable. `composite_n1050` antisymmetrises axes
(0,1,2,3,4), so those five slots are interchangeable up to sign, and only axis
5 is distinct. Each block therefore has two SLOT CLASSES:

    class A -- axes 0..4, multiplicity 5, antisymmetric
    class B -- axis 5,    multiplicity 1

A contraction is then fully described by how many edges join each ordered pair
of (block, class) endpoints. That is a small integer enumeration with degree
constraints, and it is what `enumerate_topologies` walks.

Pruning, in the order applied
-----------------------------
1. **Forbidden traces.** An edge joining two class-A slots of the SAME block
   contracts two indices inside an antisymmetric group, which vanishes
   identically. Rejected before any arithmetic.
2. **Degree feasibility.** Every block must use exactly its slot multiplicity.
3. **Young symmetry / canonicalisation.** Topologies related by relabelling
   the blocks describe the same scalar. A canonical form under the block
   permutation group removes those duplicates.
4. **Connectivity.** A disconnected topology factorises into a product of
   lower-degree invariants, which lives in the product subspace and cannot
   contribute a primitive direction.

Index placement follows the same rule the forward work established: every
contracted edge carries exactly one raised end, or it contracts with delta
instead of eta and is not a Lorentz scalar.
"""

import hashlib
import itertools
from collections import Counter

import numpy as np

from .modp import P, mod_einsum
from .stress import (
    _raise_axes,
    composite_n1050,
    composite_n4125,
    five_form_moment,
)

# name -> (number of slots, tuple of class sizes)
# class 0 is the antisymmetric group where present.
BLOCK_KINDS = {
    "M": (2, (2,)),
    "N1050": (6, (5, 1)),
    "N4125": (6, (6,)),
}

ALPHABET = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"


def block_slot_classes(kind):
    """Return a list of (class_index, slot_indices) for one block kind."""
    _, classes = BLOCK_KINDS[kind]
    out, start = [], 0
    for ci, size in enumerate(classes):
        out.append((ci, tuple(range(start, start + size))))
        start += size
    return out


def endpoints(kinds):
    """All (block, class) endpoints with their capacities."""
    eps = []
    for bi, kind in enumerate(kinds):
        for ci, slots in block_slot_classes(kind):
            eps.append(((bi, ci), len(slots)))
    return eps


def _is_forbidden_trace(kinds, a, b):
    """An edge inside one block's antisymmetric class vanishes identically."""
    (ba, ca), (bb, cb) = a, b
    if ba != bb:
        return False
    if ca != cb:
        return False                      # A-B within a block: allowed
    _, classes = BLOCK_KINDS[kinds[ba]]
    return classes[ca] > 1               # antisymmetric group -> zero


def enumerate_topologies(kinds, max_results=None):
    """Enumerate edge-count assignments between endpoints.

    Yields dicts {(endpoint_a, endpoint_b): count}, with endpoints sorted, that
    saturate every endpoint's capacity and contain no forbidden trace.
    """
    eps = endpoints(kinds)
    keys = [k for k, _ in eps]
    cap = dict(eps)
    pairs = [(a, b) for i, a in enumerate(keys) for b in keys[i:]]
    pairs = [(a, b) for a, b in pairs if not _is_forbidden_trace(kinds, a, b)]

    results = []

    def recurse(idx, remaining, chosen):
        if max_results is not None and len(results) >= max_results:
            return
        if idx == len(pairs):
            if all(v == 0 for v in remaining.values()):
                results.append(dict(chosen))
            return
        a, b = pairs[idx]
        # an edge consumes one slot at each end (two at the same end if a == b)
        limit = min(remaining[a], remaining[b])
        if a == b:
            limit = remaining[a] // 2
        # prune: if the total remaining capacity is odd the branch is dead
        for count in range(limit, -1, -1):
            nxt = dict(remaining)
            if a == b:
                nxt[a] -= 2 * count
            else:
                nxt[a] -= count
                nxt[b] -= count
            if any(v < 0 for v in nxt.values()):
                continue
            if sum(nxt.values()) % 2:
                continue
            if count:
                chosen[(a, b)] = count
            recurse(idx + 1, nxt, chosen)
            chosen.pop((a, b), None)

    recurse(0, dict(cap), {})
    return results


def iter_topologies(kinds, cap=None):
    """Generator form of `enumerate_topologies`.

    The list form materialises every topology before anything is filtered,
    which is what drove enumeration RSS to ~2.6 GB across 21 sectors: up to
    30 000 dicts per sector, all live at once, none released until the sector
    finished. Yielding lets the caller canonicalise and discard immediately, so
    only the compact canonical KEYS accumulate.
    """
    eps = endpoints(kinds)
    keys = [k for k, _ in eps]
    cap_by = dict(eps)
    pairs = [(a, b) for i, a in enumerate(keys) for b in keys[i:]]
    pairs = [(a, b) for a, b in pairs if not _is_forbidden_trace(kinds, a, b)]
    produced = 0

    def recurse(idx, remaining, chosen):
        nonlocal produced
        if cap is not None and produced >= cap:
            return
        if idx == len(pairs):
            if all(v == 0 for v in remaining.values()):
                produced += 1
                yield dict(chosen)
            return
        a, b = pairs[idx]
        limit = remaining[a] // 2 if a == b else min(remaining[a], remaining[b])
        for count in range(limit, -1, -1):
            nxt = dict(remaining)
            if a == b:
                nxt[a] -= 2 * count
            else:
                nxt[a] -= count
                nxt[b] -= count
            if any(v < 0 for v in nxt.values()) or sum(nxt.values()) % 2:
                continue
            if count:
                chosen[(a, b)] = count
            yield from recurse(idx + 1, nxt, chosen)
            chosen.pop((a, b), None)
            if cap is not None and produced >= cap:
                return

    yield from recurse(0, dict(cap_by), {})


def stream_candidates(kinds, cap=None, shard_residue=0, shard_modulus=1,
                      require_connected=True, max_candidates=None):
    """Stream canonical candidates, retaining only compact canonical keys.

    Memory is bounded by the number of DISTINCT canonical keys, not by the raw
    topology count, and a key is a tuple of ints rather than a dict of tuples.

    Sharding is deterministic: a topology belongs to this shard when its
    canonical key hashes to `shard_residue` mod `shard_modulus`. Because the
    test is on the CANONICAL key, two shards can never both claim the same
    scalar, and the union over residues is exactly the unsharded set.
    """
    seen = set()
    stats = {"raw": 0, "disconnected": 0, "duplicate": 0, "yielded": 0,
             "truncated": False}
    for topo in iter_topologies(kinds, cap=cap):
        stats["raw"] += 1
        if require_connected and not is_connected(kinds, topo):
            stats["disconnected"] += 1
            continue
        key = canonical_form(kinds, topo)
        if key in seen:
            stats["duplicate"] += 1
            continue
        seen.add(key)
        if shard_modulus > 1:
            digest = int(hashlib.sha256(repr(key).encode()).hexdigest()[:8], 16)
            if digest % shard_modulus != shard_residue:
                continue
        stats["yielded"] += 1
        yield topo, stats
        if max_candidates is not None and stats["yielded"] >= max_candidates:
            break
    if cap is not None and stats["raw"] >= cap:
        stats["truncated"] = True


def canonical_form(kinds, topology):
    """Canonical key under relabelling blocks of the same kind.

    Two topologies related by permuting like blocks describe the same scalar.
    Taking the lexicographic minimum over that group is the canonicalisation
    step; it is what collapses the raw enumeration to distinct candidates.
    """
    n = len(kinds)
    best = None
    for perm in itertools.permutations(range(n)):
        if any(kinds[perm[i]] != kinds[i] for i in range(n)):
            continue                      # only permute like kinds
        relabelled = []
        for (a, b), count in topology.items():
            aa = ((perm.index(a[0])), a[1])
            bb = ((perm.index(b[0])), b[1])
            relabelled.append((min(aa, bb), max(aa, bb), count))
        key = tuple(sorted(relabelled))
        if best is None or key < best:
            best = key
    return best


def is_connected(kinds, topology):
    """Disconnected topologies factorise into products of lower invariants."""
    n = len(kinds)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for (a, b), count in topology.items():
        if count and a[0] != b[0]:
            ra, rb = find(a[0]), find(b[0])
            if ra != rb:
                parent[ra] = rb
    return len({find(i) for i in range(n)}) == 1


def build_einsum(kinds, topology):
    """Materialise a topology as an einsum spec plus a raise-set per block.

    Slots inside a class are interchangeable, so they are handed out in order.
    Each edge gets one letter; the END ON THE HIGHER-NUMBERED BLOCK is raised,
    which is a deterministic choice of the "exactly one raised end" rule. Which
    end carries the metric cannot change the value of a genuine contraction.
    """
    free = {}
    for bi, kind in enumerate(kinds):
        for ci, slots in block_slot_classes(kind):
            free[(bi, ci)] = list(slots)

    subscripts = [[None] * BLOCK_KINDS[k][0] for k in kinds]
    raises = [set() for _ in kinds]
    letter = 0
    for (a, b), count in sorted(topology.items()):
        for _ in range(count):
            if letter >= len(ALPHABET):
                raise ValueError("einsum letter budget exceeded")
            ch = ALPHABET[letter]
            letter += 1
            sa = free[a].pop(0)
            sb = free[b].pop(0)
            subscripts[a[0]][sa] = ch
            subscripts[b[0]][sb] = ch
            # raise exactly one end, deterministically
            hi = b if b[0] >= a[0] else a
            raises[hi[0]].add(sb if hi is b else sa)
    if any(None in s for s in subscripts):
        raise ValueError("topology left a slot unassigned")
    spec = ",".join("".join(s) for s in subscripts) + "->"
    return spec, [tuple(sorted(r)) for r in raises]


GREEK = ["mu", "nu", "rho", "sigma", "tau", "alpha", "beta", "gamma",
         "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa",
         "lambda", "chi", "psi", "omega", "phi"]


def render_formula(kinds, topology, index_names=None):
    """Render a topology as a readable indexed contraction.

    A canonical identifier is not mathematics. This produces the actual
    contraction: every tensor factor, every index, and its variance, with `^`
    marking a raised (contravariant) slot and `_` a lowered one. Exactly one
    end of every contracted pair is raised, which is what makes the expression
    a Lorentz scalar.

    Round-trips through `parse_formula`.
    """
    spec, raises = build_einsum(kinds, topology)
    subs = spec.split("->")[0].split(",")
    letters = sorted({ch for sub in subs for ch in sub})
    if index_names is None:
        index_names = {ch: GREEK[i % len(GREEK)] + ("'" * (i // len(GREEK)))
                       for i, ch in enumerate(letters)}
    parts = []
    for bi, (kind, sub) in enumerate(zip(kinds, subs)):
        slots = []
        for slot, ch in enumerate(sub):
            up = slot in raises[bi]
            slots.append(("^" if up else "_") + index_names[ch])
        parts.append(f"{kind}[{bi}]{{{' '.join(slots)}}}")
    return " * ".join(parts)


def parse_formula(text):
    """Parse a rendered formula back into (kinds, topology).

    Deliberately reconstructs the topology from the INDEX PAIRINGS rather than
    trusting any stored structure, so a round-trip failure means the rendered
    formula does not actually describe the topology it came from.
    """
    kinds, per_block = [], []
    for factor in [f.strip() for f in text.split("*")]:
        head, _, rest = factor.partition("[")
        _, _, body = rest.partition("{")
        body = body.rstrip("}").strip()
        kinds.append(head.strip())
        slots = []
        for token in body.split():
            slots.append((token[0] == "^", token[1:]))
        per_block.append(slots)

    # locate the two ends of each index, then map slots back to slot classes
    where = {}
    for bi, slots in enumerate(per_block):
        for slot, (_up, name) in enumerate(slots):
            where.setdefault(name, []).append((bi, slot))

    def slot_class(kind, slot):
        for ci, s in block_slot_classes(kind):
            if slot in s:
                return ci
        raise ValueError(f"slot {slot} not in any class of {kind}")

    topology = {}
    for name, ends in where.items():
        if len(ends) != 2:
            raise ValueError(f"index {name!r} appears {len(ends)} times, not 2")
        (b1, s1), (b2, s2) = ends
        a = (b1, slot_class(kinds[b1], s1))
        b = (b2, slot_class(kinds[b2], s2))
        key = (min(a, b), max(a, b))
        topology[key] = topology.get(key, 0) + 1
    return kinds, topology


def make_blocks(five_form, mod=P, backend="optimized"):
    """The three verified quadratic-block channels, ALL-LOWER.

    `build_einsum` assumes every operand is fully covariant and raises exactly
    one end of each contracted edge itself. `five_form_moment` returns `mixed`
    as M_{a}{}^{b} -- slot 0 DOWN, slot 1 UP -- so handing it over unchanged
    made every M-block edge carry either two raised ends or none.

    This shipped, and a boost test over generated candidates is what caught it:
    pure-N sectors were 100% boost invariant while M-containing sectors were
    not, and 129 of 579 pilot candidates failed to lie in the degree-10 atlas
    span -- which a genuine Lorentz scalar cannot do. Lowering the second index
    (an involution for a diagonal metric) makes M all-lower like the others.
    """
    _, mixed = five_form_moment(five_form, mod, backend)
    return {
        "M": _raise_axes(mixed, (1,), mod),        # M_{a}{}^{b} -> M_{ab}
        "N1050": composite_n1050(five_form, mod, backend),
        "N4125": composite_n4125(five_form, mod, backend),
    }


def evaluate(kinds, topology, blocks, mod=P):
    """Exact modular evaluation of one candidate."""
    spec, raises = build_einsum(kinds, topology)
    operands = []
    for kind, raise_set in zip(kinds, raises):
        base = blocks[kind]
        operands.append(_raise_axes(base, raise_set, mod) if raise_set
                        else base % mod)
    return int(mod_einsum(spec, operands, mod) % mod)


def multiset_slot_count(kinds):
    return sum(BLOCK_KINDS[k][0] for k in kinds)


def block_multisets(size=5):
    """All block multisets of the given size, CHEAPEST FIRST.

    Ordered by total slot count ascending, so `M`-heavy multisets (10 slots)
    come before `N`-heavy ones (30 slots). The enumeration cost grows steeply
    with slot count -- 30 slots means 15 edges distributed over ~50 endpoint
    pairs -- so ordering by cost lets a bounded pilot cover whole sectors
    instead of stalling inside the most expensive one.

    Enumerated exhaustively rather than chosen: restricting to a favoured
    multiset would import knowledge of the published answer.
    """
    kinds = sorted(BLOCK_KINDS)
    out = {tuple(sorted(c))
           for c in itertools.combinations_with_replacement(kinds, size)}
    return sorted(out, key=lambda c: (multiset_slot_count(c), c))


def generate_candidates(kinds, max_topologies=None, require_connected=True):
    """Full pruned pipeline for one block multiset.

    Returns (candidates, stats) where candidates are canonical topologies.
    """
    raw = enumerate_topologies(kinds, max_results=max_topologies)
    stats = {"raw_topologies": len(raw)}

    connected = [t for t in raw if (not require_connected)
                 or is_connected(kinds, t)]
    stats["after_connectivity"] = len(connected)

    seen, canonical = set(), []
    for topo in connected:
        key = canonical_form(kinds, topo)
        if key in seen:
            continue
        seen.add(key)
        canonical.append(topo)
    stats["canonical"] = len(canonical)
    return canonical, stats
