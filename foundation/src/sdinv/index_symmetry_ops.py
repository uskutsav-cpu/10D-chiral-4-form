"""Exact nested index (anti)symmetrisation with explicit execution order.

The source convention that makes this necessary, quoted from
"Some remarks on invariants", J. Phys. A 59 (2026) 065203, eq (4.24):

    "the symmetrization and/or anti-symmetrization of the indices within the
     red brackets is made upon the anti-symmetrization within the black
     brackets"

So bracket groups are NOT interchangeable: black-stage antisymmetrisations run
first, red-stage operations second. Flattening them, or applying red before
black, gives a different tensor. `test_index_symmetry_ops.py` contains an
explicit example where the two orders disagree, so the ordering cannot be
silently dropped.

Arithmetic is modular throughout. Nothing here uses a bare ``np.einsum`` on
unbounded operands -- permutation sums are accumulated with an explicit ``%``
after every term, which keeps every intermediate below ``p`` and therefore far
inside int64 regardless of tensor rank.
"""

from dataclasses import dataclass, field
from itertools import permutations
from typing import Tuple

import numpy as np

from .modp import P, inv

BLACK = "black"
RED = "red"


def permutation_sign(perm):
    """Sign of a permutation given as a sequence of positions."""
    perm = list(perm)
    seen = [False] * len(perm)
    sign = 1
    for i in range(len(perm)):
        if seen[i]:
            continue
        j, length = i, 0
        while not seen[j]:
            seen[j] = True
            j = perm[j]
            length += 1
        if length % 2 == 0:
            sign = -sign
    return sign


@dataclass(frozen=True)
class BracketOp:
    """One bracket operation.

    kind        : "antisym" or "sym"
    slots       : tensor axes participating, as a tuple
    stage       : BLACK or RED -- BLACK executes first
    normalized  : divide by len(slots)! if True (the paper's convention for
                  the displayed brackets); False leaves the raw sum
    source      : provenance string, e.g. "eq (4.24) I^(5)_10 black group 1"
    """

    kind: str
    slots: Tuple[int, ...]
    stage: str = BLACK
    normalized: bool = True
    source: str = ""

    def __post_init__(self):
        if self.kind not in ("antisym", "sym"):
            raise ValueError(f"unknown bracket kind {self.kind!r}")
        if self.stage not in (BLACK, RED):
            raise ValueError(f"unknown stage {self.stage!r}")
        if len(set(self.slots)) != len(self.slots):
            raise ValueError(f"repeated slot in {self.slots}")

    def serialize(self):
        return (f"{self.stage}:{self.kind}{list(self.slots)}"
                f"{'/norm' if self.normalized else '/raw'}")


def apply_bracket(tensor, op, mod=P):
    """Apply a single bracket operation exactly, modulo `mod`."""
    slots = list(op.slots)
    if len(slots) < 2:
        return tensor % mod
    total = np.zeros_like(tensor)
    for perm in permutations(range(len(slots))):
        axes = list(range(tensor.ndim))
        for target, source in zip(slots, perm):
            axes[target] = slots[source]
        term = np.transpose(tensor, axes)
        if op.kind == "antisym":
            s = permutation_sign(perm)
            total = (total + s * term) % mod
        else:
            total = (total + term) % mod
    if op.normalized:
        total = (total * inv(_factorial(len(slots)) % mod, mod)) % mod
    return total % mod


def _factorial(n):
    out = 1
    for k in range(2, n + 1):
        out *= k
    return out


@dataclass
class BracketProgram:
    """An ordered list of bracket operations with staged execution.

    Execution order is BLACK stage first, in listed order, then RED stage, in
    listed order. That mirrors the source statement exactly and is asserted
    rather than assumed: `stages()` exposes the order for inspection and the
    tests compare it against the reversed order.
    """

    ops: list = field(default_factory=list)
    source: str = ""

    def stages(self):
        black = [o for o in self.ops if o.stage == BLACK]
        red = [o for o in self.ops if o.stage == RED]
        return black, red

    def serialize(self):
        black, red = self.stages()
        return " ; ".join(o.serialize() for o in black + red)

    def apply(self, tensor, mod=P, reverse_stages=False):
        black, red = self.stages()
        order = (red + black) if reverse_stages else (black + red)
        out = tensor % mod
        for op in order:
            out = apply_bracket(out, op, mod)
        return out % mod
