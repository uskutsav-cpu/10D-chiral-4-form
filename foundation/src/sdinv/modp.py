"""Exact linear algebra over F_p.

Why: invariant values span a huge dynamic range at high order, so float SVD
forces you to guess a rank tolerance. Over F_p, zero means zero.

p is chosen ~2^15 so that products fit comfortably in int64 during einsum
accumulation (product < 2^30, summed over <2^17 terms stays < 2^47).
Rank over F_p at a random point equals the true rank with probability
1 - O(degree/p). Re-run with ALT_P to confirm.
"""

import numpy as np

P = 32749  # prime, ~2^15 -- see the overflow guard in mod_einsum
ALT_P = 32719  # second prime for confirmation runs
FLOAT_BLAS_MIN_WORK = 10_000_000


def inv(a, p=P):
    """Modular inverse via Fermat."""
    return pow(int(a) % p, p - 2, p)


class RankSieve:
    """Streaming row-reduced basis over F_p.

    Feed it one Jacobian row at a time. Returns True if the row was
    independent of everything seen so far (i.e. a genuinely new invariant),
    False if it reduced to zero (i.e. you found a syzygy).

    Memory is O(rank * ncols), not O(ncandidates * ncols).
    """

    def __init__(self, ncols, p=P):
        self.p = p
        self.ncols = ncols
        self.rows = []   # reduced basis rows
        self.pivots = []  # pivot column of each basis row

    def add(self, row):
        row = np.asarray(row, dtype=np.int64) % self.p
        for r, pc in zip(self.rows, self.pivots):
            if row[pc]:
                row = (row - int(row[pc]) * r) % self.p
        nz = np.nonzero(row)[0]
        if len(nz) == 0:
            return False
        pc = int(nz[0])
        row = (row * inv(row[pc], self.p)) % self.p
        # back-substitute into existing rows to keep basis fully reduced
        for i, (r, q) in enumerate(zip(self.rows, self.pivots)):
            if r[pc]:
                self.rows[i] = (r - int(r[pc]) * row) % self.p
        self.rows.append(row)
        self.pivots.append(pc)
        return True

    @property
    def rank(self):
        return len(self.rows)

    def to_state(self):
        """Return a JSON-serializable, exact checkpoint state."""
        return {
            "p": int(self.p),
            "ncols": int(self.ncols),
            "pivots": [int(x) for x in self.pivots],
            "rows": [
                [int(x) for x in np.asarray(row, dtype=np.int64)]
                for row in self.rows
            ],
        }

    @classmethod
    def from_state(cls, state):
        """Restore and validate a sieve produced by to_state()."""
        sieve = cls(int(state["ncols"]), int(state["p"]))
        pivots = [int(x) for x in state["pivots"]]
        rows = [
            np.asarray(row, dtype=np.int64) % sieve.p
            for row in state["rows"]
        ]
        if len(rows) != len(pivots):
            raise ValueError("rank-sieve row/pivot count mismatch")
        if any(row.shape != (sieve.ncols,) for row in rows):
            raise ValueError("rank-sieve row has wrong column count")
        if len(pivots) != len(set(pivots)):
            raise ValueError("rank-sieve pivots must be unique")
        for row, pivot in zip(rows, pivots):
            if not 0 <= pivot < sieve.ncols or int(row[pivot]) != 1:
                raise ValueError("rank-sieve pivot is not normalized")
            if np.any(row[:pivot]):
                raise ValueError("rank-sieve row is not in echelon form")
        for i, pivot in enumerate(pivots):
            if any(int(rows[j][pivot]) for j in range(len(rows)) if j != i):
                raise ValueError("rank-sieve basis is not fully reduced")
        sieve.rows = rows
        sieve.pivots = pivots
        return sieve


def mod_einsum(subscripts, operands, p=P):
    """einsum over F_p, contracting STRICTLY TWO AT A TIME with reduction
    after every step.

    Two traps this avoids, both of which produce wrong answers silently:

    1. OVERFLOW. np.einsum multiplies all operands before summing, so at
       order n the products reach p^n. Even asking numpy for an optimised
       path does not help -- for a symmetric contraction it will happily
       return a single step contracting all n operands at once. int64
       wraps with no warning and the result is garbage. We therefore never
       hand einsum more than two operands.

    2. MEMORY. Contracting in a bad order builds enormous intermediates
       (order 8 in 6D wanted 16 GiB). We greedily pick the pair whose
       result has the fewest free indices.

    An explicit guard asserts p^2 * (terms summed) stays inside int64.
    """
    ins, out = subscripts.split("->")
    terms = ins.split(",")
    ops = [np.asarray(o, dtype=np.int64) % p for o in operands]

    while len(ops) > 1:
        best = None
        for i in range(len(ops)):
            for j in range(i + 1, len(ops)):
                rest = set("".join(terms[k] for k in range(len(terms))
                                   if k not in (i, j))) | set(out)
                keep = [c for c in dict.fromkeys(terms[i] + terms[j])
                        if c in rest]
                # Cost of this pair is output size TIMES summed size, i.e.
                # the product over every index appearing in either operand.
                # Ranking by len(keep) alone minimises the intermediate but
                # not the work: a step can have a tiny output and still sum
                # ~1e9 terms, which passes the overflow guard below and then
                # runs for hours on numpy's scalar integer loop (there is no
                # BLAS for int64). Rank by work, break ties on output size so
                # the original memory behaviour is preserved.
                d = {}
                for t, o in ((terms[i], ops[i]), (terms[j], ops[j])):
                    for c, n in zip(t, o.shape):
                        d[c] = n
                work = 1
                for c in set(terms[i] + terms[j]):
                    work *= d[c]
                outsz = 1
                for c in keep:
                    outsz *= d[c]
                if best is None or (work, outsz) < (best[0], best[1]):
                    best = (work, outsz, i, j, keep)
        _, _, i, j, keep = best
        ti, tj, oi, oj = terms[i], terms[j], ops[i], ops[j]

        dims = {}
        for t, o in ((ti, oi), (tj, oj)):
            for c, n in zip(t, o.shape):
                dims[c] = n
        nterms = 1
        for c in set(ti + tj) - set(keep):
            nterms *= dims[c]
        max_sum = (p - 1) ** 2 * nterms
        if max_sum >= 2 ** 63:
            raise OverflowError(
                f"int64 would overflow: p^2 * {nterms} summed terms. "
                f"Lower P in modp.py.")

        pair_subscripts = f"{ti},{tj}->{''.join(keep)}"
        work = best[0]
        if work >= FLOAT_BLAS_MIN_WORK and max_sum < 2 ** 53:
            # NumPy's integer einsum does not use BLAS and can be hundreds of
            # times slower on wide contractions. Float64 is still exact here:
            # inputs and products are integers, and the guard proves every
            # unreduced dot-product sum is below the 53-bit integer limit.
            raw = np.einsum(
                pair_subscripts,
                oi.astype(np.float64),
                oj.astype(np.float64),
                optimize=True,
            )
            res = np.rint(raw).astype(np.int64) % p
        else:
            res = np.einsum(pair_subscripts, oi, oj) % p
        for k in sorted((i, j), reverse=True):
            terms.pop(k)
            ops.pop(k)
        terms.append("".join(keep))
        ops.append(res)

    if terms[0] != out:
        ops[0] = np.einsum(f"{terms[0]}->{out}", ops[0]) % p
    return ops[0]
