"""Robust nauty/genbg enumeration of loopless 5-regular contraction multigraphs.

For n1+n2 > 32, use genbgL. Brendan McKay's documented limits for genbgL are
n1 <= 30 and n1+n2 <= 64, which cover the degree-14 (14+35=49) and
degree-16 (16+40=56) incidence graphs used here.

Every multiedge occurrence is encoded by a degree-2 edge vertex in a bipartite
incidence graph; collapsing the right colour class recovers the multigraph.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


def _executable_for_sizes(n1: int, n2: int) -> str:
    total = n1 + n2
    if total <= 32:
        exe = shutil.which("genbg") or shutil.which("genbgL")
    else:
        exe = shutil.which("genbgL")
        if exe is None:
            raise RuntimeError(
                f"degree incidence graph has {total} vertices; nauty genbgL is required "
                "(standard genbg may be limited to 32 vertices). "
                "Install nauty 2.9+ with Homebrew: brew install nauty"
            )
    if exe is None:
        raise RuntimeError(
            "nauty genbg/genbgL not found. Install with: brew install nauty"
        )
    return exe


def require_genbg(total_vertices: int | None = None) -> str:
    if total_vertices is None:
        return shutil.which("genbg") or shutil.which("genbgL") or _raise_missing()
    # n1/n2 are not known here, so for >32 simply require genbgL.
    if total_vertices > 32:
        exe = shutil.which("genbgL")
        if exe is None:
            raise RuntimeError(
                f"genbgL required for {total_vertices} vertices; install nauty 2.9+"
            )
        return exe
    exe = shutil.which("genbg") or shutil.which("genbgL")
    if exe is None:
        raise RuntimeError("genbg/genbgL not found")
    return exe


def _raise_missing():
    raise RuntimeError("genbg/genbgL not found")


def decode_graph6(line):
    data = line.strip().encode("ascii") if isinstance(line, str) else line.strip()
    if not data:
        raise ValueError("empty graph6 record")
    if data.startswith(b">>graph6<<"):
        data = data[len(b">>graph6<<"):]
    if not data:
        raise ValueError("empty graph6 payload")
    if data[0] == 126:
        # Support standard extended graph6 even though current degree14/16
        # incidence graphs are below 63 vertices.
        if len(data) < 4:
            raise ValueError("truncated extended graph6 header")
        if data[1] != 126:
            n = ((data[1]-63) << 12) | ((data[2]-63) << 6) | (data[3]-63)
            payload = data[4:]
        else:
            raise ValueError("graph6 order >= 258048 is unsupported here")
    else:
        n = data[0] - 63
        payload = data[1:]
    if not 0 <= n <= 258047:
        raise ValueError("invalid graph6 order")

    bits = []
    for c in payload:
        v = c - 63
        if not 0 <= v < 64:
            raise ValueError("invalid graph6 character")
        bits.extend((v >> k) & 1 for k in (5,4,3,2,1,0))

    needed = n * (n - 1) // 2
    if len(bits) < needed:
        raise ValueError(
            f"truncated graph6 record: have {len(bits)} bits, need {needed}"
        )

    A = [[0] * n for _ in range(n)]
    q = 0
    for j in range(1, n):
        for i in range(j):
            if bits[q]:
                A[i][j] = A[j][i] = 1
            q += 1
    return A


def collapse_incidence_graph(A, left_n):
    n = len(A)
    if not 1 <= left_n < n:
        raise ValueError("invalid left colour-class size")
    right = range(left_n, n)
    M = [[0] * left_n for _ in range(left_n)]

    if any(sum(A[i]) != 5 for i in range(left_n)):
        raise ValueError("left incidence degrees are not all five")

    for r in right:
        nbr = [i for i in range(left_n) if A[r][i]]
        if len(nbr) != 2:
            raise ValueError("edge-vertex degree is not two")
        if any(A[r][s] for s in right):
            raise ValueError("right-right edge in supposedly bipartite graph")
        a, b = sorted(nbr)
        if a == b:
            raise ValueError("loop encountered in incidence graph")
        M[a][b] += 1
        M[b][a] += 1

    if any(sum(row) != 5 for row in M):
        raise ValueError("collapsed graph is not 5-regular")
    return tuple(tuple(x for x in row) for row in M)


def graph_label(M):
    parts = []
    n = len(M)
    for i in range(n):
        for j in range(i + 1, n):
            if M[i][j]:
                parts.append(f"{i}-{j}^{M[i][j]}")
    return f"n{n}[" + ",".join(parts) + "]"


def genbg_command(degree):
    if degree < 2 or degree % 2:
        raise ValueError("5-regular contraction graph needs even degree >=2")
    edge_vertices = 5 * degree // 2
    total = degree + edge_vertices

    # genbgL is the L1/large-size flavour required beyond the usual 32-vertex
    # small executable range.
    exe = _executable_for_sizes(degree, edge_vertices)

    # First colour class: exact degree 5.
    # Second colour class: exact degree 2.
    # Total incidence edges = 5*degree = 2*edge_vertices.
    incidence_edges = 5 * degree
    return [
        exe,
        "-q",
        "-l",
        "-g",
        "-d5:2",
        "-D5:2",
        str(degree),
        str(edge_vertices),
        f"{incidence_edges}:{incidence_edges}",
    ]


def _construct_known_incidence(degree):
    """Construct one explicit valid incidence graph for smoke testing.

    Uses the 5-regular circulant multigraph-free graph on even n:
    connect i to i±1, i±2, and i+n/2.
    """
    if degree < 6 or degree % 2:
        raise ValueError("smoke construction needs even degree >=6")

    edges = set()
    half = degree // 2
    for i in range(degree):
        for step in (1, 2):
            j = (i + step) % degree
            edges.add(tuple(sorted((i, j))))
            j = (i - step) % degree
            edges.add(tuple(sorted((i, j))))
        edges.add(tuple(sorted((i, (i + half) % degree))))

    edges = sorted(edges)
    if len(edges) != 5 * degree // 2:
        raise AssertionError(
            f"constructed graph has {len(edges)} edges, expected {5*degree//2}"
        )

    total = degree + len(edges)
    A = [[0] * total for _ in range(total)]
    for q, (a, b) in enumerate(edges, start=degree):
        A[a][q] = A[q][a] = 1
        A[b][q] = A[q][b] = 1
    return A


def local_incidence_smoke_test(degree):
    A = _construct_known_incidence(degree)
    M = collapse_incidence_graph(A, degree)
    return {
        "degree": degree,
        "incidence_vertices": len(A),
        "collapsed_label": graph_label(M),
        "passed": True,
    }


def stream_multigraphs(degree):
    cmd = genbg_command(degree)
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )
    stderr_text = ""
    produced = 0
    try:
        assert proc.stdout is not None
        assert proc.stderr is not None
        for line in proc.stdout:
            if not line.strip():
                continue
            produced += 1
            yield collapse_incidence_graph(decode_graph6(line), degree)

        stderr_text = proc.stderr.read()
        rc = proc.wait()
        if rc != 0:
            raise RuntimeError(
                "genbg process failed\n"
                f"command: {' '.join(cmd)}\n"
                f"return code: {rc}\n"
                f"stderr:\n{stderr_text}"
            )
        if produced == 0:
            raise RuntimeError(
                "genbg completed successfully but produced zero graphs\n"
                f"command: {' '.join(cmd)}\n"
                f"stderr:\n{stderr_text}\n"
                "This is treated as an enumerator failure, not as an empty "
                "invariant space."
            )
    finally:
        if proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()


def enumerator_smoke_test(degree, count=1):
    from itertools import islice
    local = local_incidence_smoke_test(degree)
    graphs = list(islice(stream_multigraphs(degree), count))
    if len(graphs) != count:
        raise RuntimeError(
            f"expected {count} streamed degree-{degree} graphs, got {len(graphs)}"
        )
    return {
        "local_constructed": local,
        "command": genbg_command(degree),
        "streamed_count": len(graphs),
        "first_labels": [graph_label(g) for g in graphs],
        "passed": True,
    }
