"""Readable formulas with exactly one raised endpoint per graph edge."""

from .graphs import validate_graph
from .serialize import to_record


def graph_to_latex(matrix):
    validate_graph(matrix, valence=5, max_mult=4)
    lower, upper = [[] for _ in matrix], [[] for _ in matrix]
    for i, j, multiplicity in to_record(matrix)["edges"]:
        for k in range(1, multiplicity + 1):
            name = rf"\mu_{{{i},{j},{k}}}"
            lower[i].append(name)
            upper[j].append(name)
    # Lexicographic edge order places every incoming (raised) slot before
    # every outgoing (lowered) slot. This is the evaluator's exact slot order.
    return r"\,".join("F" + ("^{" + " ".join(upper[v]) + "}" if upper[v] else "")
                      + ("_{" + " ".join(lower[v]) + "}" if lower[v] else "")
                      for v in range(len(matrix)))


def factor_lines(matrix, factors_per_line=2):
    """A multiline aligned expression suitable for a paper appendix."""
    factors = graph_to_latex(matrix).split(r"\,")
    return " \\\\\n &\\quad ".join(r"\,".join(factors[k:k + factors_per_line])
                                for k in range(0, len(factors), factors_per_line))
