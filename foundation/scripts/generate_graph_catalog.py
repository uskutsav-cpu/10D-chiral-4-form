"""Generate the exact 10D contraction-graph catalog with nauty gtools.

Install nauty (https://pallini.di.uniroma1.it/) so that `geng` and `multig`
are on PATH, or pass their paths explicitly.
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sdinv.graphs import generate_graphs_nauty, graph_label


D, VALENCE = 10, 5
MAX_MULTIPLICITY = VALENCE - 1


def build_catalog(orders, geng="geng", multig="multig"):
    catalog = {
        "schema": 1,
        "description": (
            "Connected loop-free valence-5 multigraphs, maximum edge "
            "multiplicity 4, deduplicated exactly up to isomorphism."
        ),
        "generator": {
            "software": "nauty gtools 2.9.3",
            "pipeline": (
                "geng -cq ORDER | multig -q -T -m4 -r5"
            ),
            "geng": os.path.basename(geng),
            "multig": os.path.basename(multig),
        },
        "dimension": D,
        "valence": VALENCE,
        "max_multiplicity": MAX_MULTIPLICITY,
        "orders": {},
    }
    for order in orders:
        graphs = generate_graphs_nauty(
            order,
            VALENCE,
            MAX_MULTIPLICITY,
            geng=geng,
            multig=multig,
        )
        labels = [graph_label(M) for M in graphs]
        catalog["orders"][str(order)] = {
            "count": len(labels),
            "graphs": labels,
        }
        print(f"order {order}: {len(labels)} exact isomorphism classes")
    return catalog


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders", type=int, nargs="+", default=[4, 6, 8])
    parser.add_argument("--geng", default="geng")
    parser.add_argument("--multig", default="multig")
    parser.add_argument("--out", default="results/10d_graph_catalog.json")
    args = parser.parse_args()

    result = build_catalog(args.orders, args.geng, args.multig)
    parent = os.path.dirname(args.out)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(args.out, "w") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
