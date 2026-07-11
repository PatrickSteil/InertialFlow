#!/usr/bin/env python3
"""Export a real InertialFlow max-flow subproblem as a DIMACS .max file.

InertialFlow never solves a single big max-flow problem on the whole graph.
Each recursion step instead: sorts the current node set by projecting
coordinates onto some direction, grows a source set from the low end and a
sink set from the high end until each accumulates >= `fraction` of the total
vertex weight, and asks for max-flow between "all sources" and "all sinks"
(implemented in the C++ code via infinite-capacity terminal edges, see
Partitioner::evaluate_cut). Contracting every source-set vertex into one
physical node and every sink-set vertex into another is mathematically
equivalent to that terminal-edge formulation, and produces an ordinary
two-terminal max-flow instance -- the kind maxflow_algorithms' `demo`/`bench`
tools expect.

Usage:
    python3 export_dimacs_maxflow.py -g graph.gr -c coords.co -o problem.max
    python3 export_dimacs_maxflow.py -g graph.gr -c coords.co -o problem.max \
        --fraction 0.2 --angle-deg 0

Convert the result to the repo's binary format for benchmarking:
    bench_io dimacs_to_bbk problem.max
"""

import argparse
import math
import sys
from collections import defaultdict


def read_coordinates(coord_path):
    """Mirrors Partitioner::loadCoordinates: 'v <id> <lon> <lat> [weight]'."""
    lon, lat, weight = {}, {}, {}
    next_id = 1
    with open(coord_path) as f:
        for line in f:
            if not line.startswith("v"):
                continue
            parts = line.split()
            vid = int(parts[1])
            if vid != next_id:
                sys.exit(f"Malformed coordinate file: expected id {next_id}, got {vid}")
            lon[vid] = float(parts[2])
            lat[vid] = float(parts[3])
            weight[vid] = int(parts[4]) if len(parts) > 4 else 1
            next_id += 1
    if not lon:
        sys.exit(f"Coordinate file has no vertices: {coord_path}")
    return lon, lat, weight


def read_edges(graph_path, num_nodes, weight):
    """Mirrors Partitioner::loadDimacsEdges: 'a <u> <v> <cap>', 'n <id> <w>'."""
    edges = []
    with open(graph_path) as f:
        for line in f:
            if line.startswith("n"):
                _, vid, w = line.split()
                weight[int(vid)] = int(w)
            elif line.startswith("a"):
                _, u, v, cap = line.split()
                u, v, cap = int(u), int(v), int(cap)
                if u == v:
                    continue  # self loops never affect a cut
                edges.append((u, v, cap))
    if not edges:
        sys.exit(f"Graph file has no edges: {graph_path}")
    return edges


def pick_source_sink_sets(node_ids, lon, lat, weight, fraction, angle_deg):
    """Mirrors the two-pointer selection in Partitioner::run_projection_task."""
    dx, dy = math.cos(math.radians(angle_deg)), math.sin(math.radians(angle_deg))
    ordered = sorted(node_ids, key=lambda i: lon[i] * dx + lat[i] * dy)

    total_weight = sum(weight[i] for i in ordered)
    # Partitioner::submit_bisect_task does
    # `target_weight = static_cast<long long>(total_weight * fraction)`,
    # which truncates towards zero -- not rounds. int() on a positive float
    # does the same, but keep this explicit since getting it wrong shifts
    # the source/sink boundary by a node and changes the resulting flow.
    target = int(total_weight * fraction)

    src_end, lo_weight = 0, 0
    while src_end < len(ordered) and lo_weight < target:
        lo_weight += weight[ordered[src_end]]
        src_end += 1

    snk_start, hi_weight = len(ordered), 0
    while snk_start > src_end and hi_weight < target:
        snk_start -= 1
        hi_weight += weight[ordered[snk_start]]

    sources = set(ordered[:src_end])
    sinks = set(ordered[snk_start:])
    if not sources or not sinks:
        sys.exit("fraction too small for this graph: source or sink set is empty")
    return sources, sinks


def contract_and_write(out_path, node_ids, edges, sources, sinks):
    S, T = 1, 2
    middle_ids = [i for i in node_ids if i not in sources and i not in sinks]
    remap = {i: S for i in sources}
    remap.update((i, T) for i in sinks)
    remap.update((i, 3 + k) for k, i in enumerate(middle_ids))

    # Partitioner::evaluate_cut() calls graph.add_edge(u, v, cap, cap) --
    # i.e. it gives every DIMACS 'a' arc capacity cap in *both* directions,
    # not just forward. Mirror that here, or the exported instance's max
    # flow won't match what InertialFlow actually solves.
    #
    # Sum parallel arcs created by contraction (flow-equivalent to keeping
    # them as separate arcs, but keeps the output file small); self loops
    # on S or T (both endpoints were in the same contracted set) are simply
    # dropped, same as in the original graph loader.
    merged = defaultdict(int)
    for u, v, cap in edges:
        if u not in remap or v not in remap:
            continue  # edge touches a vertex outside this subproblem
        ru, rv = remap[u], remap[v]
        if ru != rv:
            merged[(ru, rv)] += cap
            merged[(rv, ru)] += cap

    num_nodes = 2 + len(middle_ids)
    with open(out_path, "w") as f:
        f.write(f"c InertialFlow subproblem: {len(sources)} sources, "
                f"{len(sinks)} sinks, {len(middle_ids)} middle vertices\n")
        f.write(f"p max {num_nodes} {len(merged)}\n")
        f.write(f"n {S} s\n")
        f.write(f"n {T} t\n")
        for (u, v), cap in merged.items():
            f.write(f"a {u} {v} {cap}\n")

    print(f"wrote {out_path}: {num_nodes} nodes, {len(merged)} arcs "
          f"({len(sources)} sources + {len(sinks)} sinks contracted)",
          file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-g", "--graph", required=True, help="DIMACS .gr edge file")
    p.add_argument("-c", "--coord", required=True, help="DIMACS .co coordinate file")
    p.add_argument("-o", "--output", required=True, help="Output .max file")
    p.add_argument("--fraction", type=float, default=0.2,
                    help="Fraction of total vertex weight per side (default: 0.2)")
    p.add_argument("--angle-deg", type=float, default=0.0,
                    help="Projection direction in degrees; 0 matches InertialFlow's "
                         "first projection (default: 0)")
    args = p.parse_args()

    if not (0 < args.fraction < 0.5):
        sys.exit("--fraction must be in (0, 0.5)")

    lon, lat, weight = read_coordinates(args.coord)
    edges = read_edges(args.graph, len(lon), weight)
    node_ids = list(lon.keys())

    sources, sinks = pick_source_sink_sets(
        node_ids, lon, lat, weight, args.fraction, args.angle_deg)
    contract_and_write(args.output, node_ids, edges, sources, sinks)


if __name__ == "__main__":
    main()
