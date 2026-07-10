#!/usr/bin/env python3
"""Plot a DIMACS graph, optionally colored by an InertialFlow partition.

Reads the same `.co` (coordinates) and `.gr` (edges) files the
`partitioner` CLI tool reads, plus (optionally) the partition file it
writes via `-o`, and renders a scatter/line plot: vertices positioned by
their (lon, lat) coordinates, edges drawn between them, and vertices
colored by cell id if a partition file is given.

Usage:
    python3 plot_partition.py -c coords.co -g graph.gr [-p partition.txt] [-o plot.png]

Examples:
    # Just the raw graph, single color.
    python3 plot_partition.py -c examples/example2.co -g examples/example2.gr

    # Graph colored by partition, after running the partitioner.
    ./build/partitioner -g examples/example2.gr -c examples/example2.co \\
        -k 4 -o /tmp/partition.txt
    python3 plot_partition.py -c examples/example2.co -g examples/example2.gr \\
        -p /tmp/partition.txt -o /tmp/partition.png
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, List, Tuple


def read_coords(path: str) -> Dict[int, Tuple[float, float]]:
    """Parses a `.co` file into {vertex_id: (lon, lat)}.

    Mirrors Partitioner::loadGraph's parsing: each `v <id> <lon> <lat>`
    line is read in that field order.
    """
    coords: Dict[int, Tuple[float, float]] = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line[0] != "v":
                continue
            parts = line.split()
            if len(parts) != 4:
                continue
            _, vid, lon, lat = parts
            coords[int(vid)] = (float(lon), float(lat))
    if not coords:
        raise ValueError(f"No vertices found in coordinate file: {path}")
    return coords


def read_edges(path: str) -> List[Tuple[int, int]]:
    """Parses a `.gr` file into a list of (u, v) edges (1-based ids).

    Self loops (u == v) are dropped since they render as nothing useful.
    """
    edges: List[Tuple[int, int]] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line[0] != "a":
                continue
            parts = line.split()
            if len(parts) != 4:
                continue
            _, u, v, _cap = parts
            u, v = int(u), int(v)
            if u != v:
                edges.append((u, v))
    return edges


def read_partition(path: str) -> Dict[int, int]:
    """Parses a partition file (Partitioner::saveResults output) into
    {vertex_id: partition_id}."""
    partition: Dict[int, int] = {}
    with open(path) as f:
        first = f.readline()
        if not first.startswith("N "):
            raise ValueError(
                f"Unexpected partition file header (expected 'N <count>'): {first!r}"
            )
        for line in f:
            line = line.strip()
            if not line:
                continue
            vid, pid = line.split()
            partition[int(vid)] = int(pid)
    return partition


def plot(
    coords: Dict[int, Tuple[float, float]],
    edges: List[Tuple[int, int]],
    partition: Dict[int, int] | None,
    output_path: str | None,
    show: bool,
    draw_edges: bool,
    node_size: float,
    edge_alpha: float,
    title: str | None,
    dpi: int,
) -> None:
    import matplotlib

    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    fig, ax = plt.subplots(figsize=(10, 10))

    if draw_edges and edges:
        segments = [
            (coords[u], coords[v]) for u, v in edges if u in coords and v in coords
        ]
        ax.add_collection(
            LineCollection(segments, colors="0.6", linewidths=0.4, alpha=edge_alpha, zorder=1)
        )

    ids = sorted(coords)
    xs = [coords[i][0] for i in ids]
    ys = [coords[i][1] for i in ids]

    if partition is not None:
        missing = [i for i in ids if i not in partition]
        if missing:
            print(
                f"Warning: {len(missing)} vertex(es) missing from the partition "
                "file; plotting them in gray.",
                file=sys.stderr,
            )
        colors = [partition.get(i, -1) for i in ids]
        num_cells = max(colors) + 1 if colors else 0
        cmap = plt.get_cmap("tab20" if num_cells <= 20 else "hsv")
        scatter = ax.scatter(
            xs, ys, c=colors, cmap=cmap, s=node_size, zorder=2, linewidths=0
        )
        cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Cell id")
    else:
        ax.scatter(xs, ys, c="steelblue", s=node_size, zorder=2, linewidths=0)

    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("lon")
    ax.set_ylabel("lat")
    ax.set_title(
        title
        or (
            f"{len(coords)} vertices, {len(edges)} edges"
            + (f", {max(partition.values()) + 1} cells" if partition else "")
        )
    )
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=dpi)
        print(f"Wrote {output_path}")
    if show:
        plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot a DIMACS graph, optionally colored by an InertialFlow partition.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-c", "--coords", required=True, help="Path to the .co coordinate file.")
    parser.add_argument("-g", "--graph", required=True, help="Path to the .gr graph file.")
    parser.add_argument(
        "-p",
        "--partition",
        default=None,
        help="Path to a partition file written by `partitioner -o ...`. "
        "If omitted, all vertices are plotted in a single color.",
    )
    parser.add_argument(
        "-o", "--output", default=None, help="Path to save the plot to (e.g. plot.png)."
    )
    parser.add_argument(
        "--show", action="store_true", help="Also open an interactive window."
    )
    parser.add_argument(
        "--no-edges", action="store_true", help="Skip drawing edges (faster for large graphs)."
    )
    parser.add_argument(
        "--edge-alpha", type=float, default=0.3, help="Edge line opacity (default: 0.3)."
    )
    parser.add_argument(
        "--node-size", type=float, default=8.0, help="Vertex marker size (default: 8)."
    )
    parser.add_argument("--dpi", type=int, default=150, help="Output image DPI (default: 150).")
    parser.add_argument("--title", default=None, help="Custom plot title.")
    args = parser.parse_args()

    if not args.output and not args.show:
        parser.error("Nothing to do: pass -o/--output, --show, or both.")

    coords = read_coords(args.coords)
    edges = read_edges(args.graph)
    partition = read_partition(args.partition) if args.partition else None

    # Large graphs: default to skipping edges unless explicitly requested,
    # since drawing tens of thousands of line segments is slow and tends
    # to just produce a gray smear anyway.
    draw_edges = not args.no_edges
    if draw_edges and len(edges) > 200_000:
        print(
            f"Graph has {len(edges)} edges; skipping edge drawing for speed "
            "(pass nothing to change this, edges are auto-disabled above "
            "200k; there is currently no override).",
            file=sys.stderr,
        )
        draw_edges = False

    plot(
        coords,
        edges,
        partition,
        args.output,
        args.show,
        draw_edges,
        args.node_size,
        args.edge_alpha,
        args.title,
        args.dpi,
    )


if __name__ == "__main__":
    main()
