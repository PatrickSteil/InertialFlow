# InertialFlow

A multithreaded implementation of Inertial Flow graph partitioning:
recursive bisection where, at each step, several candidate cuts are
proposed by projecting vertex coordinates onto evenly-spaced directions,
and the cheapest cut (by max flow / min cut between the two extremes of
each projection) is kept. Recursion continues until the graph is split
into the requested number of cells.

Both axes of parallelism — the several candidate projections at one
recursion node, and independent left/right subtrees — are scheduled as
tasks on a shared persistent thread pool (see `thread_pool.hpp`), rather
than spawning a fixed number of threads per node.

## Project layout

```
include/         Public headers for this project (partitioner.hpp, thread_pool.hpp)
src/             Implementation (partitioner.cpp, thread_pool.cpp) and the CLI entry point (main.cpp)
external/        Third-party dependencies (maxflow, cmdparser, generation_checker, status_log)
tests/           Unit tests (doctest)
examples/        Small sample graphs in DIMACS format
plots/           Plotting helper script
```

## Building

```bash
./compile.sh
```

This produces both a Release and a Debug build (`build/` and
`build-debug/`). To drive CMake directly instead:

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

Unit tests are built by default (`inertialflow_tests`) and can be run with
`ctest` or by invoking the test binary directly. Disable them with
`-DINERTIALFLOW_BUILD_TESTS=OFF`.

## Usage

```bash
./build/partitioner -g <graph_file> -c <coord_file> -k <num_cells> [options]
```

| Flag | Long form         | Description                                                                                          | Default                  |
|------|--------------------|-------------------------------------------------------------------------------------------------------|---------------------------|
| `-g` | `--input_graph`    | Input graph (DIMACS `a`/`n` lines, or METIS format with `-m`). **Required.**                          | —                         |
| `-c` | `--input_coord`    | Input vertex coordinates, DIMACS-style `v <id> <lon> <lat> [weight]`. **Required.**                   | —                         |
| `-k` | `--num_cells`      | Number of cells to partition the graph into. **Required.**                                            | —                         |
| `-o` | `--output_file`    | Where to write the resulting partition.                                                               | `dump.txt`                |
| `-f` | `--fraction`       | Fraction of total vertex weight used to grow sources/sinks at each step; must be in `(0, 0.5)`.        | `0.25`                    |
| `-m` | `--metis_format`   | Treat `-g` as a METIS graph file instead of DIMACS. `-c` is still read for coordinates.                | `false`                   |
| `-t` | `--threads`        | Worker threads used to compute the partition.                                                         | hardware concurrency      |
| `-s` | `--show_stats`     | Print partition statistics (cut size, balance, disconnected count) after computing.                    | `false`                   |
| `-l` | `--verbose_log`    | Log one CSV row per recursion step to `stderr` — see below.                                            | `false`                   |

Example, using the small graph under `examples/`:

```bash
./build/partitioner -g examples/example.gr -c examples/example.co -k 4 -s
```

### Detailed step logging (`-l` / `--verbose_log`)

With `-l`, every recursion step (i.e. every node of the bisection tree)
writes one CSV row to `stderr` describing the cut it picked and how long
it took to compute it. `stdout` is left with only the usual progress/stats
output, so the two can be separated with a normal shell redirect:

```bash
./build/partitioner -g examples/example2.gr -c examples/example2.co -k 8 -l \
    2> steps.csv
```

Columns:

| Column                | Meaning                                                                 |
|------------------------|--------------------------------------------------------------------------|
| `level`                | Recursion depth of this step (root is `0`).                             |
| `lo`, `hi`             | Cell-id range `[lo, hi)` this step is bisecting.                        |
| `num_vertices`         | Number of vertices considered at this step.                             |
| `best_projection_deg`  | Angle (degrees) of the projection direction whose cut was chosen.       |
| `best_flow`            | Max-flow value of the chosen cut (i.e. the cut's edge weight).          |
| `time_ms`              | Wall-clock time to evaluate this step's candidate cuts and pick a winner.|

Steps are logged as soon as they complete, so rows from different
recursion levels/subtrees may be interleaved depending on which pool
worker finishes first; sort by `level` (and then by `lo`) if a strict
top-down ordering is needed.

## License

The `external/maxflow` library (Boykov-Kolmogorov max-flow implementation)
is distributed under the GPL; see `external/maxflow/GPL.TXT`.
