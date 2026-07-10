# InertialFlow

```bash
./compile.sh
```

Or drive CMake directly:

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Usage

```
./partitioner -g <graph.gr> -c <coords.co> -k <num_cells> [-f fraction] [-o output_file] [-s]
```

Example, using the graphs in `examples/`:

```bash
./build/partitioner -g examples/example2.gr -c examples/example2.co -k 2 -s -o /tmp/partition.txt
```

### Input format

Both input files use a DIMACS-style, whitespace-separated line format.
Vertex ids are 1-based and must be dense (i.e. the coordinate file must
list ids `1..n` in order, with no gaps).

`coords.co` — one `v <id> <lon> <lat>` line per vertex:

```
p aux sp 4
v 1 -73.98 40.75
v 2 -73.99 40.75
v 3 -73.98 40.76
v 4 -73.99 40.76
```

`graph.gr` — one `a <from> <to> <capacity>` line per (directed) edge; the
partitioner treats edges as undirected with the given capacity, so a
single `a u v cap` line is enough to connect `u` and `v` (no need to also
list `a v u cap`):

```
p sp 4 4
a 1 2 10
a 2 4 5
a 1 3 5
a 3 4 10
```

Lines starting with `c` are comments and are ignored, as is the `p` header
line.

### Output format

The output file lists, for every vertex, which of the `k` cells (numbered
`0..k-1`) it was assigned to:

```
N <num_vertices>
<vertex_id> <partition_id>
<vertex_id> <partition_id>
...
```

## License

The `external/maxflow` library (Boykov-Kolmogorov max-flow implementation)
is distributed under the GPL; see `external/maxflow/GPL.TXT`. 
