#pragma once

#include <maxflow/graph.h>  // From maxflow library

#include <array>
#include <atomic>
#include <memory>
#include <span>
#include <string>
#include <thread>
#include <vector>

#include "generation_checker.hpp"
#include "thread_pool.hpp"

typedef maxflow::Graph<int, int, int> MaxGraph;

struct Node {
  int id;
  double lat, lon;
  // Vertex weight, used for balancing instead of raw vertex counts. Defaults
  // to 1, which reproduces the old unweighted behavior. Populated either
  // from an optional 4th column in the DIMACS coordinate file, or from the
  // METIS graph file's own vertex weights when format == kMetis (see
  // loadGraph).
  long long weight = 1;
  int partition_id = 0;
};

struct Edge {
  int to;
  int capacity;
};

// Selects which format loadGraph() expects the graph file (the -g argument)
// to be in. The coordinate file (-c) is always the DIMACS-style `v <id>
// <lon> <lat> [weight]` format described in the README, regardless of
// graphPath's format, since METIS graph files carry no geometric
// coordinates.
enum class GraphFormat {
  kDimacs,  // `a <from> <to> <capacity>` lines, see README.
  kMetis,   // Standard METIS graph format (header line `n m [fmt] [ncon]`
            // followed by one adjacency line per vertex).
};

class Partitioner {
 public:
  // Throws std::invalid_argument if k < 1.
  //
  // num_threads sizes the persistent worker pool used to evaluate cuts in
  // parallel (see the "Persistent worker pool" section below); it defaults
  // to the hardware's reported concurrency (falling back to 1 if that's
  // unknown) and is clamped to at least 1 either way.
  explicit Partitioner(int k, unsigned num_threads =
                                  std::thread::hardware_concurrency());

  // Stops and joins the persistent worker pool.
  ~Partitioner();

  // The persistent worker pool (a ThreadPool member, which itself owns
  // std::thread workers) makes this type inherently non-copyable and
  // non-movable; say so explicitly rather than relying on the implicit
  // deletion, since the reason is no longer obvious just from the member
  // list.
  Partitioner(const Partitioner&) = delete;
  Partitioner& operator=(const Partitioner&) = delete;
  Partitioner(Partitioner&&) = delete;
  Partitioner& operator=(Partitioner&&) = delete;

  // Throws std::runtime_error if either file cannot be opened, or if the
  // graph/coordinate files are malformed (e.g. non-contiguous vertex ids,
  // or edges referencing unknown vertices).
  //
  // The loaded graph is cleaned before use: self loops are dropped, and
  // duplicate (parallel) edges between the same ordered pair of vertices
  // are merged into one, keeping the larger capacity.
  //
  // Vertex weights: for format == kDimacs, an optional 4th column on each
  // coordinate line (`v <id> <lon> <lat> <weight>`) sets that vertex's
  // weight; if omitted, weight defaults to 1. For format == kMetis, vertex
  // weights come from the METIS file itself when its header flags indicate
  // they are present (if the METIS file declares multiple weight
  // constraints, they are summed into one scalar weight); otherwise they
  // fall back to the coordinate file's 4th column (or 1).
  MaxGraph loadGraph(const std::string& graphPath, const std::string& coordPath,
                     GraphFormat format = GraphFormat::kDimacs);

  // Throws std::invalid_argument if fraction is not in (0, 0.5).
  //
  // fraction is a fraction of total *vertex weight* (not vertex count): at
  // each recursion step, sources/sinks are grown from the extremes of the
  // current projection until their accumulated weight reaches
  // fraction * (total weight of the current node set). With unit weights
  // (the default) this reproduces the original count-based behavior.
  void run(MaxGraph& graph, const double fraction);
  void saveResults(const std::string& outputPath);
  std::size_t numVertices() const;
  void printStats() const;

  // Exposed mainly for testing.
  int numCells() const { return num_cells; }
  const std::vector<Node>& getNodes() const { return nodes; }
  const std::vector<std::vector<Edge>>& getAdjacency() const { return adj; }
  std::size_t numDisconnected() const { return disconnected_count; }

 private:
  // Parses the DIMACS-style coordinate file into `nodes`, defaulting every
  // node's weight to 1 (overwritten by loadMetisEdges if the METIS file
  // supplies its own vertex weights). Shared by both graph formats.
  void loadCoordinates(const std::string& coordPath);

  // Parses a DIMACS-style `a <from> <to> <capacity>` graph file into `adj`.
  // Returns the number of edges kept (i.e. after dropping self loops, before
  // clean_adjacency's dedup pass).
  std::size_t loadDimacsEdges(const std::string& graphPath);

  // Parses a METIS-format graph file into `adj`, overwriting `nodes[i].weight`
  // in place when the file's header flags indicate vertex weights are
  // present. METIS conventionally lists each undirected edge on both
  // endpoints' lines, but files listing it only once (from either endpoint)
  // are also accepted: edges are collected keyed by their unordered vertex
  // pair, merging any duplicate occurrence by keeping the larger capacity,
  // then emitted as exactly one directed entry per pair. Returns the number
  // of (deduplicated) edges kept.
  std::size_t loadMetisEdges(const std::string& graphPath);

  // Each recursion node evaluates one candidate cut per projection
  // direction and keeps whichever has the smallest flow. Directions are
  // evenly spaced over 180 degrees (see run_projection_task) -- 180, not
  // 360, because projecting onto a direction and its opposite produces the
  // same ordering of vertices (just reversed), so anything past 180
  // degrees is redundant. Higher values give the root of the recursion
  // (by far the most expensive single node, since it's the whole graph)
  // more independent work to parallelize, and typically also find a
  // somewhat better cut by sampling more candidate directions; the
  // tradeoff is more max-flow evaluations overall. 4 was the original
  // value; 8 is a reasonable default for machines with more than a
  // handful of cores.
  static constexpr int kNumProjections = 8;

  int num_cells;
  std::vector<Node> nodes;
  // Loading-time / reference adjacency: built directly by
  // loadDimacsEdges/loadMetisEdges, deduplicated in place by
  // clean_adjacency(), and consulted by filter_disconnected(), printStats(),
  // and getAdjacency() (tests). Not used by evaluate_cut() itself; see
  // adj_edges/adj_offset below for the copy that the hot path actually
  // walks.
  std::vector<std::vector<Edge>> adj;

  // CSR-style copy of `adj`, built once by build_csr() right after `adj` is
  // finalized (after filter_disconnected() in loadGraph()): edges for
  // vertex u live in adj_edges[adj_offset[u] .. adj_offset[u + 1]).
  // evaluate_cut() is called on the order of thousands of times over one
  // partitioning run and re-walks the neighbor lists of every vertex in its
  // subproblem each time, so having those lists live in one contiguous
  // array (rather than n separate std::vector<Edge> heap allocations
  // scattered across the heap, as `adj` has) keeps that repeated traversal
  // sequential in memory instead of pointer-chasing.
  std::vector<Edge> adj_edges;
  std::vector<int> adj_offset;  // size nodes.size() + 1

  // Builds adj_edges/adj_offset from the current contents of `adj`. Call
  // once, after `adj` is fully loaded and cleaned.
  void build_csr();

  // Indices (into nodes/adj) of vertices with at least one incident edge
  // (incoming or outgoing) after cleaning. Populated once by loadGraph and
  // consumed by run() as the working set for recursive bisection: vertices
  // not in this list are fully disconnected and are left at partition id 0
  // without being counted towards cell balance.
  std::vector<int> connected_indices;
  std::size_t disconnected_count = 0;

  // Per-worker max-flow scratch state. Every candidate cut (one per
  // projection direction, per recursion node -- see run_projection_task)
  // is evaluated on whichever pool worker happens to pick up its task, and
  // the maxflow library's Graph is stateful (reset()/add_node()/
  // add_edge()/maxflow() all mutate shared internal buffers), so each
  // worker needs its own Graph instance plus its own GenerationChecker/
  // global_to_local/local_indices scratch space rather than sharing one
  // across workers. Sized to pool.size(); indexed by worker id, not by
  // projection index, since with a task queue any worker can end up
  // evaluating any projection. Worker 0 is unused (it uses the MaxGraph
  // passed into run() by the caller, see external_graph); workers 1..N-1
  // each own a dedicated Graph sized identically to the caller's.
  struct FlowWorkspace {
    std::unique_ptr<MaxGraph> owned_graph;  // null for worker 0
    GenerationChecker<uint32_t> active;
    std::vector<int> global_to_local;
    // Reusable scratch copy of the current task's node_indices, sorted by
    // this task's projection direction. Reassigned every task (see
    // run_projection_task); kept as a per-worker member instead of a fresh
    // local vector so its buffer's capacity carries over between tasks
    // instead of being freed and reallocated every time.
    std::vector<int> local_indices;
  };
  std::vector<FlowWorkspace> flow_workspaces;

  void init_buffers(size_t n, size_t numEdges);

  // Removes duplicate (parallel) edges from adj, merging each duplicate
  // pair by keeping the larger capacity. Self loops are filtered out
  // earlier, while the file is still being parsed. Returns the number of
  // duplicate edges removed.
  std::size_t clean_adjacency();

  // Finds vertices with zero in-degree and zero out-degree, assigns them
  // partition id 0 directly, and populates connected_indices with every
  // other vertex. Disconnected vertices are identified by their original
  // (dense, 1-based-minus-1) index, so no renumbering of nodes/adj is
  // needed: run() simply starts the bisection from connected_indices
  // instead of the full [0, n) range, and every existing global-to-local
  // mapping (see evaluate_cut) keeps working unchanged.
  void filter_disconnected();

  // ---- Parallel recursive bisection ----
  //
  // The original implementation recursed on the calling thread, spawning
  // kNumProjections=4 threads per recursion node just to evaluate that
  // node's 4 candidate cuts, then joining them before recursing into the
  // left and right halves *sequentially*. That caps parallelism at 4
  // regardless of how many cells (and therefore how much independent
  // subtree work) are available, and wastes it entirely near the leaves,
  // where a whole subtree's cut is often cheap but the recursion is still
  // strictly serial.
  //
  // Instead, both axes of parallelism -- the 4 projections of one node,
  // and independent left/right subtrees -- are expressed as tasks on one
  // shared `pool` (see thread_pool.hpp), sized to `num_threads` workers
  // (typically all available cores) rather than hardcoded to 4:
  //
  //   1. submit_bisect_task(node_indices, lo, hi, fraction) is the unit of
  //      recursion. It handles the base case directly; otherwise it builds
  //      a BisectContext (shared_ptr, so it outlives this call and is
  //      shared by its own 4 projection tasks) and submits 4
  //      run_projection_task calls to `pool` -- one per projection -- then
  //      returns immediately without waiting for them.
  //   2. run_projection_task(worker_id, ctx, p) does exactly what the old
  //      per-projection thread lambda did (sort, pick source/sink,
  //      evaluate_cut), writing into ctx->flows[p]/lefts[p]/rights[p], and
  //      finally atomically decrements ctx->remaining. Whichever of the 4
  //      tasks is the *last* to finish (remaining reaches 0 on its
  //      decrement) picks the best of the 4 flows and calls
  //      submit_bisect_task() again for the left and right halves -- all
  //      still inline, on that same worker, no extra task hop needed.
  //
  // Because child subtree tasks are submitted to the same pool as
  // everything else, an idle worker that finishes early naturally picks up
  // whatever's next in the queue, whether that's another projection for
  // the current recursion level or a subtree several levels deeper -- the
  // pool self-balances instead of pairing work to threads up front.
  //
  // run() kicks off the root call and then calls pool.wait_idle(), which
  // blocks until the entire task graph (including every task any task
  // submitted while it was waiting) has finished.
  struct BisectContext {
    std::vector<int> node_indices;
    int lo, hi;
    double fraction;
    long long target_weight;
    // Starts at kNumProjections; each run_projection_task call decrements
    // it after writing its slot below. The task whose decrement brings it
    // to 0 is guaranteed (via the fetch_sub's acq_rel ordering) to see
    // every other task's writes, and is the one that proceeds to pick the
    // best cut and spawn the child subtree tasks.
    std::atomic<int> remaining{kNumProjections};
    std::array<long long, kNumProjections> flows;
    std::array<std::vector<int>, kNumProjections> lefts, rights;
  };

  unsigned num_threads;
  ThreadPool pool;
  // Set at the start of each run() call; worker 0's FlowWorkspace uses
  // this instead of an owned Graph (see FlowWorkspace). Only ever read
  // while `pool` has outstanding work from the run() call that set it, so
  // there's no risk of it being read after the pointed-to MaxGraph goes
  // away.
  MaxGraph* external_graph = nullptr;

  // Handles the base case directly (assigns partition id `lo` to every
  // vertex in node_indices), or otherwise builds a BisectContext and
  // submits its 4 projection tasks to `pool`. Returns immediately in
  // either case; does not block on the projections it just submitted.
  void submit_bisect_task(std::vector<int> node_indices, int lo, int hi,
                          double fraction);

  // Evaluates the candidate cut for projection p of ctx, exactly like the
  // old per-projection thread lambda did; if it's the last of the 4 to
  // finish, also picks the winning cut and recurses (see BisectContext and
  // submit_bisect_task above).
  void run_projection_task(std::size_t worker_id,
                           std::shared_ptr<BisectContext> ctx, int p);

  // Computes the max flow (and, if it beats current_best_flow, the
  // induced left/right vertex sets) for one candidate source/sink split.
  // active and global_to_local are scratch space owned by the caller
  // (see FlowWorkspace) rather than Partitioner members, so that
  // concurrent workers each use their own scratch space.
  long long evaluate_cut(MaxGraph& graph, GenerationChecker<uint32_t>& active,
                         std::vector<int>& global_to_local,
                         const std::vector<int>& node_indices,
                         std::span<const int> sources,
                         std::span<const int> sinks, std::vector<int>& out_left,
                         std::vector<int>& out_right,
                         long long current_best_flow);
};
