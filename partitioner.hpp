#pragma once

#include <maxflow/graph.h>  // From maxflow library

#include <array>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "generation_checker.hpp"

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
  explicit Partitioner(int k);

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

  // recursive_bisect evaluates one candidate cut per projection direction
  // and keeps whichever has the smallest flow. There are exactly four
  // fixed directions (see recursive_bisect), so that's how many candidate
  // cuts are evaluated per call, one per thread.
  static constexpr int kNumProjections = 4;

  int num_cells;
  std::vector<Node> nodes;
  std::vector<std::vector<Edge>> adj;

  // Indices (into nodes/adj) of vertices with at least one incident edge
  // (incoming or outgoing) after cleaning. Populated once by loadGraph and
  // consumed by run() as the working set for recursive bisection: vertices
  // not in this list are fully disconnected and are left at partition id 0
  // without being counted towards cell balance.
  std::vector<int> connected_indices;
  std::size_t disconnected_count = 0;

  // Per-projection max-flow scratch state. Each projection's candidate cut
  // is evaluated on its own thread, and the maxflow library's Graph is
  // stateful (reset()/add_node()/add_edge()/maxflow() all mutate shared
  // internal buffers), so each thread needs its own Graph instance plus
  // its own GenerationChecker/global_to_local scratch space rather than
  // sharing one across threads. Slot 0 is unused (it uses the MaxGraph
  // passed into run()/recursive_bisect by the caller); slots 1..N-1 own a
  // dedicated Graph sized identically to the caller's.
  struct FlowWorkspace {
    std::unique_ptr<MaxGraph> owned_graph;  // null for slot 0
    GenerationChecker<uint32_t> active;
    std::vector<int> global_to_local;
  };
  std::array<FlowWorkspace, kNumProjections> flow_workspaces;

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
  // needed: run() simply starts recursive_bisect from connected_indices
  // instead of the full [0, n) range, and every existing global-to-local
  // mapping (see evaluate_cut) keeps working unchanged.
  void filter_disconnected();

  // Recursively splits node_indices into the partition id range [lo, hi).
  // Using an explicit id range (rather than bit tricks on a power-of-two
  // cell count) means num_cells does not need to be a power of two.
  void recursive_bisect(MaxGraph& graph, std::vector<int>& node_indices, int lo,
                        int hi, const double fraction);

  // Computes the max flow (and, if it beats current_best_flow, the
  // induced left/right vertex sets) for one candidate source/sink split.
  // active and global_to_local are scratch space owned by the caller
  // (see FlowWorkspace) rather than Partitioner members, so that
  // recursive_bisect can run several calls concurrently, each with its
  // own scratch space.
  long long evaluate_cut(MaxGraph& graph, GenerationChecker<uint32_t>& active,
                         std::vector<int>& global_to_local,
                         const std::vector<int>& node_indices,
                         std::span<const int> sources,
                         std::span<const int> sinks, std::vector<int>& out_left,
                         std::vector<int>& out_right,
                         long long current_best_flow);
};
