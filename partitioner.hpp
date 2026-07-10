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
  int partition_id = 0;
};

struct Edge {
  int to;
  int capacity;
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
  MaxGraph loadGraph(const std::string& graphPath,
                     const std::string& coordPath);

  // Throws std::invalid_argument if fraction is not in (0, 0.5).
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
