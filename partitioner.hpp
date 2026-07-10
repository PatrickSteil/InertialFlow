#pragma once

#include <maxflow/graph.h>  // From maxflow library

#include <array>
#include <barrier>
#include <memory>
#include <span>
#include <string>
#include <thread>
#include <vector>

#include "generation_checker.hpp"

typedef maxflow::Graph<int, int, int> MaxGraph;

struct Node {
  int id;
  double lat, lon;
  long long weight = 1;
  int partition_id = 0;
};

struct Edge {
  int to;
  int capacity;
};

enum class GraphFormat {
  kDimacs,
  kMetis,
};

class Partitioner {
 public:
  explicit Partitioner(int k);

  ~Partitioner();

  Partitioner(const Partitioner&) = delete;
  Partitioner& operator=(const Partitioner&) = delete;
  Partitioner(Partitioner&&) = delete;
  Partitioner& operator=(Partitioner&&) = delete;

  MaxGraph loadGraph(const std::string& graphPath, const std::string& coordPath,
                     GraphFormat format = GraphFormat::kDimacs);

  void run(MaxGraph& graph, const double fraction);
  void saveResults(const std::string& outputPath);
  std::size_t numVertices() const;
  void printStats() const;

  int numCells() const { return num_cells; }
  const std::vector<Node>& getNodes() const { return nodes; }
  const std::vector<std::vector<Edge>>& getAdjacency() const { return adj; }
  std::size_t numDisconnected() const { return disconnected_count; }

 private:
  void loadCoordinates(const std::string& coordPath);

  std::size_t loadDimacsEdges(const std::string& graphPath);

  std::size_t loadMetisEdges(const std::string& graphPath);

  static constexpr int kNumProjections = 4;

  int num_cells;
  std::vector<Node> nodes;
  std::vector<std::vector<Edge>> adj;

  std::vector<Edge> adj_edges;
  std::vector<int> adj_offset;

  void build_csr();

  std::vector<int> connected_indices;
  std::size_t disconnected_count = 0;

  struct FlowWorkspace {
    std::unique_ptr<MaxGraph> owned_graph;
    GenerationChecker<uint32_t> active;
    std::vector<int> global_to_local;
    std::vector<int> local_indices;
  };
  std::array<FlowWorkspace, kNumProjections> flow_workspaces;

  void init_buffers(size_t n, size_t numEdges);
  std::size_t clean_adjacency();

  void filter_disconnected();

  void recursive_bisect(MaxGraph& graph, std::vector<int>& node_indices, int lo,
                        int hi, const double fraction);

  struct RoundData {
    const std::vector<int>* node_indices = nullptr;
    long long target_weight = 0;
    std::array<MaxGraph*, kNumProjections>* graphs = nullptr;
    std::array<long long, kNumProjections>* flows = nullptr;
    std::array<std::vector<int>, kNumProjections>* lefts = nullptr;
    std::array<std::vector<int>, kNumProjections>* rights = nullptr;
  };
  RoundData round_data;
  bool pool_stop = false;
  std::barrier<> start_barrier{kNumProjections + 1};
  std::barrier<> done_barrier{kNumProjections + 1};
  std::array<std::thread, kNumProjections> workers;

  void start_pool();
  void stop_pool();
  void worker_loop(int p);

  void run_projection(int p);

  long long evaluate_cut(MaxGraph& graph, GenerationChecker<uint32_t>& active,
                         std::vector<int>& global_to_local,
                         const std::vector<int>& node_indices,
                         std::span<const int> sources,
                         std::span<const int> sinks, std::vector<int>& out_left,
                         std::vector<int>& out_right,
                         long long current_best_flow);
};
