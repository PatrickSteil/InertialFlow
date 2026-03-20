#pragma once

#include "generation_checker.hpp"

#include <maxflow/graph.h> // From maxflow library
#include <span>
#include <string>
#include <vector>

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
  Partitioner(int k) : num_cells(k) {}

  MaxGraph loadGraph(const std::string &graphPath,
                     const std::string &coordPath);
  void run(MaxGraph &graph, const double fraction);
  void saveResults(const std::string &outputPath);
  std::size_t numVertices() const;
  void printStats() const;

private:
  int num_cells;
  std::vector<Node> nodes;
  std::vector<std::vector<Edge>> adj;

  std::vector<int> buf_src;
  std::vector<int> buf_snk;
  std::vector<int> buf_L;
  std::vector<int> buf_R;

  GenerationChecker<uint32_t> active;
  std::vector<int> global_to_local;

  void init_buffers(size_t n);
  void recursive_bisect(MaxGraph &graph, std::vector<int> &node_indices,
                        int current_k, const double fraction);
  long long evaluate_cut(MaxGraph &graph, const std::vector<int> &node_indices,
                         std::span<const int> sources,
                         std::span<const int> sinks, std::vector<int> &out_left,
                         std::vector<int> &out_right,
                         long long current_best_flow);
};
