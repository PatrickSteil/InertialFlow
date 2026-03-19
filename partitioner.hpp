#pragma once

#include <maxflow/graph.h> // From maxflow library
#include <string>
#include <vector>

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

  void loadGraph(const std::string &graphPath, const std::string &coordPath);
  void run();
  void saveResults(const std::string &outputPath);
  std::size_t numVertices() const;

private:
  int num_cells;
  std::vector<Node> nodes;
  std::vector<std::vector<Edge>> adj;
  std::vector<int> global_to_local;

  void recursive_bisect(std::vector<int> &node_indices, int current_k);
  long long evaluate_cut(const std::vector<int> &node_indices,
                         const std::vector<int> &sources,
                         const std::vector<int> &sinks,
                         std::vector<int> &out_left,
                         std::vector<int> &out_right, long long current_best_flow);
};
