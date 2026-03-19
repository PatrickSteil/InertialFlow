#include "partitioner.hpp"
#include "status_log.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

typedef maxflow::Graph<int, int, int> MaxGraph;

std::size_t Partitioner::numVertices() const { return nodes.size(); }

void Partitioner::loadGraph(const std::string &graphPath,
                            const std::string &coordPath) {
  StatusLog log("Reading graph and coordinates");
  std::ifstream gFile(graphPath), cFile(coordPath);
  std::string line;

  while (std::getline(cFile, line)) {
    if (line.empty() || line[0] != 'v')
      continue;
    int id;
    double lat, lon;
    if (sscanf(line.c_str(), "v %d %lf %lf", &id, &lon, &lat) == 3) {
      nodes.push_back({id, lat, lon, 0});
      assert((size_t)id == nodes.size());
    }
  }

  adj.resize(nodes.size());
  global_to_local.assign(nodes.size(), -1);

  while (std::getline(gFile, line)) {
    if (line.empty() || line[0] != 'a')
      continue;
    int u, v, cap;
    if (sscanf(line.c_str(), "a %d %d %d", &u, &v, &cap) == 3) {
      assert((size_t)(u - 1) < nodes.size());
      assert((size_t)(v - 1) < nodes.size());
      adj[u - 1].push_back({v - 1, 1});
    }
  }
}

long long Partitioner::evaluate_cut(const std::vector<int> &node_indices,
                                    const std::vector<int> &sources,
                                    const std::vector<int> &sinks,
                                    std::vector<int> &out_left,
                                    std::vector<int> &out_right,
                                    long long current_best_flow) {
  MaxGraph g(node_indices.size(), node_indices.size() * 4);
  g.add_node(node_indices.size());

  for (int i = 0; i < node_indices.size(); ++i) {
    assert(node_indices[i] < global_to_local.size());
    global_to_local[node_indices[i]] = i;
  }

  const int INF = 1e9;
  for (int s : sources) {
    g.add_tweights(global_to_local[s], INF, 0);
  }
  for (int t : sinks) {
    g.add_tweights(global_to_local[t], 0, INF);
  }

  for (int u_global : node_indices) {
    int u_local = global_to_local[u_global];
    for (const auto &e : adj[u_global]) {
      int v_local = global_to_local[e.to];
      assert(u_local < nodes.size());
      assert(v_local < nodes.size());
      if (v_local != -1) {
        g.add_edge(u_local, v_local, e.capacity, e.capacity);
      }
    }
  }

  long long flow = g.maxflow();

  if (flow < current_best_flow) {
    out_left.clear();
    out_right.clear();
    out_left.reserve(node_indices.size() / 2);
    out_right.reserve(node_indices.size() / 2);

    for (int i = 0; i < node_indices.size(); ++i) {
      int g_idx = node_indices[i];
      if (g.what_segment(i) == MaxGraph::SOURCE) {
        out_left.push_back(g_idx);
      } else {
        out_right.push_back(g_idx);
      }
    }
  }

  for (int i = 0; i < node_indices.size(); ++i) {
    global_to_local[node_indices[i]] = -1;
  }
  return flow;
}

void Partitioner::recursive_bisect(std::vector<int> &node_indices,
                                   int k_remaining) {
  if (k_remaining <= 1 || node_indices.size() < 4)
    return;

  long long best_flow = std::numeric_limits<long long>::max();
  std::vector<int> best_left, best_right;

  struct {
    double dx, dy;
  } projs[] = {{1, 0}, {0, 1}, {1, 1}, {1, -1}};

  for (auto &p : projs) {
    std::vector<int> sorted_indices = node_indices;
    std::sort(sorted_indices.begin(), sorted_indices.end(), [&](int a, int b) {
      return (nodes[a].lon * p.dx + nodes[a].lat * p.dy) <
             (nodes[b].lon * p.dx + nodes[b].lat * p.dy);
    });

    int q = sorted_indices.size() / 4;
    std::vector<int> src(sorted_indices.begin(), sorted_indices.begin() + q);
    std::vector<int> snk(sorted_indices.end() - q, sorted_indices.end());

    std::vector<int> L, R;
    long long flow = evaluate_cut(node_indices, src, snk, L, R, best_flow);

    if (flow < best_flow) {
      best_flow = flow;
      best_left = std::move(L);
      best_right = std::move(R);
    }
  }

  int bit = (int)log2(num_cells / k_remaining);
  for (int idx : best_right) {
    nodes[idx].partition_id |= (1 << bit);
  }

  recursive_bisect(best_left, k_remaining / 2);
  recursive_bisect(best_right, k_remaining / 2);
}

void Partitioner::run() {
  StatusLog log("Computing Partition");
  std::vector<int> all_indices(nodes.size());
  for (int i = 0; i < nodes.size(); ++i)
    all_indices[i] = i;

  recursive_bisect(all_indices, num_cells);
}

void Partitioner::saveResults(const std::string &outputPath) {
  StatusLog log("Saving partition to file " + outputPath);
  std::ofstream outFile(outputPath);
  if (!outFile.is_open()) {
    std::cerr << "Error: Could not open output file " << outputPath
              << std::endl;
    return;
  }

  outFile << "N " << nodes.size() << "\n";

  for (const auto &node : nodes) {
    outFile << node.id << " " << node.partition_id << "\n";
  }

  outFile.close();
}
