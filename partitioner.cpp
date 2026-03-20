#include "partitioner.hpp"
#include "status_log.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <span>

std::size_t Partitioner::numVertices() const { return nodes.size(); }

MaxGraph Partitioner::loadGraph(const std::string &graphPath,
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

  active.resize(nodes.size());
  global_to_local.resize(nodes.size());

  std::size_t numEdges = 0;

  while (std::getline(gFile, line)) {
    if (line.empty() || line[0] != 'a')
      continue;
    int u, v, cap;
    if (sscanf(line.c_str(), "a %d %d %d", &u, &v, &cap) == 3) {
      assert((size_t)(u - 1) < nodes.size());
      assert((size_t)(v - 1) < nodes.size());
      adj[u - 1].push_back({v - 1, 1});

      ++numEdges;
    }
  }

  init_buffers(nodes.size());

  return MaxGraph{(int)nodes.size(), (int)numEdges};
}

long long Partitioner::evaluate_cut(MaxGraph &graph,
                                    const std::vector<int> &node_indices,
                                    std::span<const int> sources,
                                    std::span<const int> sinks,
                                    std::vector<int> &out_left,
                                    std::vector<int> &out_right,
                                    long long current_best_flow) {
  graph.reset();
  graph.add_node(node_indices.size());

  active.reset();

  for (int i = 0; i < node_indices.size(); ++i) {
    int g = node_indices[i];
    active.mark(g);
    global_to_local[g] = i;
  }

  const int INF = 1e9;
  for (int s : sources) {
    graph.add_tweights(global_to_local[s], INF, 0);
  }
  for (int t : sinks) {
    graph.add_tweights(global_to_local[t], 0, INF);
  }

  for (int u_global : node_indices) {
    int u_local = global_to_local[u_global];
    assert(u_local < nodes.size());

    for (const auto &e : adj[u_global]) {
      if (active.isMarked(e.to)) {
        int v_local = global_to_local[e.to];
        assert(v_local < nodes.size());

        graph.add_edge(u_local, v_local, e.capacity, e.capacity);
      }
    }
  }

  long long flow = graph.maxflow();

  if (flow < current_best_flow) {
    out_left.clear();
    out_right.clear();
    out_left.reserve(node_indices.size() / 2);
    out_right.reserve(node_indices.size() / 2);

    for (int i = 0; i < node_indices.size(); ++i) {
      int g_idx = node_indices[i];
      if (graph.what_segment(i) == MaxGraph::SOURCE) {
        out_left.push_back(g_idx);
      } else {
        out_right.push_back(g_idx);
      }
    }
  }
  return flow;
}

void Partitioner::recursive_bisect(MaxGraph &graph,
                                   std::vector<int> &node_indices,
                                   int k_remaining,
                                   const double fraction = 0.25) {
  if (k_remaining <= 1 || node_indices.size() < 4)
    return;

  long long best_flow = std::numeric_limits<long long>::max();
  std::vector<int> best_left, best_right;

  struct {
    double dx, dy;
  } projs[] = {{1, 0}, {0, 1}, {1, 1}, {1, -1}};

  for (auto &p : projs) {
    buf_sorted = node_indices;

    std::sort(buf_sorted.begin(), buf_sorted.end(), [&](int a, int b) {
      return (nodes[a].lon * p.dx + nodes[a].lat * p.dy) <
             (nodes[b].lon * p.dx + nodes[b].lat * p.dy);
    });

    int q = (int)(buf_sorted.size() * fraction);

    std::span<const int> src(buf_sorted.begin(), q);
    std::span<const int> snk(buf_sorted.end() - q, q);

    buf_L.clear();
    buf_R.clear();

    long long flow =
        evaluate_cut(graph, node_indices, src, snk, buf_L, buf_R, best_flow);

    if (flow < best_flow) {
      best_flow = flow;
      best_left = buf_L;
      best_right = buf_R;
    }
  }

  int bit = (int)log2(num_cells / k_remaining);
  for (int idx : best_right) {
    nodes[idx].partition_id |= (1 << bit);
  }

  recursive_bisect(graph, best_left, k_remaining / 2);
  recursive_bisect(graph, best_right, k_remaining / 2);
}

void Partitioner::run(MaxGraph &graph, const double fraction = 0.25) {
  StatusLog log("Computing Partition");
  std::vector<int> all_indices(nodes.size());
  for (int i = 0; i < nodes.size(); ++i)
    all_indices[i] = i;

  recursive_bisect(graph, all_indices, num_cells, fraction);
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

void Partitioner::printStats() const {
  if (nodes.empty())
    return;

  std::vector<int> cell_sizes(num_cells, 0);
  for (const auto &node : nodes) {
    if (node.partition_id < num_cells) {
      cell_sizes[node.partition_id]++;
    }
  }

  long long total_cut_size = 0;
  for (int u = 0; u < nodes.size(); ++u) {
    for (const auto &edge : adj[u]) {
      int v = edge.to;
      total_cut_size += (nodes[u].partition_id != nodes[v].partition_id);
    }
  }

  int min_size = std::numeric_limits<int>::max();
  int max_size = 0;
  double avg_size = static_cast<double>(nodes.size()) / num_cells;

  for (int size : cell_sizes) {
    if (size < min_size)
      min_size = size;
    if (size > max_size)
      max_size = size;
  }

  double imbalance = (max_size / avg_size) - 1.0;

  std::cout << std::string(30, '=') << "\n";
  std::cout << "   PARTITION STATISTICS\n";
  std::cout << std::string(30, '=') << "\n";
  std::cout << "Total Nodes:      " << nodes.size() << "\n";
  std::cout << "Total Cells (k):  " << num_cells << "\n";
  std::cout << "Total Cut Edges:  " << total_cut_size << "\n";
  std::cout << "------------------------------\n";
  std::cout << "Max Cell Size:    " << max_size << "\n";
  std::cout << "Min Cell Size:    " << min_size << "\n";
  std::cout << "Avg Cell Size:    " << avg_size << "\n";
  printf("Imbalance:        %.2f%%\n", imbalance * 100.0);
  std::cout << "------------------------------\n";
}

void Partitioner::init_buffers(size_t n) {
  buf_sorted.reserve(n);
  buf_src.reserve(n);
  buf_snk.reserve(n);
  buf_L.reserve(n);
  buf_R.reserve(n);
}
