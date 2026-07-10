#include "partitioner.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <thread>

#include "status_log.hpp"

Partitioner::Partitioner(int k) : num_cells(k) {
  if (k < 1) {
    throw std::invalid_argument("Number of cells k must be >= 1, was " +
                                std::to_string(k));
  }
}

std::size_t Partitioner::numVertices() const { return nodes.size(); }

MaxGraph Partitioner::loadGraph(const std::string& graphPath,
                                const std::string& coordPath) {
  StatusLog log("Reading graph and coordinates");
  std::ifstream cFile(coordPath);
  if (!cFile.is_open()) {
    throw std::runtime_error("Could not open coordinate file: " + coordPath);
  }
  std::ifstream gFile(graphPath);
  if (!gFile.is_open()) {
    throw std::runtime_error("Could not open graph file: " + graphPath);
  }

  std::string line;

  while (std::getline(cFile, line)) {
    if (line.empty() || line[0] != 'v') continue;
    int id;
    double lat, lon;
    if (sscanf(line.c_str(), "v %d %lf %lf", &id, &lon, &lat) == 3) {
      if (static_cast<std::size_t>(id) != nodes.size() + 1) {
        throw std::runtime_error(
            "Malformed coordinate file: vertex ids must be dense and "
            "1-based, expected id " +
            std::to_string(nodes.size() + 1) + " but got " +
            std::to_string(id));
      }
      nodes.push_back({id, lat, lon, 0});
    }
  }

  if (nodes.empty()) {
    throw std::runtime_error("Coordinate file contains no vertices: " +
                             coordPath);
  }

  adj.resize(nodes.size());

  std::size_t numEdges = 0;
  std::size_t selfLoops = 0;

  while (std::getline(gFile, line)) {
    if (line.empty() || line[0] != 'a') continue;
    int u, v, cap;
    if (sscanf(line.c_str(), "a %d %d %d", &u, &v, &cap) == 3) {
      if (static_cast<std::size_t>(u - 1) >= nodes.size() ||
          static_cast<std::size_t>(v - 1) >= nodes.size() || u < 1 || v < 1) {
        throw std::runtime_error("Malformed graph file: edge (" +
                                 std::to_string(u) + ", " + std::to_string(v) +
                                 ") references a vertex id outside [1, " +
                                 std::to_string(nodes.size()) + "]");
      }

      // Self loops never affect a cut (a node is always on the same side
      // as itself), so drop them here rather than carrying them through
      // the rest of the pipeline.
      if (u == v) {
        ++selfLoops;
        continue;
      }

      adj[u - 1].push_back({v - 1, cap});

      ++numEdges;
    }
  }

  std::size_t duplicates = clean_adjacency();
  numEdges -= duplicates;

  if (selfLoops > 0 || duplicates > 0) {
    std::cout << "Cleaned graph: removed " << selfLoops << " self-loop(s) and "
              << duplicates
              << " duplicate edge(s) (kept max capacity per pair)\n";
  }

  filter_disconnected();
  if (disconnected_count > 0) {
    std::cout << "Found " << disconnected_count
              << " disconnected vertex(es) (no incident edges); assigned "
                 "to cell 0 and excluded from balancing\n";
  }

  init_buffers(nodes.size(), numEdges);

  return MaxGraph{(int)nodes.size(), (int)numEdges};
}

std::size_t Partitioner::clean_adjacency() {
  std::size_t duplicates = 0;

  for (auto& edges : adj) {
    if (edges.size() < 2) continue;

    std::sort(edges.begin(), edges.end(),
              [](const Edge& a, const Edge& b) { return a.to < b.to; });

    std::size_t write = 0;
    for (std::size_t read = 0; read < edges.size(); ++read) {
      if (write > 0 && edges[write - 1].to == edges[read].to) {
        // Same (u, v) pair seen again: merge by keeping the larger
        // capacity rather than silently dropping information.
        edges[write - 1].capacity =
            std::max(edges[write - 1].capacity, edges[read].capacity);
        ++duplicates;
      } else {
        edges[write++] = edges[read];
      }
    }
    edges.resize(write);
  }

  return duplicates;
}

void Partitioner::filter_disconnected() {
  std::vector<char> has_edge(nodes.size(), 0);

  for (std::size_t u = 0; u < adj.size(); ++u) {
    if (!adj[u].empty()) has_edge[u] = 1;
    for (const auto& e : adj[u]) {
      has_edge[e.to] = 1;
    }
  }

  connected_indices.clear();
  connected_indices.reserve(nodes.size());
  disconnected_count = 0;

  for (std::size_t i = 0; i < nodes.size(); ++i) {
    if (has_edge[i]) {
      connected_indices.push_back(static_cast<int>(i));
    } else {
      nodes[i].partition_id = 0;
      ++disconnected_count;
    }
  }
}

long long Partitioner::evaluate_cut(
    MaxGraph& graph, GenerationChecker<uint32_t>& active,
    std::vector<int>& global_to_local, const std::vector<int>& node_indices,
    std::span<const int> sources, std::span<const int> sinks,
    std::vector<int>& out_left, std::vector<int>& out_right,
    long long current_best_flow) {
  graph.reset();
  graph.add_node(node_indices.size());

  active.reset();

  for (std::size_t i = 0; i < node_indices.size(); ++i) {
    int g = node_indices[i];
    active.mark(g);
    global_to_local[g] = static_cast<int>(i);
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

    for (const auto& e : adj[u_global]) {
      if (active.isMarked(e.to)) {
        int v_local = global_to_local[e.to];
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

    for (std::size_t i = 0; i < node_indices.size(); ++i) {
      int g_idx = node_indices[i];
      if (graph.what_segment(static_cast<int>(i)) == MaxGraph::SOURCE) {
        out_left.push_back(g_idx);
      } else {
        out_right.push_back(g_idx);
      }
    }
  }
  return flow;
}

void Partitioner::recursive_bisect(MaxGraph& graph,
                                   std::vector<int>& node_indices, int lo,
                                   int hi, const double fraction) {
  const int k_remaining = hi - lo;

  // Base case: only one cell id left for this group of nodes (or too few
  // nodes to meaningfully bisect further), so every node here gets id `lo`.
  // Previously this fell through leaving whatever partition_id an ancestor
  // bit-OR had accumulated, which only produced correct results when
  // num_cells was a power of two; assigning `lo` explicitly here works for
  // any k.
  if (k_remaining <= 1 || node_indices.size() < 4) {
    for (int idx : node_indices) {
      nodes[idx].partition_id = lo;
    }
    return;
  }

  struct {
    double dx, dy;
  } projs[] = {{1, 0}, {0, 1}, {1, 1}, {1, -1}};

  const auto n = node_indices.size();

  // One candidate cut per projection direction, evaluated concurrently.
  // Each projection needs its own copy of node_indices (nth_element
  // reorders in place) and its own MaxGraph + scratch space (see
  // FlowWorkspace) since the maxflow library's Graph is not safe to share
  // across threads. Slot 0 reuses the caller-supplied graph; slots 1..N-1
  // use the Partitioner-owned graphs allocated in init_buffers().
  std::array<MaxGraph*, kNumProjections> graphs = {
      &graph, flow_workspaces[1].owned_graph.get(),
      flow_workspaces[2].owned_graph.get(),
      flow_workspaces[3].owned_graph.get()};

  std::array<long long, kNumProjections> flows;
  flows.fill(std::numeric_limits<long long>::max());
  std::array<std::vector<int>, kNumProjections> lefts, rights;

  std::vector<std::thread> threads;
  threads.reserve(kNumProjections);

  for (int p = 0; p < kNumProjections; ++p) {
    int q = static_cast<int>(n * fraction);
    if (q < 1) {
      // Too few nodes to pick a non-empty source/sink set for this
      // fraction; skip this projection rather than evaluating a
      // degenerate (empty source or sink) cut.
      continue;
    }

    threads.emplace_back([&, p, q]() {
      const auto& proj = projs[p];
      auto cmp = [&](int a, int b) {
        return (nodes[a].lon * proj.dx + nodes[a].lat * proj.dy) <
               (nodes[b].lon * proj.dx + nodes[b].lat * proj.dy);
      };

      // Private copy: nth_element must not mutate node_indices, which is
      // shared by every projection's thread.
      std::vector<int> local_indices = node_indices;

      std::nth_element(local_indices.begin(), local_indices.begin() + q,
                       local_indices.end(), cmp);
      std::nth_element(local_indices.begin() + q,
                       local_indices.begin() + (n - q), local_indices.end(),
                       cmp);

      std::span<const int> src(local_indices.begin(), q);
      std::span<const int> snk(local_indices.end() - q, q);

      flows[p] = evaluate_cut(*graphs[p], flow_workspaces[p].active,
                              flow_workspaces[p].global_to_local, local_indices,
                              src, snk, lefts[p], rights[p],
                              std::numeric_limits<long long>::max());
    });
  }

  for (auto& t : threads) t.join();

  long long best_flow = std::numeric_limits<long long>::max();
  int best_p = -1;
  for (int p = 0; p < kNumProjections; ++p) {
    if (flows[p] < best_flow) {
      best_flow = flows[p];
      best_p = p;
    }
  }

  if (best_p < 0) {
    // Every projection was skipped (graph too small for this fraction);
    // fall back to the base case behavior.
    for (int idx : node_indices) {
      nodes[idx].partition_id = lo;
    }
    return;
  }

  std::vector<int>& best_left = lefts[best_p];
  std::vector<int>& best_right = rights[best_p];

  // Split the id range proportionally to how many cells each side should
  // get, rather than assuming an even power-of-two split.
  const int mid = lo + k_remaining / 2;

  recursive_bisect(graph, best_left, lo, mid, fraction);
  recursive_bisect(graph, best_right, mid, hi, fraction);
}

void Partitioner::run(MaxGraph& graph, const double fraction) {
  if (!(fraction > 0.0 && fraction < 0.5)) {
    throw std::invalid_argument(
        "Given fraction should be between (0, 0.5), was " +
        std::to_string(fraction));
  }

  StatusLog log("Computing Partition");
  // Only connected vertices participate in bisection/balancing; vertices
  // with no incident edges were already assigned to cell 0 by
  // filter_disconnected() in loadGraph().
  std::vector<int> all_indices = connected_indices;

  recursive_bisect(graph, all_indices, 0, num_cells, fraction);
}

void Partitioner::saveResults(const std::string& outputPath) {
  StatusLog log("Saving partition to file " + outputPath);
  std::ofstream outFile(outputPath);
  if (!outFile.is_open()) {
    std::cerr << "Error: Could not open output file " << outputPath
              << std::endl;
    return;
  }

  outFile << "N " << nodes.size() << "\n";

  for (const auto& node : nodes) {
    outFile << node.id << " " << node.partition_id << "\n";
  }

  outFile.close();
}

void Partitioner::printStats() const {
  if (nodes.empty()) return;

  std::vector<int> cell_sizes(num_cells, 0);
  for (const auto& node : nodes) {
    if (node.partition_id < num_cells) {
      cell_sizes[node.partition_id]++;
    }
  }

  long long total_cut_size = 0;
  for (std::size_t u = 0; u < nodes.size(); ++u) {
    for (const auto& edge : adj[u]) {
      int v = edge.to;
      total_cut_size += (nodes[u].partition_id != nodes[v].partition_id);
    }
  }

  int min_size = std::numeric_limits<int>::max();
  int max_size = 0;
  double avg_size = static_cast<double>(nodes.size()) / num_cells;

  for (int size : cell_sizes) {
    if (size < min_size) min_size = size;
    if (size > max_size) max_size = size;
  }

  double imbalance = (max_size / avg_size) - 1.0;

  std::cout << std::string(30, '=') << "\n";
  std::cout << "   PARTITION STATISTICS\n";
  std::cout << std::string(30, '=') << "\n";
  std::cout << "Total Nodes:      " << nodes.size() << "\n";
  std::cout << "Total Cells (k):  " << num_cells << "\n";
  std::cout << "Disconnected:     " << disconnected_count
            << " (assigned to cell 0)\n";
  std::cout << "Total Cut Edges:  " << total_cut_size << "\n";
  std::cout << "------------------------------\n";
  std::cout << "Max Cell Size:    " << max_size << "\n";
  std::cout << "Min Cell Size:    " << min_size << "\n";
  std::cout << "Avg Cell Size:    " << avg_size << "\n";
  printf("Imbalance:        %.2f%%\n", imbalance * 100.0);
  std::cout << "------------------------------\n";
}

void Partitioner::init_buffers(size_t n, size_t numEdges) {
  for (int p = 0; p < kNumProjections; ++p) {
    auto& ws = flow_workspaces[p];
    ws.active.resize(n);
    ws.global_to_local.resize(n);
    // Slot 0 uses the MaxGraph the caller passes into run()/
    // recursive_bisect; slots 1..N-1 need their own, sized the same way
    // the caller's was in loadGraph().
    if (p != 0) {
      ws.owned_graph = std::make_unique<MaxGraph>(static_cast<int>(n),
                                                  static_cast<int>(numEdges));
    }
  }
}
