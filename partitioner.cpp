#include "partitioner.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <span>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>

#include "status_log.hpp"

Partitioner::Partitioner(int k) : num_cells(k) {
  if (k < 1) {
    throw std::invalid_argument("Number of cells k must be >= 1, was " +
                                std::to_string(k));
  }
  start_pool();
}

Partitioner::~Partitioner() { stop_pool(); }

void Partitioner::start_pool() {
  for (int p = 0; p < kNumProjections; ++p) {
    workers[p] = std::thread([this, p]() { worker_loop(p); });
  }
}

void Partitioner::stop_pool() {
  pool_stop = true;
  start_barrier.arrive_and_wait();
  for (auto& t : workers) {
    if (t.joinable()) t.join();
  }
}

void Partitioner::worker_loop(int p) {
  while (true) {
    start_barrier.arrive_and_wait();
    if (pool_stop) return;
    run_projection(p);
    done_barrier.arrive_and_wait();
  }
}

std::size_t Partitioner::numVertices() const { return nodes.size(); }

MaxGraph Partitioner::loadGraph(const std::string& graphPath,
                                const std::string& coordPath,
                                GraphFormat format) {
  StatusLog log("Reading graph and coordinates");

  loadCoordinates(coordPath);

  std::size_t numEdges = (format == GraphFormat::kMetis)
                             ? loadMetisEdges(graphPath)
                             : loadDimacsEdges(graphPath);

  std::size_t duplicates = clean_adjacency();
  numEdges -= duplicates;

  if (duplicates > 0) {
    std::cout << "Cleaned graph: removed " << duplicates
              << " duplicate edge(s) (kept max capacity per pair)\n";
  }

  filter_disconnected();
  build_csr();
  init_buffers(nodes.size(), numEdges);

  return MaxGraph{(int)nodes.size(), (int)numEdges};
}

void Partitioner::loadCoordinates(const std::string& coordPath) {
  std::ifstream cFile(coordPath);
  if (!cFile.is_open()) {
    throw std::runtime_error("Could not open coordinate file: " + coordPath);
  }

  std::string line;
  while (std::getline(cFile, line)) {
    if (line.empty() || line[0] != 'v') continue;
    int id;
    double lat, lon;
    long long weight = 1;
    int matched =
        sscanf(line.c_str(), "v %d %lf %lf %lld", &id, &lon, &lat, &weight);
    if (matched >= 3) {
      if (static_cast<std::size_t>(id) != nodes.size() + 1) {
        throw std::runtime_error(
            "Malformed coordinate file: vertex ids must be dense and "
            "1-based, expected id " +
            std::to_string(nodes.size() + 1) + " but got " +
            std::to_string(id));
      }
      nodes.push_back({id, lat, lon, weight, 0});
    }
  }

  if (nodes.empty()) {
    throw std::runtime_error("Coordinate file contains no vertices: " +
                             coordPath);
  }

  adj.resize(nodes.size());
}

std::size_t Partitioner::loadDimacsEdges(const std::string& graphPath) {
  std::ifstream gFile(graphPath);
  if (!gFile.is_open()) {
    throw std::runtime_error("Could not open graph file: " + graphPath);
  }

  std::string line;
  std::size_t numEdges = 0;
  std::size_t selfLoops = 0;

  while (std::getline(gFile, line)) {
    if (line.empty()) continue;

    if (line[0] == 'n') {
      int id;
      long long weight;
      if (sscanf(line.c_str(), "n %d %lld", &id, &weight) == 2) {
        if (static_cast<std::size_t>(id - 1) >= nodes.size() || id < 1) {
          throw std::runtime_error(
              "Malformed graph file: vertex weight line references vertex "
              "id " +
              std::to_string(id) + " outside [1, " +
              std::to_string(nodes.size()) + "]");
        }
        nodes[id - 1].weight = weight;
      }
      continue;
    }

    if (line[0] != 'a') continue;
    int u, v, cap;
    if (sscanf(line.c_str(), "a %d %d %d", &u, &v, &cap) == 3) {
      if (static_cast<std::size_t>(u - 1) >= nodes.size() ||
          static_cast<std::size_t>(v - 1) >= nodes.size() || u < 1 || v < 1) {
        throw std::runtime_error("Malformed graph file: edge (" +
                                 std::to_string(u) + ", " + std::to_string(v) +
                                 ") references a vertex id outside [1, " +
                                 std::to_string(nodes.size()) + "]");
      }

      if (u == v) {
        ++selfLoops;
        continue;
      }

      adj[u - 1].push_back({v - 1, cap});

      ++numEdges;
    }
  }

  if (selfLoops > 0) {
    std::cout << "Dropped " << selfLoops << " self-loop(s) from graph file\n";
  }

  return numEdges;
}

std::size_t Partitioner::loadMetisEdges(const std::string& graphPath) {
  std::ifstream gFile(graphPath);
  if (!gFile.is_open()) {
    throw std::runtime_error("Could not open graph file: " + graphPath);
  }

  std::string line;
  auto nextRealLine = [&]() -> bool {
    while (std::getline(gFile, line)) {
      std::size_t start = line.find_first_not_of(" \t\r");
      if (start == std::string::npos) continue;
      if (line[start] == '%') continue;
      return true;
    }
    return false;
  };

  if (!nextRealLine()) {
    throw std::runtime_error("METIS graph file has no header line: " +
                             graphPath);
  }

  long long declaredN = 0, declaredM = 0, fmt = 0, ncon = 1;
  std::istringstream header(line);
  if (!(header >> declaredN >> declaredM)) {
    throw std::runtime_error("Malformed METIS header in file: " + graphPath);
  }
  bool hasFmt = static_cast<bool>(header >> fmt);
  bool hasNcon = hasFmt && static_cast<bool>(header >> ncon);
  if (hasFmt && !hasNcon) ncon = 1;

  if (static_cast<std::size_t>(declaredN) != nodes.size()) {
    throw std::runtime_error("METIS header vertex count (" +
                             std::to_string(declaredN) +
                             ") does not match coordinate file vertex count (" +
                             std::to_string(nodes.size()) + ")");
  }

  const bool hasVSizes = (fmt / 100) % 10 == 1;
  const bool hasVWeights = (fmt / 10) % 10 == 1;
  const bool hasEdgeWeights = fmt % 10 == 1;
  if (ncon < 1) ncon = 1;

  std::size_t selfLoops = 0;

  // METIS conventionally lists every undirected edge twice (once on each
  // endpoint's line), but well-formed files that instead list it only once
  // (from either endpoint) are also valid input. Rather than assume one
  // convention, collect edges keyed by their unordered vertex pair and
  // merge duplicates (taking the larger capacity, matching clean_adjacency's
  // convention), then emit exactly one directed entry per pair afterwards.
  std::unordered_map<std::uint64_t, int> edgeCapByPair;

  for (std::size_t u1 = 1; u1 <= nodes.size(); ++u1) {
    if (!nextRealLine()) {
      throw std::runtime_error("METIS graph file ends before the declared " +
                               std::to_string(nodes.size()) +
                               " vertex line(s) were read");
    }
    std::istringstream tok(line);

    if (hasVSizes) {
      long long vsize;
      if (!(tok >> vsize)) {
        throw std::runtime_error("Malformed METIS vertex line " +
                                 std::to_string(u1) +
                                 " (missing declared vertex size)");
      }
    }

    if (hasVWeights) {
      long long summedWeight = 0;
      for (long long c = 0; c < ncon; ++c) {
        long long w;
        if (!(tok >> w)) {
          throw std::runtime_error(
              "Malformed METIS vertex line " + std::to_string(u1) +
              " (missing one or more of the " + std::to_string(ncon) +
              " declared vertex weight(s))");
        }
        summedWeight += w;
      }
      // Overwrites whatever loadCoordinates() defaulted/parsed: the METIS
      // file's own weights take precedence when it declares any.
      nodes[u1 - 1].weight = summedWeight;
    }

    long long v1;
    while (tok >> v1) {
      long long cap = 1;
      if (hasEdgeWeights && !(tok >> cap)) {
        throw std::runtime_error("Malformed METIS vertex line " +
                                 std::to_string(u1) + " (neighbor " +
                                 std::to_string(v1) +
                                 " is missing its declared edge weight)");
      }

      if (v1 < 1 || static_cast<std::size_t>(v1) > nodes.size()) {
        throw std::runtime_error("Malformed METIS graph file: vertex " +
                                 std::to_string(u1) + " references neighbor " +
                                 std::to_string(v1) + " outside [1, " +
                                 std::to_string(nodes.size()) + "]");
      }

      if (v1 == static_cast<long long>(u1)) {
        ++selfLoops;
        continue;
      }

      std::uint64_t lo =
          static_cast<std::uint64_t>(std::min<long long>(u1, v1)) - 1;
      std::uint64_t hi =
          static_cast<std::uint64_t>(std::max<long long>(u1, v1)) - 1;
      std::uint64_t key = (lo << 32) | hi;

      auto it = edgeCapByPair.find(key);
      if (it == edgeCapByPair.end()) {
        edgeCapByPair.emplace(key, static_cast<int>(cap));
      } else {
        it->second = std::max(it->second, static_cast<int>(cap));
      }
    }
  }

  if (selfLoops > 0) {
    std::cout << "Dropped " << selfLoops << " self-loop(s) from graph file\n";
  }

  for (const auto& [key, cap] : edgeCapByPair) {
    int lo = static_cast<int>(key >> 32);
    int hi = static_cast<int>(key & 0xffffffffu);
    adj[lo].push_back({hi, cap});
  }

  return edgeCapByPair.size();
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

void Partitioner::build_csr() {
  adj_offset.assign(nodes.size() + 1, 0);
  std::size_t total = 0;
  for (std::size_t u = 0; u < nodes.size(); ++u) {
    adj_offset[u] = static_cast<int>(total);
    total += adj[u].size();
  }
  adj_offset[nodes.size()] = static_cast<int>(total);

  adj_edges.resize(total);
  for (std::size_t u = 0; u < nodes.size(); ++u) {
    std::copy(adj[u].begin(), adj[u].end(), adj_edges.begin() + adj_offset[u]);
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

    for (int ei = adj_offset[u_global]; ei < adj_offset[u_global + 1]; ++ei) {
      const Edge& e = adj_edges[ei];
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

  // One candidate cut per projection direction, evaluated concurrently by
  // the persistent worker pool (see run_projection). Each worker needs its
  // own MaxGraph + scratch space (see FlowWorkspace) since the maxflow
  // library's Graph is not safe to share across threads. Slot 0 reuses the
  // caller-supplied graph; slots 1..N-1 use the Partitioner-owned graphs
  // allocated in init_buffers().
  std::array<MaxGraph*, kNumProjections> graphs = {
      &graph, flow_workspaces[1].owned_graph.get(),
      flow_workspaces[2].owned_graph.get(),
      flow_workspaces[3].owned_graph.get()};

  std::array<long long, kNumProjections> flows;
  flows.fill(std::numeric_limits<long long>::max());
  std::array<std::vector<int>, kNumProjections> lefts, rights;

  // Total vertex weight of this node set, used to turn `fraction` into a
  // weight target rather than a vertex-count target. Same across all
  // projections (order-independent sum), so computed once here.
  long long total_weight = 0;
  for (int idx : node_indices) total_weight += nodes[idx].weight;
  const long long target_weight =
      static_cast<long long>(total_weight * fraction);

  // Dispatch one round to the persistent worker pool (see the "Persistent
  // worker pool" comment in partitioner.hpp) instead of spawning and
  // joining kNumProjections std::threads here. round_data is only safe to
  // write at this point because every worker is currently parked at
  // start_barrier from the end of the *previous* round (or pool startup,
  // for the very first round); arrive_and_wait() below is what lets them
  // see it.
  round_data.node_indices = &node_indices;
  round_data.target_weight = target_weight;
  round_data.graphs = &graphs;
  round_data.flows = &flows;
  round_data.lefts = &lefts;
  round_data.rights = &rights;

  start_barrier.arrive_and_wait();
  done_barrier.arrive_and_wait();

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

void Partitioner::run_projection(int p) {
  const std::vector<int>& node_indices = *round_data.node_indices;
  const long long target_weight = round_data.target_weight;
  auto& graphs = *round_data.graphs;
  auto& flows = *round_data.flows;
  auto& lefts = *round_data.lefts;
  auto& rights = *round_data.rights;

  struct {
    double dx, dy;
  } projs[] = {{1, 0}, {0, 1}, {1, 1}, {1, -1}};
  const auto& proj = projs[p];
  const auto n = node_indices.size();

  auto cmp = [&](int a, int b) {
    return (nodes[a].lon * proj.dx + nodes[a].lat * proj.dy) <
           (nodes[b].lon * proj.dx + nodes[b].lat * proj.dy);
  };

  std::vector<int>& local_indices = flow_workspaces[p].local_indices;
  local_indices = node_indices;
  std::sort(local_indices.begin(), local_indices.end(), cmp);

  std::size_t lo = 0;
  long long lo_weight = 0;
  while (lo < n && lo_weight < target_weight) {
    lo_weight += nodes[local_indices[lo]].weight;
    ++lo;
  }

  std::size_t hi = n;
  std::size_t sink_count = 0;
  long long hi_weight = 0;
  while (hi > lo && hi_weight < target_weight) {
    --hi;
    hi_weight += nodes[local_indices[hi]].weight;
    ++sink_count;
  }

  if (lo < 1 || sink_count < 1) {
    return;
  }

  std::span<const int> src(local_indices.data(), lo);
  std::span<const int> snk(local_indices.data() + hi, sink_count);

  flows[p] =
      evaluate_cut(*graphs[p], flow_workspaces[p].active,
                   flow_workspaces[p].global_to_local, local_indices, src, snk,
                   lefts[p], rights[p], std::numeric_limits<long long>::max());
}

void Partitioner::run(MaxGraph& graph, const double fraction) {
  if (!(fraction > 0.0 && fraction < 0.5)) {
    throw std::invalid_argument(
        "Given fraction should be between (0, 0.5), was " +
        std::to_string(fraction));
  }

  StatusLog log("Computing Partition");
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
  std::vector<long long> cell_weights(num_cells, 0);
  long long total_weight = 0;
  for (const auto& node : nodes) {
    total_weight += node.weight;
    if (node.partition_id < num_cells) {
      cell_sizes[node.partition_id]++;
      cell_weights[node.partition_id] += node.weight;
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

  long long min_weight = std::numeric_limits<long long>::max();
  long long max_weight = 0;
  double avg_weight = static_cast<double>(total_weight) / num_cells;
  for (long long w : cell_weights) {
    if (w < min_weight) min_weight = w;
    if (w > max_weight) max_weight = w;
  }
  double weight_imbalance =
      avg_weight > 0 ? (max_weight / avg_weight) - 1.0 : 0.0;

  bool anyWeighted = (total_weight != static_cast<long long>(nodes.size()));

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
  if (anyWeighted) {
    std::cout << "------------------------------\n";
    std::cout << "Total Weight:     " << total_weight << "\n";
    std::cout << "Max Cell Weight:  " << max_weight << "\n";
    std::cout << "Min Cell Weight:  " << min_weight << "\n";
    std::cout << "Avg Cell Weight:  " << avg_weight << "\n";
    printf("Weight Imbalance: %.2f%%\n", weight_imbalance * 100.0);
  }
  std::cout << "------------------------------\n";
}

void Partitioner::init_buffers(size_t n, size_t numEdges) {
  for (int p = 0; p < kNumProjections; ++p) {
    auto& ws = flow_workspaces[p];
    ws.active.resize(n);
    ws.global_to_local.resize(n);
    if (p != 0) {
      ws.owned_graph = std::make_unique<MaxGraph>(static_cast<int>(n),
                                                  static_cast<int>(numEdges));
    }
  }
}
