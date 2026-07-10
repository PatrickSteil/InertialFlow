#include "external/doctest.h"
#include "partitioner.hpp"

#include <atomic>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>

namespace {

// Writes `content` to a fresh temp file with the given suffix and returns
// its path. Files are created under the system temp directory and are
// intentionally left behind (they are tiny and the OS cleans temp dirs
// periodically); this keeps the test helper simple.
std::string writeTempFile(const std::string &content,
                          const std::string &suffix) {
  static std::atomic<int> counter{0};
  auto path = std::filesystem::temp_directory_path() /
             ("inertialflow_test_" + std::to_string(counter++) + suffix);
  std::ofstream f(path);
  f << content;
  f.close();
  return path.string();
}

// A single 3x3 grid of 9 nodes (ids 1-9), used as a building block for
// constructing test graphs with one or more disconnected clusters.
std::string gridCoordLines(int idOffset, double xOffset, double yOffset) {
  std::ostringstream out;
  int id = idOffset;
  for (int y = 0; y < 3; ++y) {
    for (int x = 0; x < 3; ++x) {
      out << "v " << id << " " << (xOffset + x) << " " << (yOffset + y)
          << "\n";
      ++id;
    }
  }
  return out.str();
}

std::string gridEdgeLines(int idOffset, int capacity) {
  // 3x3 grid node ids, laid out row-major starting at idOffset:
  //   idOffset+0 idOffset+1 idOffset+2
  //   idOffset+3 idOffset+4 idOffset+5
  //   idOffset+6 idOffset+7 idOffset+8
  std::ostringstream out;
  auto id = [&](int x, int y) { return idOffset + y * 3 + x; };
  for (int y = 0; y < 3; ++y) {
    for (int x = 0; x < 3; ++x) {
      if (x + 1 < 3)
        out << "a " << id(x, y) << " " << id(x + 1, y) << " " << capacity
            << "\n";
      if (y + 1 < 3)
        out << "a " << id(x, y) << " " << id(x, y + 1) << " " << capacity
            << "\n";
    }
  }
  return out.str();
}

// Builds a coordinate/graph file pair containing `numClusters` disconnected
// 3x3-grid clusters (9 nodes each), spaced far apart along the x axis so
// that a coordinate-based cut trivially separates them. Returns
// {coordPath, graphPath}.
std::pair<std::string, std::string> makeClusteredGraph(int numClusters,
                                                        int capacity = 1) {
  std::ostringstream coords, edges;
  int numNodes = numClusters * 9;
  coords << "p aux sp " << numNodes << "\n";
  edges << "p sp " << numNodes << " " << (numClusters * 12) << "\n";
  for (int c = 0; c < numClusters; ++c) {
    int idOffset = c * 9 + 1;
    coords << gridCoordLines(idOffset, c * 100.0, 0.0);
    edges << gridEdgeLines(idOffset, capacity);
  }
  std::string coordPath = writeTempFile(coords.str(), ".co");
  std::string graphPath = writeTempFile(edges.str(), ".gr");
  return {coordPath, graphPath};
}

} // namespace

TEST_CASE("Partitioner: constructor rejects k < 1") {
  CHECK_THROWS_AS(Partitioner(0), std::invalid_argument);
  CHECK_THROWS_AS(Partitioner(-1), std::invalid_argument);
  CHECK_NOTHROW(Partitioner(1));
}

TEST_CASE("Partitioner: loadGraph throws when files are missing") {
  Partitioner p(2);
  CHECK_THROWS_AS(
      p.loadGraph("/nonexistent/path.gr", "/nonexistent/path.co"),
      std::runtime_error);
}

TEST_CASE("Partitioner: loadGraph throws on non-dense/non-1-based vertex ids") {
  std::string coordPath = writeTempFile("p aux sp 2\nv 1 0 0\nv 3 1 1\n", ".co");
  std::string graphPath = writeTempFile("p sp 2 0\n", ".gr");

  Partitioner p(2);
  CHECK_THROWS_AS(p.loadGraph(graphPath, coordPath), std::runtime_error);
}

TEST_CASE("Partitioner: loadGraph throws on edges referencing unknown vertices") {
  std::string coordPath =
      writeTempFile("p aux sp 2\nv 1 0 0\nv 2 1 1\n", ".co");
  std::string graphPath = writeTempFile("p sp 2 1\na 1 5 1\n", ".gr");

  Partitioner p(2);
  CHECK_THROWS_AS(p.loadGraph(graphPath, coordPath), std::runtime_error);
}

TEST_CASE("Partitioner: loadGraph correctly reports vertex count") {
  auto [coordPath, graphPath] = makeClusteredGraph(2);
  Partitioner p(2);
  p.loadGraph(graphPath, coordPath);
  CHECK(p.numVertices() == 18);
}

TEST_CASE("Partitioner: loadGraph parses edge capacities from the file "
         "(regression test: capacities used to be silently hardcoded to 1)") {
  std::string coordPath =
      writeTempFile("p aux sp 2\nv 1 0 0\nv 2 1 0\n", ".co");
  std::string graphPath = writeTempFile("p sp 2 1\na 1 2 42\n", ".gr");

  Partitioner p(2);
  p.loadGraph(graphPath, coordPath);

  const auto &adj = p.getAdjacency();
  REQUIRE(adj.size() == 2);
  REQUIRE(adj[0].size() == 1);
  CHECK(adj[0][0].to == 1);
  CHECK(adj[0][0].capacity == 42);
}

TEST_CASE("Partitioner: loadGraph drops self loops") {
  std::string coordPath =
      writeTempFile("p aux sp 2\nv 1 0 0\nv 2 1 0\n", ".co");
  std::string graphPath =
      writeTempFile("p sp 2 2\na 1 1 99\na 1 2 5\n", ".gr");

  Partitioner p(2);
  p.loadGraph(graphPath, coordPath);

  const auto &adj = p.getAdjacency();
  REQUIRE(adj.size() == 2);
  REQUIRE(adj[0].size() == 1);
  CHECK(adj[0][0].to == 1);
  CHECK(adj[0][0].capacity == 5);
}

TEST_CASE("Partitioner: loadGraph merges duplicate edges, keeping the max "
         "capacity") {
  std::string coordPath =
      writeTempFile("p aux sp 2\nv 1 0 0\nv 2 1 0\n", ".co");
  std::string graphPath =
      writeTempFile("p sp 2 3\na 1 2 3\na 1 2 7\na 1 2 5\n", ".gr");

  Partitioner p(2);
  p.loadGraph(graphPath, coordPath);

  const auto &adj = p.getAdjacency();
  REQUIRE(adj.size() == 2);
  REQUIRE(adj[0].size() == 1);
  CHECK(adj[0][0].to == 1);
  CHECK(adj[0][0].capacity == 7);
}

TEST_CASE("Partitioner: run() rejects out-of-range fraction") {
  auto [coordPath, graphPath] = makeClusteredGraph(2);
  Partitioner p(2);
  MaxGraph g = p.loadGraph(graphPath, coordPath);

  CHECK_THROWS_AS(p.run(g, 0.0), std::invalid_argument);
  CHECK_THROWS_AS(p.run(g, 0.5), std::invalid_argument);
  CHECK_THROWS_AS(p.run(g, 0.6), std::invalid_argument);
  CHECK_THROWS_AS(p.run(g, -0.1), std::invalid_argument);
}

TEST_CASE("Partitioner: k=2 cleanly separates two disconnected clusters "
         "with zero cut edges") {
  auto [coordPath, graphPath] = makeClusteredGraph(2);
  Partitioner p(2);
  MaxGraph g = p.loadGraph(graphPath, coordPath);
  p.run(g, 0.25);

  const auto &nodes = p.getNodes();
  REQUIRE(nodes.size() == 18);

  // All of cluster A (first 9 nodes) must share one partition id, all of
  // cluster B (last 9 nodes) must share a different one.
  int clusterAId = nodes[0].partition_id;
  int clusterBId = nodes[9].partition_id;
  CHECK(clusterAId != clusterBId);
  for (int i = 0; i < 9; ++i) {
    CHECK(nodes[i].partition_id == clusterAId);
  }
  for (int i = 9; i < 18; ++i) {
    CHECK(nodes[i].partition_id == clusterBId);
  }

  // Every partition id used must be within [0, k).
  for (const auto &n : nodes) {
    CHECK(n.partition_id >= 0);
    CHECK(n.partition_id < p.numCells());
  }
}

TEST_CASE("Partitioner: k=3 (non-power-of-two) uses all three partition ids "
         "(regression test: previously only powers of two worked)") {
  auto [coordPath, graphPath] = makeClusteredGraph(3);
  Partitioner p(3);
  MaxGraph g = p.loadGraph(graphPath, coordPath);
  p.run(g, 0.25);

  const auto &nodes = p.getNodes();
  std::set<int> usedIds;
  for (const auto &n : nodes) {
    usedIds.insert(n.partition_id);
    CHECK(n.partition_id >= 0);
    CHECK(n.partition_id < p.numCells());
  }
  CHECK(usedIds.size() == 3);
}

TEST_CASE("Partitioner: k=5 (non-power-of-two, larger) uses all five "
         "non-empty partition ids") {
  // Note: inertial flow always bisects the *current* node set roughly in
  // half at every recursion level, regardless of how many cells remain, so
  // pre-existing clusters are only guaranteed to end up as single
  // partitions when their sizes happen to align with those recursive
  // halves. With 5 equal-sized clusters and k=5, a cluster can legitimately
  // get split across two of the halves (e.g. the first split pulls in a
  // whole neighboring cluster once even one of its nodes is required to
  // fill the quota). So we only assert the invariants that must hold for
  // *any* k, not exact cluster-to-partition alignment.
  auto [coordPath, graphPath] = makeClusteredGraph(5);
  Partitioner p(5);
  MaxGraph g = p.loadGraph(graphPath, coordPath);
  p.run(g, 0.25);

  const auto &nodes = p.getNodes();
  REQUIRE(nodes.size() == 45);

  std::vector<int> cellSizes(5, 0);
  for (const auto &n : nodes) {
    CHECK(n.partition_id >= 0);
    CHECK(n.partition_id < p.numCells());
    cellSizes[n.partition_id]++;
  }

  // Every one of the 5 partition ids should have been used at least once.
  for (int c = 0; c < 5; ++c) {
    CHECK(cellSizes[c] > 0);
  }
}

TEST_CASE("Partitioner: k=4 uses all four non-empty partition ids") {
  auto [coordPath, graphPath] = makeClusteredGraph(4);
  Partitioner p(4);
  MaxGraph g = p.loadGraph(graphPath, coordPath);
  p.run(g, 0.25);

  const auto &nodes = p.getNodes();
  REQUIRE(nodes.size() == 36);

  std::vector<int> cellSizes(4, 0);
  for (const auto &n : nodes) {
    CHECK(n.partition_id >= 0);
    CHECK(n.partition_id < p.numCells());
    cellSizes[n.partition_id]++;
  }
  for (int c = 0; c < 4; ++c) {
    CHECK(cellSizes[c] > 0);
  }
}

TEST_CASE("Partitioner: k=1 assigns every node to partition 0") {
  auto [coordPath, graphPath] = makeClusteredGraph(2);
  Partitioner p(1);
  MaxGraph g = p.loadGraph(graphPath, coordPath);
  p.run(g, 0.25);

  for (const auto &n : p.getNodes()) {
    CHECK(n.partition_id == 0);
  }
}

TEST_CASE("Partitioner: saveResults writes the documented file format") {
  auto [coordPath, graphPath] = makeClusteredGraph(2);
  Partitioner p(2);
  MaxGraph g = p.loadGraph(graphPath, coordPath);
  p.run(g, 0.25);

  std::string outPath = writeTempFile("", ".out");
  p.saveResults(outPath);

  std::ifstream in(outPath);
  std::string firstLine;
  std::getline(in, firstLine);
  CHECK(firstLine == "N 18");

  int lineCount = 0;
  std::string line;
  while (std::getline(in, line)) {
    ++lineCount;
    std::istringstream ls(line);
    int id, partitionId;
    ls >> id >> partitionId;
    CHECK_FALSE(ls.fail());
    CHECK(id == lineCount);
  }
  CHECK(lineCount == 18);
}
