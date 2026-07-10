#include <iostream>
#include <thread>

#include "cmdparser.hpp"
#include "partitioner.hpp"

void configure_parser(cli::Parser& parser) {
  parser.set_required<std::string>("g", "input_graph",
                                   "Input graph in DIMACs format.");
  parser.set_required<std::string>("c", "input_coord",
                                   "Input coords in DIMACs format.");
  parser.set_required<int>("k", "num_cells",
                           "The number of cells to partition the graph into.");
  parser.set_optional<bool>("s", "show_stats", false, "Show statistics.");
  parser.set_optional<std::string>("o", "output_file", "dump.txt",
                                   "Write the partition to this file.");
  parser.set_optional<double>("f", "fraction", 0.25,
                              "Fraction of vertices to pick, must be < 0.5.");
};

int main(int argc, char* argv[]) {
  cli::Parser parser(argc, argv, "Inertial-Flow");
  configure_parser(parser);
  parser.run_and_exit_if_error();

  const std::string gFile = parser.get<std::string>("g");
  const std::string cFile = parser.get<std::string>("c");
  const int k = parser.get<int>("k");
  const bool showStats = parser.get<bool>("s");
  const std::string oFile = parser.get<std::string>("o");
  const double fraction = parser.get<double>("f");

  if (fraction <= 0.001 || fraction >= 0.49) {
    std::cerr << "Given fraction should be between (0, 0.5), was " << fraction
              << "!\n";
    return -1;
  }

  if (k < 1) {
    std::cerr << "Number of cells k must be >= 1, was " << k << "!\n";
    return -1;
  }

  try {
    Partitioner p(k);
    MaxGraph graph = p.loadGraph(gFile, cFile);

    if (showStats) {
      std::cout << "** Loaded a graph with " << p.numVertices()
                << " many nodes **\n";
    }

    p.run(graph, fraction);
    p.saveResults(oFile);

    if (showStats) {
      p.printStats();
    }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return -1;
  }

  return 0;
}
