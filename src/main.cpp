/// @file
/// Command-line driver: load the graphs, filter (FilterEngine), order (OrderEngine), compute
/// symmetry-breaking constraints when counting (Symmetry), search (SearchEngine), and report the
/// match count and the time of each phase.

#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <algorithm>

#include "AppGraph.h"
#include "FilterEngine.h"
#include "OrderEngine.h"
#include "SearchEngine.h"
#include "Symmetry.h"

void print_usage() {
    std::cerr << "Usage: ./build/iterorbit --data <path> --pattern <path> [options]\n\n"
              << "Options:\n"
              << "  --graphlet-size <3|4|5>  Set the graphlet size for orbit counting (default: 4)\n"
              << "  --iterate <N>            Set max filter iterations (0 = until convergence, default: 1)\n"
              << "  --induced                Perform an induced subgraph isomorphism search\n"
              << "  --use-full-graph         Use the full data graph for orbit filtering instead of a subgraph\n"
              << "  --cbj                    Use conflict-directed backjumping during enumeration\n"
              << "  --filter-only            Run filtering only; skip ordering and search\n"
              << "  --verbose                Print all found matches to the console\n"
              << "  --orbit-throttle         Skip the per-iteration orbit recount once candidates\n"
              << "                           have nearly stabilized (opt-in; changes what a round\n"
              << "                           costs, not the fixpoint it converges to)\n"
              << "  --orbit-throttle-threshold <f>  Min candidate-shrink fraction to keep running the\n"
              << "                                  orbit recount every iteration (default: 0.05)\n"
              << "  --disable-semi-join      Skip the edge-consistency (semi-join) filter stage\n"
              << "                           (enabled by default)\n"
              << "  --no-orbit-filter        Skip the orbit filter stage in every round (ablation);\n"
              << "                           contraction and NLF/semi-join rounds still run\n"
              << "  --single-induction       Induce G_d^(t) once per round, before the orbit stage;\n"
              << "                           by default it is induced again after the orbit stage,\n"
              << "                           for the next round's NLF/semi-join\n"
              << "  --no-arc-consistency     Run the semi-join as one pass per filter round instead of\n"
              << "                           to a fixpoint (arc consistency, the default)\n"
              << "  --count-only             Count matches without materializing each embedding\n"
              << "                           (avoids O(matches x pattern-size) memory for patterns\n"
              << "                           with huge result counts; incompatible with --verbose)\n"
              << "  --order <u,u,...>        Search in this matching order (pattern vertex IDs, a\n"
              << "                           permutation) instead of OrderEngine's (diagnostic)\n"
              << "  --order-mode <m>         Matching order: core (default; 2-core first, fewest\n"
              << "                           candidates), orbit (by orbit strength), core-orbit (core,\n"
              << "                           dense-orbit score before candidates) -- see OrderEngine.h\n"
              << "  --search-stats           Print per-depth search-tree node and local-candidate counts\n"
              << "  --no-edge-orbit-filter   Keep candidate-space index entries whose data edge lies in\n"
              << "                           fewer triangles than the pattern edge (ablation)\n"
              << "  --no-symmetry            With --count-only: count every embedding instead of one per\n"
              << "                           class of pattern automorphisms (ablation)\n";
}

void print_ms(const char* label, double ms) {
    std::ios::fmtflags flags = std::cout.flags();
    std::streamsize precision = std::cout.precision();
    std::cout << label << std::fixed << std::setprecision(3) << ms << " ms." << std::endl;
    std::cout.flags(flags);
    std::cout.precision(precision);
}

int main(int argc, char* argv[]) {
    std::vector<std::string> args(argv + 1, argv + argc);

    if (args.empty()) {
        print_usage();
        return 1;
    }

    std::string data_path, pattern_path;
    int graphlet_size = 4;
    int iterations = 1;
    bool induced_search = false;
    bool use_full_graph = false;
    bool cbj_search = false;
    bool verbose = false;
    bool filter_only = false;
    bool orbit_throttle = false;
    double orbit_throttle_threshold = 0.05;
    bool enable_semi_join = true;
    bool arc_consistency = true;
    bool enable_orbit_filter = true;
    bool post_orbit_induction = true;
    bool count_only = false;
    bool search_stats = false;
    bool edge_orbit_filter = true;
    bool symmetry_breaking = true;
    std::string order_arg;
    OrderMode order_mode = OrderMode::Core;

    bool iterate_specified = false;

    for (size_t i = 0; i < args.size(); ++i) {
        if (args[i] == "--data" && i + 1 < args.size()) {
            data_path = args[++i];
        } else if (args[i] == "--pattern" && i + 1 < args.size()) {
            pattern_path = args[++i];
        } else if (args[i] == "--graphlet-size" && i + 1 < args.size()) {
            try {
                graphlet_size = std::stoi(args[++i]);
            } catch(...) {
                std::cerr << "Invalid graphlet size.\n";
                print_usage();
                return 1;
            }
        } else if (args[i] == "--iterate" && i + 1 < args.size()) {
            iterate_specified = true;
            try {
                iterations = std::stoi(args[++i]);
                if (iterations < 0) throw std::invalid_argument("negative");
            } catch(...) {
                std::cerr << "Invalid iteration count. --iterate must be an integer >= 0.\n";
                print_usage();
                return 1;
            }
        } else if (args[i] == "--induced") {
            induced_search = true;
        } else if (args[i] == "--use-full-graph") {
            use_full_graph = true;
        } else if (args[i] == "--cbj") {
            cbj_search = true;
        } else if (args[i] == "--filter-only") {
            filter_only = true;
        } else if (args[i] == "--verbose") {
            verbose = true;
        } else if (args[i] == "--orbit-throttle") {
            orbit_throttle = true;
        } else if (args[i] == "--orbit-throttle-threshold" && i + 1 < args.size()) {
            try {
                orbit_throttle_threshold = std::stod(args[++i]);
            } catch(...) {
                std::cerr << "Invalid --orbit-throttle-threshold value.\n";
                print_usage();
                return 1;
            }
        } else if (args[i] == "--disable-semi-join") {
            enable_semi_join = false;
        } else if (args[i] == "--arc-consistency") {
            arc_consistency = true;   // already the default: the flag is accepted and changes nothing
        } else if (args[i] == "--no-arc-consistency") {
            arc_consistency = false;
        } else if (args[i] == "--no-orbit-filter") {
            enable_orbit_filter = false;
        } else if (args[i] == "--single-induction") {
            post_orbit_induction = false;
        } else if (args[i] == "--count-only") {
            count_only = true;
        } else if (args[i] == "--no-edge-orbit-filter") {
            edge_orbit_filter = false;
        } else if (args[i] == "--no-symmetry") {
            symmetry_breaking = false;
        } else if (args[i] == "--search-stats") {
            search_stats = true;
        } else if (args[i] == "--order" && i + 1 < args.size()) {
            order_arg = args[++i];
        } else if (args[i] == "--order-mode" && i + 1 < args.size()) {
            const std::string& m = args[++i];
            if (m == "orbit") order_mode = OrderMode::Orbit;
            else if (m == "core") order_mode = OrderMode::Core;
            else if (m == "core-orbit") order_mode = OrderMode::CoreOrbit;
            else {
                std::cerr << "Invalid --order-mode (orbit, core or core-orbit).\n";
                print_usage();
                return 1;
            }
        }
    }

    if (data_path.empty() || pattern_path.empty()) {
        std::cerr << "Error: --data and --pattern paths are required.\n";
        print_usage();
        return 1;
    }
    if (graphlet_size != 3 && graphlet_size != 4 && graphlet_size != 5) {
        std::cerr << "Error: --graphlet-size must be 3, 4, or 5.\n";
        print_usage();
        return 1;
    }
    if (use_full_graph && iterate_specified) {
        std::cerr << "Error: --use-full-graph cannot be used together with --iterate.\n";
        print_usage();
        return 1;
    }
    if (count_only && verbose) {
        std::cerr << "Error: --count-only cannot be used together with --verbose "
                     "(count-only mode never materializes matches to print).\n";
        print_usage();
        return 1;
    }

    auto load_start = std::chrono::high_resolution_clock::now();
    AppGraph data_graph, pattern_graph;
    if (!data_graph.readFromFile(data_path)) {
        std::cerr << "Error reading data graph from " << data_path << std::endl;
        return 1;
    }
    if (!pattern_graph.readFromFile(pattern_path)) {
        std::cerr << "Error reading pattern graph from " << pattern_path << std::endl;
        return 1;
    }
    std::chrono::duration<double, std::milli> load_elapsed = std::chrono::high_resolution_clock::now() - load_start;
    std::cout << "Graphs loaded successfully." << std::endl;
    print_ms("Total loading time: ", load_elapsed.count());
    std::cout << "Search type: " << (induced_search ? "Induced" : "Non-Induced") << std::endl;
    std::cout << "Backjumping: " << (cbj_search ? "Enabled (CBJ)" : "Disabled") << std::endl;
    std::cout << "Graphlet size: " << graphlet_size << std::endl;
    std::cout << "Filter mode: " << (use_full_graph ? "Full Graph" : "Subgraph") << std::endl;
    std::cout << "Max iterations: " << (iterations == 0 ? "Until Convergence (0)" : std::to_string(iterations)) << std::endl;
    if (orbit_throttle) {
        std::cout << "Orbit recount throttle: Enabled (threshold " << orbit_throttle_threshold << ")" << std::endl;
    }
    if (!enable_semi_join) {
        std::cout << "Semi-Join filter: Disabled" << std::endl;
    } else if (!arc_consistency) {
        std::cout << "Semi-Join filter: Single pass (arc consistency disabled)" << std::endl;
    }
    if (!enable_orbit_filter) {
        std::cout << "Orbit filter: Disabled" << std::endl;
    }
    if (!post_orbit_induction) {
        std::cout << "Induction: once per round (single)" << std::endl;
    }

    auto start = std::chrono::high_resolution_clock::now();

    bool print_filter_stats = filter_only || verbose;
    FilterEngine filter_engine(data_graph, pattern_graph, graphlet_size, use_full_graph, iterations,
                                print_filter_stats, orbit_throttle, orbit_throttle_threshold,
                                enable_semi_join, arc_consistency, enable_orbit_filter, post_orbit_induction);
    const bool filter_ok = filter_engine.run();

    auto filter_end = std::chrono::high_resolution_clock::now();
    // Fractional ms, printed fixed-point so large values never switch to exponent notation.
    // Also printed when filtering refutes the instance, so a refuted run's cost is known too.
    std::chrono::duration<double, std::milli> filter_elapsed = filter_end - start;
    print_ms("\nTotal filtering time: ", filter_elapsed.count());
    if (!filter_ok) {
        std::cout << "Matches found: 0" << std::endl;
        return 0;
    }

    if (filter_only) {
        return 0;
    }

    auto order_start = std::chrono::high_resolution_clock::now();
    OrderEngine order_engine(pattern_graph, filter_engine.getPatternOrbits(),
                             &filter_engine.getCandidateSets(), order_mode);
    order_engine.run();
    std::vector<int> order = order_engine.getOrder();
    if (!order_arg.empty()) {
        std::vector<int> given;
        std::stringstream ss(order_arg);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            if (!tok.empty()) given.push_back(std::stoi(tok));
        }
        std::vector<int> sorted_given = given, ids = pattern_graph.original_node_ids;
        std::sort(sorted_given.begin(), sorted_given.end());
        std::sort(ids.begin(), ids.end());
        if (sorted_given != ids) {
            std::cerr << "Error: --order must be a permutation of the pattern's vertex IDs.\n";
            return 1;
        }
        order = given;
    }
    auto order_end = std::chrono::high_resolution_clock::now();
    auto order_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(order_end - order_start);
    std::cout << "Total ordering time: " << order_elapsed.count() << " us." << std::endl;
    if (search_stats) {
        std::cout << "Matching order:";
        for (int u : order) std::cout << ' ' << u;
        std::cout << std::endl;
    }

    // Symmetry breaking only applies when counting: listing matches needs every one of them.
    SymmetryBreaking symmetry;
    const bool use_symmetry = symmetry_breaking && count_only;
    if (use_symmetry) {
        auto sym_start = std::chrono::high_resolution_clock::now();
        symmetry = computeSymmetryBreaking(pattern_graph, filter_engine.getPatternOrbits(), order);
        std::chrono::duration<double, std::milli> sym_elapsed = std::chrono::high_resolution_clock::now() - sym_start;
        std::cout << "Symmetry breaking: " << symmetry.levels << " levels, multiplier "
                  << symmetry.multiplier << (symmetry.truncated ? " (truncated)" : "") << std::endl;
        print_ms("Total symmetry time: ", sym_elapsed.count());
    }

    auto search_start = std::chrono::high_resolution_clock::now();
    SearchEngine search_engine(filter_engine.getCandidateSubgraph(), pattern_graph,
                               filter_engine.getCandidateSets(), order,
                               order_engine.getPivot(), induced_search, cbj_search, count_only,
                               edge_orbit_filter, use_symmetry ? &symmetry : nullptr);
    std::cout << "Index entries: " << search_engine.getIndexEntriesAfter() << " of "
              << search_engine.getIndexEntriesBefore() << " after the edge-orbit filter" << std::endl;
    search_engine.run();
    if (search_engine.countOverflowed()) {
        std::cerr << "Warning: the match count exceeds 64 bits; reported as INT64_MAX." << std::endl;
    }
    auto search_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> search_elapsed = search_end - search_start;
    print_ms("Total searching time: ", search_elapsed.count());

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> total_elapsed = end - start;

    std::cout << "\nMatches found: " << search_engine.getMatchCount() << std::endl;

    if (search_stats) {
        const auto& nodes = search_engine.getNodesPerDepth();
        const auto& local = search_engine.getLocalCandidatesPerDepth();
        long long total_nodes = 0;
        std::cout << "Search stats (depth: pattern vertex, nodes, local candidates):" << std::endl;
        for (size_t d = 0; d < nodes.size(); ++d) {
            std::cout << "  depth " << d << ": u" << order[d] << ' ' << nodes[d] << ' ' << local[d] << std::endl;
            total_nodes += nodes[d];
        }
        std::cout << "Search nodes: " << total_nodes << std::endl;
    }

    if (verbose) {
        const auto& matches = search_engine.getMatches();
        std::cout << "\n--- Found Matches ---" << std::endl;
        for (const auto& match : matches) {
            std::cout << "{";
            bool first = true;
            for (size_t u = 0; u < match.size(); ++u) {
                if (!first) std::cout << ", ";
                std::cout << u << ": " << match[u];
                first = false;
            }
            std::cout << "}" << std::endl;
        }
        std::cout << "---------------------" << std::endl;
    }

    print_ms("\nTotal execution time: ", total_elapsed.count());

    return 0;
}