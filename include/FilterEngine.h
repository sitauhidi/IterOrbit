#ifndef FILTER_ENGINE_H
#define FILTER_ENGINE_H

#include "AppGraph.h"
#include "OrbitCounter.h"
#include <string>

using CandidateSets = std::unordered_map<int, std::vector<int>>;

class FilterEngine {
public:
    FilterEngine(const AppGraph& data, const AppGraph& pattern, int g_size, bool use_full_graph_mode,
                 int max_iters = 1, bool print_stats = true,
                 bool enable_orbit_throttle = false, double orbit_throttle_threshold = 0.05,
                 bool enable_semi_join = true);

    bool run();
    const CandidateSets& getCandidateSets() const { return candidate_sets; }
    const OrbitCounts& getPatternOrbits() const { return pattern_orbits; }
    const AppGraph& getCandidateSubgraph() const { return candidate_subgraph; }

    void countCandidates(const std::string& filter_step_name) const;

private:
    bool ldfFilter();
    bool nlfFilter(const AppGraph& neighbor_graph);
    bool semiJoinFilter(const AppGraph& neighbor_graph);
    bool orbitFilter();
    size_t getCandidateCount() const;

    const AppGraph& data_graph;
    const AppGraph& pattern_graph;
    int graphlet_size;
    bool use_subgraph;
    int max_iterations;
    bool print_stats;
    bool enable_orbit_throttle;
    double orbit_throttle_threshold;
    bool enable_semi_join;

    CandidateSets candidate_sets;
    OrbitCounts pattern_orbits;
    std::unordered_map<int, std::unordered_map<int, int>> pattern_nlf;
    int max_data_label = 0;
    AppGraph candidate_subgraph;
};

#endif // FILTER_ENGINE_H

