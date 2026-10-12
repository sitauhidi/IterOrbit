#ifndef FILTER_ENGINE_H
#define FILTER_ENGINE_H

#include "AppGraph.h"
#include "OrbitCounter.h"
#include <cstdint>
#include <string>

using CandidateSets = std::unordered_map<int, std::vector<int>>;

class FilterEngine {
public:
    FilterEngine(const AppGraph& data, const AppGraph& pattern, int g_size, bool use_full_graph_mode,
                 int max_iters = 1, bool print_stats = true,
                 bool enable_orbit_throttle = false, double orbit_throttle_threshold = 0.05,
                 bool enable_semi_join = true, bool arc_consistency = true,
                 bool enable_orbit_filter = true, bool post_orbit_induction = true);

    bool run();
    const CandidateSets& getCandidateSets() const { return candidate_sets; }
    const OrbitCounts& getPatternOrbits() const { return pattern_orbits; }
    const AppGraph& getCandidateSubgraph() const { return candidate_subgraph; }

    void countCandidates(const std::string& filter_step_name) const;

private:
    bool ldfNlfFilter(size_t& ldf_only_count);
    bool nlfFilter(const AppGraph& neighbor_graph);
    bool semiJoinFilter(const AppGraph& neighbor_graph);
    bool arcConsistencyFilter(const AppGraph& neighbor_graph);
    bool orbitFilter();
    void induceCandidateSubgraph();
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
    bool arc_consistency;
    bool enable_orbit_filter;
    bool post_orbit_induction;

    CandidateSets candidate_sets;
    OrbitCounts pattern_orbits;
    std::vector<bool> needed_orbits;
    std::unordered_map<int, std::unordered_map<int, int>> pattern_nlf;
    int max_data_label = 0;
    AppGraph candidate_subgraph;
    size_t subgraph_candidate_count = SIZE_MAX;
};

#endif