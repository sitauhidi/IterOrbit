#ifndef ORBIT_COUNTER_H
#define ORBIT_COUNTER_H

#include "AppGraph.h"
#include <vector>
#include <unordered_map>

using OrbitCounts = std::vector<std::vector<long long>>;

struct OrbitCountTimings {
    double relabel_dag_ms = 0;
    double three_ms = 0;
    double four_closed_ms = 0;
    double four_other_ms = 0;
    double four_cycle_ms = 0;
    double four_clique_ms = 0;
    double copy_ms = 0;
    double five_setup_ms = 0;
    double five_tree_ms = 0;
    double five_parallel_ms = 0;
    double five_section_ms[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
};

class EVOKEOrbitCounter {
public:
    EVOKEOrbitCounter(const AppGraph& g, int size, const std::vector<bool>* needed_orbits = nullptr);
    OrbitCounts count();
    const OrbitCountTimings& timings() const { return last_timings; }
    void setCliqueVertexCap(int orbit, const std::vector<double>* cap) {
        (orbit == 14 ? clique4_cap : clique5_cap) = cap;
    }

private:
    OrbitCountTimings last_timings;
    const AppGraph& graph;
    int graphlet_size;
    const std::vector<bool>* needed_orbits;
    const std::vector<double>* clique4_cap = nullptr;
    const std::vector<double>* clique5_cap = nullptr;
    bool needed(int orbit) const { return !needed_orbits || (*needed_orbits)[orbit]; }
};

#endif