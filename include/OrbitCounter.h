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
};

class EVOKEOrbitCounter {
public:
    EVOKEOrbitCounter(const AppGraph& g, int size);
    OrbitCounts count();
    const OrbitCountTimings& timings() const { return last_timings; }

private:
    OrbitCountTimings last_timings;
    const AppGraph& graph;
    int graphlet_size;
};

#endif // ORBIT_COUNTER_H
