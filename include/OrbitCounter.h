#ifndef ORBIT_COUNTER_H
#define ORBIT_COUNTER_H

#include "AppGraph.h"
#include <vector>
#include <unordered_map>

using OrbitCounts = std::vector<std::vector<long long>>;

class EVOKEOrbitCounter {
public:
    EVOKEOrbitCounter(const AppGraph& g, int size);
    OrbitCounts count();

private:
    const AppGraph& graph;
    int graphlet_size;
};

#endif // ORBIT_COUNTER_H
