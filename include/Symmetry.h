#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <cstdint>
#include <vector>
#include "AppGraph.h"
#include "OrbitCounter.h"

struct SymmetryBreaking {
    std::vector<std::vector<int>> smaller_than;
    uint64_t multiplier = 1;
    int levels = 0;
    bool truncated = false;
    long long work = 0;
};

SymmetryBreaking computeSymmetryBreaking(const AppGraph& pattern, const OrbitCounts& orbits,
                                         const std::vector<int>& order, long long budget = 20000000);

#endif
