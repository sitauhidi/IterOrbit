#ifndef ORDER_ENGINE_H
#define ORDER_ENGINE_H

#include "AppGraph.h"
#include "FilterEngine.h"
#include "OrbitCounter.h"

enum class OrderMode { Orbit, Core, CoreOrbit };

class OrderEngine {
public:
    OrderEngine(const AppGraph& pattern, const OrbitCounts& orbits,
                const CandidateSets* candidates = nullptr, OrderMode mode = OrderMode::Core);
    void run();

    const std::vector<int>& getOrder() const { return order; }
    const std::unordered_map<int, int>& getPivot() const { return pivot; }

private:
    void runOrbit();
    void runCore(bool use_orbits);

    const AppGraph& pattern_graph;
    const OrbitCounts& pattern_orbits;
    const CandidateSets* candidates;
    OrderMode mode;
    std::vector<int> order;
    std::unordered_map<int, int> pivot;
};

#endif
