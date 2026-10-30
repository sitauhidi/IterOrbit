/// @file
/// Matching orders for the search.

#ifndef ORDER_ENGINE_H
#define ORDER_ENGINE_H

#include "AppGraph.h"
#include "FilterEngine.h"
#include "OrbitCounter.h"

/// How the matching order is chosen. A vertex's "backward" neighbours are its neighbours that
/// are already in the order.
enum class OrderMode {
    /// Start at the pattern vertex with the largest orbit strength (sum of squared orbit
    /// counts); then repeatedly take the vertex with the most backward neighbours, ties by
    /// orbit strength.
    Orbit,
    /// The pattern's 2-core before its trees. Start at the core vertex with the
    /// fewest candidates per neighbour; then most backward neighbours, ties by fewest
    /// candidates.
    Core,
    /// The default. As Core, but the dense-orbit score (triangle, 4-cycle, paw-triangle, diamond and
    /// 4-clique orbit counts) is compared before the candidate count, and the start is the
    /// vertex with the highest such score.
    CoreOrbit
};

/// Chooses the order in which the search maps the pattern vertices.
class OrderEngine {
public:
    /// @param pattern     the pattern graph
    /// @param orbits      its orbit counts (FilterEngine::getPatternOrbits())
    /// @param candidates  the filtered candidate sets; needed by Core and CoreOrbit
    /// @param mode        the ordering rule
    OrderEngine(const AppGraph& pattern, const OrbitCounts& orbits,
                const CandidateSets* candidates = nullptr, OrderMode mode = OrderMode::CoreOrbit);

    /// Computes the order and the pivots.
    void run();

    /// Pattern node IDs in matching order.
    const std::vector<int>& getOrder() const { return order; }
    /// For each pattern node ID after the first, an earlier neighbour in the order (absent
    /// when it has none, e.g. the first vertex of another connected component).
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
