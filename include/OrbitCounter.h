/// @file
/// Graphlet orbit counting on top of the Escape/EVOKE library.

#ifndef ORBIT_COUNTER_H
#define ORBIT_COUNTER_H

#include "AppGraph.h"
#include <vector>
#include <unordered_map>

/// Graphlet orbit counts of a graph: counts[i][k] is the number of non-induced occurrences of
/// orbit k rooted at the vertex with internal index i. Orbits are numbered as in Escape/EVOKE
/// (0-3 for graphlets on up to three vertices, up to 14 with four, up to 72 with five).
using OrbitCounts = std::vector<std::vector<long long>>;

/// Wall-clock time of each step of the last EVOKEOrbitCounter::count() call, in milliseconds.
struct OrbitCountTimings {
    double relabel_dag_ms = 0;      ///< degree-ordered relabelling and DAG construction
    double three_ms = 0;            ///< orbits of graphlets on up to three vertices
    double four_closed_ms = 0;      ///< 4-vertex orbits with a closed form
    double four_other_ms = 0;       ///< paw and diamond orbits
    double four_cycle_ms = 0;       ///< 4-cycle orbit
    double four_clique_ms = 0;      ///< 4-clique orbit
    double copy_ms = 0;             ///< copying the counts back to internal-index order
    double five_setup_ms = 0;       ///< per-edge triangle lists for the 5-vertex routines
    double five_tree_ms = 0;        ///< 5-vertex tree orbits
    double five_parallel_ms = 0;    ///< the nine independent 5-vertex routines, run concurrently
    double five_section_ms[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};   ///< each of those nine
};

/// Counts graphlet orbits with the Escape/EVOKE library.
///
/// Several of Escape's routines are replaced by reimplementations in OrbitCounter.cpp that
/// return the same counts (checked by test/orbit_equivalence.cpp). The counter is used on the
/// pattern (all orbits) and, in every filtering round, on the data graph or its
/// candidate-induced subgraph.
class EVOKEOrbitCounter {
public:
    /// @param g              the graph to count on
    /// @param size           graphlet size: 3, 4 or 5 vertices
    /// @param needed_orbits  if given, orbits whose entry is false may be left at zero. The
    ///                       filter passes the orbits some pattern vertex has: an orbit no
    ///                       pattern vertex has cannot reject a candidate.
    EVOKEOrbitCounter(const AppGraph& g, int size, const std::vector<bool>* needed_orbits = nullptr);

    /// Counts the orbits. The result is indexed by internal index, then orbit.
    OrbitCounts count();

    /// Timings of the last count() call.
    const OrbitCountTimings& timings() const { return last_timings; }

    /// Caps the clique count of each vertex (orbit 14: 4-cliques; any other value: orbit 72,
    /// 5-cliques). `cap` is indexed by internal index. A vertex's count is then either exact
    /// or at least its cap, so every test "count >= requirement" with requirement <= cap is
    /// unchanged, while cliques among vertices that have all reached their cap are skipped.
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
