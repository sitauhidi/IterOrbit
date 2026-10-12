/// @file
/// Symmetry-breaking constraints for counting embeddings.

#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <cstdint>
#include <vector>
#include "AppGraph.h"
#include "OrbitCounter.h"

/// Symmetry-breaking constraints for counting embeddings (Grochow and Kellis).
///
/// Embeddings that differ by an automorphism of the pattern are counted once: the search
/// enforces the constraints below, which keep exactly one embedding out of every `multiplier`,
/// and multiplies its count by `multiplier`.
struct SymmetryBreaking {
    /// smaller_than[w]: the pattern vertices (internal indices) whose image must have a
    /// smaller data node ID than the image of w.
    std::vector<std::vector<int>> smaller_than;
    /// Number of embeddings each counted one stands for: the product of the orbit sizes used,
    /// which is |Aut(G_p)| when the computation was not truncated.
    uint64_t multiplier = 1;
    /// Number of vertices whose orbit contributed constraints.
    int levels = 0;
    /// The work budget ran out, or the multiplier would overflow. The constraints found so far
    /// are still valid: they break a subgroup of the automorphisms, with the matching multiplier.
    bool truncated = false;
    /// Work units spent, against the budget.
    long long work = 0;
};

/// Computes the constraints for `pattern` under the matching order `order` (node IDs).
///
/// Vertices are fixed one at a time in matching order; at each step the orbit of the vertex
/// under the automorphisms that fix the earlier ones is found by individualisation-refinement,
/// starting from a colouring by label and orbit-count vector (`orbits`, by internal index).
/// `budget` bounds the total refinement work.
SymmetryBreaking computeSymmetryBreaking(const AppGraph& pattern, const OrbitCounts& orbits,
                                         const std::vector<int>& order, long long budget = 20000000);

#endif
