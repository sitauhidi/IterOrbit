/// @file
/// Backtracking enumeration over a candidate-space index.

#ifndef SEARCH_ENGINE_H
#define SEARCH_ENGINE_H

#include <cstdint>
#include "AppGraph.h"
#include "FilterEngine.h"
#include "Symmetry.h"

/// One embedding: the data node ID each pattern vertex is mapped to, indexed by pattern node ID
/// (pattern node IDs are 0..n-1).
using Mapping = std::vector<int>;

/// Backtracking enumeration of the embeddings of the pattern, within the candidate sets.
///
/// The constructor builds a candidate-space index: for every pattern edge {u', u} with u' before
/// u in the matching order ("backward edge" of u) and every candidate v' of u', the candidates
/// of u adjacent to v'. A pattern vertex's local candidates are then the intersection of one
/// such list per already-mapped backward neighbour. Every result is exact: the optional parts
/// (backjumping, the edge-orbit filter, symmetry breaking, counting the last level) change the
/// work done, never the embeddings found or their number.
class SearchEngine {
public:
    /// @param data        the graph to search: the data graph or an induced subgraph of it that
    ///                    contains every candidate (FilterEngine::getCandidateSubgraph())
    /// @param pattern     the pattern graph
    /// @param candidates  the candidate sets
    /// @param order       the matching order (pattern node IDs)
    /// @param pivot       OrderEngine::getPivot()
    /// @param induced     induced subgraph isomorphism: non-adjacent pattern vertices must map
    ///                    to non-adjacent data vertices
    /// @param cbj         conflict-directed backjumping
    /// @param count_only  count the embeddings without storing them; the last pattern vertex's
    ///                    matches are then counted without being visited one by one
    /// @param edge_orbit_filter  drop an index entry (v', v) along pattern edge {u', u} when
    ///                    the data edge lies in fewer triangles than the pattern edge
    /// @param symmetry    with count_only: count one embedding per class of pattern
    ///                    automorphisms and multiply by the class size
    SearchEngine(const AppGraph& data, const AppGraph& pattern,
                 const CandidateSets& candidates, const std::vector<int>& order,
                 const std::unordered_map<int, int>& pivot, bool induced, bool cbj,
                 bool count_only = false, bool edge_orbit_filter = false,
                 const SymmetryBreaking* symmetry = nullptr);

    /// Enumeration with a match limit (the protocol of Sun and Luo): without count_only, visit the
    /// embeddings one by one and stop after `limit` of them (0: no limit); `store` keeps them for
    /// getMatches(). Call before run().
    void setMatchLimit(long long limit, bool store) { match_limit = limit; store_matches = store; }

    /// Runs the search.
    void run();

    /// The embeddings found (empty with count_only, or when not stored).
    const std::vector<Mapping>& getMatches() const { return matches; }
    /// The number of embeddings, saturated at INT64_MAX (see countOverflowed()); with a match
    /// limit, at most the limit.
    long long getMatchCount() const { return match_count; }
    /// Whether the search stopped at the match limit.
    bool limitReached() const { return stopped; }
    /// Search-tree nodes visited at each order depth.
    const std::vector<long long>& getNodesPerDepth() const { return nodes_per_depth; }
    /// Local candidates computed at each order depth, summed over the visits.
    const std::vector<long long>& getLocalCandidatesPerDepth() const { return local_per_depth; }
    /// Entries of the candidate-space index before the edge-orbit filter.
    long long getIndexEntriesBefore() const { return index_entries_before; }
    /// Entries of the candidate-space index after the edge-orbit filter.
    long long getIndexEntriesAfter() const { return index_entries_after; }
    /// Whether the count exceeded 64 bits and was saturated.
    bool countOverflowed() const { return count_overflow; }

private:
    /// The index along one backward edge of the pattern vertex at some depth: for the i-th
    /// candidate of `from_u`, targets[offsets[i]..offsets[i+1]) are the adjacent candidates of
    /// that vertex, as ascending positions in its candidate list.
    struct BackwardEdge {
        int from_u;
        std::vector<uint32_t> offsets;
        std::vector<uint32_t> targets;
    };

    /// Chronological backtracking from `depth`.
    void backtrack(int depth);
    /// Backtracking with conflict-directed backjumping. Returns the depth to return to.
    int backtrack_cbj(int depth);
    /// Computes the local candidates of order[depth] from the index (local_ptr, local_len).
    /// Returns whether the symmetry constraints removed any of them, which backjumping must
    /// then treat as a conflict with the earlier vertices of those constraints.
    bool computeLocalCandidates(int depth);
    /// Counts the embeddings completed at the last pattern vertex without visiting each one;
    /// the pattern vertices that block a candidate are added to the conflict set `used_by`.
    long long countLastLevel(int depth, uint64_t* used_by) const;
    /// Under induced semantics: whether mapping u to v puts an edge between the images of two
    /// non-adjacent pattern vertices; `conflict_u` receives the other pattern vertex.
    bool inducedConflict(int u, int v, int& conflict_u) const;
    uint64_t* conflicts(int depth) { return conflict_bits.data() + static_cast<size_t>(depth) * conflict_words; }
    void setConflict(uint64_t* set, int d) const { set[d >> 6] |= uint64_t{1} << (d & 63); }
    void mergeConflicts(int target_depth, int source_depth);

    const AppGraph& data_graph;
    const AppGraph& pattern_graph;
    const std::vector<int>& order;
    bool is_induced;
    bool enable_cbj;
    bool count_only;
    long long match_count = 0;
    long long match_limit = 0;
    bool store_matches = true;
    bool stopped = false;

    /// The partial embedding: pattern node ID -> data node ID (-1 if unmapped). Like the three
    /// arrays below, a flat array, so that assigning and undoing take O(1).
    std::vector<int> fast_mapping;
    /// Pattern node ID -> position of its image in its candidate list.
    std::vector<int> fast_mapping_local_index;
    /// Data node ID -> the pattern node ID mapped to it (-1 if none).
    std::vector<int> fast_inverse_mapping;
    /// Data node ID -> whether some pattern vertex is mapped to it.
    std::vector<uint8_t> in_use;

    /// Per pattern vertex, its candidates sorted by data node ID.
    std::vector<std::vector<int>> candidate_lists;
    /// Per depth, the index along each backward edge of order[depth].
    std::vector<std::vector<BackwardEdge>> backward_edges;
    /// Per depth, the earlier pattern vertices not adjacent to order[depth] (induced search).
    std::vector<std::vector<int>> backward_non_neighbors;
    /// Per depth, a buffer for the intersection; swapped with intersect_buffer instead of
    /// copied (see computeLocalCandidates).
    std::vector<std::vector<int>> local_candidates;
    /// The second buffer of that swap.
    std::vector<int> intersect_buffer;
    /// 0, 1, 2, ...: "every candidate", for a vertex with no backward neighbour.
    std::vector<int> identity;
    /// Per depth, the local candidates, as positions in the candidate list.
    std::vector<const int*> local_ptr;
    /// Per depth, the number of local candidates.
    std::vector<size_t> local_len;

    /// Symmetry breaking: per depth, the earlier depths whose image must have a smaller data
    /// node ID than this depth's.
    std::vector<std::vector<int>> smaller_depths;
    /// The factor the final count is multiplied by (SymmetryBreaking::multiplier).
    uint64_t symmetry_multiplier = 1;
    bool count_overflow = false;
    long long index_entries_before = 0;
    long long index_entries_after = 0;

    std::vector<long long> nodes_per_depth;
    std::vector<long long> local_per_depth;

    /// Pattern node ID -> depth in the matching order.
    std::vector<int> order_position;
    /// Backjumping: the number of 64-bit words in one conflict set.
    size_t conflict_words = 0;
    /// Backjumping: one conflict set per depth, a bitset in which bit d stands for the pattern
    /// vertex at depth d.
    std::vector<uint64_t> conflict_bits;
    /// Backjumping: per depth, the bitset of the depths of order[depth]'s backward neighbours.
    std::vector<uint64_t> backward_masks;

    std::vector<Mapping> matches;
};
#endif
