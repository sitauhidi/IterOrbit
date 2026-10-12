#ifndef SEARCH_ENGINE_H
#define SEARCH_ENGINE_H

#include <cstdint>
#include "AppGraph.h"
#include "FilterEngine.h"
#include "Symmetry.h"

using Mapping = std::vector<int>;

class SearchEngine {
public:
    SearchEngine(const AppGraph& data, const AppGraph& pattern,
                 const CandidateSets& candidates, const std::vector<int>& order,
                 const std::unordered_map<int, int>& pivot, bool induced, bool cbj,
                 bool count_only = false, bool edge_orbit_filter = false,
                 const SymmetryBreaking* symmetry = nullptr);

    void run();
    const std::vector<Mapping>& getMatches() const { return matches; }
    long long getMatchCount() const { return count_only ? match_count : static_cast<long long>(matches.size()); }
    const std::vector<long long>& getNodesPerDepth() const { return nodes_per_depth; }
    const std::vector<long long>& getLocalCandidatesPerDepth() const { return local_per_depth; }
    long long getIndexEntriesBefore() const { return index_entries_before; }
    long long getIndexEntriesAfter() const { return index_entries_after; }
    bool countOverflowed() const { return count_overflow; }

private:
    struct BackwardEdge {
        int from_u;
        std::vector<uint32_t> offsets;
        std::vector<uint32_t> targets;
    };

    void backtrack(int depth);
    int backtrack_cbj(int depth);
    bool computeLocalCandidates(int depth);
    long long countLastLevel(int depth, uint64_t* used_by) const;
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

    std::vector<int> fast_mapping;
    std::vector<int> fast_mapping_local_index;
    std::vector<int> fast_inverse_mapping;
    std::vector<uint8_t> in_use;

    std::vector<std::vector<int>> candidate_lists;
    std::vector<std::vector<BackwardEdge>> backward_edges;
    std::vector<std::vector<int>> backward_non_neighbors;
    std::vector<std::vector<int>> local_candidates;
    std::vector<int> intersect_buffer;
    std::vector<int> identity;
    std::vector<const int*> local_ptr;
    std::vector<size_t> local_len;

    std::vector<std::vector<int>> smaller_depths;
    uint64_t symmetry_multiplier = 1;
    bool count_overflow = false;
    long long index_entries_before = 0;
    long long index_entries_after = 0;

    std::vector<long long> nodes_per_depth;
    std::vector<long long> local_per_depth;

    std::vector<int> order_position;
    size_t conflict_words = 0;
    std::vector<uint64_t> conflict_bits;
    std::vector<uint64_t> backward_masks;

    std::vector<Mapping> matches;
};
#endif
