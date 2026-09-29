#ifndef SEARCH_ENGINE_H
#define SEARCH_ENGINE_H

#include <cstdint>
#include "AppGraph.h"
#include "FilterEngine.h"

using Mapping = std::vector<int>;

class SearchEngine {
public:
    SearchEngine(const AppGraph& data, const AppGraph& pattern,
                 const CandidateSets& candidates, const std::vector<int>& order,
                 const std::unordered_map<int, int>& pivot, bool induced, bool cbj,
                 bool count_only = false);

    void run();
    const std::vector<Mapping>& getMatches() const { return matches; }
    long long getMatchCount() const { return count_only ? match_count : static_cast<long long>(matches.size()); }

private:
    struct BackwardEdge {
        int from_u;
        std::vector<uint32_t> offsets;
        std::vector<uint32_t> targets;
    };

    void backtrack(int depth);
    int backtrack_cbj(int depth);
    void computeLocalCandidates(int depth);
    bool inducedConflict(int u, int v, int& conflict_u) const;
    static void mergeConflicts(std::vector<bool>& target, const std::vector<bool>& source);

    const AppGraph& data_graph;
    const AppGraph& pattern_graph;
    const std::vector<int>& order;
    bool is_induced;
    bool enable_cbj;
    bool count_only;
    long long match_count = 0;

    std::vector<int> fast_mapping;             // pattern vertex -> data vertex ID, -1 if unmatched
    std::vector<int> fast_mapping_local_index; // pattern vertex -> index into its candidate list
    std::vector<int> fast_inverse_mapping;     // data vertex ID -> pattern vertex, -1 if unused

    std::vector<std::vector<int>> candidate_lists;       // sorted by data vertex ID
    std::vector<std::vector<BackwardEdge>> backward_edges; // per depth
    std::vector<std::vector<int>> backward_non_neighbors;  // per depth, for induced search
    std::vector<std::vector<int>> local_candidates;        // per depth, reused buffers
    std::vector<int> intersect_buffer;

    std::vector<int> order_position;
    std::vector<std::vector<bool>> conflict_set;

    std::vector<Mapping> matches;
};
#endif
