#include "SearchEngine.h"
#include <algorithm>
#include <unordered_map>

namespace {

size_t intersectSorted(const int* a, size_t na, const int* b, size_t nb, int* out) {
    if (na > nb) {
        std::swap(a, b);
        std::swap(na, nb);
    }
    size_t count = 0;
    if (na * 32 < nb) {
        size_t lo = 0;
        for (size_t i = 0; i < na && lo < nb; ++i) {
            size_t step = 1, hi = lo;
            while (hi < nb && b[hi] < a[i]) {
                lo = hi;
                hi += step;
                step <<= 1;
            }
            lo = std::lower_bound(b + lo, b + std::min(hi + 1, nb), a[i]) - b;
            if (lo < nb && b[lo] == a[i]) out[count++] = a[i];
        }
        return count;
    }
    size_t i = 0, j = 0;
    while (i < na && j < nb) {
        if (a[i] < b[j]) {
            ++i;
        } else if (b[j] < a[i]) {
            ++j;
        } else {
            out[count++] = a[i];
            ++i;
            ++j;
        }
    }
    return count;
}

}

SearchEngine::SearchEngine(const AppGraph& data, const AppGraph& pattern,
                           const CandidateSets& candidates, const std::vector<int>& o,
                           const std::unordered_map<int, int>&,
                           bool induced, bool cbj, bool count_only)
    : data_graph(data), pattern_graph(pattern), order(o), is_induced(induced),
      enable_cbj(cbj), count_only(count_only) {

    const size_t num_p_nodes = pattern_graph.original_node_ids.size();
    int max_v_id = data_graph.original_node_ids.empty() ? 0 : data_graph.original_node_ids.back();
    for (const auto& pair : candidates) {
        for (int v : pair.second) max_v_id = std::max(max_v_id, v);
    }
    fast_mapping.assign(num_p_nodes, -1);
    fast_mapping_local_index.assign(num_p_nodes, -1);
    fast_inverse_mapping.assign(max_v_id + 1, -1);

    order_position.assign(num_p_nodes, -1);
    for (size_t d = 0; d < order.size(); ++d) {
        order_position[order[d]] = static_cast<int>(d);
    }
    conflict_set.assign(order.size(), std::vector<bool>(num_p_nodes, false));

    candidate_lists.assign(num_p_nodes, {});
    for (const auto& pair : candidates) {
        int u = pair.first;
        if (u >= 0 && static_cast<size_t>(u) < num_p_nodes) {
            candidate_lists[u] = pair.second;
            std::sort(candidate_lists[u].begin(), candidate_lists[u].end());
        }
    }

    std::vector<int> index_of(max_v_id + 1, -1);
    size_t max_candidates = 0;
    backward_edges.assign(order.size(), {});
    backward_non_neighbors.assign(order.size(), {});
    for (size_t d = 0; d < order.size(); ++d) {
        int u = order[d];
        const auto& cands_u = candidate_lists[u];
        max_candidates = std::max(max_candidates, cands_u.size());
        for (size_t j = 0; j < cands_u.size(); ++j) index_of[cands_u[j]] = static_cast<int>(j);

        int u_id = pattern_graph.original_node_ids[u];
        for (size_t e = 0; e < d; ++e) {
            int from_u = order[e];
            if (!pattern_graph.hasEdge(pattern_graph.original_node_ids[from_u], u_id)) {
                backward_non_neighbors[d].push_back(from_u);
                continue;
            }
            BackwardEdge edge;
            edge.from_u = from_u;
            const auto& cands_from = candidate_lists[from_u];
            edge.offsets.assign(cands_from.size() + 1, 0);
            for (size_t i = 0; i < cands_from.size(); ++i) {
                int v_idx = data_graph.get_node_idx(cands_from[i]);
                if (v_idx != -1) {
                    for (Escape::EdgeIdx k = data_graph.c_graph->offsets[v_idx];
                         k < data_graph.c_graph->offsets[v_idx + 1]; ++k) {
                        int j = index_of[data_graph.original_node_ids[data_graph.c_graph->nbors[k]]];
                        if (j != -1) edge.targets.push_back(static_cast<uint32_t>(j));
                    }
                }
                edge.offsets[i + 1] = static_cast<uint32_t>(edge.targets.size());
            }
            backward_edges[d].push_back(std::move(edge));
        }

        for (int v : cands_u) index_of[v] = -1;
    }

    local_candidates.assign(order.size() + 1, std::vector<int>(max_candidates));
    intersect_buffer.assign(max_candidates, 0);
}

void SearchEngine::mergeConflicts(std::vector<bool>& target, const std::vector<bool>& source) {
    for (size_t x = 0; x < source.size(); ++x) {
        if (source[x]) target[x] = true;
    }
}

void SearchEngine::run() {
    matches.clear();
    match_count = 0;
    std::fill(fast_mapping.begin(), fast_mapping.end(), -1);
    std::fill(fast_mapping_local_index.begin(), fast_mapping_local_index.end(), -1);
    std::fill(fast_inverse_mapping.begin(), fast_inverse_mapping.end(), -1);
    if (order.empty()) return;
    if (enable_cbj) {
        backtrack_cbj(0);
    } else {
        backtrack(0);
    }
}

void SearchEngine::computeLocalCandidates(int depth) {
    auto& out = local_candidates[depth];
    const auto& edges = backward_edges[depth];
    const size_t num_candidates = candidate_lists[order[depth]].size();

    if (edges.empty()) {
        out.resize(num_candidates);
        for (size_t j = 0; j < num_candidates; ++j) out[j] = static_cast<int>(j);
        return;
    }

    int shortest = 0;
    size_t shortest_len = SIZE_MAX;
    for (size_t k = 0; k < edges.size(); ++k) {
        int i = fast_mapping_local_index[edges[k].from_u];
        size_t len = edges[k].offsets[i + 1] - edges[k].offsets[i];
        if (len < shortest_len) {
            shortest_len = len;
            shortest = static_cast<int>(k);
        }
    }
    const BackwardEdge& first = edges[shortest];
    int first_i = fast_mapping_local_index[first.from_u];
    out.assign(first.targets.begin() + first.offsets[first_i],
               first.targets.begin() + first.offsets[first_i + 1]);

    size_t size = out.size();
    for (size_t k = 0; k < edges.size() && size > 0; ++k) {
        if (static_cast<int>(k) == shortest) continue;
        int i = fast_mapping_local_index[edges[k].from_u];
        const int* list = reinterpret_cast<const int*>(edges[k].targets.data()) + edges[k].offsets[i];
        size_t len = edges[k].offsets[i + 1] - edges[k].offsets[i];
        size = intersectSorted(out.data(), size, list, len, intersect_buffer.data());
        std::copy(intersect_buffer.begin(), intersect_buffer.begin() + size, out.begin());
    }
    out.resize(size);
}

bool SearchEngine::inducedConflict(int u, int v, int& conflict_u) const {
    for (int u_non_nb : backward_non_neighbors[order_position[u]]) {
        if (data_graph.hasEdge(v, fast_mapping[u_non_nb])) {
            conflict_u = u_non_nb;
            return true;
        }
    }
    return false;
}

void SearchEngine::backtrack(int depth) {
    if (static_cast<size_t>(depth) == order.size()) {
        ++match_count;
        if (!count_only) {
            matches.push_back(fast_mapping);
        }
        return;
    }

    int u = order[depth];
    computeLocalCandidates(depth);
    const auto& cands_u = candidate_lists[u];
    const auto& local = local_candidates[depth];
    for (size_t k = 0; k < local.size(); ++k) {
        int i = local[k];
        int v = cands_u[i];
        if (fast_inverse_mapping[v] != -1) continue;
        int unused;
        if (is_induced && inducedConflict(u, v, unused)) continue;

        fast_mapping[u] = v;
        fast_mapping_local_index[u] = i;
        fast_inverse_mapping[v] = u;

        backtrack(depth + 1);

        fast_mapping[u] = -1;
        fast_mapping_local_index[u] = -1;
        fast_inverse_mapping[v] = -1;
    }
}

int SearchEngine::backtrack_cbj(int depth) {
    if (static_cast<size_t>(depth) == order.size()) {
        ++match_count;
        if (!count_only) {
            matches.push_back(fast_mapping);
        }
        return depth - 1;
    }

    int u = order[depth];
    std::fill(conflict_set[depth].begin(), conflict_set[depth].end(), false);
    long long match_count_before = match_count;

    computeLocalCandidates(depth);
    const auto& cands_u = candidate_lists[u];
    const auto& local = local_candidates[depth];

    if (local.size() < cands_u.size()) {
        for (const auto& edge : backward_edges[depth]) {
            conflict_set[depth][edge.from_u] = true;
        }
    }

    for (size_t k = 0; k < local.size(); ++k) {
        int i = local[k];
        int v = cands_u[i];

        int conflict_u = fast_inverse_mapping[v];
        if (conflict_u == -1 && is_induced) inducedConflict(u, v, conflict_u);
        if (conflict_u != -1) {
            conflict_set[depth][conflict_u] = true;
            continue;
        }

        fast_mapping[u] = v;
        fast_mapping_local_index[u] = i;
        fast_inverse_mapping[v] = u;

        int backjump_depth = backtrack_cbj(depth + 1);

        fast_mapping[u] = -1;
        fast_mapping_local_index[u] = -1;
        fast_inverse_mapping[v] = -1;

        if (backjump_depth >= 0 && backjump_depth < depth) {
            conflict_set[depth][order[backjump_depth]] = true;
        }

        if (backjump_depth < depth - 1) {
            if (backjump_depth >= 0) {
                mergeConflicts(conflict_set[backjump_depth], conflict_set[depth]);
            }
            return backjump_depth;
        }
    }

    int max_conflict_depth = -1;
    for (size_t conf_u = 0; conf_u < conflict_set[depth].size(); ++conf_u) {
        if (!conflict_set[depth][conf_u]) continue;
        int d_conf = order_position[conf_u];
        if (d_conf > max_conflict_depth && d_conf < depth) {
            max_conflict_depth = d_conf;
        }
    }

    int target_depth = (max_conflict_depth != -1) ? max_conflict_depth : (depth - 1);
    if (match_count > match_count_before) {
        target_depth = depth - 1;
    }
    if (target_depth >= 0 && target_depth < static_cast<int>(depth)) {
        mergeConflicts(conflict_set[target_depth], conflict_set[depth]);
    }

    return target_depth;
}