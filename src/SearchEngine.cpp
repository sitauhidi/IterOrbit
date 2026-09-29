#include "SearchEngine.h"
#include <unordered_map>

SearchEngine::SearchEngine(const AppGraph& data, const AppGraph& pattern,
                           const CandidateSets& candidates, const std::vector<int>& o,
                           const std::unordered_map<int, int>& p, bool induced, bool cbj,
                           bool count_only)
    : data_graph(data), pattern_graph(pattern), order(o), pivot(p), is_induced(induced),
      enable_cbj(cbj), count_only(count_only) {

    int max_v_id = 0;
    for (int v_id : data_graph.original_node_ids) {
        if (v_id > max_v_id) max_v_id = v_id;
    }
    fast_mapping.assign(pattern_graph.original_node_ids.size(), -1);
    fast_mapping_local_index.assign(pattern_graph.original_node_ids.size(), -1);
    fast_inverse_mapping.assign(max_v_id + 1, -1);

    size_t num_p_nodes = pattern_graph.original_node_ids.size();
    pattern_neighbors.resize(num_p_nodes);
    pattern_non_neighbors.resize(num_p_nodes);

    for (size_t u = 0; u < num_p_nodes; ++u) {
        int u_id = pattern_graph.original_node_ids[u];
        for (size_t u2 = 0; u2 < num_p_nodes; ++u2) {
            if (u == u2) continue;
            int u2_id = pattern_graph.original_node_ids[u2];
            if (pattern_graph.hasEdge(u_id, u2_id)) {
                pattern_neighbors[u].push_back(u2);
            } else {
                pattern_non_neighbors[u].push_back(u2);
            }
        }
    }

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
        }
    }

    std::vector<std::unordered_map<int, int>> local_index(num_p_nodes);
    for (size_t u = 0; u < num_p_nodes; ++u) {
        local_index[u].reserve(candidate_lists[u].size());
        for (size_t i = 0; i < candidate_lists[u].size(); ++i) {
            local_index[u][candidate_lists[u][i]] = static_cast<int>(i);
        }
    }

    pivot_slot.assign(num_p_nodes, -1);
    for (size_t u = 0; u < num_p_nodes; ++u) {
        auto pivot_it = pivot.find(static_cast<int>(u));
        if (pivot_it == pivot.end()) continue;
        const auto& nbrs = pattern_neighbors[u];
        for (size_t k = 0; k < nbrs.size(); ++k) {
            if (nbrs[k] == pivot_it->second) {
                pivot_slot[u] = static_cast<int>(k);
                break;
            }
        }
    }

    candidate_adj.assign(num_p_nodes, {});
    for (size_t u = 0; u < num_p_nodes; ++u) {
        const auto& nbrs = pattern_neighbors[u];
        candidate_adj[u].resize(nbrs.size());

        for (size_t k = 0; k < nbrs.size(); ++k) {
            int u_nb = nbrs[k];
            int nb_words = static_cast<int>((candidate_lists[u_nb].size() + 63) / 64);
            auto& adj_for_edge = candidate_adj[u][k];
            adj_for_edge.assign(candidate_lists[u].size(), std::vector<uint64_t>(nb_words, 0ULL));

            for (size_t i = 0; i < candidate_lists[u].size(); ++i) {
                int v = candidate_lists[u][i];
                int v_idx = data_graph.get_node_idx(v);
                if (v_idx == -1) continue;

                Escape::EdgeIdx start = data_graph.c_graph->offsets[v_idx];
                Escape::EdgeIdx end = data_graph.c_graph->offsets[v_idx + 1];
                for (Escape::EdgeIdx e = start; e < end; ++e) {
                    int nbr_v_idx = data_graph.c_graph->nbors[e];
                    int nbr_v_id = data_graph.original_node_ids[nbr_v_idx];
                    auto found = local_index[u_nb].find(nbr_v_id);
                    if (found != local_index[u_nb].end()) {
                        int j = found->second;
                        adj_for_edge[i][j / 64] |= (1ULL << (j % 64));
                    }
                }
            }
        }
    }
}

bool SearchEngine::bit_test(const std::vector<uint64_t>& bits, int pos) {
    return (bits[pos / 64] >> (pos % 64)) & 1ULL;
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
    if (enable_cbj) {
        backtrack_cbj(0);
    } else {
        backtrack(0);
    }
}

bool SearchEngine::is_valid(int u, int i, int v) {
    if (fast_inverse_mapping[v] != -1) return false;

    int slot = pivot_slot[u];
    if (slot != -1) {
        int pivot_u = pattern_neighbors[u][slot];
        if (fast_mapping[pivot_u] != -1) {
            if (!bit_test(candidate_adj[u][slot][i], fast_mapping_local_index[pivot_u])) {
                return false;
            }
        }
    }

    const auto& nbrs = pattern_neighbors[u];
    for (size_t k = 0; k < nbrs.size(); ++k) {
        int u_nb = nbrs[k];
        if (fast_mapping[u_nb] != -1) {
            if (!bit_test(candidate_adj[u][k][i], fast_mapping_local_index[u_nb])) {
                return false;
            }
        }
    }

    if (is_induced) {
        for (int u_non_nb : pattern_non_neighbors[u]) {
            int v_non_nb = fast_mapping[u_non_nb];
            if (v_non_nb != -1) {
                if (data_graph.hasEdge(v, v_non_nb)) {
                    return false;
                }
            }
        }
    }

    return true;
}

bool SearchEngine::is_valid_cbj(int u, int i, int v, int depth, int& conflict_u) {
    int current_mapped_u = fast_inverse_mapping[v];
    if (current_mapped_u != -1) {
        conflict_u = current_mapped_u;
        return false;
    }

    int slot = pivot_slot[u];
    if (slot != -1) {
        int pivot_u = pattern_neighbors[u][slot];
        if (fast_mapping[pivot_u] != -1) {
            if (!bit_test(candidate_adj[u][slot][i], fast_mapping_local_index[pivot_u])) {
                conflict_u = pivot_u;
                return false;
            }
        }
    }

    const auto& nbrs = pattern_neighbors[u];
    for (size_t k = 0; k < nbrs.size(); ++k) {
        int u_nb = nbrs[k];
        if (fast_mapping[u_nb] != -1) {
            if (!bit_test(candidate_adj[u][k][i], fast_mapping_local_index[u_nb])) {
                conflict_u = u_nb;
                return false;
            }
        }
    }

    if (is_induced) {
        for (int u_non_nb : pattern_non_neighbors[u]) {
            int v_non_nb = fast_mapping[u_non_nb];
            if (v_non_nb != -1) {
                if (data_graph.hasEdge(v, v_non_nb)) {
                    conflict_u = u_non_nb;
                    return false;
                }
            }
        }
    }

    return true;
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
    const auto& candidates_u = candidate_lists[u];
    for (size_t i = 0; i < candidates_u.size(); ++i) {
        int v = candidates_u[i];

        if (is_valid(u, static_cast<int>(i), v)) {
            fast_mapping[u] = v;
            fast_mapping_local_index[u] = static_cast<int>(i);
            fast_inverse_mapping[v] = u;

            backtrack(depth + 1);

            fast_mapping[u] = -1;
            fast_mapping_local_index[u] = -1;
            fast_inverse_mapping[v] = -1;
        }
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

    const auto& candidates_u = candidate_lists[u];
    for (size_t i = 0; i < candidates_u.size(); ++i) {
        int v = candidates_u[i];

        int conflict_u = -1;
        if (is_valid_cbj(u, static_cast<int>(i), v, depth, conflict_u)) {
            fast_mapping[u] = v;
            fast_mapping_local_index[u] = static_cast<int>(i);
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
        } else {
            if (conflict_u != -1) {
                conflict_set[depth][conflict_u] = true;
            }
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