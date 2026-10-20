/// @file
/// Iterative orbit filtering: Algorithm 1 (IterOrbitFilter) and the filter operators of
/// Section "Methodology" of the paper. See FilterEngine.h.

#include "FilterEngine.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <chrono>
#include <cstdint>
#include <limits>

FilterEngine::FilterEngine(const AppGraph& data, const AppGraph& pattern, int g_size, bool use_full_graph_mode,
                            int max_iters, bool print_stats_,
                            bool enable_orbit_throttle_, double orbit_throttle_threshold_,
                            bool enable_semi_join_, bool arc_consistency_,
                            bool enable_orbit_filter_, bool post_orbit_induction_,
                            bool edge_thinning_)
    : data_graph(data),
      pattern_graph(pattern),
      graphlet_size(g_size),
      use_subgraph(!use_full_graph_mode),
      max_iterations(!use_full_graph_mode ? max_iters : 1),
      print_stats(print_stats_),
      enable_orbit_throttle(enable_orbit_throttle_),
      orbit_throttle_threshold(orbit_throttle_threshold_),
      enable_semi_join(enable_semi_join_),
      arc_consistency(arc_consistency_),
      enable_orbit_filter(enable_orbit_filter_),
      post_orbit_induction(post_orbit_induction_),
      edge_thinning(edge_thinning_)
{}

// G_d^(t): the subgraph of the data graph induced by the union of the candidate sets. Skipped when
// no candidate has been removed since the last induction (the same union gives the same graph).
// Induced from the previous G_d^(t) when there is one: it contains every candidate (sets only
// shrink), and the subgraph it induces on them is the one the data graph induces. The union is
// marked over the source graph's internal indices, which are in node-ID order, so it comes out
// sorted without sorting (candidate sets overlap heavily on unlabelled graphs).
void FilterEngine::induceCandidateSubgraph() {
    const size_t count = getCandidateCount();
    if (count == subgraph_candidate_count) return;
    const bool from_previous = subgraph_candidate_count != SIZE_MAX && candidate_subgraph.c_graph;
    const AppGraph& source = from_previous ? candidate_subgraph : data_graph;
    std::vector<char> in_union(source.original_node_ids.size(), 0);
    for (const auto& pair : candidate_sets) {
        for (int v : pair.second) {
            const int idx = source.get_node_idx(v);
            if (idx != -1) in_union[idx] = 1;
        }
    }
    std::vector<int> candidate_nodes;
    for (size_t idx = 0; idx < in_union.size(); ++idx) {
        if (in_union[idx]) candidate_nodes.push_back(source.original_node_ids[idx]);
    }
    candidate_subgraph = AppGraph::createSubgraph(source, candidate_nodes);
    subgraph_candidate_count = count;
}

size_t FilterEngine::getCandidateCount() const {
    size_t total = 0;
    for (const auto& pair : candidate_sets) {
        total += pair.second.size();
    }
    return total;
}

bool FilterEngine::run() {
    auto start = std::chrono::high_resolution_clock::now();

    pattern_nlf.clear();
    for (int u : pattern_graph.original_node_ids) {
        pattern_nlf[u];  // isolated pattern vertices get an empty (no-requirement) entry too
        for (int nbr : pattern_graph.getNeighbors(u)) {
            pattern_nlf[u][pattern_graph.get_label_fast(nbr)]++;
        }
    }

    max_data_label = 0;
    if (data_graph.label_index) {
        max_data_label = std::max(0, data_graph.label_index->max_label);
    } else {
        for (int label : data_graph.fast_labels) {
            if (label > max_data_label) max_data_label = label;
        }
    }

    // Round 1's NLF runs on the full data graph, so it is applied in the same pass as LDF.
    // ldf_only_count is LDF's count alone: round 1 is compared against it, so that a round 1 in
    // which only NLF removes candidates still counts as having changed something.
    size_t ldf_only_count = 0;
    if (!ldfNlfFilter(ldf_only_count)) {
        countCandidates("After LDF Filter (Failed)");
        return false;
    }
    auto ldf_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ldf_duration = ldf_end - start;
    if (print_stats) {
        std::cout << "\nLDF filtering took " << ldf_duration.count() << " ms." << std::endl;
    }
    countCandidates("After LDF Filter");

    EVOKEOrbitCounter pattern_counter(pattern_graph, graphlet_size);
    pattern_orbits = pattern_counter.count();

    // An orbit no pattern vertex has can't reject any candidate, so the data graph need not count it.
    needed_orbits.clear();
    for (const auto& orbits_u : pattern_orbits) {
        if (needed_orbits.empty()) needed_orbits.assign(orbits_u.size(), false);
        for (size_t o = 0; o < orbits_u.size(); ++o) {
            if (orbits_u[o] > 0) needed_orbits[o] = true;
        }
    }

    bool converged = false;
    bool last_round_ran_orbit = true;
    size_t candidates_at_last_orbit_run = ldf_only_count;

    int iter = 0;
    while (max_iterations == 0 || iter < max_iterations) {
        iter++;
        if (print_stats && max_iterations != 1) {
            std::cout << "\n--- Iteration " << iter << " ---" << std::endl;
        }

        size_t prev_count = (iter == 1) ? ldf_only_count : getCandidateCount();

        const AppGraph& nlf_graph = (use_subgraph && iter > 1) ? candidate_subgraph : data_graph;

        auto nlf_start = std::chrono::high_resolution_clock::now();
        if (iter > 1 && !nlfFilter(nlf_graph)) {
            countCandidates("After NLF Filter (Failed)");
            return false;
        }
        auto nlf_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> nlf_duration = nlf_end - nlf_start;
        if (print_stats) {
            std::cout << "NLF filtering took " << nlf_duration.count() << " ms." << std::endl;
        }
        countCandidates("After NLF Filter");

        if (enable_semi_join) {
            auto sj_start = std::chrono::high_resolution_clock::now();
            if (!semiJoinFilter(nlf_graph)) {
                countCandidates("After Semi-Join Filter (Failed)");
                return false;
            }
            auto sj_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> sj_duration = sj_end - sj_start;
            if (print_stats) {
                std::cout << "Semi-Join filtering took " << sj_duration.count() << " ms." << std::endl;
            }
            countCandidates("After Semi-Join Filter");
        }

        // Without the orbit stage (ablation) the branch below still induces G_d^(t), so later
        // rounds' NLF and semi-join run on the contracted graph.
        bool run_orbit_this_round = enable_orbit_filter;
        if (enable_orbit_filter && enable_orbit_throttle && iter > 1) {
            size_t count_before_orbit = getCandidateCount();
            bool cheap_filters_still_progressing = (count_before_orbit < prev_count);
            double shrink_since_last_orbit =
                (candidates_at_last_orbit_run > 0)
                    ? 1.0 - static_cast<double>(count_before_orbit) /
                                static_cast<double>(candidates_at_last_orbit_run)
                    : 1.0;
            if (cheap_filters_still_progressing && shrink_since_last_orbit < orbit_throttle_threshold) {
                run_orbit_this_round = false;
            }
        }

        if (run_orbit_this_round) {
            auto orbit_start = std::chrono::high_resolution_clock::now();
            if (!orbitFilter()) {
                countCandidates("After Orbit Filter (Failed)");
                return false;
            }
            auto orbit_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> orbit_duration = orbit_end - orbit_start;
            if (print_stats) {
                std::cout << "Orbit filtering took " << orbit_duration.count() << " ms." << std::endl;
            }
            countCandidates("After Orbit Filter");
            candidates_at_last_orbit_run = getCandidateCount();
            // Second induction: the next round's NLF and arc consistency run on the graph the
            // orbit stage has just contracted, not on the one it counted orbits on.
            if (use_subgraph && post_orbit_induction) induceCandidateSubgraph();
        } else if (use_subgraph) {
            induceCandidateSubgraph();
        }
        last_round_ran_orbit = run_orbit_this_round;

        size_t curr_count = getCandidateCount();
        if ((run_orbit_this_round || !enable_orbit_filter) && curr_count == prev_count) {
            converged = true;
            if (print_stats && max_iterations != 1) {
                std::cout << "\nFiltering converged after " << iter << " iterations." << std::endl;
            }
            break;
        }
    }

    if (enable_orbit_filter && enable_orbit_throttle && !last_round_ran_orbit) {
        auto orbit_start = std::chrono::high_resolution_clock::now();
        if (!orbitFilter()) {
            countCandidates("After Orbit Filter (Failed)");
            return false;
        }
        auto orbit_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> orbit_duration = orbit_end - orbit_start;
        if (print_stats) {
            std::cout << "Orbit filtering took " << orbit_duration.count() << " ms." << std::endl;
        }
        countCandidates("After Orbit Filter");
    }

    if (use_subgraph && !converged) induceCandidateSubgraph();

    return true;
}

// LDF (label, degree) and NLF on the full data graph: each data vertex's neighbour-label counts
// (only labels the pattern requires) are read once and shared by every same-label pattern vertex.
// One pass over the data vertices in ID order; a vertex's counts come from the data graph's label
// index (built at load, as the SOTA tools' loaders build theirs) or, without it, from its adjacency.
// Either way the same candidates, in the same (node-ID) order, as LDF followed by nlfFilter(data_graph).
bool FilterEngine::ldfNlfFilter(size_t& ldf_only_count) {
    struct PatternVertex {
        int u;
        int degree;
        bool unsatisfiable;                        // requires a label the data graph lacks
        std::vector<std::pair<int, int>> needs;    // (label slot, minimum neighbour count)
        std::vector<int>* out = nullptr;           // candidate_sets[u]
    };
    std::vector<std::vector<PatternVertex>> pattern_by_label(max_data_label + 1);
    std::vector<int> slot_of_label(max_data_label + 1, -1);
    int num_slots = 0;
    int max_requirement = 0;
    for (int u : pattern_graph.original_node_ids) {
        int label_u = pattern_graph.get_label_fast(u);
        if (label_u < 0 || label_u > max_data_label) continue;
        PatternVertex pv{u, pattern_graph.degree(u), false, {}};
        auto it = pattern_nlf.find(u);
        if (it != pattern_nlf.end()) {
            for (const auto& req : it->second) {
                if (req.first < 0 || req.first > max_data_label) {
                    pv.unsatisfiable = true;
                    continue;
                }
                if (slot_of_label[req.first] == -1) slot_of_label[req.first] = num_slots++;
                pv.needs.emplace_back(slot_of_label[req.first], req.second);
                max_requirement = std::max(max_requirement, req.second);
            }
        }
        pattern_by_label[label_u].push_back(std::move(pv));
    }
    const int cap = std::min(max_requirement, 65535);

    candidate_sets.clear();
    for (int u : pattern_graph.original_node_ids) candidate_sets[u];
    for (auto& group : pattern_by_label) {
        for (PatternVertex& pv : group) pv.out = &candidate_sets[pv.u];
    }
    ldf_only_count = 0;

    const Escape::CGraph* g = data_graph.c_graph;
    const size_t n = g ? g->nVertices : 0;
    const LabelIndex* index = data_graph.label_index.get();
    std::vector<int16_t> slot_of_vertex(index ? 0 : n, -1);   // only needed without the index
    for (size_t idx = 0; idx < slot_of_vertex.size(); ++idx) {
        int lbl = data_graph.get_label_fast(data_graph.original_node_ids[idx]);
        if (lbl >= 0 && lbl <= max_data_label) slot_of_vertex[idx] = static_cast<int16_t>(slot_of_label[lbl]);
    }

    std::vector<uint16_t> counts(num_slots);
    for (size_t idx = 0; idx < n; ++idx) {
        int v = data_graph.original_node_ids[idx];
        int lbl = data_graph.get_label_fast(v);
        if (lbl < 0 || lbl > max_data_label || pattern_by_label[lbl].empty()) continue;
        const Escape::EdgeIdx begin = g->offsets[idx], end = g->offsets[idx + 1];
        const int deg = static_cast<int>(end - begin);

        bool counted = false;
        for (const PatternVertex& pv : pattern_by_label[lbl]) {
            if (deg < pv.degree) continue;
            ++ldf_only_count;
            if (pv.unsatisfiable) continue;
            if (!counted) {
                std::fill(counts.begin(), counts.end(), 0);
                if (index) {
                    for (size_t t = index->nlf_offsets[idx]; t < index->nlf_offsets[idx + 1]; ++t) {
                        const int l = index->nlf_label[t];
                        const int slot = (l >= 0 && l <= max_data_label) ? slot_of_label[l] : -1;
                        if (slot >= 0) counts[slot] = static_cast<uint16_t>(std::min(index->nlf_count[t], cap));
                    }
                } else {
                    for (Escape::EdgeIdx e = begin; e < end; ++e) {
                        int slot = slot_of_vertex[g->nbors[e]];
                        if (slot >= 0 && counts[slot] < cap) ++counts[slot];
                    }
                }
                counted = true;
            }
            bool is_valid = true;
            for (const auto& req : pv.needs) {
                if (counts[req.first] < req.second) {
                    is_valid = false;
                    break;
                }
            }
            if (is_valid) pv.out->push_back(v);
        }
    }

    for (const auto& pair : candidate_sets) {
        if (pair.second.empty()) return false;
    }
    return true;
}

bool FilterEngine::nlfFilter(const AppGraph& neighbor_graph) {
    // Only labels some pattern vertex requires matter; give each a slot, and count each
    // candidate's neighbours per slot once per call (shared by every pattern vertex it's a
    // candidate of), capped at the largest requirement since counts above it never reject.
    std::vector<int> slot_of_label(max_data_label + 1, -1);
    int num_slots = 0;
    int max_requirement = 0;
    for (const auto& u_entry : pattern_nlf) {
        for (const auto& req : u_entry.second) {
            if (req.first <= max_data_label && slot_of_label[req.first] == -1) {
                slot_of_label[req.first] = num_slots++;
            }
            max_requirement = std::max(max_requirement, req.second);
        }
    }
    const int cap = std::min(max_requirement, 65535);

    const Escape::CGraph* g = neighbor_graph.c_graph;
    const size_t n = g ? g->nVertices : 0;
    std::vector<int16_t> slot_of_vertex(n, -1);
    if (num_slots > 0) {
        for (size_t idx = 0; idx < n; ++idx) {
            int lbl = neighbor_graph.get_label_fast(neighbor_graph.original_node_ids[idx]);
            if (lbl >= 0 && lbl <= max_data_label) slot_of_vertex[idx] = static_cast<int16_t>(slot_of_label[lbl]);
        }
    }

    std::vector<int> position(n, -1);
    std::vector<uint16_t> counts;
    auto countsOf = [&](int v_idx) -> const uint16_t* {
        if (position[v_idx] == -1) {
            position[v_idx] = static_cast<int>(counts.size() / std::max(num_slots, 1));
            size_t base = counts.size();
            counts.resize(base + num_slots, 0);
            for (Escape::EdgeIdx e = g->offsets[v_idx]; e < g->offsets[v_idx + 1]; ++e) {
                int slot = slot_of_vertex[g->nbors[e]];
                if (slot >= 0 && counts[base + slot] < cap) ++counts[base + slot];
            }
        }
        return counts.data() + static_cast<size_t>(position[v_idx]) * num_slots;
    };

    for (auto& pair : candidate_sets) {
        auto& candidates_v = pair.second;
        std::vector<std::pair<int, int>> requirements;
        bool unsatisfiable = false;
        for (const auto& req : pattern_nlf.at(pair.first)) {
            if (req.first > max_data_label) {
                unsatisfiable = true;
            } else {
                requirements.emplace_back(slot_of_label[req.first], req.second);
            }
        }

        size_t write = 0;
        if (!unsatisfiable) {
            for (size_t read = 0; read < candidates_v.size(); ++read) {
                int v = candidates_v[read];
                int v_idx = neighbor_graph.get_node_idx(v);
                if (v_idx == -1) continue;
                const uint16_t* v_counts = countsOf(v_idx);
                bool is_valid = true;
                for (const auto& req : requirements) {
                    if (v_counts[req.first] < req.second) {
                        is_valid = false;
                        break;
                    }
                }
                if (is_valid) candidates_v[write++] = v;
            }
        }

        candidates_v.resize(write);
        if (candidates_v.empty()) return false;
    }
    return true;
}

// Semi-join (edge consistency): v stays in C(u) only if, for every pattern neighbour u' of u, some
// neighbour of v is in C(u'). By default (arc_consistency) it runs to a fixpoint (AC-3 over
// pattern vertices): removals update the candidate bits at once, and every pattern neighbour of
// a pattern vertex whose set shrank is re-checked, until no set changes. Without it
// (--no-arc-consistency) it is one pass against the candidate sets as they were at its start.
bool FilterEngine::semiJoinFilter(const AppGraph& neighbor_graph) {
    if (arc_consistency) return arcConsistencyFilter(neighbor_graph);
    int num_words = (neighbor_graph.c_graph && neighbor_graph.c_graph->nVertices > 0)
                         ? static_cast<int>((neighbor_graph.c_graph->nVertices - 1) / 64) + 1
                         : 1;

    std::unordered_map<int, std::vector<uint64_t>> candidate_bits;
    candidate_bits.reserve(candidate_sets.size());
    for (const auto& pair : candidate_sets) {
        std::vector<uint64_t> bits(num_words, 0ULL);
        for (int v : pair.second) {
            int v_idx = neighbor_graph.get_node_idx(v);
            if (v_idx != -1) {
                bits[v_idx / 64] |= (1ULL << (v_idx % 64));
            }
        }
        candidate_bits[pair.first] = std::move(bits);
    }

    for (auto& pair : candidate_sets) {
        int u = pair.first;
        auto& candidates_v = pair.second;
        const auto u_neighbors = pattern_graph.getNeighbors(u);

        size_t write = 0;
        for (size_t read = 0; read < candidates_v.size(); ++read) {
            int v = candidates_v[read];
            int v_idx = neighbor_graph.get_node_idx(v);
            bool keeps_all_edges = (v_idx != -1);

            if (keeps_all_edges) {
                Escape::EdgeIdx start = neighbor_graph.c_graph->offsets[v_idx];
                Escape::EdgeIdx end = neighbor_graph.c_graph->offsets[v_idx + 1];

                for (int u_prime : u_neighbors) {
                    auto lookup_it = candidate_bits.find(u_prime);
                    if (lookup_it == candidate_bits.end()) {
                        keeps_all_edges = false;
                        break;
                    }
                    const auto& target_bits = lookup_it->second;

                    bool has_supporting_neighbor = false;
                    for (Escape::EdgeIdx i = start; i < end; ++i) {
                        int nbr_v_idx = neighbor_graph.c_graph->nbors[i];
                        if ((target_bits[nbr_v_idx / 64] >> (nbr_v_idx % 64)) & 1ULL) {
                            has_supporting_neighbor = true;
                            break;
                        }
                    }
                    if (!has_supporting_neighbor) {
                        keeps_all_edges = false;
                        break;
                    }
                }
            }

            if (keeps_all_edges) {
                candidates_v[write++] = v;
            }
        }

        candidates_v.resize(write);
        if (candidates_v.empty()) return false;
    }
    return true;
}

// AC-3 over pattern vertices. A re-check of u tests only the pattern neighbours u' whose set has
// shrunk since u was last checked: a support in an unchanged set is still there. The same fixpoint
// as re-checking every neighbour.
bool FilterEngine::arcConsistencyFilter(const AppGraph& neighbor_graph) {
    const Escape::CGraph* g = neighbor_graph.c_graph;
    const size_t num_words = (g && g->nVertices > 0) ? static_cast<size_t>((g->nVertices - 1) / 64) + 1 : 1;

    // Pattern vertices by position, their candidate bits and their pattern neighbours' positions.
    std::vector<int> pattern_ids;
    std::unordered_map<int, int> position;
    for (const auto& pair : candidate_sets) {
        position[pair.first] = static_cast<int>(pattern_ids.size());
        pattern_ids.push_back(pair.first);
    }
    const size_t k = pattern_ids.size();
    std::vector<std::vector<uint64_t>> bits(k, std::vector<uint64_t>(num_words, 0ULL));
    std::vector<std::vector<int>> pattern_nbrs(k);
    std::vector<std::vector<int>*> cands(k);
    for (size_t p = 0; p < k; ++p) {
        cands[p] = &candidate_sets[pattern_ids[p]];
        for (int v : *cands[p]) {
            const int v_idx = neighbor_graph.get_node_idx(v);
            if (v_idx != -1) bits[p][v_idx / 64] |= 1ULL << (v_idx % 64);
        }
        for (int u_prime : pattern_graph.getNeighbors(pattern_ids[p])) {
            auto it = position.find(u_prime);
            // A pattern neighbour without a candidate set supports nothing.
            pattern_nbrs[p].push_back(it == position.end() ? -1 : it->second);
        }
    }
    // back[p][n]: the position of p among the pattern neighbours of pattern_nbrs[p][n]
    std::vector<std::vector<int>> back(k);
    for (size_t p = 0; p < k; ++p) {
        for (int q : pattern_nbrs[p]) {
            int pos = -1;
            if (q >= 0) {
                for (size_t m = 0; m < pattern_nbrs[q].size(); ++m) {
                    if (pattern_nbrs[q][m] == static_cast<int>(p)) { pos = static_cast<int>(m); break; }
                }
            }
            back[p].push_back(pos);
        }
    }
    // dirty[p][n]: C(pattern_nbrs[p][n]) shrank since p was last checked.
    std::vector<std::vector<char>> dirty(k);
    for (size_t p = 0; p < k; ++p) dirty[p].assign(pattern_nbrs[p].size(), 1);

    std::vector<int> queue(k);
    std::iota(queue.begin(), queue.end(), 0);
    std::vector<char> queued(k, 1);
    std::vector<int> check;
    size_t head = 0;
    while (head < queue.size()) {
        const int p = queue[head++];
        queued[p] = 0;
        const size_t d = pattern_nbrs[p].size();
        check.clear();
        for (size_t n = 0; n < d; ++n) {
            if (dirty[p][n]) { check.push_back(static_cast<int>(n)); dirty[p][n] = 0; }
        }
        auto& candidates_v = *cands[p];
        size_t write = 0;
        for (size_t read = 0; read < candidates_v.size(); ++read) {
            const int v = candidates_v[read];
            const int v_idx = neighbor_graph.get_node_idx(v);
            bool supported = v_idx != -1;
            for (size_t c = 0; supported && c < check.size(); ++c) {
                const int n = check[c];
                const int q = pattern_nbrs[p][n];
                supported = false;
                if (q < 0) break;
                const auto& target = bits[q];
                for (Escape::EdgeIdx e = g->offsets[v_idx]; e < g->offsets[v_idx + 1]; ++e) {
                    const Escape::VertexIdx w = g->nbors[e];
                    if ((target[w / 64] >> (w % 64)) & 1ULL) { supported = true; break; }
                }
            }
            if (supported) {
                candidates_v[write++] = v;
            } else if (v_idx != -1) {
                bits[p][v_idx / 64] &= ~(1ULL << (v_idx % 64));
            }
        }
        const bool shrank = write < candidates_v.size();
        candidates_v.resize(write);
        if (candidates_v.empty()) return false;
        if (!shrank) continue;
        for (size_t n = 0; n < d; ++n) {
            const int q = pattern_nbrs[p][n];
            if (q < 0) continue;
            if (back[p][n] >= 0) dirty[q][back[p][n]] = 1;
            if (!queued[q]) { queued[q] = 1; queue.push_back(q); }
        }
    }
    return true;
}

bool FilterEngine::orbitFilter() {
    using Clock = std::chrono::high_resolution_clock;
    auto ms_since = [](Clock::time_point start) {
        return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    };
    const AppGraph* graph_for_orbit_counting = nullptr;

    auto subgraph_start = Clock::now();
    if (this->use_subgraph) {
        if (print_stats) {
            std::cout << "\nCreating subgraph for orbit counting..." << std::endl;
        }
        induceCandidateSubgraph();
        graph_for_orbit_counting = &this->candidate_subgraph;
    } else {
        if (print_stats) {
            std::cout << "\nUsing full data graph for orbit counting..." << std::endl;
        }
        graph_for_orbit_counting = &this->data_graph;
        this->candidate_subgraph = this->data_graph;
        subgraph_candidate_count = SIZE_MAX;
    }

    double subgraph_ms = ms_since(subgraph_start);

    auto thin_start = Clock::now();
    OwnedCGraph thinned;
    const Escape::CGraph* counted = graph_for_orbit_counting->c_graph;
    if (edge_thinning && counted) {
        thinEdges(*graph_for_orbit_counting, thinned);
        counted = &thinned.view;
    }
    double thin_ms = ms_since(thin_start);

    auto count_start = Clock::now();
    EVOKEOrbitCounter data_counter(counted, graphlet_size, &needed_orbits);

    // A vertex's clique count (orbit 14: 4-clique, 72: 5-clique) is only ever compared against
    // that orbit's count for the pattern vertices it's a candidate of, so counting it past the
    // largest of those changes nothing.
    std::vector<double> clique_caps[2];
    const int clique_orbits[2] = {14, 72};
    for (int c = 0; c < 2; ++c) {
        const int orbit = clique_orbits[c];
        if (static_cast<int>(needed_orbits.size()) <= orbit || !needed_orbits[orbit] ||
            !graph_for_orbit_counting->c_graph) {
            continue;
        }
        clique_caps[c].assign(graph_for_orbit_counting->c_graph->nVertices, 0);
        for (const auto& pair : candidate_sets) {
            double requirement = pattern_orbits.at(pattern_graph.get_node_idx(pair.first))[orbit];
            if (requirement <= 0) continue;
            for (int v : pair.second) {
                int v_idx = graph_for_orbit_counting->get_node_idx(v);
                if (v_idx != -1) clique_caps[c][v_idx] = std::max(clique_caps[c][v_idx], requirement);
            }
        }
        data_counter.setCliqueVertexCap(orbit, &clique_caps[c]);
    }
    // Only the orbits some pattern vertex has are read back, one column per orbit; a candidate of u
    // is then tested only on u's non-zero orbits (counts are never negative).
    std::vector<int> columns_of;
    std::vector<int> column(needed_orbits.size(), -1);
    for (size_t o = 0; o < needed_orbits.size(); ++o) {
        if (needed_orbits[o]) {
            column[o] = static_cast<int>(columns_of.size());
            columns_of.push_back(static_cast<int>(o));
        }
    }
    OrbitCounts data_columns = data_counter.countColumns(columns_of);
    double count_ms = ms_since(count_start);

    auto compare_start = Clock::now();
    const size_t num_counted = graph_for_orbit_counting->c_graph ? graph_for_orbit_counting->c_graph->nVertices : 0;
    std::vector<std::pair<const long long*, long long>> requirements;
    for (auto& pair : candidate_sets) {
        int u = pair.first;
        const auto& orbit_u = pattern_orbits.at(pattern_graph.get_node_idx(u));
        auto& candidates_v = pair.second;
        requirements.clear();
        for (size_t o = 0; o < orbit_u.size(); ++o) {
            if (orbit_u[o] > 0) requirements.emplace_back(data_columns[column[o]].data(), orbit_u[o]);
        }

        size_t write = 0;
        for (size_t read = 0; read < candidates_v.size(); ++read) {
            int v = candidates_v[read];
            int v_idx = graph_for_orbit_counting->get_node_idx(v);
            if (v_idx == -1 || static_cast<size_t>(v_idx) >= num_counted) continue;

            bool compatible = true;
            for (const auto& req : requirements) {
                if (req.first[v_idx] < req.second) {
                    compatible = false;
                    break;
                }
            }
            if (compatible) {
                candidates_v[write++] = v;
            }
        }
        candidates_v.resize(write);
        if (candidates_v.empty()) return false;
    }

    if (print_stats) {
        const OrbitCountTimings& ot = data_counter.timings();
        std::cout << "Orbit stage breakdown (ms): subgraph " << subgraph_ms << ", thin " << thin_ms
                  << ", relabel+DAG " << ot.relabel_dag_ms << ", 3-orbits " << ot.three_ms
                  << ", 4-closed-form " << ot.four_closed_ms << ", paw/diamond " << ot.four_other_ms
                  << ", 4-cycle " << ot.four_cycle_ms << ", 4-clique " << ot.four_clique_ms
                  << ", copy " << ot.copy_ms;
        if (graphlet_size == 5) {
            std::cout << ", 5-setup " << ot.five_setup_ms << ", 5-tree " << ot.five_tree_ms
                      << ", 5-parallel " << ot.five_parallel_ms;
            for (int s = 0; s < 9; ++s) std::cout << ", 5-sec" << s << " " << ot.five_section_ms[s];
        }
        std::cout << ", count-total " << count_ms
                  << ", compare " << ms_since(compare_start)
                  << " | subgraph |V|=" << graph_for_orbit_counting->original_node_ids.size()
                  << " |E|=" << (graph_for_orbit_counting->c_graph ? graph_for_orbit_counting->c_graph->nEdges / 2 : 0)
                  << " counted |E|=" << (counted ? counted->nEdges / 2 : 0) << std::endl;
    }
    return true;
}

// An edge {a, b} is kept iff some pattern edge {u, u'} has a in C(u) and b in C(u'): with mask(v) the
// pattern vertices v is a candidate of and adj(v) the union of their pattern neighbourhoods, iff
// mask(a) & adj(b) is non-empty (the relation is symmetric).
void FilterEngine::thinEdges(const AppGraph& g, OwnedCGraph& out) const {
    const Escape::CGraph* cg = g.c_graph;
    const size_t n = cg->nVertices;
    std::vector<int> position(pattern_graph.original_node_ids.size(), -1);
    std::vector<int> ids;
    for (const auto& pair : candidate_sets) {
        position[pattern_graph.get_node_idx(pair.first)] = static_cast<int>(ids.size());
        ids.push_back(pair.first);
    }
    const size_t words = (ids.size() + 63) / 64;
    std::vector<uint64_t> pattern_adj(ids.size() * words, 0);
    for (size_t p = 0; p < ids.size(); ++p) {
        for (int u_prime : pattern_graph.getNeighbors(ids[p])) {
            const int q = position[pattern_graph.get_node_idx(u_prime)];
            if (q >= 0) pattern_adj[p * words + q / 64] |= uint64_t{1} << (q % 64);
        }
    }
    std::vector<uint64_t> mask(n * words, 0), adj(n * words, 0);
    for (size_t p = 0; p < ids.size(); ++p) {
        for (int v : candidate_sets.at(ids[p])) {
            const int v_idx = g.get_node_idx(v);
            if (v_idx == -1) continue;
            mask[v_idx * words + p / 64] |= uint64_t{1} << (p % 64);
            for (size_t w = 0; w < words; ++w) adj[v_idx * words + w] |= pattern_adj[p * words + w];
        }
    }
    out.offsets.assign(n + 1, 0);
    out.nbors.clear();
    out.nbors.reserve(cg->nEdges);
    for (size_t a = 0; a < n; ++a) {
        const uint64_t* ma = mask.data() + a * words;
        for (Escape::EdgeIdx e = cg->offsets[a]; e < cg->offsets[a + 1]; ++e) {
            const Escape::VertexIdx b = cg->nbors[e];
            const uint64_t* ab = adj.data() + b * words;
            bool keep = false;
            for (size_t w = 0; w < words && !keep; ++w) keep = (ma[w] & ab[w]) != 0;
            if (keep) out.nbors.push_back(b);
        }
        out.offsets[a + 1] = static_cast<Escape::EdgeIdx>(out.nbors.size());
    }
    out.view.nVertices = static_cast<Escape::VertexIdx>(n);
    out.view.nEdges = static_cast<Escape::EdgeIdx>(out.nbors.size());
    out.view.offsets = out.offsets.data();
    out.view.nbors = out.nbors.data();
}

void FilterEngine::countCandidates(const std::string& filter_step_name) const {
    if (!print_stats) return;

    std::cout << "\n--- Candidate Set Statistics (" << filter_step_name << ") ---" << std::endl;
    if (candidate_sets.empty()) {
        std::cout << "  Status: Empty (Filter likely failed)" << std::endl;
        std::cout << "-------------------------------------------------" << std::endl;
        return;
    }

    std::vector<double> sizes;
    sizes.reserve(candidate_sets.size());
    double min_size = std::numeric_limits<double>::max();
    double max_size = 0.0;

    for (const auto& pair : candidate_sets) {
        double current_size = static_cast<double>(pair.second.size());
        sizes.push_back(current_size);
        if (current_size < min_size) min_size = current_size;
        if (current_size > max_size) max_size = current_size;
    }

    double sum = std::accumulate(sizes.begin(), sizes.end(), 0.0);
    double mean = sum / sizes.size();

    double squared_diff_sum = 0.0;
    for (const double size : sizes) {
        squared_diff_sum += (size - mean) * (size - mean);
    }
    double std_dev = std::sqrt(squared_diff_sum / sizes.size());

    std::cout << "  Avg candidates per vertex: " << mean << std::endl;
    std::cout << "  Std. dev. of candidates: " << std_dev << std::endl;
    std::cout << "  Min candidates for a vertex: " << min_size << std::endl;
    std::cout << "  Max candidates for a vertex: " << max_size << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
}
