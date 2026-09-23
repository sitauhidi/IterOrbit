#include "FilterEngine.h"
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <chrono>
#include <limits>

FilterEngine::FilterEngine(const AppGraph& data, const AppGraph& pattern, int g_size, bool use_full_graph_mode,
                            int max_iters, bool print_stats_,
                            bool enable_orbit_throttle_, double orbit_throttle_threshold_,
                            bool enable_semi_join_)
    : data_graph(data),
      pattern_graph(pattern),
      graphlet_size(g_size),
      use_subgraph(!use_full_graph_mode),
      max_iterations(!use_full_graph_mode ? max_iters : 1),
      print_stats(print_stats_),
      enable_orbit_throttle(enable_orbit_throttle_),
      orbit_throttle_threshold(orbit_throttle_threshold_),
      enable_semi_join(enable_semi_join_)
{}

size_t FilterEngine::getCandidateCount() const {
    size_t total = 0;
    for (const auto& pair : candidate_sets) {
        total += pair.second.size();
    }
    return total;
}

bool FilterEngine::run() {
    auto start = std::chrono::high_resolution_clock::now();

    if (!ldfFilter()) {
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

    pattern_nlf.clear();
    for (int u : pattern_graph.original_node_ids) {
        for (int nbr : pattern_graph.getNeighbors(u)) {
            pattern_nlf[u][pattern_graph.get_label_fast(nbr)]++;
        }
    }

    max_data_label = 0;
    for (const auto& kv : data_graph.labels) {
        if (kv.second > max_data_label) max_data_label = kv.second;
    }

    bool converged = false;
    bool last_round_ran_orbit = true;
    size_t candidates_at_last_orbit_run = getCandidateCount();

    int iter = 0;
    while (max_iterations == 0 || iter < max_iterations) {
        iter++;
        if (print_stats && max_iterations != 1) {
            std::cout << "\n--- Iteration " << iter << " ---" << std::endl;
        }

        size_t prev_count = getCandidateCount();

        const AppGraph& nlf_graph = (use_subgraph && iter > 1) ? candidate_subgraph : data_graph;

        auto nlf_start = std::chrono::high_resolution_clock::now();
        if (!nlfFilter(nlf_graph)) {
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

        bool run_orbit_this_round = true;
        if (enable_orbit_throttle && iter > 1) {
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
        } else if (use_subgraph) {
            std::unordered_set<int> candidate_nodes;
            for (const auto& pair : candidate_sets) {
                candidate_nodes.insert(pair.second.begin(), pair.second.end());
            }
            candidate_subgraph = AppGraph::createSubgraph(data_graph, candidate_nodes);
        }
        last_round_ran_orbit = run_orbit_this_round;

        size_t curr_count = getCandidateCount();
        if (run_orbit_this_round && curr_count == prev_count) {
            converged = true;
            if (print_stats && max_iterations != 1) {
                std::cout << "\nFiltering converged after " << iter << " iterations." << std::endl;
            }
            break;
        }
    }

    if (enable_orbit_throttle && !last_round_ran_orbit) {
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

    if (use_subgraph && !converged) {
        std::unordered_set<int> final_candidate_nodes;
        for (const auto& pair : candidate_sets) {
            final_candidate_nodes.insert(pair.second.begin(), pair.second.end());
        }
        candidate_subgraph = AppGraph::createSubgraph(data_graph, final_candidate_nodes);
    }

    return true;
}

bool FilterEngine::ldfFilter() {
    std::unordered_map<int, std::vector<int>> label_index;
    for (int node : data_graph.original_node_ids) {
        label_index[data_graph.get_label_fast(node)].push_back(node);
    }

    for (int u : pattern_graph.original_node_ids) {
        int label_u = pattern_graph.get_label_fast(u);
        int deg_u = pattern_graph.degree(u);

        if (label_index.find(label_u) == label_index.end()) return false;

        std::vector<int> candidates;
        for (int v : label_index.at(label_u)) {
            if (data_graph.degree(v) >= deg_u) {
                candidates.push_back(v);
            }
        }

        if (candidates.empty()) return false;
        candidate_sets[u] = candidates;
    }
    return true;
}

bool FilterEngine::nlfFilter(const AppGraph& neighbor_graph) {
    std::vector<int> label_counts(max_data_label + 1, 0);

    for (auto& pair : candidate_sets) {
        int u = pair.first;
        auto& candidates_v = pair.second;
        const auto& u_nlf = pattern_nlf.at(u);

        size_t write = 0;
        for (size_t read = 0; read < candidates_v.size(); ++read) {
            int v = candidates_v[read];
            int v_idx = neighbor_graph.get_node_idx(v);
            if (v_idx == -1) continue;

            Escape::EdgeIdx start = neighbor_graph.c_graph->offsets[v_idx];
            Escape::EdgeIdx end = neighbor_graph.c_graph->offsets[v_idx + 1];

            for (Escape::EdgeIdx i = start; i < end; ++i) {
                int nbr_v_idx = neighbor_graph.c_graph->nbors[i];
                int nbr_v_id = neighbor_graph.original_node_ids[nbr_v_idx];
                int lbl = neighbor_graph.get_label_fast(nbr_v_id);
                if (lbl <= max_data_label) {
                    label_counts[lbl]++;
                }
            }

            bool is_valid = true;
            for (const auto& nlf_pair : u_nlf) {
                int target_lbl = nlf_pair.first;
                int req_count = nlf_pair.second;
                int actual_count = (target_lbl <= max_data_label) ? label_counts[target_lbl] : 0;
                if (actual_count < req_count) {
                    is_valid = false;
                    break;
                }
            }

            for (Escape::EdgeIdx i = start; i < end; ++i) {
                int nbr_v_idx = neighbor_graph.c_graph->nbors[i];
                int nbr_v_id = neighbor_graph.original_node_ids[nbr_v_idx];
                int lbl = neighbor_graph.get_label_fast(nbr_v_id);
                if (lbl <= max_data_label) {
                    label_counts[lbl] = 0;
                }
            }

            if (is_valid) {
                candidates_v[write++] = v;
            }
        }

        candidates_v.resize(write);
        if (candidates_v.empty()) return false;
    }
    return true;
}

bool FilterEngine::semiJoinFilter(const AppGraph& neighbor_graph) {
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

bool FilterEngine::orbitFilter() {
    const AppGraph* graph_for_orbit_counting = nullptr;

    if (this->use_subgraph) {
        if (print_stats) {
            std::cout << "\nCreating subgraph for orbit counting..." << std::endl;
        }
        std::unordered_set<int> candidate_nodes;
        for (const auto& pair : candidate_sets) {
            candidate_nodes.insert(pair.second.begin(), pair.second.end());
        }
        this->candidate_subgraph = AppGraph::createSubgraph(data_graph, candidate_nodes);
        graph_for_orbit_counting = &this->candidate_subgraph;
    } else {
        if (print_stats) {
            std::cout << "\nUsing full data graph for orbit counting..." << std::endl;
        }
        graph_for_orbit_counting = &this->data_graph;
        this->candidate_subgraph = this->data_graph;
    }

    EVOKEOrbitCounter data_counter(*graph_for_orbit_counting, graphlet_size);
    OrbitCounts data_orbits = data_counter.count();

    for (auto& pair : candidate_sets) {
        int u = pair.first;
        const auto& orbit_u = pattern_orbits.at(pattern_graph.get_node_idx(u));
        auto& candidates_v = pair.second;

        size_t write = 0;
        for (size_t read = 0; read < candidates_v.size(); ++read) {
            int v = candidates_v[read];
            int v_idx = graph_for_orbit_counting->get_node_idx(v);
            if (v_idx == -1 || static_cast<size_t>(v_idx) >= data_orbits.size()) continue;
            const auto& orbit_v = data_orbits[v_idx];

            bool compatible = true;
            for (size_t i = 0; i < orbit_u.size(); ++i) {
                if (orbit_v[i] < orbit_u[i]) {
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

    return true;
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
