/// @file
/// Matching orders (Section "Ordering Strategy" of the paper). See OrderEngine.h.

#include "OrderEngine.h"
#include <algorithm>
#include <unordered_set>
#include <limits>

OrderEngine::OrderEngine(const AppGraph& pattern, const OrbitCounts& orbits,
                         const CandidateSets* candidates, OrderMode mode)
    : pattern_graph(pattern), pattern_orbits(orbits), candidates(candidates), mode(mode) {}

void OrderEngine::run() {
    order.clear();
    pivot.clear();
    if (mode == OrderMode::Orbit) {
        runOrbit();
    } else {
        runCore(mode == OrderMode::CoreOrbit);
    }
}

void OrderEngine::runCore(bool use_orbits) {
    const std::vector<int>& ids = pattern_graph.original_node_ids;
    const size_t n = ids.size();
    if (n == 0) return;

    std::vector<std::vector<int>> nbrs(n);
    for (size_t i = 0; i < n; ++i) {
        for (int w : pattern_graph.getNeighbors(ids[i])) nbrs[i].push_back(pattern_graph.get_node_idx(w));
    }

    // 2-core: peel degree <= 1 vertices until none is left. Empty for a tree, in which case
    // every vertex counts as core (nothing to put first).
    std::vector<int> deg(n);
    std::vector<int> peel;
    std::vector<char> in_core(n, 1);
    for (size_t i = 0; i < n; ++i) {
        deg[i] = static_cast<int>(nbrs[i].size());
        if (deg[i] <= 1) peel.push_back(static_cast<int>(i));
    }
    while (!peel.empty()) {
        int x = peel.back();
        peel.pop_back();
        if (!in_core[x]) continue;
        in_core[x] = 0;
        for (int w : nbrs[x]) {
            if (in_core[w] && --deg[w] <= 1) peel.push_back(w);
        }
    }
    bool core_empty = std::find(in_core.begin(), in_core.end(), 1) == in_core.end();
    if (core_empty) std::fill(in_core.begin(), in_core.end(), 1);

    static const int kDenseOrbits[] = {3, 8, 10, 11, 12, 13, 14};
    std::vector<double> dense(n, 0.0), cand(n, std::numeric_limits<double>::max());
    for (size_t i = 0; i < n; ++i) {
        const auto& counts = pattern_orbits.at(i);
        for (int o : kDenseOrbits) {
            if (static_cast<size_t>(o) < counts.size()) dense[i] += static_cast<double>(counts[o]);
        }
        if (candidates) {
            auto it = candidates->find(ids[i]);
            if (it != candidates->end()) cand[i] = static_cast<double>(it->second.size());
        }
    }

    // Start: highest dense-orbit score (use_orbits), then fewest candidates per neighbour.
    auto better_start = [&](int a, int b) {
        if (in_core[a] != in_core[b]) return in_core[a] > in_core[b];
        if (use_orbits && dense[a] != dense[b]) return dense[a] > dense[b];
        double ra = cand[a] / std::max<size_t>(1, nbrs[a].size());
        double rb = cand[b] / std::max<size_t>(1, nbrs[b].size());
        if (ra != rb) return ra < rb;
        return a < b;
    };

    std::vector<char> visited(n, 0);
    std::vector<int> backward(n, 0);
    while (order.size() < n) {
        int best = -1;
        bool frontier = false;
        for (size_t i = 0; i < n; ++i) {
            if (!visited[i] && backward[i] > 0) frontier = true;
        }
        for (size_t i = 0; i < n; ++i) {
            int c = static_cast<int>(i);
            if (visited[c]) continue;
            if (frontier && backward[c] == 0) continue;
            if (best == -1) { best = c; continue; }
            if (!frontier) {
                if (better_start(c, best)) best = c;
                continue;
            }
            // core before trees, then most backward neighbours, then (use_orbits) highest
            // dense-orbit score, then fewest candidates
            if (in_core[c] != in_core[best]) { if (in_core[c] > in_core[best]) best = c; continue; }
            if (backward[c] != backward[best]) { if (backward[c] > backward[best]) best = c; continue; }
            if (use_orbits && dense[c] != dense[best]) { if (dense[c] > dense[best]) best = c; continue; }
            if (cand[c] != cand[best]) { if (cand[c] < cand[best]) best = c; continue; }
        }
        visited[best] = 1;
        order.push_back(ids[best]);
        for (int w : nbrs[best]) ++backward[w];
    }

    for (size_t d = 1; d < order.size(); ++d) {
        for (size_t e = 0; e < d; ++e) {
            if (pattern_graph.hasEdge(order[d], order[e])) {
                pivot[order[d]] = order[e];
                break;
            }
        }
    }
}

void OrderEngine::runOrbit() {
    std::unordered_map<int, double> orbit_strength;
    for (int u : pattern_graph.original_node_ids) {
        double score = 0.0;
        for (long long val : pattern_orbits.at(pattern_graph.get_node_idx(u))) {
            score += static_cast<double>(val) * val;
        }
        orbit_strength[u] = score;
    }

    int start_node = -1;
    double max_score = -1.0;
    for (const auto& pair : orbit_strength) {
        if (pair.second > max_score) {
            max_score = pair.second;
            start_node = pair.first;
        }
    }

    if (start_node == -1) {
        if (!pattern_graph.original_node_ids.empty()) {
            order = pattern_graph.original_node_ids;
        }
        return;
    }


    order.push_back(start_node);
    std::unordered_set<int> visited = {start_node};

    while (order.size() < pattern_graph.original_node_ids.size()) {
        int selected_node = -1;
        int max_bn = -1;
        double best_score = -1.0;

        for (int u : pattern_graph.original_node_ids) {
            if (visited.count(u)) continue;

            int backward_neighbors = 0;
            for (int v_ordered : order) {
                if (pattern_graph.hasEdge(u, v_ordered)) {
                    backward_neighbors++;
                }
            }

            if (backward_neighbors > max_bn) {
                max_bn = backward_neighbors;
                selected_node = u;
                best_score = orbit_strength.at(u);
            } else if (backward_neighbors == max_bn) {
                if (orbit_strength.count(u) && (selected_node == -1 || orbit_strength.at(u) > best_score)) {
                    selected_node = u;
                    best_score = orbit_strength.at(u);
                }
            }
        }

        if (selected_node == -1) {
            max_score = -1.0;
            for (int u : pattern_graph.original_node_ids) {
                if (!visited.count(u)) {
                    if (orbit_strength.at(u) > max_score) {
                        max_score = orbit_strength.at(u);
                        selected_node = u;
                    }
                }
            }
        }

        for (int ordered_node : order) {
            if (pattern_graph.hasEdge(selected_node, ordered_node)) {
                pivot[selected_node] = ordered_node;
                break;
            }
        }
        
        order.push_back(selected_node);
        visited.insert(selected_node);
    }
}