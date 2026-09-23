#include "OrderEngine.h"
#include <unordered_set>
#include <limits>

OrderEngine::OrderEngine(const AppGraph& pattern, const OrbitCounts& orbits)
    : pattern_graph(pattern), pattern_orbits(orbits) {}

void OrderEngine::run() {
    order.clear();
    pivot.clear();

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