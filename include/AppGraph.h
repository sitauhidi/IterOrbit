#ifndef APP_GRAPH_H
#define APP_GRAPH_H

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "Escape/Graph.h"

class AppGraph {
public:
    Escape::CGraph* c_graph = nullptr;
    std::vector<int> original_node_ids;
    std::unordered_map<int, int> node_to_idx;
    std::unordered_map<int, int> labels;
    std::vector<int> fast_node_to_idx;
    std::vector<int> fast_labels;

    AppGraph() = default;
    ~AppGraph();
    AppGraph(const AppGraph& other);
    AppGraph& operator=(const AppGraph& other);

    inline int get_node_idx(int node_id) const {
        if (node_id >= 0 && static_cast<size_t>(node_id) < fast_node_to_idx.size()) {
            return fast_node_to_idx[node_id];
        }
        return -1;
    }

    inline int get_label_fast(int node_id) const {
        if (node_id >= 0 && static_cast<size_t>(node_id) < fast_labels.size()) {
            return fast_labels[node_id];
        }
        return -1;
    }

    bool readFromFile(const std::string& filepath);
    int degree(int node_id) const;
    bool hasEdge(int u, int v) const;
    std::vector<int> getNeighbors(int node_id) const;
    
    static AppGraph createSubgraph(const AppGraph& original, const std::unordered_set<int>& node_subset);

private:
    void clear();
    void copy_from(const AppGraph& other);
};

#endif