#ifndef APP_GRAPH_H
#define APP_GRAPH_H

#include <string>
#include <vector>
#include <unordered_set>
#include "Escape/Graph.h"

class AppGraph {
public:
    Escape::CGraph* c_graph = nullptr;
    std::vector<int> original_node_ids;     // sorted; position = internal index
    std::vector<int> fast_node_to_idx;      // node ID -> internal index, -1 if absent
    std::vector<int> fast_labels;           // node ID -> label, -1 if absent

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
    // Build from an in-memory graph whose vertex IDs are 0..vertex_labels.size()-1
    // (each undirected edge listed once), e.g. one another tool has already loaded.
    bool buildFromLabelsAndEdges(const std::vector<int>& vertex_labels,
                                 const std::vector<std::pair<int, int>>& edges);
    int degree(int node_id) const;
    bool hasEdge(int u, int v) const;
    std::vector<int> getNeighbors(int node_id) const;
    
    static AppGraph createSubgraph(const AppGraph& original, const std::unordered_set<int>& node_subset);

private:
    void build(const std::vector<int>& vertex_ids, const std::vector<int>& vertex_labels,
               const std::vector<std::pair<int, int>>& edges);
    void clear();
    void copy_from(const AppGraph& other);
};

#endif