#include "AppGraph.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

AppGraph::~AppGraph() {
    clear();
}

void AppGraph::clear() {
    if (c_graph) {
        delete[] c_graph->offsets;
        delete[] c_graph->nbors;
        delete c_graph;
        c_graph = nullptr;
    }
}

AppGraph::AppGraph(const AppGraph& other) {
    copy_from(other);
}

AppGraph& AppGraph::operator=(const AppGraph& other) {
    if (this != &other) {
        clear();
        copy_from(other);
    }
    return *this;
}

void AppGraph::copy_from(const AppGraph& other) {
    this->labels = other.labels;
    this->original_node_ids = other.original_node_ids;
    this->node_to_idx = other.node_to_idx;
    this->fast_node_to_idx = other.fast_node_to_idx;
    this->fast_labels = other.fast_labels;
    
    if (other.c_graph) {
        this->c_graph = new Escape::CGraph();
        this->c_graph->nVertices = other.c_graph->nVertices;
        this->c_graph->nEdges = other.c_graph->nEdges;

        this->c_graph->offsets = new Escape::EdgeIdx[c_graph->nVertices + 1];
        std::copy(other.c_graph->offsets, other.c_graph->offsets + c_graph->nVertices + 1, this->c_graph->offsets);

        this->c_graph->nbors = new Escape::VertexIdx[c_graph->nEdges];
        std::copy(other.c_graph->nbors, other.c_graph->nbors + c_graph->nEdges, this->c_graph->nbors);
    }
}

bool AppGraph::readFromFile(const std::string& filepath) {
    clear();
    std::ifstream file(filepath);
    if (!file.is_open()) return false;

    std::string line;
    std::vector<std::pair<int, int>> edges;
    std::unordered_set<int> node_set;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        ss >> token;
        if (token == "v") {
            int node_id, label;
            ss >> node_id >> label;
            labels[node_id] = label;
            node_set.insert(node_id);
        } else if (token == "e") {
            int u, v;
            ss >> u >> v;
            edges.emplace_back(u, v);
        }
    }

    original_node_ids.assign(node_set.begin(), node_set.end());
    std::sort(original_node_ids.begin(), original_node_ids.end());

    int max_id = original_node_ids.empty() ? -1 : original_node_ids.back();
    fast_node_to_idx.assign(max_id + 1, -1);
    fast_labels.assign(max_id + 1, -1);

    for (size_t i = 0; i < original_node_ids.size(); ++i) {
        node_to_idx[original_node_ids[i]] = i;
        fast_node_to_idx[original_node_ids[i]] = i;
        fast_labels[original_node_ids[i]] = labels.at(original_node_ids[i]);
    }

    std::vector<std::vector<int>> adj_temp(original_node_ids.size());
    for (const auto& edge : edges) {
        int u_idx = fast_node_to_idx[edge.first];
        int v_idx = fast_node_to_idx[edge.second];
        adj_temp[u_idx].push_back(v_idx);
        adj_temp[v_idx].push_back(u_idx);
    }
    
    c_graph = new Escape::CGraph();
    c_graph->nVertices = original_node_ids.size();
    c_graph->offsets = new Escape::EdgeIdx[c_graph->nVertices + 1];
    
    std::vector<Escape::VertexIdx> nbors_flat;
    c_graph->offsets[0] = 0;
    for (size_t i = 0; i < (size_t)c_graph->nVertices; ++i) {
        std::sort(adj_temp[i].begin(), adj_temp[i].end());
        nbors_flat.insert(nbors_flat.end(), adj_temp[i].begin(), adj_temp[i].end());
        c_graph->offsets[i+1] = nbors_flat.size();
    }

    c_graph->nEdges = nbors_flat.size();
    c_graph->nbors = new Escape::VertexIdx[c_graph->nEdges];
    std::copy(nbors_flat.begin(), nbors_flat.end(), c_graph->nbors);

    return true;
}

int AppGraph::degree(int node_id) const {
    int idx = get_node_idx(node_id);
    if (idx == -1) return 0;
    return c_graph->offsets[idx + 1] - c_graph->offsets[idx];
}

bool AppGraph::hasEdge(int u_id, int v_id) const {
    int u_idx = get_node_idx(u_id);
    int v_idx = get_node_idx(v_id);
    if (u_idx == -1 || v_idx == -1) return false;

    Escape::EdgeIdx start = c_graph->offsets[u_idx];
    Escape::EdgeIdx end = c_graph->offsets[u_idx + 1];
    
    return std::binary_search(c_graph->nbors + start, c_graph->nbors + end, v_idx);
}

std::vector<int> AppGraph::getNeighbors(int node_id) const {
    int u_idx = get_node_idx(node_id);
    if (u_idx == -1) return {};

    std::vector<int> neighbors;
    Escape::EdgeIdx start = c_graph->offsets[u_idx];
    Escape::EdgeIdx end = c_graph->offsets[u_idx + 1];
    neighbors.reserve(end - start);
    for (Escape::EdgeIdx i = start; i < end; ++i) {
        neighbors.push_back(original_node_ids[c_graph->nbors[i]]);
    }
    return neighbors;
}

AppGraph AppGraph::createSubgraph(const AppGraph& original, const std::unordered_set<int>& node_subset) {
    AppGraph subgraph;
    
    for (int node_id : node_subset) {
        int lbl = original.get_label_fast(node_id);
        if (lbl != -1) {
            subgraph.labels[node_id] = lbl;
            subgraph.original_node_ids.push_back(node_id);
        }
    }
    std::sort(subgraph.original_node_ids.begin(), subgraph.original_node_ids.end());

    int max_id = subgraph.original_node_ids.empty() ? -1 : subgraph.original_node_ids.back();
    subgraph.fast_node_to_idx.assign(max_id + 1, -1);
    subgraph.fast_labels.assign(max_id + 1, -1);

    for (size_t i = 0; i < subgraph.original_node_ids.size(); ++i) {
        int node_id = subgraph.original_node_ids[i];
        subgraph.node_to_idx[node_id] = i;
        subgraph.fast_node_to_idx[node_id] = i;
        subgraph.fast_labels[node_id] = subgraph.labels.at(node_id);
    }

    std::vector<std::vector<int>> adj_temp(subgraph.original_node_ids.size());
    for (size_t u_sub_idx = 0; u_sub_idx < subgraph.original_node_ids.size(); ++u_sub_idx) {
        int u_orig_id = subgraph.original_node_ids[u_sub_idx];
        int u_orig_idx = original.get_node_idx(u_orig_id);

        Escape::EdgeIdx start = original.c_graph->offsets[u_orig_idx];
        Escape::EdgeIdx end = original.c_graph->offsets[u_orig_idx + 1];
        
        for (Escape::EdgeIdx i = start; i < end; ++i) {
            int v_orig_idx = original.c_graph->nbors[i];
            int v_orig_id = original.original_node_ids[v_orig_idx];
            int v_sub_idx = subgraph.get_node_idx(v_orig_id);
            if (v_sub_idx != -1) {
                adj_temp[u_sub_idx].push_back(v_sub_idx);
            }
        }
    }
    
    subgraph.c_graph = new Escape::CGraph();
    subgraph.c_graph->nVertices = subgraph.original_node_ids.size();
    subgraph.c_graph->offsets = new Escape::EdgeIdx[subgraph.c_graph->nVertices + 1];
    
    std::vector<Escape::VertexIdx> nbors_flat;
    subgraph.c_graph->offsets[0] = 0;
    for (size_t i = 0; i < (size_t)subgraph.c_graph->nVertices; ++i) {
        std::sort(adj_temp[i].begin(), adj_temp[i].end());
        nbors_flat.insert(nbors_flat.end(), adj_temp[i].begin(), adj_temp[i].end());
        subgraph.c_graph->offsets[i+1] = nbors_flat.size();
    }

    subgraph.c_graph->nEdges = nbors_flat.size();
    subgraph.c_graph->nbors = new Escape::VertexIdx[subgraph.c_graph->nEdges];
    std::copy(nbors_flat.begin(), nbors_flat.end(), subgraph.c_graph->nbors);
    
    return subgraph;
}
