#include "AppGraph.h"
#include <fstream>
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
    this->original_node_ids = other.original_node_ids;
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

namespace {

inline void skipBlanks(const char*& p, const char* end) {
    while (p < end && (*p == ' ' || *p == '\t' || *p == '\r')) ++p;
}

inline void skipLine(const char*& p, const char* end) {
    while (p < end && *p != '\n') ++p;
    if (p < end) ++p;
}

inline int parseInt(const char*& p, const char* end) {
    skipBlanks(p, end);
    bool negative = (p < end && *p == '-');
    if (negative) ++p;
    int value = 0;
    while (p < end && *p >= '0' && *p <= '9') {
        value = value * 10 + (*p - '0');
        ++p;
    }
    return negative ? -value : value;
}

}

bool AppGraph::readFromFile(const std::string& filepath) {
    clear();
    std::ifstream file(filepath, std::ios::binary);
    if (!file.is_open()) return false;
    file.seekg(0, std::ios::end);
    std::string buffer(static_cast<size_t>(file.tellg()), '\0');
    file.seekg(0, std::ios::beg);
    file.read(&buffer[0], buffer.size());

    std::vector<int> vertex_ids, vertex_labels;
    std::vector<std::pair<int, int>> edges;
    const char* p = buffer.data();
    const char* end = p + buffer.size();
    while (p < end) {
        skipBlanks(p, end);
        char type = (p + 1 < end && (p[1] == ' ' || p[1] == '\t')) ? *p : '\0';
        if (type == 'v') {
            ++p;
            int node_id = parseInt(p, end);
            int label = parseInt(p, end);
            vertex_ids.push_back(node_id);
            vertex_labels.push_back(label);
        } else if (type == 'e') {
            ++p;
            int u = parseInt(p, end);
            int v = parseInt(p, end);
            edges.emplace_back(u, v);
        }
        skipLine(p, end);
    }

    build(vertex_ids, vertex_labels, edges);
    return true;
}

bool AppGraph::buildFromLabelsAndEdges(const std::vector<int>& vertex_labels,
                                       const std::vector<std::pair<int, int>>& edges) {
    clear();
    std::vector<int> vertex_ids(vertex_labels.size());
    for (size_t id = 0; id < vertex_ids.size(); ++id) {
        vertex_ids[id] = static_cast<int>(id);
    }
    build(vertex_ids, vertex_labels, edges);
    return true;
}

void AppGraph::build(const std::vector<int>& vertex_ids, const std::vector<int>& vertex_labels,
                     const std::vector<std::pair<int, int>>& edges) {
    original_node_ids = vertex_ids;
    std::sort(original_node_ids.begin(), original_node_ids.end());
    original_node_ids.erase(std::unique(original_node_ids.begin(), original_node_ids.end()),
                            original_node_ids.end());

    int max_id = original_node_ids.empty() ? -1 : original_node_ids.back();
    fast_node_to_idx.assign(max_id + 1, -1);
    fast_labels.assign(max_id + 1, -1);
    for (size_t i = 0; i < original_node_ids.size(); ++i) {
        fast_node_to_idx[original_node_ids[i]] = static_cast<int>(i);
    }
    for (size_t k = 0; k < vertex_ids.size(); ++k) {
        fast_labels[vertex_ids[k]] = vertex_labels[k];
    }

    const size_t n = original_node_ids.size();
    c_graph = new Escape::CGraph();
    c_graph->nVertices = n;
    c_graph->offsets = new Escape::EdgeIdx[n + 1]();
    for (const auto& edge : edges) {
        int u_idx = get_node_idx(edge.first);
        int v_idx = get_node_idx(edge.second);
        if (u_idx == -1 || v_idx == -1) continue;
        ++c_graph->offsets[u_idx + 1];
        ++c_graph->offsets[v_idx + 1];
    }
    for (size_t i = 0; i < n; ++i) {
        c_graph->offsets[i + 1] += c_graph->offsets[i];
    }

    c_graph->nEdges = c_graph->offsets[n];
    c_graph->nbors = new Escape::VertexIdx[c_graph->nEdges];
    std::vector<Escape::EdgeIdx> next(c_graph->offsets, c_graph->offsets + n);
    for (const auto& edge : edges) {
        int u_idx = get_node_idx(edge.first);
        int v_idx = get_node_idx(edge.second);
        if (u_idx == -1 || v_idx == -1) continue;
        c_graph->nbors[next[u_idx]++] = v_idx;
        c_graph->nbors[next[v_idx]++] = u_idx;
    }
    for (size_t i = 0; i < n; ++i) {
        std::sort(c_graph->nbors + c_graph->offsets[i], c_graph->nbors + c_graph->offsets[i + 1]);
    }
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

    subgraph.original_node_ids.reserve(node_subset.size());
    for (int node_id : node_subset) {
        if (original.get_label_fast(node_id) != -1) {
            subgraph.original_node_ids.push_back(node_id);
        }
    }
    std::sort(subgraph.original_node_ids.begin(), subgraph.original_node_ids.end());

    int max_id = subgraph.original_node_ids.empty() ? -1 : subgraph.original_node_ids.back();
    subgraph.fast_node_to_idx.assign(max_id + 1, -1);
    subgraph.fast_labels.assign(max_id + 1, -1);
    for (size_t i = 0; i < subgraph.original_node_ids.size(); ++i) {
        int node_id = subgraph.original_node_ids[i];
        subgraph.fast_node_to_idx[node_id] = static_cast<int>(i);
        subgraph.fast_labels[node_id] = original.get_label_fast(node_id);
    }

    const size_t n = subgraph.original_node_ids.size();
    subgraph.c_graph = new Escape::CGraph();
    subgraph.c_graph->nVertices = n;
    subgraph.c_graph->offsets = new Escape::EdgeIdx[n + 1]();
    auto for_each_kept_neighbor = [&](size_t u_sub_idx, auto&& emit) {
        int u_orig_idx = original.get_node_idx(subgraph.original_node_ids[u_sub_idx]);
        for (Escape::EdgeIdx e = original.c_graph->offsets[u_orig_idx];
             e < original.c_graph->offsets[u_orig_idx + 1]; ++e) {
            int v_sub_idx = subgraph.get_node_idx(original.original_node_ids[original.c_graph->nbors[e]]);
            if (v_sub_idx != -1) emit(v_sub_idx);
        }
    };
    for (size_t u = 0; u < n; ++u) {
        for_each_kept_neighbor(u, [&](int) { ++subgraph.c_graph->offsets[u + 1]; });
    }
    for (size_t u = 0; u < n; ++u) {
        subgraph.c_graph->offsets[u + 1] += subgraph.c_graph->offsets[u];
    }
    subgraph.c_graph->nEdges = subgraph.c_graph->offsets[n];
    subgraph.c_graph->nbors = new Escape::VertexIdx[subgraph.c_graph->nEdges];
    for (size_t u = 0; u < n; ++u) {
        Escape::EdgeIdx pos = subgraph.c_graph->offsets[u];
        for_each_kept_neighbor(u, [&](int v_sub_idx) { subgraph.c_graph->nbors[pos++] = v_sub_idx; });
    }

    return subgraph;
}
