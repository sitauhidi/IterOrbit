/// @file
/// The graph type shared by the pattern, the data graph and its induced subgraphs.

#ifndef APP_GRAPH_H
#define APP_GRAPH_H

#include <memory>
#include <string>
#include <vector>
#include <unordered_set>
#include "Escape/Graph.h"

/// A vertex-labelled, undirected graph: the pattern, the data graph, or an induced subgraph of
/// the data graph.
///
/// Vertices keep the IDs of the input file ("node IDs"). Internally they are numbered
/// 0..n-1 in ascending ID order ("internal index"), which is the numbering of `c_graph` and of
/// every per-vertex array produced from it (e.g. OrbitCounts). A subgraph keeps the node IDs of
/// the graph it was taken from, so candidate sets, which hold node IDs, stay valid across
/// subgraphs.
/// Per-label lookups on a data graph, built once at load (AppGraph::buildLabelIndex), as the SOTA
/// tools' graph loader builds its neighbour-label frequencies.
struct LabelIndex {
    /// Largest label, -1 if the graph is empty.
    int max_label = -1;
    /// Vertex v has nlf_count[i] neighbours of label nlf_label[i], for i in
    /// [nlf_offsets[v], nlf_offsets[v + 1]), labels increasing; labels it has no neighbour of are absent.
    std::vector<size_t> nlf_offsets;
    std::vector<int> nlf_label, nlf_count;
};

class AppGraph {
public:
    /// Adjacency in CSR form, indexed by internal index; neighbour lists are sorted. Owned.
    Escape::CGraph* c_graph = nullptr;
    /// Node IDs in ascending order; the position of an ID is its internal index.
    std::vector<int> original_node_ids;
    /// Node ID -> internal index, or -1 if the vertex is not in this graph.
    std::vector<int> fast_node_to_idx;
    /// Node ID -> label, or -1 if the vertex is not in this graph.
    std::vector<int> fast_labels;
    /// Built by buildLabelIndex() (the data graph's); null otherwise. Shared by copies.
    std::shared_ptr<const LabelIndex> label_index;

    /// Builds label_index.
    void buildLabelIndex();
    AppGraph() = default;
    ~AppGraph();
    AppGraph(const AppGraph& other);
    AppGraph& operator=(const AppGraph& other);
    AppGraph(AppGraph&& other) noexcept;
    AppGraph& operator=(AppGraph&& other) noexcept;

    /// Internal index of a node ID, or -1 if the vertex is not in this graph.
    inline int get_node_idx(int node_id) const {
        if (node_id >= 0 && static_cast<size_t>(node_id) < fast_node_to_idx.size()) {
            return fast_node_to_idx[node_id];
        }
        return -1;
    }

    /// Label of a node ID, or -1 if the vertex is not in this graph.
    inline int get_label_fast(int node_id) const {
        if (node_id >= 0 && static_cast<size_t>(node_id) < fast_labels.size()) {
            return fast_labels[node_id];
        }
        return -1;
    }

    /// Reads a graph in the "t / v / e" text format: a line `v <id> <label> [...]` per vertex
    /// and a line `e <u> <v> [...]` per undirected edge; other lines are ignored.
    /// Returns false if the file cannot be opened.
    bool readFromFile(const std::string& filepath);

    /// Builds the graph from memory. Vertex IDs are 0..vertex_labels.size()-1 and each
    /// undirected edge is listed once.
    bool buildFromLabelsAndEdges(const std::vector<int>& vertex_labels,
                                 const std::vector<std::pair<int, int>>& edges);

    /// Degree of a node ID in this graph (0 if absent).
    int degree(int node_id) const;
    /// Whether node IDs u and v are adjacent in this graph.
    bool hasEdge(int u, int v) const;
    /// Neighbours of a node ID, as node IDs.
    std::vector<int> getNeighbors(int node_id) const;

    /// The subgraph of `original` induced by `node_subset` (node IDs): those vertices and
    /// every edge of `original` between two of them. Node IDs and labels are kept.
    static AppGraph createSubgraph(const AppGraph& original, const std::unordered_set<int>& node_subset);
    /// The same, for node IDs given sorted and without duplicates.
    static AppGraph createSubgraph(const AppGraph& original, const std::vector<int>& sorted_ids);

private:
    void build(const std::vector<int>& vertex_ids, const std::vector<int>& vertex_labels,
               const std::vector<std::pair<int, int>>& edges);
    void clear();
    void copy_from(const AppGraph& other);
};

#endif
