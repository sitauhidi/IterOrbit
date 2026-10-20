/// @file
/// Iterative orbit filtering: candidate sets and Algorithm 1 (IterOrbitFilter) of the paper.

#ifndef FILTER_ENGINE_H
#define FILTER_ENGINE_H

#include "AppGraph.h"
#include "OrbitCounter.h"
#include <cstdint>
#include <string>

/// A candidate assignment: pattern node ID -> the data node IDs it may be mapped to.
using CandidateSets = std::unordered_map<int, std::vector<int>>;

/// Iterative orbit filtering (Algorithm 1 of the paper, IterOrbitFilter).
///
/// Starting from the label-and-degree filter (LDF), each round applies the neighbourhood label
/// filter (NLF) and the semi-join (to arc consistency), induces the data graph on the
/// surviving candidates, applies the orbit filter on that graph, and induces the data graph
/// again on its survivors for the next round. Candidate sets only ever shrink, and no round
/// removes the image of an embedding.
class FilterEngine {
public:
    /// @param data                 the data graph G_d
    /// @param pattern              the pattern graph G_p
    /// @param g_size               graphlet size for orbit counting: 3, 4 or 5
    /// @param use_full_graph_mode  count orbits on G_d instead of the candidate-induced
    ///                             subgraph; forces a single round
    /// @param max_iters            I_max: the maximum number of rounds; 0 iterates until a
    ///                             round removes no candidate
    /// @param print_stats          print candidate statistics and timings after each stage
    /// @param enable_orbit_throttle     skip the orbit stage in rounds where the cheaper
    ///                             filters are still removing candidates but the sets have
    ///                             shrunk by less than `orbit_throttle_threshold` since it
    ///                             last ran; it always runs once more before returning
    /// @param orbit_throttle_threshold  that fraction
    /// @param enable_semi_join     run the semi-join stage
    /// @param arc_consistency      iterate the semi-join until it removes nothing (default);
    ///                             false: one pass per round
    /// @param enable_orbit_filter  run the orbit stage; false leaves contraction with NLF and
    ///                             the semi-join only (ablation)
    /// @param post_orbit_induction induce the data graph again after the orbit stage; false:
    ///                             the next round's NLF and semi-join read the graph the
    ///                             orbit stage counted on
    /// @param edge_thinning        count orbits on the candidate subgraph without the edges no
    ///                             pattern edge can map to (no pattern edge {u, u'} with one end
    ///                             a candidate of u and the other of u'); an embedding uses none
    ///                             of them, so each of its graphlets is still counted at its image
    FilterEngine(const AppGraph& data, const AppGraph& pattern, int g_size, bool use_full_graph_mode,
                 int max_iters = 1, bool print_stats = true,
                 bool enable_orbit_throttle = false, double orbit_throttle_threshold = 0.05,
                 bool enable_semi_join = true, bool arc_consistency = true,
                 bool enable_orbit_filter = true, bool post_orbit_induction = true,
                 bool edge_thinning = false);

    /// Runs the filter. Returns false if some candidate set becomes empty, in which case the
    /// pattern has no embedding in the data graph.
    bool run();

    /// The candidate sets after run().
    const CandidateSets& getCandidateSets() const { return candidate_sets; }
    /// Orbit counts of the pattern, indexed by the pattern's internal index.
    const OrbitCounts& getPatternOrbits() const { return pattern_orbits; }
    /// The data graph induced on the final candidates (the whole data graph in full-graph
    /// mode). It contains the image of every embedding; the search runs on it.
    const AppGraph& getCandidateSubgraph() const { return candidate_subgraph; }

    /// Prints size statistics of the current candidate sets under the given heading (only
    /// when statistics are enabled).
    void countCandidates(const std::string& filter_step_name) const;

private:
    /// LDF, and round 1's NLF on the full data graph, in one pass. Sets `ldf_only_count` to
    /// the number of candidates LDF alone leaves. Returns false if a set is empty.
    bool ldfNlfFilter(size_t& ldf_only_count);
    /// NLF against `neighbor_graph`: keeps v in C(u) only if v's neighbours there cover the
    /// label multiset of u's pattern neighbours. Returns false if a set becomes empty.
    bool nlfFilter(const AppGraph& neighbor_graph);
    /// Semi-join against `neighbor_graph`, to arc consistency unless disabled.
    bool semiJoinFilter(const AppGraph& neighbor_graph);
    bool arcConsistencyFilter(const AppGraph& neighbor_graph);
    /// Induces the candidate subgraph, counts its orbits and keeps v in C(u) only if every
    /// orbit count of v is at least that of u. Returns false if a set becomes empty.
    bool orbitFilter();
    /// Sets `candidate_subgraph` to the data graph induced on the union of the candidate sets.
    void induceCandidateSubgraph();
    /// A graph in CSR form that owns its arrays; `view` points into them.
    struct OwnedCGraph {
        std::vector<Escape::EdgeIdx> offsets;
        std::vector<Escape::VertexIdx> nbors;
        Escape::CGraph view{};
    };
    /// `g` without the edges no pattern edge can map to (see edge_thinning): same vertices and
    /// internal indices, fewer edges.
    void thinEdges(const AppGraph& g, OwnedCGraph& out) const;
    /// Total number of candidates over all pattern vertices.
    size_t getCandidateCount() const;

    const AppGraph& data_graph;
    const AppGraph& pattern_graph;
    int graphlet_size;
    bool use_subgraph;
    int max_iterations;
    bool print_stats;
    bool enable_orbit_throttle;
    double orbit_throttle_threshold;
    bool enable_semi_join;
    bool arc_consistency;
    bool enable_orbit_filter;
    bool post_orbit_induction;
    bool edge_thinning;

    CandidateSets candidate_sets;
    OrbitCounts pattern_orbits;
    /// Orbits some pattern vertex has; the data graph's other orbits are not counted.
    std::vector<bool> needed_orbits;
    /// Per pattern node ID: label -> number of pattern neighbours with that label.
    std::unordered_map<int, std::unordered_map<int, int>> pattern_nlf;
    int max_data_label = 0;
    AppGraph candidate_subgraph;
    /// Candidate count when `candidate_subgraph` was last induced (SIZE_MAX: not yet).
    /// Candidate sets only shrink, so an equal count means the same sets and the same graph.
    size_t subgraph_candidate_count = SIZE_MAX;
};

#endif
