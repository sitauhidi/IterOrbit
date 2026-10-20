/// @file
/// Enumeration (Section "Enumeration" of the paper): the candidate-space index, the edge-orbit
/// index filter, backtracking with optional conflict-directed backjumping, symmetry-broken
/// counting, and counting the last level without visiting it. See SearchEngine.h.

#include "SearchEngine.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#ifdef __AVX2__
#include <immintrin.h>
#endif

namespace {

#ifdef __AVX2__
// kCompress.idx[m]: the lanes set in the 8-bit mask m, in order, then zeros -- the permutation
// that packs an AVX2 vector's selected lanes to its front.
struct CompressTable {
    alignas(32) int idx[256][8];
    CompressTable() {
        for (int m = 0; m < 256; ++m) {
            int k = 0;
            for (int lane = 0; lane < 8; ++lane) {
                if (m >> lane & 1) idx[m][k++] = lane;
            }
            while (k < 8) idx[m][k++] = 0;
        }
    }
};
const CompressTable kCompress;

// Block-wise merge, 8 x 8 elements per step: each a-block is compared with every rotation of the
// b-block, its matching elements packed to out, and the block with the smaller maximum
// advances (both on a tie). Sorted, duplicate-free inputs give the same sorted output as a
// scalar merge. Writes up to 8 elements past the result, so out needs that much slack.
size_t intersectMergeAVX2(const int* a, size_t na, const int* b, size_t nb, int* out) {
    size_t i = 0, j = 0, count = 0;
    const size_t na8 = na & ~size_t{7}, nb8 = nb & ~size_t{7};
    const __m256i rotate = _mm256_set_epi32(0, 7, 6, 5, 4, 3, 2, 1);
    while (i < na8 && j < nb8) {
        const __m256i va = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a + i));
        __m256i vb = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(b + j));
        __m256i eq = _mm256_cmpeq_epi32(va, vb);
        for (int r = 1; r < 8; ++r) {
            vb = _mm256_permutevar8x32_epi32(vb, rotate);
            eq = _mm256_or_si256(eq, _mm256_cmpeq_epi32(va, vb));
        }
        const int mask = _mm256_movemask_ps(_mm256_castsi256_ps(eq));
        if (mask) {
            const __m256i perm = _mm256_load_si256(reinterpret_cast<const __m256i*>(kCompress.idx[mask]));
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(out + count), _mm256_permutevar8x32_epi32(va, perm));
            count += static_cast<size_t>(__builtin_popcount(mask));
        }
        const int a_max = a[i + 7], b_max = b[j + 7];
        if (a_max <= b_max) i += 8;
        if (b_max <= a_max) j += 8;
    }
    while (i < na && j < nb) {
        if (a[i] < b[j]) {
            ++i;
        } else if (b[j] < a[i]) {
            ++j;
        } else {
            out[count++] = a[i];
            ++i;
            ++j;
        }
    }
    return count;
}
#endif

size_t intersectSorted(const int* a, size_t na, const int* b, size_t nb, int* out) {
    if (na > nb) {
        std::swap(a, b);
        std::swap(na, nb);
    }
    size_t count = 0;
    if (na * 32 < nb) {
        size_t lo = 0;
        for (size_t i = 0; i < na && lo < nb; ++i) {
            size_t step = 1, hi = lo;
            while (hi < nb && b[hi] < a[i]) {
                lo = hi;
                hi += step;
                step <<= 1;
            }
            lo = std::lower_bound(b + lo, b + std::min(hi + 1, nb), a[i]) - b;
            if (lo < nb && b[lo] == a[i]) out[count++] = a[i];
        }
        return count;
    }
#ifdef __AVX2__
    if (na >= 8) return intersectMergeAVX2(a, na, b, nb, out);
#endif
    size_t i = 0, j = 0;
    while (i < na && j < nb) {
        if (a[i] < b[j]) {
            ++i;
        } else if (b[j] < a[i]) {
            ++j;
        } else {
            out[count++] = a[i];
            ++i;
            ++j;
        }
    }
    return count;
}

}

// Per adjacency position p = (x -> y) of g: the triangles containing the edge {x, y}. Each
// triangle is found once from its lowest-ranked vertex (rank = degree, then index), credited to
// its three edges at the positions met, then copied to the reverse positions. About
// sum over edges of min(deg) work: the right choice when most edges are needed.
static std::vector<int64_t> edgeTriangles(const Escape::CGraph* g) {
    const int64_t n = g->nVertices;
    const Escape::EdgeIdx* off = g->offsets;
    const Escape::VertexIdx* nb = g->nbors;
    auto deg = [&](int64_t v) { return off[v + 1] - off[v]; };
    auto higher = [&](int64_t a, int64_t b) { return deg(b) > deg(a) || (deg(b) == deg(a) && b > a); };

    std::vector<int64_t> once(g->nEdges, 0);
    std::vector<Escape::EdgeIdx> pos_from_x(n, -1);
    for (int64_t x = 0; x < n; ++x) {
        for (Escape::EdgeIdx p = off[x]; p < off[x + 1]; ++p) if (higher(x, nb[p])) pos_from_x[nb[p]] = p;
        for (Escape::EdgeIdx pxy = off[x]; pxy < off[x + 1]; ++pxy) {
            const int64_t y = nb[pxy];
            if (!higher(x, y)) continue;
            for (Escape::EdgeIdx pyz = off[y]; pyz < off[y + 1]; ++pyz) {
                const int64_t z = nb[pyz];
                if (!higher(y, z) || pos_from_x[z] < 0) continue;
                ++once[pxy];
                ++once[pyz];
                ++once[pos_from_x[z]];
            }
        }
        for (Escape::EdgeIdx p = off[x]; p < off[x + 1]; ++p) pos_from_x[nb[p]] = -1;
    }
    std::vector<int64_t> tri(g->nEdges);
    for (int64_t x = 0; x < n; ++x) {
        for (Escape::EdgeIdx p = off[x]; p < off[x + 1]; ++p) {
            const int64_t y = nb[p];
            const Escape::EdgeIdx r = std::lower_bound(nb + off[y], nb + off[y + 1], x) - nb;
            tri[p] = once[p] + once[r];
        }
    }
    return tri;
}

// |N(x) & N(y)| in g (the triangles on edge {x, y}), galloping through the longer list when the
// lengths are far apart; and the cost estimate the index filter uses to choose between this and
// edgeTriangles.
static int64_t commonNeighbours(const Escape::CGraph* g, int64_t x, int64_t y) {
    const Escape::VertexIdx* a = g->nbors + g->offsets[x];
    const Escape::VertexIdx* a_end = g->nbors + g->offsets[x + 1];
    const Escape::VertexIdx* b = g->nbors + g->offsets[y];
    const Escape::VertexIdx* b_end = g->nbors + g->offsets[y + 1];
    if (a_end - a > b_end - b) { std::swap(a, b); std::swap(a_end, b_end); }
    int64_t count = 0;
    if ((a_end - a) * 32 < (b_end - b)) {
        for (; a < a_end && b < b_end; ++a) {
            int64_t step = 1;
            while (b + step < b_end && b[step] < *a) step <<= 1;
            b = std::lower_bound(b + (step >> 1), std::min(b + step + 1, b_end), *a);
            if (b < b_end && *b == *a) ++count;
        }
        return count;
    }
    while (a < a_end && b < b_end) {
        if (*a < *b) ++a;
        else if (*b < *a) ++b;
        else { ++count; ++a; ++b; }
    }
    return count;
}

static double intersectionCost(int64_t da, int64_t db) {
    const double lo = static_cast<double>(std::min(da, db)), hi = static_cast<double>(std::max(da, db));
    return 1.0 + (lo * 32 < hi ? lo * (1.0 + std::log2(hi / std::max(lo, 1.0))) : lo + hi);
}

SearchEngine::SearchEngine(const AppGraph& data, const AppGraph& pattern,
                           const CandidateSets& candidates, const std::vector<int>& o,
                           const std::unordered_map<int, int>&,
                           bool induced, bool cbj, bool count_only, bool edge_orbit_filter,
                           const SymmetryBreaking* symmetry)
    : data_graph(data), pattern_graph(pattern), order(o), is_induced(induced),
      enable_cbj(cbj), count_only(count_only) {

    const size_t num_p_nodes = pattern_graph.original_node_ids.size();
    int max_v_id = data_graph.original_node_ids.empty() ? 0 : data_graph.original_node_ids.back();
    for (const auto& pair : candidates) {
        for (int v : pair.second) max_v_id = std::max(max_v_id, v);
    }
    fast_mapping.assign(num_p_nodes, -1);
    fast_mapping_local_index.assign(num_p_nodes, -1);
    fast_inverse_mapping.assign(max_v_id + 1, -1);
    in_use.assign(max_v_id + 1, 0);

    order_position.assign(num_p_nodes, -1);
    for (size_t d = 0; d < order.size(); ++d) {
        order_position[order[d]] = static_cast<int>(d);
    }
    conflict_words = (order.size() + 63) / 64;
    conflict_bits.assign(order.size() * conflict_words, 0);
    backward_masks.assign(order.size() * conflict_words, 0);

    candidate_lists.assign(num_p_nodes, {});
    for (const auto& pair : candidates) {
        int u = pair.first;
        if (u >= 0 && static_cast<size_t>(u) < num_p_nodes) {
            candidate_lists[u] = pair.second;
            std::sort(candidate_lists[u].begin(), candidate_lists[u].end());
        }
    }

    // Edge-orbit filter: an index entry (v', v) along pattern edge {u', u} survives only if the
    // data edge {v', v} lies in at least as many triangles as {u', u} does. An embedding maps
    // each triangle through {u', u} to a distinct one through its image, inside any induced
    // subgraph containing the embedding's image -- such as the candidate subgraph searched here
    // -- so no entry an embedding uses is removed (non-induced counts: sound for induced search
    // too). The data counts come either from one pass over the whole subgraph (edgeTriangles) or
    // per entry by intersecting neighbour lists. Per entry is used while it costs at most ~4x
    // the adjacency scan that built that pattern edge's entries; beyond that the whole-graph
    // pass, computed once for all pattern edges, if it costs at most ~4x the scan for the whole
    // index; otherwise the entries stay unfiltered (always safe). This caps the overhead on easy
    // searches.
    std::vector<int64_t> pattern_tri;
    if (edge_orbit_filter && pattern_graph.c_graph && data_graph.c_graph) {
        pattern_tri = edgeTriangles(pattern_graph.c_graph);
        if (std::none_of(pattern_tri.begin(), pattern_tri.end(), [](int64_t t) { return t > 0; })) pattern_tri.clear();
    }
    const Escape::CGraph* dg = data_graph.c_graph;
    std::vector<int64_t> data_tri;          // whole-graph counts, once computed
    double whole_graph_cost = -1;
    auto dataDeg = [&](int64_t v) { return dg->offsets[v + 1] - dg->offsets[v]; };
    double index_scan_work = 0;   // adjacency positions the whole index build scans
    if (!pattern_tri.empty()) {
        std::vector<double> scan_of(num_p_nodes, 0);
        for (size_t x = 0; x < num_p_nodes; ++x) {
            for (int v : candidate_lists[x]) {
                const int v_idx = data_graph.get_node_idx(v);
                if (v_idx != -1) scan_of[x] += static_cast<double>(dataDeg(v_idx));
            }
        }
        for (size_t d = 0; d < order.size(); ++d) {
            for (size_t e = 0; e < d; ++e) {
                if (pattern_graph.hasEdge(pattern_graph.original_node_ids[order[e]], pattern_graph.original_node_ids[order[d]]))
                    index_scan_work += scan_of[order[e]];
            }
        }
    }
    auto patternEdgePos = [&](int a, int b) -> Escape::EdgeIdx {
        const Escape::CGraph* pg = pattern_graph.c_graph;
        const int ia = pattern_graph.get_node_idx(pattern_graph.original_node_ids[a]);
        const int ib = pattern_graph.get_node_idx(pattern_graph.original_node_ids[b]);
        const Escape::VertexIdx* it = std::lower_bound(pg->nbors + pg->offsets[ia], pg->nbors + pg->offsets[ia + 1], ib);
        return it - pg->nbors;
    };

    std::vector<int> index_of(max_v_id + 1, -1);
    size_t max_candidates = 0;
    backward_edges.assign(order.size(), {});
    backward_non_neighbors.assign(order.size(), {});
    for (size_t d = 0; d < order.size(); ++d) {
        int u = order[d];
        const auto& cands_u = candidate_lists[u];
        max_candidates = std::max(max_candidates, cands_u.size());
        for (size_t j = 0; j < cands_u.size(); ++j) index_of[cands_u[j]] = static_cast<int>(j);

        int u_id = pattern_graph.original_node_ids[u];
        for (size_t e = 0; e < d; ++e) {
            int from_u = order[e];
            if (!pattern_graph.hasEdge(pattern_graph.original_node_ids[from_u], u_id)) {
                backward_non_neighbors[d].push_back(from_u);
                continue;
            }
            setConflict(backward_masks.data() + d * conflict_words, static_cast<int>(e));
            BackwardEdge edge;
            edge.from_u = from_u;
            const auto& cands_from = candidate_lists[from_u];
            edge.offsets.assign(cands_from.size() + 1, 0);
            double scan_work = 0;
            for (size_t i = 0; i < cands_from.size(); ++i) {
                int v_idx = data_graph.get_node_idx(cands_from[i]);
                if (v_idx != -1) {
                    scan_work += static_cast<double>(dataDeg(v_idx));
                    for (Escape::EdgeIdx k = dg->offsets[v_idx]; k < dg->offsets[v_idx + 1]; ++k) {
                        int j = index_of[data_graph.original_node_ids[dg->nbors[k]]];
                        if (j != -1) edge.targets.push_back(static_cast<uint32_t>(j));
                    }
                }
                edge.offsets[i + 1] = static_cast<uint32_t>(edge.targets.size());
            }
            index_entries_before += static_cast<long long>(edge.targets.size());

            const int64_t need_t = pattern_tri.empty() ? 0 : pattern_tri[patternEdgePos(from_u, u)];
            if (need_t > 0 && !edge.targets.empty()) {
                // Same walk as above, so the t-th entry is met at the same position k.
                bool per_entry = data_tri.empty();
                bool apply = true;
                if (per_entry) {
                    double cost = 0;
                    for (size_t i = 0; i < cands_from.size(); ++i) {
                        int v_idx = data_graph.get_node_idx(cands_from[i]);
                        if (v_idx == -1) continue;
                        for (Escape::EdgeIdx k = dg->offsets[v_idx]; k < dg->offsets[v_idx + 1]; ++k) {
                            if (index_of[data_graph.original_node_ids[dg->nbors[k]]] != -1)
                                cost += intersectionCost(dataDeg(v_idx), dataDeg(dg->nbors[k]));
                        }
                    }
                    const double whole_budget = 4 * index_scan_work + 1000;
                    if (cost > 4 * scan_work + 1000) {
                        // The whole-graph pass costs at least m; estimate it (in O(m)) only if
                        // that fits its budget.
                        if (whole_graph_cost < 0 && static_cast<double>(dg->nEdges) <= whole_budget) {
                            whole_graph_cost = 0;
                            for (int64_t x = 0; x < dg->nVertices; ++x)
                                for (Escape::EdgeIdx k = dg->offsets[x]; k < dg->offsets[x + 1]; ++k)
                                    whole_graph_cost += 1.0 + std::min(dataDeg(x), dataDeg(dg->nbors[k]));
                        }
                        if (whole_graph_cost >= 0 && whole_graph_cost <= whole_budget) {
                            data_tri = edgeTriangles(dg);
                            per_entry = false;
                        } else {
                            apply = false;
                        }
                    }
                }
                size_t t = 0, kept = 0;
                std::vector<uint32_t> new_offsets(cands_from.size() + 1, 0);
                if (apply) {
                for (size_t i = 0; i < cands_from.size(); ++i) {
                    int v_idx = data_graph.get_node_idx(cands_from[i]);
                    if (v_idx != -1) {
                        for (Escape::EdgeIdx k = dg->offsets[v_idx]; k < dg->offsets[v_idx + 1]; ++k) {
                            if (index_of[data_graph.original_node_ids[dg->nbors[k]]] == -1) continue;
                            int64_t tri;
                            if (!per_entry) {
                                tri = data_tri[k];
                            } else {
                                tri = commonNeighbours(dg, v_idx, dg->nbors[k]);
                            }
                            const uint32_t target = edge.targets[t++];
                            if (tri >= need_t) edge.targets[kept++] = target;
                        }
                    }
                    new_offsets[i + 1] = static_cast<uint32_t>(kept);
                }
                edge.targets.resize(kept);
                edge.offsets = std::move(new_offsets);
                }
            }
            backward_edges[d].push_back(std::move(edge));
        }

        for (int v : cands_u) index_of[v] = -1;
    }
    for (const auto& edges : backward_edges)
        for (const auto& e : edges) index_entries_after += static_cast<long long>(e.targets.size());

    smaller_depths.assign(order.size(), {});
    if (symmetry && count_only) {
        symmetry_multiplier = symmetry->multiplier;
        for (size_t d = 0; d < order.size(); ++d) {
            const int u_idx = pattern_graph.get_node_idx(pattern_graph.original_node_ids[order[d]]);
            for (int x_idx : symmetry->smaller_than[u_idx]) {
                smaller_depths[d].push_back(order_position[pattern_graph.original_node_ids[x_idx]]);
            }
        }
    }

    // +8: the AVX2 intersection writes up to 8 elements past its result
    local_candidates.assign(order.size() + 1, std::vector<int>(max_candidates + 8));
    intersect_buffer.assign(max_candidates + 8, 0);
    identity.resize(max_candidates);
    for (size_t j = 0; j < max_candidates; ++j) identity[j] = static_cast<int>(j);
    local_ptr.assign(order.size() + 1, nullptr);
    local_len.assign(order.size() + 1, 0);
}

void SearchEngine::mergeConflicts(int target_depth, int source_depth) {
    uint64_t* target = conflicts(target_depth);
    const uint64_t* source = conflicts(source_depth);
    for (size_t w = 0; w < conflict_words; ++w) target[w] |= source[w];
}

void SearchEngine::run() {
    matches.clear();
    match_count = 0;
    stopped = false;
    std::fill(fast_mapping.begin(), fast_mapping.end(), -1);
    std::fill(fast_mapping_local_index.begin(), fast_mapping_local_index.end(), -1);
    std::fill(fast_inverse_mapping.begin(), fast_inverse_mapping.end(), -1);
    std::fill(in_use.begin(), in_use.end(), 0);
    nodes_per_depth.assign(order.size(), 0);
    local_per_depth.assign(order.size(), 0);
    if (order.empty()) return;
    if (enable_cbj) {
        backtrack_cbj(0);
    } else {
        backtrack(0);
    }
    if (symmetry_multiplier > 1) {
        const __int128 total = static_cast<__int128>(match_count) * symmetry_multiplier;
        count_overflow = total > static_cast<__int128>(INT64_MAX);
        match_count = count_overflow ? INT64_MAX : static_cast<long long>(total);
    }
}

// Points local_ptr[depth] at order[depth]'s candidates (indices into its candidate list) that are
// adjacent to every matched backward neighbour: straight into the candidate-space index when
// there is one backward neighbour, else the intersection, ping-ponged between
// local_candidates[depth] and intersect_buffer without copying.
bool SearchEngine::computeLocalCandidates(int depth) {
    const auto& edges = backward_edges[depth];

    if (edges.empty()) {
        local_ptr[depth] = identity.data();
        local_len[depth] = candidate_lists[order[depth]].size();
    } else {
        int shortest = 0;
        size_t shortest_len = SIZE_MAX;
        for (size_t k = 0; k < edges.size(); ++k) {
            int i = fast_mapping_local_index[edges[k].from_u];
            size_t len = edges[k].offsets[i + 1] - edges[k].offsets[i];
            if (len < shortest_len) {
                shortest_len = len;
                shortest = static_cast<int>(k);
            }
        }
        const BackwardEdge& first = edges[shortest];
        const int* cur = reinterpret_cast<const int*>(first.targets.data()) +
                         first.offsets[fast_mapping_local_index[first.from_u]];
        size_t size = shortest_len;
        bool into_local = true;
        for (size_t k = 0; k < edges.size() && size > 0; ++k) {
            if (static_cast<int>(k) == shortest) continue;
            int i = fast_mapping_local_index[edges[k].from_u];
            const int* list = reinterpret_cast<const int*>(edges[k].targets.data()) + edges[k].offsets[i];
            size_t len = edges[k].offsets[i + 1] - edges[k].offsets[i];
            if (into_local) {
                size = intersectSorted(cur, size, list, len, local_candidates[depth].data());
                into_local = false;
            } else {
                size = intersectSorted(cur, size, list, len, intersect_buffer.data());
                std::swap(local_candidates[depth], intersect_buffer);
            }
            cur = local_candidates[depth].data();
        }
        local_ptr[depth] = cur;
        local_len[depth] = size;
    }
    // Symmetry breaking: this vertex's image must exceed the images of smaller_depths[depth].
    // Local candidates are ascending candidate-list positions, and the list is sorted by data
    // vertex ID, so the allowed ones are a suffix.
    bool trimmed = false;
    if (!smaller_depths[depth].empty() && local_len[depth] > 0) {
        int bound = -1;
        for (int d : smaller_depths[depth]) bound = std::max(bound, fast_mapping[order[d]]);
        const auto& cands_u = candidate_lists[order[depth]];
        const int* begin = local_ptr[depth];
        const int* first = std::partition_point(begin, begin + local_len[depth],
                                                [&](int i) { return cands_u[i] <= bound; });
        trimmed = first != begin;
        local_len[depth] -= static_cast<size_t>(first - begin);
        local_ptr[depth] = first;
    }
    ++nodes_per_depth[depth];
    local_per_depth[depth] += static_cast<long long>(local_len[depth]);
    return trimmed;
}

bool SearchEngine::inducedConflict(int u, int v, int& conflict_u) const {
    for (int u_non_nb : backward_non_neighbors[order_position[u]]) {
        if (data_graph.hasEdge(v, fast_mapping[u_non_nb])) {
            conflict_u = u_non_nb;
            return true;
        }
    }
    return false;
}

// Matches completed at the last pattern vertex, without visiting each: every local candidate
// completes exactly one match unless a (distinct) already-matched pattern vertex uses it.
// Those pattern vertices -- the injectivity conflicts of the per-candidate loop -- are added
// to used_by (a conflict bitset, by order depth) when given.
long long SearchEngine::countLastLevel(int depth, uint64_t* used_by) const {
    const auto& cands_u = candidate_lists[order[depth]];
    const int* local = local_ptr[depth];
    const size_t len = local_len[depth];
    long long count = static_cast<long long>(len);

    // Short lists: an O(1) check per candidate beats two binary searches per matched vertex.
    if (len <= 8 * static_cast<size_t>(depth)) {
        for (size_t k = 0; k < len; ++k) {
            const int v = cands_u[local[k]];
            if (in_use[v]) {
                --count;
                if (used_by) setConflict(used_by, order_position[fast_inverse_mapping[v]]);
            }
        }
        return count;
    }

    for (int d = 0; d < depth; ++d) {
        const int x = order[d];
        auto c = std::lower_bound(cands_u.begin(), cands_u.end(), fast_mapping[x]);
        if (c == cands_u.end() || *c != fast_mapping[x]) continue;
        if (std::binary_search(local, local + len, static_cast<int>(c - cands_u.begin()))) {
            --count;
            if (used_by) setConflict(used_by, d);
        }
    }
    return count;
}

void SearchEngine::backtrack(int depth) {
    if (static_cast<size_t>(depth) == order.size()) {
        ++match_count;
        if (!count_only && store_matches) {
            matches.push_back(fast_mapping);
        }
        if (match_limit > 0 && match_count >= match_limit) stopped = true;
        return;
    }

    int u = order[depth];
    computeLocalCandidates(depth);
    if (count_only && !is_induced && static_cast<size_t>(depth) + 1 == order.size()) {
        match_count += countLastLevel(depth, nullptr);
        return;
    }
    const auto& cands_u = candidate_lists[u];
    const int* local = local_ptr[depth];
    const size_t len = local_len[depth];
    for (size_t k = 0; k < len; ++k) {
        int i = local[k];
        int v = cands_u[i];
        if (in_use[v]) continue;
        int unused;
        if (is_induced && inducedConflict(u, v, unused)) continue;

        fast_mapping[u] = v;
        fast_mapping_local_index[u] = i;
        fast_inverse_mapping[v] = u;
        in_use[v] = 1;

        backtrack(depth + 1);

        fast_mapping[u] = -1;
        fast_mapping_local_index[u] = -1;
        fast_inverse_mapping[v] = -1;
        in_use[v] = 0;
        if (stopped) return;
    }
}

int SearchEngine::backtrack_cbj(int depth) {
    if (static_cast<size_t>(depth) == order.size()) {
        ++match_count;
        if (!count_only && store_matches) {
            matches.push_back(fast_mapping);
        }
        if (match_limit > 0 && match_count >= match_limit) stopped = true;
        return depth - 1;
    }

    int u = order[depth];
    uint64_t* conflict = conflicts(depth);
    std::fill(conflict, conflict + conflict_words, 0);
    long long match_count_before = match_count;

    const bool trimmed = computeLocalCandidates(depth);
    const auto& cands_u = candidate_lists[u];
    const int* local = local_ptr[depth];
    const size_t len = local_len[depth];

    if (len < cands_u.size()) {
        const uint64_t* backward = backward_masks.data() + static_cast<size_t>(depth) * conflict_words;
        for (size_t w = 0; w < conflict_words; ++w) conflict[w] |= backward[w];
    }
    if (trimmed) {
        for (int d : smaller_depths[depth]) setConflict(conflict, d);
    }

    const bool count_last = count_only && !is_induced && static_cast<size_t>(depth) + 1 == order.size();
    if (count_last) match_count += countLastLevel(depth, conflict);

    for (size_t k = 0; !count_last && k < len; ++k) {
        int i = local[k];
        int v = cands_u[i];

        int conflict_u = in_use[v] ? fast_inverse_mapping[v] : -1;
        if (conflict_u == -1 && is_induced) inducedConflict(u, v, conflict_u);
        if (conflict_u != -1) {
            setConflict(conflict, order_position[conflict_u]);
            continue;
        }

        fast_mapping[u] = v;
        fast_mapping_local_index[u] = i;
        fast_inverse_mapping[v] = u;
        in_use[v] = 1;

        int backjump_depth = backtrack_cbj(depth + 1);

        fast_mapping[u] = -1;
        fast_mapping_local_index[u] = -1;
        fast_inverse_mapping[v] = -1;
        in_use[v] = 0;
        if (stopped) return -1;   // the match limit unwinds the whole search

        if (backjump_depth >= 0 && backjump_depth < depth) {
            setConflict(conflict, backjump_depth);
        }

        if (backjump_depth < depth - 1) {
            if (backjump_depth >= 0) {
                mergeConflicts(backjump_depth, depth);
            }
            return backjump_depth;
        }
    }

    int target_depth = depth - 1;
    if (match_count == match_count_before) {
        // the deepest conflict above this depth, if any
        for (int w = static_cast<int>((depth - 1) >> 6); depth > 0 && w >= 0; --w) {
            uint64_t bits = conflict[w];
            if (w == ((depth - 1) >> 6) && ((depth & 63) != 0)) bits &= (uint64_t{1} << (depth & 63)) - 1;
            if (bits) {
                target_depth = w * 64 + 63 - __builtin_clzll(bits);
                break;
            }
        }
    }
    if (target_depth >= 0) {
        mergeConflicts(target_depth, depth);
    }

    return target_depth;
}
