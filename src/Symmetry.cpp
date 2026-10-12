/// @file
/// Symmetry-breaking constraints for counting: the pattern's automorphism orbits by
/// individualisation-refinement, turned into Grochow-Kellis ordering constraints. See Symmetry.h.

#include "Symmetry.h"
#include <algorithm>
#include <map>
#include <numeric>

namespace {

// Colour refinement (1-dimensional Weisfeiler-Leman) of two colourings of the same graph at once,
// with one shared renaming of the signatures, so equal colours mean the same thing on both sides.
// Colours are always 0..k-1. Returns false when the two sides' colour histograms differ (no
// automorphism can map one colouring to the other) or when the work budget is exhausted.
class Refiner {
public:
    Refiner(const std::vector<std::vector<int>>& adj, long long budget) : adj(adj), budget(budget) {
        for (const auto& a : adj) edges += static_cast<long long>(a.size());
    }

    bool refine(std::vector<int>& c1, std::vector<int>& c2) {
        const size_t n = adj.size();
        int classes = numClasses(c1);
        std::vector<std::vector<int>> s1(n), s2(n);
        while (true) {
            work += 2 * (static_cast<long long>(n) + edges);
            if (work > budget) { out_of_budget = true; return false; }
            auto signatures = [&](const std::vector<int>& c, std::vector<std::vector<int>>& s) {
                for (size_t v = 0; v < n; ++v) {
                    s[v].clear();
                    s[v].push_back(c[v]);
                    for (int w : adj[v]) s[v].push_back(c[w]);
                    std::sort(s[v].begin() + 1, s[v].end());
                }
            };
            signatures(c1, s1);
            signatures(c2, s2);
            std::map<std::vector<int>, int> id;
            for (size_t v = 0; v < n; ++v) { id.emplace(s1[v], 0); id.emplace(s2[v], 0); }
            int next = 0;
            for (auto& kv : id) kv.second = next++;
            std::vector<int> hist(next, 0);
            for (size_t v = 0; v < n; ++v) { c1[v] = id[s1[v]]; ++hist[c1[v]]; }
            for (size_t v = 0; v < n; ++v) { c2[v] = id[s2[v]]; --hist[c2[v]]; }
            for (int h : hist) if (h != 0) return false;
            if (next == classes) return true;
            classes = next;
        }
    }

    static int numClasses(const std::vector<int>& c) {
        return c.empty() ? 0 : *std::max_element(c.begin(), c.end()) + 1;
    }

    const std::vector<std::vector<int>>& adj;
    long long budget;
    long long edges = 0;
    long long work = 0;
    bool out_of_budget = false;
};

bool hasEdge(const std::vector<std::vector<int>>& adj, int a, int b) {
    return std::binary_search(adj[a].begin(), adj[a].end(), b);
}

// An automorphism mapping colouring c1 onto c2 (both equitable-refinable), by individualisation-
// refinement: refine, and if not yet discrete, individualise one vertex of the smallest
// non-singleton cell on the left against each vertex of that cell on the right.
bool findIsomorphism(Refiner& r, const std::vector<int>& labels, std::vector<int> c1,
                     std::vector<int> c2, std::vector<int>& sigma) {
    if (!r.refine(c1, c2)) return false;
    const size_t n = c1.size();
    const int k = Refiner::numClasses(c1);
    if (static_cast<size_t>(k) == n) {
        std::vector<int> at(n);
        for (size_t v = 0; v < n; ++v) at[c2[v]] = static_cast<int>(v);
        sigma.assign(n, -1);
        for (size_t v = 0; v < n; ++v) sigma[v] = at[c1[v]];
        for (size_t v = 0; v < n; ++v) {
            if (labels[v] != labels[sigma[v]] || r.adj[v].size() != r.adj[sigma[v]].size()) return false;
            for (int w : r.adj[v]) if (!hasEdge(r.adj, sigma[v], sigma[w])) return false;
        }
        return true;
    }
    std::vector<int> size(k, 0);
    for (int col : c1) ++size[col];
    int cell = -1;
    for (int col = 0; col < k; ++col) {
        if (size[col] > 1 && (cell == -1 || size[col] < size[cell])) cell = col;
    }
    int a = -1;
    for (size_t v = 0; v < n && a == -1; ++v) if (c1[v] == cell) a = static_cast<int>(v);
    for (size_t b = 0; b < n; ++b) {
        if (c2[b] != cell) continue;
        std::vector<int> d1 = c1, d2 = c2;
        d1[a] = k;
        d2[b] = k;
        if (findIsomorphism(r, labels, std::move(d1), std::move(d2), sigma)) return true;
        if (r.out_of_budget) return false;
    }
    return false;
}

int find(std::vector<int>& parent, int x) {
    while (parent[x] != x) x = parent[x] = parent[parent[x]];
    return x;
}

}  // namespace

SymmetryBreaking computeSymmetryBreaking(const AppGraph& pattern, const OrbitCounts& orbits,
                                         const std::vector<int>& order, long long budget) {
    SymmetryBreaking result;
    const size_t n = pattern.original_node_ids.size();
    result.smaller_than.assign(n, {});
    if (n < 2 || !pattern.c_graph) return result;

    std::vector<std::vector<int>> adj(n);
    std::vector<int> labels(n);
    for (size_t i = 0; i < n; ++i) {
        labels[i] = pattern.get_label_fast(pattern.original_node_ids[i]);
        for (Escape::EdgeIdx p = pattern.c_graph->offsets[i]; p < pattern.c_graph->offsets[i + 1]; ++p) {
            adj[i].push_back(static_cast<int>(pattern.c_graph->nbors[p]));
        }
        std::sort(adj[i].begin(), adj[i].end());
    }

    // Initial colours: (label, graphlet orbit-count vector) -- both preserved by automorphisms.
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    auto key_less = [&](int a, int b) {
        if (labels[a] != labels[b]) return labels[a] < labels[b];
        return orbits[a] < orbits[b];
    };
    std::sort(idx.begin(), idx.end(), key_less);
    std::vector<int> colour(n, 0);
    for (size_t r = 1; r < n; ++r) {
        colour[idx[r]] = colour[idx[r - 1]] + (key_less(idx[r - 1], idx[r]) ? 1 : 0);
    }

    Refiner refiner(adj, budget);
    auto refineOne = [&](std::vector<int>& c) {
        std::vector<int> copy = c;
        return refiner.refine(c, copy);
    };
    if (!refineOne(colour)) { result.truncated = refiner.out_of_budget; result.work = refiner.work; return result; }

    std::vector<int> parent(n);
    std::vector<int> sigma;
    for (int x_id : order) {
        const int x = pattern.get_node_idx(x_id);
        const int k = Refiner::numClasses(colour);
        if (static_cast<size_t>(k) == n) break;   // discrete: only the identity fixes the vertices so far

        std::iota(parent.begin(), parent.end(), 0);
        for (size_t w = 0; w < n; ++w) {
            if (static_cast<int>(w) == x || colour[w] != colour[x] || find(parent, w) == find(parent, x)) continue;
            // Twins (same label, same neighbours apart from each other): swapping them is an
            // automorphism fixing everything else.
            std::vector<int> nx, nw;
            for (int y : adj[x]) if (y != static_cast<int>(w)) nx.push_back(y);
            for (int y : adj[w]) if (y != x) nw.push_back(y);
            if (labels[x] == labels[w] && nx == nw) {
                parent[find(parent, w)] = find(parent, x);
                continue;
            }
            std::vector<int> c1 = colour, c2 = colour;
            c1[x] = k;
            c2[w] = k;
            if (findIsomorphism(refiner, labels, std::move(c1), std::move(c2), sigma)) {
                for (size_t v = 0; v < n; ++v) {
                    int a = find(parent, static_cast<int>(v)), b = find(parent, sigma[v]);
                    if (a != b) parent[a] = b;
                }
            } else if (refiner.out_of_budget) {
                result.truncated = true;
                result.work = refiner.work;
                return result;   // this level's orbit is unknown: keep only the levels before it
            }
        }

        std::vector<int> orbit;
        for (size_t w = 0; w < n; ++w) {
            if (static_cast<int>(w) != x && find(parent, static_cast<int>(w)) == find(parent, x)) orbit.push_back(static_cast<int>(w));
        }
        if (!orbit.empty()) {
            const uint64_t size = orbit.size() + 1;
            if (result.multiplier > (uint64_t{1} << 62) / size) { result.truncated = true; break; }
            result.multiplier *= size;
            ++result.levels;
            for (int w : orbit) result.smaller_than[w].push_back(x);
        }

        colour[x] = k;   // fix x for the levels after it
        if (!refineOne(colour)) { result.truncated = refiner.out_of_budget; break; }
    }
    result.work = refiner.work;
    return result;
}
