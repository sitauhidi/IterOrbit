/// @file
/// Graphlet orbit counting on top of the Escape/EVOKE library (extern/Escape). The routines in
/// the anonymous namespace are reimplementations of Escape routines that return the same counts
/// with less work, mostly by replacing binary searches for edges with marking; each one's comment
/// names the Escape routine it replaces and what differs. test/orbit_equivalence.cpp checks them
/// against the originals. See OrbitCounter.h.

#include "OrbitCounter.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <utility>
#include <stdexcept>
#include <iostream>

#include "Escape/OrbitStructure.h"
#include "Escape/GetAllOrbitCounts.h"
#include "Escape/FiveVertexOrbit.h"
#include "Escape/FourVertexOrbit.h"

using namespace Escape;

EVOKEOrbitCounter::EVOKEOrbitCounter(const AppGraph& g, int s, const std::vector<bool>* needed)
    : cgraph(g.c_graph), graphlet_size(s), needed_orbits(needed) {}

EVOKEOrbitCounter::EVOKEOrbitCounter(const Escape::CGraph* g, int s, const std::vector<bool>* needed)
    : cgraph(g), graphlet_size(s), needed_orbits(needed) {}

namespace {
using Clock = std::chrono::high_resolution_clock;
double msSince(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

// Same per-vertex 4-clique counts (orbit 14) as Escape's fourCliqueOrbitCounterNodeOnly, which
// binary-searches every pair of out-neighbours. Here each clique i<j<k<l (degree order) is
// found by marking: T = N+(i) ∩ N+(j) by scanning N+(j) against marks of N+(i), then every
// l in N+(k) marked in T.
// cap (by vertex, may be null): a vertex still needs cliques while its count is below its cap.
// A clique is enumerated from its two lowest vertices (i, j) with the rest in T, so a (i, j)
// where none of i, j, T still needs any is skipped. Counts only ever grow by real cliques, so a
// vertex ends either with its exact count or with a count >= its cap -- either way every
// "count >= requirement" test (requirement <= cap) comes out as with the exact count.
void fourCliqueOrbitsByMarking(const CGraph* gout, double* clique_counts, const std::vector<double>* cap) {
    auto need = [&](VertexIdx x) { return !cap || clique_counts[x] < (*cap)[x]; };
    const VertexIdx n = gout->nVertices;
    std::vector<char> in_out_i(n, 0), in_t(n, 0);
    std::vector<VertexIdx> t;
    for (VertexIdx i = 0; i < n; ++i) {
        const EdgeIdx i_begin = gout->offsets[i], i_end = gout->offsets[i + 1];
        if (i_end - i_begin < 3) continue;
        for (EdgeIdx p = i_begin; p < i_end; ++p) in_out_i[gout->nbors[p]] = 1;

        for (EdgeIdx pj = i_begin; pj < i_end; ++pj) {
            const VertexIdx j = gout->nbors[pj];
            t.clear();
            bool touches_needed = need(i) || need(j);
            for (EdgeIdx p = gout->offsets[j]; p < gout->offsets[j + 1]; ++p) {
                const VertexIdx k = gout->nbors[p];
                if (in_out_i[k]) {
                    t.push_back(k);
                    if (!touches_needed && need(k)) touches_needed = true;
                }
            }
            if (t.size() < 2 || !touches_needed) continue;
            for (VertexIdx k : t) in_t[k] = 1;
            double cliques_ij = 0;
            for (VertexIdx k : t) {
                double cliques_ijk = 0;
                for (EdgeIdx p = gout->offsets[k]; p < gout->offsets[k + 1]; ++p) {
                    const VertexIdx l = gout->nbors[p];
                    if (in_t[l]) {
                        ++cliques_ijk;
                        ++clique_counts[l];
                    }
                }
                clique_counts[k] += cliques_ijk;
                cliques_ij += cliques_ijk;
            }
            clique_counts[i] += cliques_ij;
            clique_counts[j] += cliques_ij;
            for (VertexIdx k : t) in_t[k] = 0;
        }

        for (EdgeIdx p = i_begin; p < i_end; ++p) in_out_i[gout->nbors[p]] = 0;
    }
}

// Same per-vertex 5-clique counts (orbit 72) as Escape's FiveCliqueOrbitCounter, which tests
// each pair of 4-clique ends with a linear scan (isEdge). Each clique i<j<k<l<o (degree order)
// is found once by marking: T1 = N+(i) ∩ N+(j), T2 = N+(k) ∩ T1, then o in N+(l) ∩ T2.
// cap as in fourCliqueOrbitsByMarking: a branch whose possible members all reached their cap
// is skipped; every vertex ends with its exact count or a count >= its cap.
void fiveCliqueOrbitsByMarking(const CGraph* gout, double* clique_counts, const std::vector<double>* cap) {
    auto need = [&](VertexIdx x) { return !cap || clique_counts[x] < (*cap)[x]; };
    const VertexIdx n = gout->nVertices;
    std::vector<char> in_out_i(n, 0), in_t1(n, 0), in_t2(n, 0);
    std::vector<VertexIdx> t1, t2;
    for (VertexIdx i = 0; i < n; ++i) {
        const EdgeIdx i_begin = gout->offsets[i], i_end = gout->offsets[i + 1];
        if (i_end - i_begin < 4) continue;
        for (EdgeIdx p = i_begin; p < i_end; ++p) in_out_i[gout->nbors[p]] = 1;

        for (EdgeIdx pj = i_begin; pj < i_end; ++pj) {
            const VertexIdx j = gout->nbors[pj];
            t1.clear();
            bool touches_needed = need(i) || need(j);
            for (EdgeIdx p = gout->offsets[j]; p < gout->offsets[j + 1]; ++p) {
                const VertexIdx k = gout->nbors[p];
                if (in_out_i[k]) {
                    t1.push_back(k);
                    if (!touches_needed && need(k)) touches_needed = true;
                }
            }
            if (t1.size() < 3 || !touches_needed) continue;
            for (VertexIdx k : t1) in_t1[k] = 1;

            double cliques_ij = 0;
            for (VertexIdx k : t1) {
                t2.clear();
                bool k_touches_needed = need(i) || need(j) || need(k);
                for (EdgeIdx p = gout->offsets[k]; p < gout->offsets[k + 1]; ++p) {
                    const VertexIdx l = gout->nbors[p];
                    if (in_t1[l]) {
                        t2.push_back(l);
                        if (!k_touches_needed && need(l)) k_touches_needed = true;
                    }
                }
                if (t2.size() < 2 || !k_touches_needed) continue;
                for (VertexIdx l : t2) in_t2[l] = 1;
                double cliques_ijk = 0;
                for (VertexIdx l : t2) {
                    double cliques_ijkl = 0;
                    for (EdgeIdx p = gout->offsets[l]; p < gout->offsets[l + 1]; ++p) {
                        const VertexIdx o = gout->nbors[p];
                        if (in_t2[o]) {
                            ++cliques_ijkl;
                            ++clique_counts[o];
                        }
                    }
                    clique_counts[l] += cliques_ijkl;
                    cliques_ijk += cliques_ijkl;
                }
                clique_counts[k] += cliques_ijk;
                cliques_ij += cliques_ijk;
                for (VertexIdx l : t2) in_t2[l] = 0;
            }
            clique_counts[i] += cliques_ij;
            clique_counts[j] += cliques_ij;
            for (VertexIdx k : t1) in_t1[k] = 0;
        }

        for (EdgeIdx p = i_begin; p < i_end; ++p) in_out_i[gout->nbors[p]] = 0;
    }
}

// Same outputs as Escape's get_all_three_orbits (orbits 0-3, plus the per-edge values its
// three_edge_orbit_counter writes: per_edge_[0], per_edge_rev_[0], per_edge_[1]), which
// enumerates every triangle twice with a binary search per out-neighbour pair. Here each
// triangle i<j<k (degree order) is found once: mark N+(i) with edge positions, then every
// marked k in N+(j). Triangle enumeration is skipped when count_triangles is false.
void threeVertexOrbitsByMarking(const CGraph* cg, const CGraph* gout, OrbitInfo& oc, bool count_triangles) {
    const VertexIdx n = cg->nVertices;
    for (VertexIdx i = 0; i < n; ++i) {
        const VertexIdx deg = cg->offsets[i + 1] - cg->offsets[i];
        oc.per_vertex_[0][i] = deg;
        for (EdgeIdx p = cg->offsets[i]; p < cg->offsets[i + 1]; ++p) {
            const VertexIdx j = cg->nbors[p];
            oc.per_vertex_[1][i] += (cg->offsets[j + 1] - cg->offsets[j]) - 1;
        }
        oc.per_vertex_[2][i] = choosec(deg, 2);
    }
    for (VertexIdx i = 0; i < n; ++i) {
        for (EdgeIdx p = gout->offsets[i]; p < gout->offsets[i + 1]; ++p) {
            oc.per_edge_[0][p] += oc.per_vertex_[0][i] - 1;
            oc.per_edge_rev_[0][p] += oc.per_vertex_[0][gout->nbors[p]] - 1;
        }
    }
    if (!count_triangles) return;

    std::vector<EdgeIdx> pos_in_out_i(n, -1);
    double* triangles = oc.per_vertex_[3];
    double* edge_triangles = oc.per_edge_[1];
    for (VertexIdx i = 0; i < n; ++i) {
        const EdgeIdx i_begin = gout->offsets[i], i_end = gout->offsets[i + 1];
        if (i_end - i_begin < 2) continue;
        for (EdgeIdx p = i_begin; p < i_end; ++p) pos_in_out_i[gout->nbors[p]] = p;
        for (EdgeIdx pij = i_begin; pij < i_end; ++pij) {
            const VertexIdx j = gout->nbors[pij];
            for (EdgeIdx pjk = gout->offsets[j]; pjk < gout->offsets[j + 1]; ++pjk) {
                const VertexIdx k = gout->nbors[pjk];
                const EdgeIdx pik = pos_in_out_i[k];
                if (pik == -1) continue;
                ++triangles[i];
                ++triangles[j];
                ++triangles[k];
                ++edge_triangles[pij];
                ++edge_triangles[pik];
                ++edge_triangles[pjk];
            }
        }
        for (EdgeIdx p = i_begin; p < i_end; ++p) pos_in_out_i[gout->nbors[p]] = -1;
    }
}

// Same orbits 9, 10, 12, 13 as Escape's otherFourPatternOrbitCounterNodeOnly, which found the
// out-edge position of every in-edge by binary search and re-enumerated triangles with a
// binary search per out-neighbour pair. Its in-edge loop is folded into the out-edge loop
// (every in-edge j->i is an out-edge of j, credited to its head), and triangles come from
// marking. Needs orbits 0 and 3 and the per-edge triangle counts per_edge_[1].
void pawDiamondOrbits(const CGraph* gout, OrbitInfo& oc) {
    const VertexIdx n = gout->nVertices;
    const double* deg = oc.per_vertex_[0];
    const double* tri = oc.per_vertex_[3];
    const double* edge_tri = oc.per_edge_[1];
    double* o9 = oc.per_vertex_[9];
    double* o10 = oc.per_vertex_[10];
    double* o12 = oc.per_vertex_[12];
    double* o13 = oc.per_vertex_[13];

    for (VertexIdx i = 0; i < n; ++i) {
        for (EdgeIdx p = gout->offsets[i]; p < gout->offsets[i + 1]; ++p) {
            const VertexIdx j = gout->nbors[p];
            const double et = edge_tri[p];
            o9[i] += tri[j] - et;
            o10[i] += et * (deg[j] - 2);
            o13[j] += choosec(et, 2);
            o9[j] += tri[i] - et;
            o10[j] += et * (deg[i] - 2);
            o13[i] += (et * (et - 1)) / 2;
        }
    }

    std::vector<EdgeIdx> pos_in_out_i(n, -1);
    for (VertexIdx i = 0; i < n; ++i) {
        const EdgeIdx i_begin = gout->offsets[i], i_end = gout->offsets[i + 1];
        if (i_end - i_begin < 2) continue;
        for (EdgeIdx p = i_begin; p < i_end; ++p) pos_in_out_i[gout->nbors[p]] = p;
        for (EdgeIdx pij = i_begin; pij < i_end; ++pij) {
            const VertexIdx j = gout->nbors[pij];
            for (EdgeIdx pjk = gout->offsets[j]; pjk < gout->offsets[j + 1]; ++pjk) {
                const VertexIdx k = gout->nbors[pjk];
                const EdgeIdx pik = pos_in_out_i[k];
                if (pik == -1) continue;
                o12[i] += edge_tri[pjk] - 1;
                o12[j] += edge_tri[pik] - 1;
                o12[k] += edge_tri[pij] - 1;
            }
        }
        for (EdgeIdx p = i_begin; p < i_end; ++p) pos_in_out_i[gout->nbors[p]] = -1;
    }
}

// Same orbit-8 counts as Escape's four_cycle_orbit_counter_node_only, minus two of its
// costs: it binary-searched two edge positions per wedge that were never used, and its
// third pass re-walked every wedge although only the first visit to each end k adds
// anything (later ones add choosec(0,2)); here the distinct ends are recorded in pass one.
void fourCycleOrbits(const CGraph* gout, const CGraph* gin, double* cycle_counts) {
    const VertexIdx n = gout->nVertices;
    std::vector<VertexIdx> wedge_count(n, 0);
    std::vector<VertexIdx> ends;
    for (VertexIdx i = 0; i < n; ++i) {
        ends.clear();
        for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i + 1]; ++pos) {
            const VertexIdx j = gin->nbors[pos];
            for (EdgeIdx e = gout->offsets[j]; e < gout->offsets[j + 1]; ++e) {
                const VertexIdx k = gout->nbors[e];
                if (k >= i) continue;
                if (wedge_count[k]++ == 0) ends.push_back(k);
            }
            for (EdgeIdx e = gin->offsets[j]; e < gin->offsets[j + 1]; ++e) {
                const VertexIdx k = gin->nbors[e];
                if (wedge_count[k]++ == 0) ends.push_back(k);
            }
        }
        for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i + 1]; ++pos) {
            const VertexIdx j = gin->nbors[pos];
            for (EdgeIdx e = gout->offsets[j]; e < gout->offsets[j + 1]; ++e) {
                const VertexIdx k = gout->nbors[e];
                if (k >= i) continue;
                cycle_counts[j] += std::max<VertexIdx>(0, wedge_count[k] - 1);
            }
            for (EdgeIdx e = gin->offsets[j]; e < gin->offsets[j + 1]; ++e) {
                cycle_counts[j] += std::max<VertexIdx>(0, wedge_count[gin->nbors[e]] - 1);
            }
        }
        for (VertexIdx k : ends) {
            const double pairs = choosec(wedge_count[k], 2);
            cycle_counts[i] += pairs;
            cycle_counts[k] += pairs;
            wedge_count[k] = 0;
        }
    }
}

// Calls f(i, j, k, l, ij, ik, il, jk, jl, kl) once per 4-clique i < j < k < l (degree order),
// with the gout positions of its six edges -- the cliques, in the order, that Escape's
// fourCliqueOrbitCounter and Pattern66OrbitCounter find by binary-searching every pair of
// triangle ends. Here: mark N+(i) with positions; T = the marked k in N+(j); then every l in
// N+(k) that is in T.
template <typename F>
void forEachFourClique(const CGraph* gout, F&& f) {
    const VertexIdx n = gout->nVertices;
    std::vector<EdgeIdx> pos_from_i(n, -1);
    std::vector<VertexIdx> index_in_t(n, -1);
    struct TriangleEnd { VertexIdx k; EdgeIdx ik, jk; };
    std::vector<TriangleEnd> t;
    for (VertexIdx i = 0; i < n; ++i) {
        const EdgeIdx i_begin = gout->offsets[i], i_end = gout->offsets[i + 1];
        if (i_end - i_begin < 3) continue;
        for (EdgeIdx p = i_begin; p < i_end; ++p) pos_from_i[gout->nbors[p]] = p;
        for (EdgeIdx ij = i_begin; ij < i_end; ++ij) {
            const VertexIdx j = gout->nbors[ij];
            t.clear();
            for (EdgeIdx p = gout->offsets[j]; p < gout->offsets[j + 1]; ++p) {
                const VertexIdx k = gout->nbors[p];
                if (pos_from_i[k] >= 0) t.push_back({k, pos_from_i[k], p});
            }
            if (t.size() < 2) continue;
            for (size_t x = 0; x < t.size(); ++x) index_in_t[t[x].k] = x;
            for (const TriangleEnd& a : t) {
                for (EdgeIdx kl = gout->offsets[a.k]; kl < gout->offsets[a.k + 1]; ++kl) {
                    const VertexIdx x = index_in_t[gout->nbors[kl]];
                    if (x < 0) continue;
                    const TriangleEnd& b = t[x];
                    f(i, j, a.k, b.k, ij, a.ik, b.ik, a.jk, b.jk, kl);
                }
            }
            for (const TriangleEnd& a : t) index_in_t[a.k] = -1;
        }
        for (EdgeIdx p = i_begin; p < i_end; ++p) pos_from_i[gout->nbors[p]] = -1;
    }
}

// Same orbit 14 and per-edge 4-clique counts (per_edge_[11]) as Escape's fourCliqueOrbitCounter.
void fourCliqueOrbitsWithEdges(const CGraph* gout, OrbitInfo& oc) {
    double* v14 = oc.per_vertex_[14];
    double* e11 = oc.per_edge_[11];
    forEachFourClique(gout, [&](VertexIdx i, VertexIdx j, VertexIdx k, VertexIdx l, EdgeIdx ij, EdgeIdx ik,
                                EdgeIdx il, EdgeIdx jk, EdgeIdx jl, EdgeIdx kl) {
        ++v14[i]; ++v14[j]; ++v14[k]; ++v14[l];
        ++e11[ij]; ++e11[ik]; ++e11[il]; ++e11[jk]; ++e11[jl]; ++e11[kl];
    });
}

// Same orbit-66 counts as Escape's Pattern66OrbitCounter: per 4-clique, each vertex gets
// per_edge_[1] - 2 of the three clique edges not at it.
void pattern66Orbits(const CGraph* gout, OrbitInfo& oc) {
    double* v66 = oc.per_vertex_[66];
    const double* e1 = oc.per_edge_[1];
    forEachFourClique(gout, [&](VertexIdx i, VertexIdx j, VertexIdx k, VertexIdx l, EdgeIdx ij, EdgeIdx ik,
                                EdgeIdx il, EdgeIdx jk, EdgeIdx jl, EdgeIdx kl) {
        v66[i] += (e1[jk] - 2) + (e1[jl] - 2) + (e1[kl] - 2);
        v66[j] += (e1[ik] - 2) + (e1[il] - 2) + (e1[kl] - 2);
        v66[k] += (e1[ij] - 2) + (e1[il] - 2) + (e1[jl] - 2);
        v66[l] += (e1[ij] - 2) + (e1[ik] - 2) + (e1[jk] - 2);
    });
}

// Appends a ∩ b (both ascending) to out, galloping through the longer list when the lengths
// are far apart.
void intersectSortedInto(const VertexIdx* a, EdgeIdx na, const VertexIdx* b, EdgeIdx nb,
                         std::vector<VertexIdx>& out) {
    if (na > nb) { std::swap(a, b); std::swap(na, nb); }
    if (na == 0) return;
    const VertexIdx* a_end = a + na;
    const VertexIdx* b_end = b + nb;
    if (na * 32 < nb) {
        for (; a < a_end && b < b_end; ++a) {
            EdgeIdx step = 1;
            while (b + step < b_end && b[step] < *a) step <<= 1;
            b = std::lower_bound(b + (step >> 1), std::min(b + step + 1, b_end), *a);
            if (b < b_end && *b == *a) out.push_back(*a);
        }
        return;
    }
    while (a < a_end && b < b_end) {
        if (*a < *b) ++a;
        else if (*b < *a) ++b;
        else { out.push_back(*a); ++a; ++b; }
    }
}

// Escape's storeAllTriangles lists, built without binary searches, plus the edge positions the
// 5-vertex routines would otherwise binary-search for. Positions are into cg (sorted by ID, and
// ID order is the DAG's degree order, so cg[x] is gin[x] followed by gout[x]); an edge {x, y},
// x < y, is "canonical" at x's copy, as in Escape.
//   canon[p]:           canonical position of the edge at position p.
//   offsets[e]:         for canonical e = {x, y}, third[offsets[e]..offsets[e+1]) are the common
//                       neighbours of x and y, ascending (Escape's triangle list for e); empty
//                       for non-canonical positions.
//   low_edge/high_edge: for the entry c on {x, y}, the canonical positions of {x, c} and {y, c}.
struct EdgeTriangles {
    std::vector<EdgeIdx> canon;
    std::vector<EdgeIdx> offsets;
    std::vector<VertexIdx> third;
    std::vector<EdgeIdx> low_edge, high_edge;
};

EdgeTriangles buildEdgeTriangles(const CGraph* cg, const CGraph* gout, const CGraph* gin) {
    const VertexIdx n = cg->nVertices;
    EdgeTriangles et;
    et.canon.resize(cg->nEdges);
    {
        // A vertex's lower neighbours come first in its list, in the order this loop meets them.
        std::vector<EdgeIdx> cursor(cg->offsets, cg->offsets + n);
        for (VertexIdx x = 0; x < n; ++x) {
            for (EdgeIdx p = cg->offsets[x]; p < cg->offsets[x + 1]; ++p) {
                const VertexIdx y = cg->nbors[p];
                if (y > x) { et.canon[p] = p; et.canon[cursor[y]++] = p; }
            }
        }
    }
    auto canon_of = [&](VertexIdx a, EdgeIdx gout_pos) {
        return cg->offsets[a] + (gin->offsets[a + 1] - gin->offsets[a]) + (gout_pos - gout->offsets[a]);
    };
    // Every triangle a < b < c once, with its three canonical edges. Lowest vertex ascending, then
    // b, then c: each edge's list is filled in ascending order of its third vertex.
    std::vector<EdgeIdx> pos_in_a(n, -1);
    auto for_each_triangle = [&](auto&& f) {
        for (VertexIdx a = 0; a < n; ++a) {
            const EdgeIdx a_begin = gout->offsets[a], a_end = gout->offsets[a + 1];
            if (a_end - a_begin < 2) continue;
            for (EdgeIdx p = a_begin; p < a_end; ++p) pos_in_a[gout->nbors[p]] = p;
            for (EdgeIdx pb = a_begin; pb < a_end; ++pb) {
                const VertexIdx b = gout->nbors[pb];
                for (EdgeIdx pc = gout->offsets[b]; pc < gout->offsets[b + 1]; ++pc) {
                    const VertexIdx c = gout->nbors[pc];
                    if (pos_in_a[c] < 0) continue;
                    f(a, b, c, canon_of(a, pb), canon_of(a, pos_in_a[c]), canon_of(b, pc));
                }
            }
            for (EdgeIdx p = a_begin; p < a_end; ++p) pos_in_a[gout->nbors[p]] = -1;
        }
    };
    et.offsets.assign(cg->nEdges + 1, 0);
    for_each_triangle([&](VertexIdx, VertexIdx, VertexIdx, EdgeIdx ab, EdgeIdx ac, EdgeIdx bc) {
        ++et.offsets[ab + 1]; ++et.offsets[ac + 1]; ++et.offsets[bc + 1];
    });
    for (EdgeIdx e = 0; e < cg->nEdges; ++e) et.offsets[e + 1] += et.offsets[e];
    const EdgeIdx total = et.offsets[cg->nEdges];
    et.third.resize(total);
    et.low_edge.resize(total);
    et.high_edge.resize(total);
    std::vector<EdgeIdx> fill(et.offsets.begin(), et.offsets.end() - 1);
    for_each_triangle([&](VertexIdx a, VertexIdx b, VertexIdx c, EdgeIdx ab, EdgeIdx ac, EdgeIdx bc) {
        EdgeIdx t = fill[ab]++;
        et.third[t] = c; et.low_edge[t] = ac; et.high_edge[t] = bc;
        t = fill[ac]++;
        et.third[t] = b; et.low_edge[t] = ab; et.high_edge[t] = bc;
        t = fill[bc]++;
        et.third[t] = a; et.low_edge[t] = ab; et.high_edge[t] = ac;
    });
    return et;
}

// For each gin position q (k in gin[j]), the gout position of the same edge (j in gout[k]).
// gin lists are ascending, so walking the k in ascending order fills each gin[j] in order.
std::vector<EdgeIdx> ginToGoutPositions(const CGraph* gout, const CGraph* gin) {
    const VertexIdx n = gout->nVertices;
    std::vector<EdgeIdx> gin_to_gout(gin->nEdges);
    std::vector<EdgeIdx> cursor(gin->offsets, gin->offsets + n);
    for (VertexIdx k = 0; k < n; ++k)
        for (EdgeIdx p = gout->offsets[k]; p < gout->offsets[k + 1]; ++p)
            gin_to_gout[cursor[gout->nbors[p]]++] = p;
    return gin_to_gout;
}

// Same orbit-8 and per-edge 4-cycle counts (per_edge_[5]) as Escape's four_cycle_orbit_counter:
// the same wedges, but the edge positions it binary-searched (j -> i per in-neighbour, k -> j
// per in-in wedge) come from ginToGoutPositions, the additions per (i, j) are summed, and its
// third pass visits each distinct wedge end once (later visits added choosec(0, 2)).
void fourCycleOrbitsWithEdges(const CGraph* gout, const CGraph* gin, OrbitInfo& oc) {
    const VertexIdx n = gout->nVertices;
    double* v8 = oc.per_vertex_[8];
    double* e5 = oc.per_edge_[5];
    const std::vector<EdgeIdx> gin_to_gout = ginToGoutPositions(gout, gin);
    std::vector<VertexIdx> wedge_count(n, 0);
    std::vector<VertexIdx> ends;
    for (VertexIdx i = 0; i < n; ++i) {
        ends.clear();
        for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i + 1]; ++pos) {
            const VertexIdx j = gin->nbors[pos];
            for (EdgeIdx e = gout->offsets[j]; e < gout->offsets[j + 1]; ++e) {
                const VertexIdx k = gout->nbors[e];
                if (k >= i) continue;
                if (wedge_count[k]++ == 0) ends.push_back(k);
            }
            for (EdgeIdx e = gin->offsets[j]; e < gin->offsets[j + 1]; ++e) {
                const VertexIdx k = gin->nbors[e];
                if (wedge_count[k]++ == 0) ends.push_back(k);
            }
        }
        for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i + 1]; ++pos) {
            const VertexIdx j = gin->nbors[pos];
            VertexIdx at_j = 0;
            for (EdgeIdx e = gout->offsets[j]; e < gout->offsets[j + 1]; ++e) {
                const VertexIdx k = gout->nbors[e];
                if (k >= i) continue;
                const VertexIdx c = std::max<VertexIdx>(0, wedge_count[k] - 1);
                at_j += c;
                e5[e] += c;
            }
            for (EdgeIdx e = gin->offsets[j]; e < gin->offsets[j + 1]; ++e) {
                const VertexIdx c = std::max<VertexIdx>(0, wedge_count[gin->nbors[e]] - 1);
                at_j += c;
                e5[gin_to_gout[e]] += c;
            }
            v8[j] += at_j;
            e5[gin_to_gout[pos]] += at_j;
        }
        for (VertexIdx k : ends) {
            const double pairs = choosec(wedge_count[k], 2);
            v8[i] += pairs;
            v8[k] += pairs;
            wedge_count[k] = 0;
        }
    }
}

// Same orbit-34 counts as Escape's FiveCycleOrbitCounter: same wedges, three-paths and
// corrections, but its four binary searches per three-path (and one per wedge) become lookups:
// marks of N(i) and N(j), and each gin edge's gout position. Its unused type-1/type-2 wedge
// arrays are dropped, and the counts it cleared by re-walking every three-path are cleared from
// a list of the ends touched.
void fiveCycleOrbits(const CGraph* cg, const CGraph* gout, const CGraph* gin, double* c34) {
    const VertexIdx n = gout->nVertices;
    const std::vector<EdgeIdx> gin_to_gout = ginToGoutPositions(gout, gin);
    // tri_inin[e] for gout edge e = (j, k): the number of i with i -> j and i -> k.
    std::vector<VertexIdx> tri_inin(gout->nEdges, 0);
    {
        std::vector<char> in_out_i(n, 0);
        for (VertexIdx i = 0; i < n; ++i) {
            const EdgeIdx i_begin = gout->offsets[i], i_end = gout->offsets[i + 1];
            if (i_end - i_begin < 2) continue;
            for (EdgeIdx p = i_begin; p < i_end; ++p) in_out_i[gout->nbors[p]] = 1;
            for (EdgeIdx pj = i_begin; pj < i_end; ++pj) {
                const VertexIdx j = gout->nbors[pj];
                for (EdgeIdx p = gout->offsets[j]; p < gout->offsets[j + 1]; ++p)
                    if (in_out_i[gout->nbors[p]]) ++tri_inin[p];
            }
            for (EdgeIdx p = i_begin; p < i_end; ++p) in_out_i[gout->nbors[p]] = 0;
        }
    }

    std::vector<VertexIdx> wedge_count(n, 0), wedge_count_type3(n, 0), three_path_counts(n, 0);
    std::vector<VertexIdx> nbr_of_i(n, -1);   // == i for neighbours of i
    std::vector<EdgeIdx> nbr_of_j(n, -1);     // == gin position of (j in gin[i]) for neighbours of j
    std::vector<VertexIdx> wedge_ends, path_ends;
    for (VertexIdx i = 0; i < n; ++i) {
        for (EdgeIdx p = cg->offsets[i]; p < cg->offsets[i + 1]; ++p) nbr_of_i[cg->nbors[p]] = i;

        wedge_ends.clear();
        for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i + 1]; ++pos) {
            const VertexIdx j = gin->nbors[pos];
            for (EdgeIdx q = gin->offsets[j]; q < gin->offsets[j + 1]; ++q) {   // i <- j <- k
                const VertexIdx k = gin->nbors[q];
                ++wedge_count_type3[k];
                if (wedge_count[k]++ == 0) wedge_ends.push_back(k);
            }
            for (EdgeIdx q = gout->offsets[j]; q < gout->offsets[j + 1]; ++q) { // i <- j -> k
                const VertexIdx k = gout->nbors[q];
                if (k == i) continue;
                if (wedge_count[k]++ == 0) wedge_ends.push_back(k);
            }
        }
        for (EdgeIdx pos = gout->offsets[i]; pos < gout->offsets[i + 1]; ++pos) {  // i -> j -> k
            const VertexIdx j = gout->nbors[pos];
            for (EdgeIdx q = gout->offsets[j]; q < gout->offsets[j + 1]; ++q) {
                const VertexIdx k = gout->nbors[q];
                if (wedge_count[k]++ == 0) wedge_ends.push_back(k);
            }
        }

        // Three-paths i <- j <- k -> ell.
        path_ends.clear();
        for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i + 1]; ++pos) {
            const VertexIdx j = gin->nbors[pos];
            if (gin->offsets[j] == gin->offsets[j + 1]) continue;
            for (EdgeIdx p = cg->offsets[j]; p < cg->offsets[j + 1]; ++p) nbr_of_j[cg->nbors[p]] = pos;
            VertexIdx sum_ij = 0;
            for (EdgeIdx q = gin->offsets[j]; q < gin->offsets[j + 1]; ++q) {
                const VertexIdx k = gin->nbors[q];
                const VertexIdx k_adj_i = nbr_of_i[k] == i;
                VertexIdx sum_k = 0;
                for (EdgeIdx r = gout->offsets[k]; r < gout->offsets[k + 1]; ++r) {
                    const VertexIdx ell = gout->nbors[r];
                    if (ell == j || ell == i) continue;
                    const VertexIdx addition = wedge_count[ell] - k_adj_i - (nbr_of_j[ell] == pos);
                    if (three_path_counts[ell]++ == 0) path_ends.push_back(ell);
                    sum_k += addition;
                    c34[ell] += addition;
                }
                c34[k] += sum_k;
                sum_ij += sum_k;
            }
            c34[i] += sum_ij;
            c34[j] += sum_ij;
        }

        // The middle node of each wedge, from the three-path counts.
        for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i + 1]; ++pos) {
            const VertexIdx j = gin->nbors[pos];
            VertexIdx add = 0;
            for (EdgeIdx q = gin->offsets[j]; q < gin->offsets[j + 1]; ++q)       // i <- j <- k
                add += three_path_counts[gin->nbors[q]] - tri_inin[gin_to_gout[q]];
            for (EdgeIdx q = gout->offsets[j]; q < gout->offsets[j + 1]; ++q) {   // i <- j -> k
                const VertexIdx k = gout->nbors[q];
                if (k == i) continue;
                const VertexIdx k_to_i = k < i && nbr_of_i[k] == i;
                add += three_path_counts[k] - (wedge_count_type3[j] - k_to_i) - tri_inin[q];
            }
            c34[j] += add;
        }
        for (EdgeIdx pos = gout->offsets[i]; pos < gout->offsets[i + 1]; ++pos) {  // i -> j -> k
            const VertexIdx j = gout->nbors[pos];
            VertexIdx add = 0;
            for (EdgeIdx q = gout->offsets[j]; q < gout->offsets[j + 1]; ++q) add += three_path_counts[gout->nbors[q]];
            c34[j] += add;
        }

        for (VertexIdx k : wedge_ends) { wedge_count[k] = 0; wedge_count_type3[k] = 0; }
        for (VertexIdx ell : path_ends) three_path_counts[ell] = 0;
    }
}

// Same orbit 62, 63, 64, 68, 69 counts as Escape's ChordalWedgeColOrbitCounter: the same
// diamonds (i, j, k, ell), ell > i, via the triangle lists, but
//  - the position of {j, k} comes from the entry (no binary search), and the ells > i of its
//    list from one upper_bound; the second pass over the diamonds reuses those ranges;
//  - Escape counted wedge_count[ell] = |N(i) ∩ N(ell)| by walking every wedge i - j - k (all of
//    N(j) for every neighbour j; three such passes) but read it only at diamond ends (the ells
//    with a j < k diamond, "ends" below). Here, per i, whichever is cheaper: one pass over the
//    wedges with k > i (also doing orbit 62's pass), or intersecting N(i) with N(ell) per end;
//  - the per-diamond additions to 63, 64 and the mid counts are summed per end and per range;
//  - parts whose orbits no pattern vertex has are skipped (do62, do6364, do6869).
void chordalWedgeOrbits(const CGraph* cg, const EdgeTriangles& et, OrbitInfo& oc,
                        bool do62, bool do6364, bool do6869) {
    const VertexIdx n = cg->nVertices;
    double* r62 = oc.per_vertex_[62];
    double* r63 = oc.per_vertex_[63];
    double* r64 = oc.per_vertex_[64];
    double* r68 = oc.per_vertex_[68];
    double* r69 = oc.per_vertex_[69];
    const VertexIdx* third = et.third.data();
    const bool need_ends = do62 || do6364;

    std::vector<VertexIdx> diamond_count(n, 0), diamond_count_final(n, 0), diamond_count_i(n, 0),
        mid_diamond_count_i(n, 0);
    std::vector<VertexIdx> end_stamp(n, -1);     // == i for the ends of i
    std::vector<VertexIdx> wedge_stamp(n, -1);   // == i once wedge_count[] is for i
    std::vector<VertexIdx> wedge_count(n, 0);
    std::vector<VertexIdx> ends;
    std::vector<EdgeIdx> above_i;                // per neighbour j: first position in N(j) past i
    std::vector<std::pair<EdgeIdx, EdgeIdx>> ranges;
    struct MidRange { VertexIdx j, k; EdgeIdx begin, end; };
    std::vector<MidRange> mid_ranges;            // the j < k diamonds of i, for orbit 64

    for (VertexIdx i = 0; i < n; ++i) {
        const VertexIdx* n_i = cg->nbors + cg->offsets[i];
        const EdgeIdx deg_i = cg->offsets[i + 1] - cg->offsets[i];
        ends.clear();
        mid_ranges.clear();
        for (EdgeIdx pos = cg->offsets[i]; pos < cg->offsets[i + 1]; ++pos) {
            const VertexIdx j = cg->nbors[pos];
            const EdgeIdx e_ij = et.canon[pos];
            const EdgeIdx t_begin = et.offsets[e_ij], t_end = et.offsets[e_ij + 1];
            if (t_begin == t_end) continue;
            const EdgeIdx* jk_edge = (j < i ? et.low_edge.data() : et.high_edge.data());
            ranges.clear();
            VertexIdx mid_j = 0;
            for (EdgeIdx t = t_begin; t < t_end; ++t) {
                const VertexIdx k = third[t];
                const EdgeIdx e_jk = jk_edge[t];
                const EdgeIdx l_end = et.offsets[e_jk + 1];
                const EdgeIdx l_begin = std::upper_bound(third + et.offsets[e_jk], third + l_end, i) - third;
                ranges.emplace_back(l_begin, l_end);
                if (do6869) {
                    for (EdgeIdx s = l_begin; s < l_end; ++s) ++diamond_count[third[s]];
                }
                if (!need_ends || j >= k || l_begin == l_end) continue;
                for (EdgeIdx s = l_begin; s < l_end; ++s) {
                    const VertexIdx ell = third[s];
                    ++diamond_count_i[ell];
                    if (end_stamp[ell] != i) { end_stamp[ell] = i; ends.push_back(ell); }
                }
                mid_j += l_end - l_begin;
                mid_diamond_count_i[k] += l_end - l_begin;
                if (do6364) mid_ranges.push_back({j, k, l_begin, l_end});
            }
            mid_diamond_count_i[j] += mid_j;
            if (!do6869) continue;
            // Escape's second pass: C(count, 2) at the first visit of each ell (it then zeroes
            // the count), and count - 1 to k for every (k, ell).
            for (EdgeIdx t = t_begin; t < t_end; ++t) {
                const auto& range = ranges[t - t_begin];
                VertexIdx to_k = 0;
                for (EdgeIdx s = range.first; s < range.second; ++s) {
                    const VertexIdx ell = third[s];
                    if (diamond_count[ell] != 0) {
                        diamond_count_final[ell] = diamond_count[ell];
                        const double pairs = choosec(diamond_count[ell], 2);
                        r69[j] += pairs;
                        r68[i] += pairs;
                        r68[ell] += pairs;
                        diamond_count[ell] = 0;
                    }
                    to_k += diamond_count_final[ell] - 1;
                }
                r68[third[t]] += to_k;
            }
        }
        if (ends.empty()) continue;

        // wedge_count[ell] for the ends; orbit 62 gets diamond_count_i[k] at j for every wedge
        // i - j - k with k > i (non-zero only at ends k).
        double wedge_cost = 0, intersect_cost = 0;
        above_i.resize(deg_i);
        for (EdgeIdx p = 0; p < deg_i; ++p) {
            const VertexIdx j = n_i[p];
            const VertexIdx* nj_begin = cg->nbors + cg->offsets[j];
            const VertexIdx* nj_end = cg->nbors + cg->offsets[j + 1];
            above_i[p] = std::upper_bound(nj_begin, nj_end, i) - cg->nbors;
            wedge_cost += cg->offsets[j + 1] - above_i[p];
        }
        for (VertexIdx ell : ends) {
            const double d_ell = cg->offsets[ell + 1] - cg->offsets[ell];
            const double lo = std::min<double>(deg_i, d_ell), hi = std::max<double>(deg_i, d_ell);
            intersect_cost += (lo * 32 < hi) ? lo * (1 + std::log2(hi / lo)) : lo + hi;
        }
        if (wedge_cost <= intersect_cost) {
            for (EdgeIdx p = 0; p < deg_i; ++p) {
                const VertexIdx j = n_i[p];
                VertexIdx to_j = 0;
                for (EdgeIdx q = above_i[p]; q < cg->offsets[j + 1]; ++q) {
                    const VertexIdx k = cg->nbors[q];
                    if (wedge_stamp[k] != i) { wedge_stamp[k] = i; wedge_count[k] = 0; }
                    ++wedge_count[k];
                    to_j += diamond_count_i[k];
                }
                if (do62) r62[j] += to_j;
            }
        } else {
            for (VertexIdx ell : ends) {
                const VertexIdx* n_ell = cg->nbors + cg->offsets[ell];
                const VertexIdx* ne = cg->nbors + cg->offsets[ell + 1];
                const VertexIdx* a = n_i;
                const VertexIdx* a_end = n_i + deg_i;
                const VertexIdx d = diamond_count_i[ell];
                VertexIdx common = 0;
                auto hit = [&](VertexIdx x) { ++common; if (do62) r62[x] += d; };
                // Galloping through the longer list when the lengths are far apart.
                const VertexIdx* b = n_ell;
                const VertexIdx* b_end = ne;
                if (a_end - a > b_end - b) { std::swap(a, b); std::swap(a_end, b_end); }
                if ((a_end - a) * 32 < (b_end - b)) {
                    for (; a < a_end && b < b_end; ++a) {
                        EdgeIdx step = 1;
                        while (b + step < b_end && b[step] < *a) step <<= 1;
                        b = std::lower_bound(b + (step >> 1), std::min(b + step + 1, b_end), *a);
                        if (b < b_end && *b == *a) hit(*a);
                    }
                } else {
                    while (a < a_end && b < b_end) {
                        if (*a < *b) ++a;
                        else if (*b < *a) ++b;
                        else { hit(*a); ++a; ++b; }
                    }
                }
                wedge_stamp[ell] = i;
                wedge_count[ell] = common;
            }
        }

        if (do6364) {
            for (VertexIdx ell : ends) {
                if (wedge_count[ell] <= 2) { wedge_count[ell] = 0; continue; }
                const VertexIdx w = wedge_count[ell] - 2;
                r63[i] += static_cast<double>(w) * diamond_count_i[ell];
                r63[ell] += static_cast<double>(w) * diamond_count_i[ell];
                wedge_count[ell] = w;   // reused below as the per-diamond amount
            }
            for (const MidRange& r : mid_ranges) {
                VertexIdx sum = 0;
                for (EdgeIdx s = r.begin; s < r.end; ++s) sum += wedge_count[third[s]];
                r64[r.j] += sum;
                r64[r.k] += sum;
            }
        }
        for (VertexIdx ell : ends) diamond_count_i[ell] = 0;
        if (do62) {
            for (EdgeIdx p = 0; p < deg_i; ++p) {
                const VertexIdx j = n_i[p];
                r62[j] -= mid_diamond_count_i[j];
                mid_diamond_count_i[j] = 0;
            }
        } else {
            for (EdgeIdx p = 0; p < deg_i; ++p) mid_diamond_count_i[n_i[p]] = 0;
        }
    }

    // Escape: orbits 68 and 69 are double counted.
    for (VertexIdx i = 0; i < n; ++i) {
        r68[i] /= 2;
        r69[i] /= 2;
    }
}

// Same orbit 70, 71 counts as Escape's AlmostFiveCliqueOrbitCounter: for each edge i < j and
// each nbr in T = N(i) ∩ N(j) with i < nbr < j, the x in T adjacent to nbr. Escape rebuilt T with
// a binary search per neighbour of i and tested each x twice by binary search; here T is the
// edge's triangle list and {x in T adjacent to nbr} = T ∩ (triangle list of {i, nbr}).
void almostFiveCliqueOrbits(const CGraph* cg, const EdgeTriangles& et, OrbitInfo& oc) {
    const VertexIdx n = cg->nVertices;
    double* r70 = oc.per_vertex_[70];
    double* r71 = oc.per_vertex_[71];
    const VertexIdx* third = et.third.data();
    std::vector<EdgeIdx> pos_from_i(n, 0);   // cg position of {i, v} in i's list, for v > i
    std::vector<VertexIdx> cliques;
    for (VertexIdx i = 0; i < n; ++i) {
        const EdgeIdx i_end = cg->offsets[i + 1];
        const EdgeIdx up_begin = std::upper_bound(cg->nbors + cg->offsets[i], cg->nbors + i_end, i) - cg->nbors;
        for (EdgeIdx p = up_begin; p < i_end; ++p) pos_from_i[cg->nbors[p]] = p;
        for (EdgeIdx p = up_begin; p < i_end; ++p) {
            const VertexIdx j = cg->nbors[p];
            const VertexIdx* t_begin = third + et.offsets[p];
            const VertexIdx* t_end = third + et.offsets[p + 1];
            if (t_end - t_begin < 2) continue;
            const VertexIdx* nbr_end = std::lower_bound(t_begin, t_end, j);
            for (const VertexIdx* it = std::upper_bound(t_begin, nbr_end, i); it < nbr_end; ++it) {
                const VertexIdx nbr = *it;
                const EdgeIdx e = pos_from_i[nbr];
                cliques.clear();
                intersectSortedInto(t_begin, t_end - t_begin, third + et.offsets[e],
                                    et.offsets[e + 1] - et.offsets[e], cliques);
                const VertexIdx fourclique = cliques.size();
                for (VertexIdx x : cliques) r70[x] += fourclique - 1;
                const double pairs = choosec(fourclique, 2);
                r71[i] += pairs;
                r71[j] += pairs;
                r71[nbr] += pairs;
            }
        }
    }
}
}

OrbitCounts EVOKEOrbitCounter::count() {
    OrbitCounts result;
    countWith([&](double* const* per_vertex, const VertexIdx* mapping, VertexIdx n, int num_orbits) {
        result.assign(n, {});
        for (VertexIdx i = 0; i < n; ++i) {
            std::vector<long long>& counts = result[i];
            counts.reserve(num_orbits);
            for (int orbit = 0; orbit < num_orbits; ++orbit) {
                counts.push_back(per_vertex[orbit][mapping[i]]);
            }
        }
    });
    return result;
}

OrbitCounts EVOKEOrbitCounter::countColumns(const std::vector<int>& orbits) {
    OrbitCounts result;
    countWith([&](double* const* per_vertex, const VertexIdx* mapping, VertexIdx n, int) {
        result.assign(orbits.size(), std::vector<long long>(n));
        for (size_t j = 0; j < orbits.size(); ++j) {
            const double* column = per_vertex[orbits[j]];
            long long* out = result[j].data();
            for (VertexIdx i = 0; i < n; ++i) out[i] = column[mapping[i]];
        }
    });
    return result;
}

void EVOKEOrbitCounter::countWith(const Extract& extract) {
    last_timings = OrbitCountTimings();
    const CGraph* g_raw = cgraph;
    if (!g_raw) return;

    auto t = Clock::now();
    VertexIdx* mapping;
    VertexIdx* inverse;
    CGraph g_relabel = g_raw->renameByDegreeOrder(mapping, inverse);
    g_relabel.sortById();

    CDAG dag = degreeOrdered(&g_relabel);
    dag.outlist.sortById();
    dag.inlist.sortById();
    last_timings.relabel_dag_ms = msSince(t);

    // Caps come indexed by the graph's internal index; the counters run on relabeled IDs.
    auto relabelCap = [&](const std::vector<double>* cap) {
        std::vector<double> relabeled;
        if (cap) {
            relabeled.assign(g_relabel.nVertices, 0);
            for (VertexIdx i = 0; i < g_relabel.nVertices; ++i) relabeled[mapping[i]] = (*cap)[i];
        }
        return relabeled;
    };

    int num_node_orbits = 0;
    int num_edge_orbits = 0;

    if (graphlet_size == 3) {
        num_node_orbits = 4;
        num_edge_orbits = 3;
    } else if (graphlet_size == 4) {
        num_node_orbits = 15;
        num_edge_orbits = 12;
    } else if (graphlet_size == 5) {
        num_node_orbits = 73;
        num_edge_orbits = 68;
    } else {
        throw std::invalid_argument("Only graphlet sizes 3, 4, and 5 are supported.");
    }
    
    OrbitInfo orbit_counts(&(dag.outlist), num_node_orbits, num_edge_orbits);
    
    t = Clock::now();
    // Triangles feed orbit 3 and everything derived from it (4, 5, 11 via closed forms,
    // paw/diamond via per-edge counts); with none of those needed they can be skipped.
    const bool count_triangles = needed(3) || (graphlet_size >= 4 && (needed(4) || needed(5) ||
        needed(9) || needed(10) || needed(11) || needed(12) || needed(13))) || graphlet_size == 5;
    threeVertexOrbitsByMarking(&g_relabel, &(dag.outlist), orbit_counts, count_triangles);
    last_timings.three_ms = msSince(t);

    if (graphlet_size == 5) {
        // The bodies of Escape's get_all_four_orbits and get_all_five_orbits (with its parallel
        // sections), expanded so that each routine and each parallel section can be timed and
        // skipped; same counts. The 5-cycle, chordal-wedge, almost-5-clique and 5-clique routines
        // are the reimplementations above.
        CGraph* cg = &g_relabel;
        CGraph* gout = &(dag.outlist);
        CGraph* gin = &(dag.inlist);

        t = Clock::now();
        for (VertexIdx i = 0; i < cg->nVertices; i++) {
            orbit_counts.per_vertex_[5][i] = orbit_counts.per_vertex_[1][i] * (orbit_counts.per_vertex_[0][i] - 1) - 2 * orbit_counts.per_vertex_[3][i];
            orbit_counts.per_vertex_[7][i] = choosec(orbit_counts.per_vertex_[0][i], 3);
            orbit_counts.per_vertex_[11][i] = orbit_counts.per_vertex_[3][i] * (orbit_counts.per_vertex_[0][i] - 2);
            for (EdgeIdx pos = cg->offsets[i]; pos < cg->offsets[i + 1]; pos++) {
                VertexIdx j = cg->nbors[pos];
                VertexIdx deg_j = orbit_counts.per_vertex_[0][j];
                orbit_counts.per_vertex_[4][i] += orbit_counts.per_vertex_[1][j];
                orbit_counts.per_vertex_[6][i] += choosec(deg_j - 1, 2);
            }
            orbit_counts.per_vertex_[4][i] -= 2 * orbit_counts.per_vertex_[2][i] + 2 * orbit_counts.per_vertex_[3][i];
        }
        last_timings.four_closed_ms = msSince(t);

        // Each 5-vertex routine writes only its own vertex orbits (and no edge orbits), and the
        // parallel ones can't read each other's, so one whose orbits no pattern vertex has is skipped.
        auto any_needed = [&](std::initializer_list<int> orbits) {
            for (int o : orbits) if (needed(o)) return true;
            return false;
        };
        const bool run_tree = any_needed({15, 16, 17, 18, 19, 20, 21, 22, 23});
        const bool run_cycle = needed(34);
        const bool run_wedge = any_needed({49, 50});
        const bool run_p51 = needed(51);
        const bool run_other = any_needed({24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 35, 37, 38, 39, 40, 41, 42,
                                           43, 44, 45, 46, 47, 48, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 65, 67});
        const bool run_p36 = needed(36);
        const bool run_chordal = any_needed({62, 63, 64, 68, 69});
        const bool run_p66 = needed(66);
        const bool run_almost = any_needed({70, 71});
        const bool run_clique = needed(72);

        t = Clock::now();
        otherFourPatternOrbitCounter(gout, gin, orbit_counts);
        last_timings.four_other_ms = msSince(t);

        // Of the 5-vertex routines only the other-five one reads the per-edge 4-cycle and
        // 4-clique counts and orbit 14 (the tree one reads orbit 8). Without it, only the vertex
        // counts a pattern vertex (or the tree routine) needs are computed, orbit 14 capped.
        t = Clock::now();
        if (run_other) fourCycleOrbitsWithEdges(gout, gin, orbit_counts);
        else if (run_tree || needed(8)) fourCycleOrbits(gout, gin, orbit_counts.per_vertex_[8]);
        last_timings.four_cycle_ms = msSince(t);

        t = Clock::now();
        if (run_other) {
            fourCliqueOrbitsWithEdges(gout, orbit_counts);
        } else if (needed(14)) {
            std::vector<double> cap_relabeled = relabelCap(clique4_cap);
            fourCliqueOrbitsByMarking(gout, orbit_counts.per_vertex_[14], clique4_cap ? &cap_relabeled : nullptr);
        }
        last_timings.four_clique_ms = msSince(t);

        // The per-edge triangle lists are only read by the chordal-wedge and almost-5-clique routines.
        t = Clock::now();
        EdgeTriangles edge_tris;
        if (run_chordal || run_almost) edge_tris = buildEdgeTriangles(cg, gout, gin);
        last_timings.five_setup_ms = msSince(t);

        t = Clock::now();
        if (run_tree) FiveTreeOrbitCounter(gout, gin, orbit_counts);
        last_timings.five_tree_ms = msSince(t);

        std::vector<double> clique5_cap_relabeled = relabelCap(clique5_cap);
        t = Clock::now();
        double* section_ms = last_timings.five_section_ms;
        #pragma omp parallel sections
        {
            #pragma omp section
            { auto s = Clock::now(); if (run_cycle) fiveCycleOrbits(cg, gout, gin, orbit_counts.per_vertex_[34]); section_ms[0] = msSince(s); }
            #pragma omp section
            { auto s = Clock::now(); if (run_wedge) WedgeCollisionOrbitCounter(cg, gout, gin, orbit_counts); section_ms[1] = msSince(s); }
            #pragma omp section
            { auto s = Clock::now(); if (run_p51) Pattern51OrbitCounter(cg, gout, gin, orbit_counts); section_ms[2] = msSince(s); }
            #pragma omp section
            { auto s = Clock::now(); if (run_other) OtherFiveOrbitsCounter(gout, gin, orbit_counts); section_ms[3] = msSince(s); }
            #pragma omp section
            { auto s = Clock::now(); if (run_p36) Pattern36OrbitCounter(gout, gin, orbit_counts); section_ms[4] = msSince(s); }
            #pragma omp section
            { auto s = Clock::now(); if (run_chordal) chordalWedgeOrbits(cg, edge_tris, orbit_counts, needed(62), needed(63) || needed(64), needed(68) || needed(69)); section_ms[5] = msSince(s); }
            #pragma omp section
            { auto s = Clock::now(); if (run_p66) pattern66Orbits(gout, orbit_counts); section_ms[6] = msSince(s); }
            #pragma omp section
            { auto s = Clock::now(); if (run_almost) almostFiveCliqueOrbits(cg, edge_tris, orbit_counts); section_ms[7] = msSince(s); }
            #pragma omp section
            {
                auto s = Clock::now();
                if (run_clique) {
                    fiveCliqueOrbitsByMarking(gout, orbit_counts.per_vertex_[72],
                                              clique5_cap ? &clique5_cap_relabeled : nullptr);
                }
                section_ms[8] = msSince(s);
            }
        }
        last_timings.five_parallel_ms = msSince(t);
    } else if (graphlet_size == 4) {
        // Body of Escape's get_all_four_orbits_node_only, expanded so each routine can be
        // timed separately; same routines, same order, same counts.
        t = Clock::now();
        CGraph* cg = &g_relabel;
        for (VertexIdx i = 0; i < cg->nVertices; i++) {
            orbit_counts.per_vertex_[5][i] = orbit_counts.per_vertex_[1][i] * (orbit_counts.per_vertex_[0][i] - 1) - 2 * orbit_counts.per_vertex_[3][i];
            orbit_counts.per_vertex_[7][i] = choosec(orbit_counts.per_vertex_[0][i], 3);
            orbit_counts.per_vertex_[11][i] = orbit_counts.per_vertex_[3][i] * (orbit_counts.per_vertex_[0][i] - 2);
            for (EdgeIdx pos = cg->offsets[i]; pos < cg->offsets[i + 1]; pos++) {
                VertexIdx j = cg->nbors[pos];
                VertexIdx deg_j = orbit_counts.per_vertex_[0][j];
                orbit_counts.per_vertex_[4][i] += orbit_counts.per_vertex_[1][j];
                orbit_counts.per_vertex_[6][i] += choosec(deg_j - 1, 2);
            }
            orbit_counts.per_vertex_[4][i] -= 2 * orbit_counts.per_vertex_[2][i] + 2 * orbit_counts.per_vertex_[3][i];
        }
        last_timings.four_closed_ms = msSince(t);

        // Each routine only writes its own orbits (paw/diamond: 9,10,12,13; cycle: 8; clique: 14),
        // so one whose orbits no pattern vertex has can be skipped: its zeros never reject.
        t = Clock::now();
        if (needed(9) || needed(10) || needed(12) || needed(13)) {
            pawDiamondOrbits(&(dag.outlist), orbit_counts);
        }
        last_timings.four_other_ms = msSince(t);

        t = Clock::now();
        if (needed(8)) {
            fourCycleOrbits(&(dag.outlist), &(dag.inlist), orbit_counts.per_vertex_[8]);
        }
        last_timings.four_cycle_ms = msSince(t);

        t = Clock::now();
        if (needed(14)) {
            std::vector<double> cap_relabeled = relabelCap(clique4_cap);
            fourCliqueOrbitsByMarking(&(dag.outlist), orbit_counts.per_vertex_[14],
                                      clique4_cap ? &cap_relabeled : nullptr);
        }
        last_timings.four_clique_ms = msSince(t);
    }

    t = Clock::now();
    extract(orbit_counts.per_vertex_, mapping, g_relabel.nVertices, num_node_orbits);

    delete[] mapping;
    delete[] inverse;
    last_timings.copy_ms = msSince(t);
}