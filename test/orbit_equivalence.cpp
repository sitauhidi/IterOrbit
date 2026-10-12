// Exactness test for the reimplemented orbit routines in src/OrbitCounter.cpp: runs
// Escape's original routines and the rewrites on the same graphs and requires identical
// per-vertex counts (orbits 8, 14, 34, 62-64, 66, 68-71), identical per-edge 4-clique counts
// and identical per-edge triangle lists.
//
// Usage: orbit_equivalence <graph file> <max vertices> <seeds>
//   Tests BFS balls of up to <max vertices> vertices grown from the highest-degree vertex and
//   from <seeds> - 1 random vertices; <max vertices> = 0 tests the whole graph once.
// Build: g++ -std=c++17 -O3 -fopenmp -Iinclude -isystem extern test/orbit_equivalence.cpp \
//            src/AppGraph.cpp extern/Escape/*.cpp -o build/orbit_equivalence
#include "../src/OrbitCounter.cpp"
#include <cstdlib>
#include <queue>
#include <random>
#include <unordered_set>

namespace {
AppGraph bfsBall(const AppGraph& g, int seed_id, size_t max_vertices) {
    std::unordered_set<int> nodes{seed_id};
    std::queue<int> frontier;
    frontier.push(seed_id);
    while (!frontier.empty() && nodes.size() < max_vertices) {
        const int u = frontier.front();
        frontier.pop();
        for (int v : g.getNeighbors(u)) {
            if (nodes.size() >= max_vertices) break;
            if (nodes.insert(v).second) frontier.push(v);
        }
    }
    return AppGraph::createSubgraph(g, nodes);
}

bool sameCounts(const char* what, const double* a, const double* b, VertexIdx n) {
    for (VertexIdx v = 0; v < n; ++v) {
        if (a[v] != b[v]) {
            std::cout << "  MISMATCH " << what << " at vertex " << v << ": escape " << a[v]
                      << " rewrite " << b[v] << "\n";
            return false;
        }
    }
    return true;
}

bool testGraph(const AppGraph& app) {
    CGraph* g_raw = app.c_graph;
    VertexIdx* mapping = nullptr;
    VertexIdx* inverse = nullptr;
    CGraph cg = g_raw->renameByDegreeOrder(mapping, inverse);
    cg.sortById();
    CDAG dag = degreeOrdered(&cg);
    dag.outlist.sortById();
    dag.inlist.sortById();
    CGraph* gout = &dag.outlist;
    CGraph* gin = &dag.inlist;
    const VertexIdx n = cg.nVertices;
    std::cout << "  n=" << n << " m=" << cg.nEdges / 2 << std::flush;

    OrbitInfo escape(gout, 73, 12), rewrite(gout, 73, 12);
    bool ok = true;

    // Triangle lists.
    auto t = Clock::now();
    TriangleInfo info = betterWedgeEnumerator(gout);
    TriangleList all_tris = storeAllTriangles(&cg, info.total);
    const double escape_tris_ms = msSince(t);
    t = Clock::now();
    EdgeTriangles et = buildEdgeTriangles(&cg, gout, gin);
    const double rewrite_tris_ms = msSince(t);
    if (static_cast<EdgeIdx>(et.third.size()) != all_tris.total) {
        std::cout << "\n  MISMATCH triangle total " << all_tris.total << " vs " << et.third.size() << "\n";
        ok = false;
    } else {
        for (EdgeIdx p = 0; p <= cg.nEdges && ok; ++p) {
            if (et.offsets[p] != all_tris.trioffsets[p]) {
                std::cout << "\n  MISMATCH trioffsets at " << p << "\n";
                ok = false;
            }
        }
        for (EdgeIdx s = 0; s < all_tris.total && ok; ++s) {
            if (et.third[s] != all_tris.triangles[s]) {
                std::cout << "\n  MISMATCH triangle entry " << s << "\n";
                ok = false;
            }
        }
    }

    t = Clock::now();
    FiveCycleOrbitCounter(gout, gin, escape);
    const double escape_cycle_ms = msSince(t);
    t = Clock::now();
    fiveCycleOrbits(&cg, gout, gin, rewrite.per_vertex_[34]);
    const double rewrite_cycle_ms = msSince(t);

    t = Clock::now();
    ChordalWedgeColOrbitCounter(&cg, gout, gin, &all_tris, escape);
    const double escape_chordal_ms = msSince(t);
    t = Clock::now();
    chordalWedgeOrbits(&cg, et, rewrite, true, true, true);
    const double rewrite_chordal_ms = msSince(t);

    t = Clock::now();
    AlmostFiveCliqueOrbitCounter(&cg, gout, gin, escape);
    const double escape_almost_ms = msSince(t);
    t = Clock::now();
    almostFiveCliqueOrbits(&cg, et, rewrite);
    const double rewrite_almost_ms = msSince(t);

    // 4-vertex routines the g=5 branch replaces: 4-clique with per-edge counts, node-only
    // 4-cycle (vs the per-edge version), and pattern 66 (reads per_edge_[1] from the
    // other-four routine, run first on both sides).
    otherFourPatternOrbitCounter(gout, gin, escape);
    otherFourPatternOrbitCounter(gout, gin, rewrite);
    t = Clock::now();
    fourCliqueOrbitCounter(gout, escape);
    const double escape_k4_ms = msSince(t);
    t = Clock::now();
    fourCliqueOrbitsWithEdges(gout, rewrite);
    const double rewrite_k4_ms = msSince(t);
    t = Clock::now();
    four_cycle_orbit_counter(gout, gin, escape);
    const double escape_c4_ms = msSince(t);
    t = Clock::now();
    fourCycleOrbits(gout, gin, rewrite.per_vertex_[8]);
    const double rewrite_c4_ms = msSince(t);
    OrbitInfo rewrite_c4e(gout, 73, 12);
    t = Clock::now();
    fourCycleOrbitsWithEdges(gout, gin, rewrite_c4e);
    const double rewrite_c4e_ms = msSince(t);
    ok &= sameCounts("orbit 8 (4-cycle with edges)", escape.per_vertex_[8], rewrite_c4e.per_vertex_[8], n);
    for (EdgeIdx e = 0; e < gout->nEdges && ok; ++e) {
        if (escape.per_edge_[5][e] != rewrite_c4e.per_edge_[5][e]) {
            std::cout << "\n  MISMATCH per_edge_[5] at " << e << "\n";
            ok = false;
        }
    }
    t = Clock::now();
    Pattern66OrbitCounter(gout, gin, escape);
    const double escape_p66_ms = msSince(t);
    t = Clock::now();
    pattern66Orbits(gout, rewrite);
    const double rewrite_p66_ms = msSince(t);
    for (EdgeIdx e = 0; e < gout->nEdges && ok; ++e) {
        if (escape.per_edge_[11][e] != rewrite.per_edge_[11][e]) {
            std::cout << "\n  MISMATCH per_edge_[11] at " << e << "\n";
            ok = false;
        }
    }

    std::cout << "  triangles=" << all_tris.total / 3 << "\n"
              << "  ms escape/rewrite: 4-clique " << escape_k4_ms << "/" << rewrite_k4_ms
              << "  4-cycle(node) " << escape_c4_ms << "/" << rewrite_c4_ms
              << "  4-cycle(edges) " << escape_c4_ms << "/" << rewrite_c4e_ms
              << "  pattern66 " << escape_p66_ms << "/" << rewrite_p66_ms << "\n"
              << "  ms escape/rewrite: triangle lists " << escape_tris_ms << "/" << rewrite_tris_ms
              << "  5-cycle " << escape_cycle_ms << "/" << rewrite_cycle_ms
              << "  chordal " << escape_chordal_ms << "/" << rewrite_chordal_ms
              << "  almost " << escape_almost_ms << "/" << rewrite_almost_ms << "\n";

    for (int o : {8, 14, 34, 62, 63, 64, 66, 68, 69, 70, 71}) {
        const std::string what = "orbit " + std::to_string(o);
        ok &= sameCounts(what.c_str(), escape.per_vertex_[o], rewrite.per_vertex_[o], n);
    }

    // The chordal parts skipped by the flags must not change the parts that run.
    const bool flags[3][3] = {{true, false, false}, {false, true, false}, {false, false, true}};
    const std::vector<int> part_orbits[3] = {{62}, {63, 64}, {68, 69}};
    for (int f = 0; f < 3; ++f) {
        OrbitInfo partial(gout, 73, 12);
        t = Clock::now();
        chordalWedgeOrbits(&cg, et, partial, flags[f][0], flags[f][1], flags[f][2]);
        std::cout << "  chordal part " << (f == 0 ? "62" : f == 1 ? "63/64" : "68/69") << " alone: " << msSince(t) << " ms\n";
        for (int o : part_orbits[f]) {
            const std::string what = "orbit " + std::to_string(o) + " (partial run " + std::to_string(f) + ")";
            ok &= sameCounts(what.c_str(), escape.per_vertex_[o], partial.per_vertex_[o], n);
        }
    }

    delete[] all_tris.triangles;
    delete[] all_tris.trioffsets;
    delete[] mapping;
    delete[] inverse;
    std::cout << (ok ? "  EQUAL\n" : "  DIFFERENT\n");
    return ok;
}
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "usage: orbit_equivalence <graph file> <max vertices> <seeds>\n";
        return 2;
    }
    AppGraph g;
    if (!g.readFromFile(argv[1])) return 2;
    const size_t max_vertices = std::strtoull(argv[2], nullptr, 10);
    const int seeds = std::atoi(argv[3]);
    std::cout << argv[1] << "\n";

    bool ok = true;
    if (max_vertices == 0) {
        std::cout << "whole graph\n";
        ok &= testGraph(g);
    } else {
        int hub = g.original_node_ids[0];
        for (int id : g.original_node_ids) if (g.degree(id) > g.degree(hub)) hub = id;
        std::mt19937 rng(12345);
        std::uniform_int_distribution<size_t> pick(0, g.original_node_ids.size() - 1);
        for (int s = 0; s < seeds; ++s) {
            const int seed_id = s == 0 ? hub : g.original_node_ids[pick(rng)];
            std::cout << "BFS ball of <= " << max_vertices << " from vertex " << seed_id << "\n";
            ok &= testGraph(bfsBall(g, seed_id, max_vertices));
        }
    }
    std::cout << (ok ? "ALL EQUAL\n" : "SOME DIFFERENT\n");
    return ok ? 0 : 1;
}
