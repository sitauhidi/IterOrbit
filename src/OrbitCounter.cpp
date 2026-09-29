#include "OrbitCounter.h"
#include <chrono>
#include <stdexcept>
#include <iostream>

#include "Escape/OrbitStructure.h"
#include "Escape/GetAllOrbitCounts.h"
#include "Escape/FiveVertexOrbit.h"
#include "Escape/FourVertexOrbit.h"

using namespace Escape;

EVOKEOrbitCounter::EVOKEOrbitCounter(const AppGraph& g, int s) : graph(g), graphlet_size(s) {}

namespace {
using Clock = std::chrono::high_resolution_clock;
double msSince(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}
}

OrbitCounts EVOKEOrbitCounter::count() {
    last_timings = OrbitCountTimings();
    CGraph* g_raw = graph.c_graph;
    if (!g_raw) return {};

    auto t = Clock::now();
    VertexIdx* mapping;
    VertexIdx* inverse;
    CGraph g_relabel = g_raw->renameByDegreeOrder(mapping, inverse);
    g_relabel.sortById();

    CDAG dag = degreeOrdered(&g_relabel);
    dag.outlist.sortById();
    dag.inlist.sortById();
    last_timings.relabel_dag_ms = msSince(t);

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
    get_all_three_orbits(&g_relabel, &dag, orbit_counts);
    last_timings.three_ms = msSince(t);

    if (graphlet_size == 5) {
        t = Clock::now();
        get_all_four_orbits(&g_relabel, &dag, orbit_counts);
        get_all_five_orbits(&g_relabel, &dag, orbit_counts, 1);
        last_timings.four_other_ms = msSince(t);
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

        t = Clock::now();
        otherFourPatternOrbitCounterNodeOnly(&(dag.outlist), &(dag.inlist), orbit_counts);
        last_timings.four_other_ms = msSince(t);

        t = Clock::now();
        four_cycle_orbit_counter_node_only(&(dag.outlist), &(dag.inlist), orbit_counts);
        last_timings.four_cycle_ms = msSince(t);

        t = Clock::now();
        fourCliqueOrbitCounterNodeOnly(&(dag.outlist), orbit_counts);
        last_timings.four_clique_ms = msSince(t);
    }

    t = Clock::now();
    OrbitCounts result(g_relabel.nVertices);
    for (VertexIdx i = 0; i < g_relabel.nVertices; ++i) {
        std::vector<long long>& counts = result[i];
        counts.reserve(num_node_orbits);
        for (int orbit = 0; orbit < num_node_orbits; ++orbit) {
            counts.push_back(orbit_counts.per_vertex_[orbit][mapping[i]]);
        }
    }

    delete[] mapping;
    delete[] inverse;
    last_timings.copy_ms = msSince(t);

    return result;
}