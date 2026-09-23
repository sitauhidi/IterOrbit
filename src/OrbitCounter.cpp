#include "OrbitCounter.h"
#include <stdexcept>
#include <iostream>

#include "Escape/OrbitStructure.h"
#include "Escape/GetAllOrbitCounts.h"
#include "Escape/FiveVertexOrbit.h"
#include "Escape/FourVertexOrbit.h"

using namespace Escape;

EVOKEOrbitCounter::EVOKEOrbitCounter(const AppGraph& g, int s) : graph(g), graphlet_size(s) {}

OrbitCounts EVOKEOrbitCounter::count() {
    CGraph* g_raw = graph.c_graph;
    if (!g_raw) return {};

    VertexIdx* mapping;
    VertexIdx* inverse;
    CGraph g_relabel = g_raw->renameByDegreeOrder(mapping, inverse);
    g_relabel.sortById();

    CDAG dag = degreeOrdered(&g_relabel);
    dag.outlist.sortById();
    dag.inlist.sortById();
    
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
    
    get_all_three_orbits(&g_relabel, &dag, orbit_counts);

    if (graphlet_size >= 4) {
        if (graphlet_size == 5) {
             get_all_four_orbits(&g_relabel, &dag, orbit_counts);
             get_all_five_orbits(&g_relabel, &dag, orbit_counts, 1);
        } else {
            get_all_four_orbits_node_only(&g_relabel, &dag, orbit_counts);
        }
    }
    
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

    return result;
}