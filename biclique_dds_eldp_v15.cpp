/*
 * biclique_dds_eldp_v15.cpp — v15: Direct v11 wrapper (stable baseline)
 *
 * Strategy: Simply call v11's eldp_biclique_dds (butterfly peeling + refinement)
 * This is the stable baseline that works well on large graphs.
 */

#include "biclique_dds_eldp_v15.h"
#include "biclique_dds_eldp.h"

#include <iostream>
#include <iomanip>

using namespace std;

// v11 baseline entry
extern DDSResult eldp_biclique_dds(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed);

DDSResultV15 eldp_v15(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed) {
    // Directly use v11's stable implementation
    DDSResult v11_result = eldp_biclique_dds(g, P, Q, epsilon, seed);
    
    DDSResultV15 out;
    out.vertex_set = v11_result.vertex_set;
    out.estimated_density = v11_result.estimated_density;
    out.real_density = v11_result.real_density;
    
    cout << "  [v15] wrapper around v11 |S|=" << out.vertex_set.size()
         << " est_d=" << fixed << setprecision(6) << out.estimated_density
         << " real_d=" << setprecision(6) << out.real_density << endl;
    
    return out;
}
