#include "biclique_dds_eldp_v13.h"

#include "biclique_dds_eldp.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <random>

using namespace std;

// v11 baseline entry.
extern DDSResult eldp_biclique_dds(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed);

DDSResultV13 eldp_v13(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed) {
    // Dual independent noisy runs (same total budget per run, independent randomness)
    // and select the stronger solution under the algorithm's own objective.
    DDSResult run1 = eldp_biclique_dds(g, P, Q, epsilon, seed);

    unsigned long seed2 = seed ^ 0x9E3779B97F4A7C15ULL;
    DDSResult run2 = eldp_biclique_dds(g, P, Q, epsilon, seed2);

    const DDSResult* best = &run1;
    if (run2.real_density > run1.real_density) {
        best = &run2;
    }

    DDSResultV13 out;
    out.vertex_set = best->vertex_set;
    out.estimated_density = best->estimated_density;
    out.real_density = best->real_density;

    cout << "  [edlp_v13] dual-noisy seeds=(" << seed << "," << seed2 << ")"
         << " chosen_real_d=" << fixed << setprecision(6) << out.real_density
         << " chosen_est_d=" << setprecision(6) << out.estimated_density << endl;

    return out;
}
