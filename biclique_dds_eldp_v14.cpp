#include "biclique_dds_eldp_v14.h"

#include "biclique_dds_eldp.h"
#include "biclique_dds_eldp_v15.h"

#include <iomanip>
#include <iostream>

using namespace std;

DDSResultV14 eldp_v14(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed) {
    // Round-14: v15 core (v11 init + v12 techniques) with multi-start
    unsigned long s1 = seed;
    unsigned long s2 = seed ^ 0x9E3779B97F4A7C15ULL;
    unsigned long s3 = seed ^ 0xC2B2AE3D27D4EB4FULL;

    DDSResultV15 r1 = eldp_v15(g, P, Q, epsilon, s1);
    DDSResultV15 r2 = eldp_v15(g, P, Q, epsilon, s2);
    DDSResultV15 r3 = eldp_v15(g, P, Q, epsilon, s3);

    const DDSResultV15* best = &r1;
    if (r2.estimated_density > best->estimated_density) best = &r2;
    if (r3.estimated_density > best->estimated_density) best = &r3;

    DDSResultV14 out;
    out.vertex_set = best->vertex_set;
    out.estimated_density = best->estimated_density;

    if (!out.vertex_set.empty()) {
        long double bc = count_bicliques_in_set(g, out.vertex_set, P, Q);
        out.real_density = bc / (long double)out.vertex_set.size();
    } else {
        out.real_density = 0.0L;
    }

    cout << "  [edlp_v14] v15-core multi-start seeds=(" << s1 << "," << s2 << "," << s3 << ")"
         << " pick_by_est=" << fixed << setprecision(6) << out.estimated_density
         << " real_d=" << setprecision(6) << out.real_density << endl;

    return out;
}
