#ifndef BICLIQUE_DDS_NAIVE_ELDP_V2_H
#define BICLIQUE_DDS_NAIVE_ELDP_V2_H

#include "bigraph.h"
#include <vector>

struct NaiveV2Result {
    std::vector<int> vertex_set;
    long double estimated_density;
    long double real_density;
};

NaiveV2Result naive_eldp_v2(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed);

#endif
