#ifndef BICLIQUE_DDS_ELDP_V13_H
#define BICLIQUE_DDS_ELDP_V13_H

#include "bigraph.h"
#include <vector>

struct DDSResultV13 {
    std::vector<int> vertex_set;
    long double estimated_density;
    long double real_density;
};

DDSResultV13 eldp_v13(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed);

#endif
