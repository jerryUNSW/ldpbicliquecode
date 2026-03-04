#ifndef BICLIQUE_DDS_ELDP_V12_H
#define BICLIQUE_DDS_ELDP_V12_H

#include "bigraph.h"
#include <vector>

struct DDSResultV12 {
    std::vector<int> vertex_set;
    long double estimated_density;
    long double real_density;
};

DDSResultV12 eldp_v12(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed);

#endif
