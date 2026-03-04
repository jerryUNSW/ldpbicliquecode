#ifndef BICLIQUE_DDS_ELDP_V14_H
#define BICLIQUE_DDS_ELDP_V14_H

#include "bigraph.h"
#include <vector>

struct DDSResultV14 {
    std::vector<int> vertex_set;
    long double estimated_density;
    long double real_density;
};

DDSResultV14 eldp_v14(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed);

#endif
