#ifndef BICLIQUE_DDS_ELDP_V15_H
#define BICLIQUE_DDS_ELDP_V15_H

#include "bigraph.h"
#include <vector>

struct DDSResultV15 {
    std::vector<int> vertex_set;
    long double estimated_density;
    long double real_density;
};

DDSResultV15 eldp_v15(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed);

#endif
