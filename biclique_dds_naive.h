#ifndef BICLIQUE_DDS_NAIVE_H
#define BICLIQUE_DDS_NAIVE_H

#include "bigraph.h"
#include "utility.h"
#include <vector>
#include <string>

struct NaiveResult {
    std::vector<int> vertex_set;
    long double estimated_density;
    long double real_density;
};

// Naive ELDP Baseline: Degree-based Filtering + Greedy Peeling on Noisy Degrees
NaiveResult eldp_biclique_dds_naive(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed);

#endif
