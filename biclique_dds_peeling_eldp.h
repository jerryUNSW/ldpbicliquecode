#ifndef BICLIQUE_DDS_PEELING_ELDP_H
#define BICLIQUE_DDS_PEELING_ELDP_H

#include "bigraph.h"
#include <vector>
#include <random>

struct DDSResult;

// Peeling-ELDP: Port of pqBicliqeDensest peeling to work on noisy graph G'.
// All operations use G' only (ELDP compliant).
DDSResult peeling_eldp(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed);

// Core peeling on given g_noisy (for reuse in peeling-expand-shrink).
std::vector<int> peeling_eldp_core(BiGraph& g, BiGraph& g_noisy, int P, int Q,
    long double p_flip, long double eps1, std::mt19937& rng);

#endif
