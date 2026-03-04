#ifndef BICLIQUE_DDS_PEELING_EXPAND_SHRINK_ELDP_H
#define BICLIQUE_DDS_PEELING_EXPAND_SHRINK_ELDP_H

#include "bigraph.h"
#include <vector>

struct DDSResult;

// Peeling-Expand-Shrink-ELDP: peeling + greedy add + greedy shrink on noisy graph G'.
// All operations use G' only. ELDP compliant.
DDSResult peeling_expand_shrink_eldp(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed);

#endif
