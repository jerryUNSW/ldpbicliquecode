/*
 * biclique_dds_eldp.h
 * (p,q)-Biclique DDS under Edge LDP — v2 Motif-Counting Approach.
 *
 * Uses direct motif counting + debiasing (like triangle DDS)
 * instead of moment-corrected estimators.
 */
#pragma once
#ifndef BICLIQUE_DDS_ELDP_H
#define BICLIQUE_DDS_ELDP_H

#include "bigraph.h"
#include <vector>
#include <string>
#include <random>

struct DDSResult {
    std::vector<int> vertex_set;
    long double estimated_density;
    long double real_density;
};

// Main ELDP DDS algorithm (greedy expand + shrink with motif-counting debiasing)
DDSResult eldp_biclique_dds(BiGraph& g, int P, int Q,
    long double epsilon, unsigned long seed);

// Count true (p,q)-bicliques in induced subgraph on vertex set S
long double count_bicliques_in_set(BiGraph& g, const std::vector<int>& S, int P, int Q);

// Lazy noisy edge query (ELDP compliant, caches in g_noisy)
bool query_noisy_edge(BiGraph& g_orig, BiGraph& g_noisy, int u, int v,
    long double p_flip, std::mt19937& rng);

// Unbiased estimator for (p,q)-bicliques using DP
long double estimate_biclique_count_dp(
    const std::vector<int>& S_upper, const std::vector<int>& S_lower,
    int P, int Q,
    BiGraph& g_noisy, BiGraph& g_orig,
    long double p_flip, long double eps1, std::mt19937& rng);

// Incremental gain: bicliques formed by adding v_new to S
long double estimate_biclique_gain(
    const std::vector<int>& S_upper, const std::vector<int>& S_lower,
    int v_new, bool is_upper_new,
    int P, int Q,
    BiGraph& g_noisy, BiGraph& g_orig,
    long double p_flip, long double eps1, std::mt19937& rng);

// Load exact optimal density from SQLite database
bool load_exact_density(const std::string& dataset, int P, int Q, long double& density);

// Set trace file for recording greedy process (for variance analysis)
void set_trace_file(std::ofstream* fp);

#endif
