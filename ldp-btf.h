#pragma once
#ifndef ABCORE_H
#define ABCORE_H
#include <sstream>

#include "bigraph.h"

void private_estimate_of_degrees(BiGraph& g);

long double BFC_EVP(BiGraph& g);

bool satisfy_bound(long double upper, long double lower, int u, int x, int v,
                   int w, BiGraph& g);

double my_genrand_real2();

void my_init_genrand(unsigned long seed);

// long double BFC_EVP_sample(BiGraph& g, long double sampling_prob);

void compute_m3_m2(long double& m4, long double& m3, long double& m2,
                   long double& m1, long double& m0, BiGraph& g2);

long double one_round_btf(BiGraph& g, unsigned long seed);

long double two_round_btf(BiGraph& g, unsigned long seed);

// naive noisy motif counts
void get_noisy_naive(BiGraph& g, BiGraph& g2, long double& local_btfs,
                     long double& local_cate, long double& res);

// process butterflies in batch
void BFC_EVP_noisy(BiGraph& g, BiGraph& g2, long double& BTF, long double& cate,
                   long double& res);

void approx_abcore(BiGraph& g, long double Eps0);

void construct_noisy_graph(BiGraph& g, BiGraph& g2, unsigned long seed);

void construct_noisy_graph_2(int upper, int lower, BiGraph& g, BiGraph& g2,
                             unsigned long seed);

void construct_noisy_graph_3(int upper, int lower, BiGraph& g, BiGraph& g2,
                             unsigned long seed);

long double get_wedges(BiGraph& g);

long double get_cate(BiGraph& g);

int findPercentile(std::vector<int>& data);

void test_random_number(unsigned long seed);

long double get_laplace(long double parameter);

void compute_m3_m2_2(long double& m4, long double& m3, long double& m2,
                     long double& m1, long double& m0, BiGraph& g2);

// Custom hash function for std::vector<int>
struct VectorHasher {
    std::size_t operator()(const std::vector<int>& v) const {
        std::size_t seed = v.size();
        for (const auto& i : v) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};
#endif