#pragma once
#ifndef ABCORE_H
#define ABCORE_H
#include <sstream>

#include "bigraph.h"


long double weighted_pair_sampling_non_DP(BiGraph& g, unsigned long seed); 
    
double locally_compute_f_given_q_and_x_ad_hoc(int q, int x, BiGraph& g, BiGraph& g2);

void randomized_response_single_bit(int u, int v, BiGraph& g, BiGraph& g2); 

// weighted sampling based approach:
long double weighted_pair_sampling(BiGraph& g, unsigned long seed); 

// (2, K)
long double wedge_based_two_round_2_K_biclique(BiGraph& g, unsigned long seed) ; 

// (3, K)
long double wedge_based_two_round_3_K_biclique(BiGraph& g, unsigned long seed) ;

// general shape: 
long double wedge_based_two_round_general_biclique(BiGraph& g, 
    unsigned long seed, int P___, int K___ );

void check_exact_result_in_DB(int P___, int K___, string dataset);

// improved new approach
long double wedge_based_btf_avg(BiGraph& g, unsigned long seed);

// VP based 
long double VP_wedge_based_two_round_btf(BiGraph& g, unsigned long seed);

double locally_compute_f_given_q_and_x_vp(int q, int x, BiGraph& g, BiGraph& g2, int& res__);

double locally_compute_f_given_q_and_x_vp_2(int q, int x, BiGraph& g, BiGraph& g2, int& res__);

double locally_compute_f_given_q_and_x(int q, int x, BiGraph& g, BiGraph& g2); 

double locally_compute_f_given_q_and_x_two_graphs(int q, int x, BiGraph& g, BiGraph& g2, BiGraph& g3);

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

void construct_noisy_graph(BiGraph& g, BiGraph& g2, unsigned long seed);

void construct_noisy_graph_2(BiGraph& g, BiGraph& g2, unsigned long seed); 

long double get_wedges(BiGraph& g);

long double get_cate(BiGraph& g);

int findPercentile(std::vector<int>& data);

void test_random_number(unsigned long seed);

long double get_laplace(long double parameter);

void compute_m3_m2_2(long double& m4, long double& m3, long double& m2,
                     long double& m1, long double& m0, BiGraph& g2);


// one-round biclique counting
long double one_round_biclique(BiGraph& g, unsigned long seed, int p__, int q__, int _switch); 


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