#include "biclique_dds_eldp.h"
#include "biclique.h"
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <random>
#include <iomanip>
#include <algorithm>

using namespace std;

// Globals required by biclique.cpp
long double Eps = 0, Eps0 = 0, Eps1 = 0, Eps2 = 0, p = 0;
long double m3__ = 0, m2__ = 0, m1__ = 0, m0__ = 0;
long double communication_cost = 0;
long double alpha0 = 0.05, alpha1 = 0.95, alpha2 = 0.0;
vector<int> priv_deg;
vector<long double> naive_estis;
int priv_dmax_1 = 0, priv_dmax_2 = 0;
vector<vector<int>> up_options, lo_options;
int iteration = 0;
int K___, P___;
long double real = 0.0, real_ld = 0.0;
long double avg_estimated_variance = 0;
long double RR_time = 0, server_side_time = 0, naive_server_side = 0,
            local_count_time = 0, deg_esti_time = 0;
int num_rounds = 1;
vector<long double> estis;

// These are extern in biclique.h
bool use_probability_filtering = false;
int post_processing_mode = 0;
long double winsorize_percentile = 95.0;
long double truncation_bound = 0.0;

// Helper to calculate mean and variance
void calculate_stats(const vector<long double>& data, double& mean, double& variance) {
    if (data.empty()) {
        mean = 0.0;
        variance = 0.0;
        return;
    }
    double sum = accumulate(data.begin(), data.end(), 0.0);
    mean = sum / data.size();
    double sq_sum = 0.0;
    for (double x : data) {
        sq_sum += (x - mean) * (x - mean);
    }
    variance = sq_sum / data.size();
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <dataset_path> <P> <Q>" << endl;
        return 1;
    }

    string dataset_path = argv[1];
    int P = atoi(argv[2]);
    int Q = atoi(argv[3]);

    cout << "Loading graph: " << dataset_path << endl;
    BiGraph g(dataset_path);
    cout << "Graph loaded: n1=" << g.num_v1 << ", n2=" << g.num_v2 << ", m=" << g.num_edges << endl;

    // Select a subgraph S
    // Strategy: Top 30 upper + Top 30 lower vertices by degree
    int n = g.num_nodes();
    cout << "Sorting degrees for " << n << " vertices..." << endl;
    vector<pair<int, int>> deg_sorted;
    for (int i = 0; i < n; i++) {
        deg_sorted.push_back({(int)g.neighbor[i].size(), i});
    }
    sort(deg_sorted.rbegin(), deg_sorted.rend());
    cout << "Degrees sorted." << endl;

    vector<int> S;
    vector<int> S_upper, S_lower;
    int count_u = 0, count_l = 0;
    
    for (auto& p : deg_sorted) {
        int v = p.second;
        if (g.is_upper(v)) {
            if (count_u < 30) {
                S.push_back(v);
                S_upper.push_back(v);
                count_u++;
            }
        } else {
            if (count_l < 30) {
                S.push_back(v);
                S_lower.push_back(v);
                count_l++;
            }
        }
        if (count_u >= 30 && count_l >= 30) break;
    }

    cout << "Selected subgraph S: |S|=" << S.size() 
         << " (|U|=" << S_upper.size() << ", |L|=" << S_lower.size() << ")" << endl;

    // Calculate TRUE count
    long double true_count = count_bicliques_in_set(g, S, P, Q);
    cout << "True (p,q)-biclique count: " << true_count << endl;

    vector<double> epsilons = {0.5, 1.0, 2.0, 3.0, 5.0};
    int rounds = 100;

    cout << "\nStarting Unbiasedness Verification (" << rounds << " rounds per epsilon)..." << endl;
    cout << "Epsilon | Mean Est. | Bias | Variance | Std Error | Rel. Bias %" << endl;
    cout << "--------|-----------|------|----------|-----------|------------" << endl;

    random_device rd;
    mt19937 rng(rd());

    for (double eps : epsilons) {
        long double eps0 = eps * 0.05; // As per algorithm
        long double eps1 = eps - eps0;
        long double p_flip = 1.0 / (exp(eps1) + 1.0);

        vector<long double> estimates;
        
        // Noisy graph structure (shared)
        BiGraph g_noisy(g);
        // Ensure edge_vector is initialized
        g_noisy.edge_vector.resize(n);

        for (int r = 0; r < rounds; r++) {
            // Clear noisy edges for new round
            for (int i = 0; i < n; i++) g_noisy.edge_vector[i].clear();

            long double est = estimate_biclique_count_dp(
                S_upper, S_lower, P, Q, g_noisy, g, p_flip, eps1, rng);
            estimates.push_back(est);
        }

        double mean, variance;
        calculate_stats(estimates, mean, variance);
        double bias = abs(mean - true_count);
        double std_err = sqrt(variance) / sqrt(rounds);
        double rel_bias = (true_count > 0) ? (bias / true_count) * 100.0 : 0.0;

        cout << fixed << setprecision(1) << eps << "     | "
             << setprecision(2) << mean << " | "
             << setprecision(2) << bias << " | "
             << scientific << setprecision(2) << variance << " | "
             << fixed << setprecision(2) << std_err << " | "
             << setprecision(2) << rel_bias << "%" << endl;
    }

    return 0;
}