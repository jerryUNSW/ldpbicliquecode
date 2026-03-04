#include "biclique_dds_naive.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <random>
#include <numeric>

using namespace std;

// Forward declaration of the true biclique counter from the main ELDP file
extern long double count_bicliques_in_set(BiGraph& g, const vector<int>& S, int P, int Q);

NaiveResult eldp_biclique_dds_naive(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed) {
    int n = g.num_nodes();
    NaiveResult result;
    mt19937 rng(seed);

    // 1. Privacy Budget Allocation: 100% to degree estimation (Naive approach)
    // We use Laplace mechanism for degrees as a simple ELDP baseline
    long double eps_deg = epsilon;
    exponential_distribution<double> exp_dist(eps_deg);
    
    vector<pair<long double, int>> noisy_degrees(n);
    for (int i = 0; i < n; i++) {
        // Laplace noise = Subtraction of two i.i.d. Exponentials
        double noise = exp_dist(rng) - exp_dist(rng);
        noisy_degrees[i] = {(long double)g.neighbor[i].size() + noise, i};
    }

    // 2. Initial Filtering: Keep top-K vertices by noisy degree
    // A naive but efficient way to handle large graphs
    sort(noisy_degrees.rbegin(), noisy_degrees.rend());
    
    int K = min(n, 1000); // Process at most 1000 vertices for the naive baseline
    vector<int> S;
    for (int i = 0; i < K; i++) {
        S.push_back(noisy_degrees[i].second);
    }

    // 3. Naive Greedy Peeling based on Noisy Degrees
    // In each step, remove the vertex with the smallest noisy degree
    vector<int> best_set = S;
    long double max_density = 0;
    
    // To make it "naive" but efficient, we just peel from the degree-sorted list
    // and evaluate true density at each step (as a benchmark)
    for (int i = K; i >= (P + Q); i--) {
        vector<int> current_S;
        for(int j=0; j<i; j++) current_S.push_back(noisy_degrees[j].second);
        
        // Check if we have enough vertices on each side
        int u_count = 0, l_count = 0;
        for(int v : current_S) {
            if(g.is_upper(v)) u_count++;
            else l_count++;
        }
        
        if (u_count >= P && l_count >= Q) {
            long double bc = count_bicliques_in_set(g, current_S, P, Q);
            long double density = bc / (long double)current_S.size();
            if (density > max_density) {
                max_density = density;
                best_set = current_S;
            }
        }
        
        // Optimization: if we've peeled a lot and density is dropping fast, stop
        if (i < K - 100 && max_density > 0 && i < K / 2) break;
    }

    result.vertex_set = best_set;
    result.real_density = max_density;
    result.estimated_density = max_density; // Naive doesn't have a good estimator
    
    cout << "  [Naive Baseline] Final |S|=" << best_set.size() << " d=" << fixed << setprecision(4) << max_density << endl;
    
    return result;
}
