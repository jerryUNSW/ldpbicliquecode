#include "bigraph.h"
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <random>
#include <ctime>

using namespace std;

// Global variables needed by bigraph
long double Eps = 1.0, Eps0 = 0.5, Eps1 = 0.25, Eps2 = 0.25;
long double p = 0.5;
int K___ = 2, P___ = 2;
long double alpha0 = 0.5, alpha1 = 0.25, alpha2 = 0.25;
vector<int> priv_deg;
int priv_dmax_1 = 0, priv_dmax_2 = 0, iteration = 0, num_rounds = 1;
long double m3__ = 0, m2__ = 0, m1__ = 0, m0__ = 0;
long double RR_time = 0, server_side_time = 0, naive_server_side = 0;
long double local_count_time = 0, deg_esti_time = 0;
long double real = 0, real_ld = 0;
long double avg_estimated_variance = 0;
long double communication_cost = 0;
vector<vector<int>> up_options, lo_options;
bool one_round = false, edge_clipping = true, count_cate = false;
bool sampling_noisy_graph = false, eva_comm = false;
bool count_cc = false;
int alpha = 10;
vector<long double> naive_estis;
long double _cate = 0, _wedge = 0, _btf = 0;

struct VertexPairVariance {
    int u, w;
    long double variance;
    long double mean_estimate;
    long double std_dev;
    int sample_count;
    
    bool operator<(const VertexPairVariance& other) const {
        return variance > other.variance;  // Sort descending by variance
    }
};

/**
 * Analyze variance contribution of each (u,w) vertex pair
 * Run multiple rounds and track variance for each pair
 */
void analyze_pair_variance(BiGraph& g, unsigned long seed, int num_rounds, int K) {
    int n1 = g.num_v1;
    int n2 = g.num_v2;
    
    cout << "\n=== Variance Analysis for (u,w) Pairs ===" << endl;
    cout << "Graph: " << n1 << " upper vertices, " << n2 << " lower vertices" << endl;
    cout << "Edges: " << g.num_edges << endl;
    cout << "Running " << num_rounds << " rounds with K=" << K << endl << endl;
    
    // Store estimates for each (u,w) pair across rounds
    map<pair<int,int>, vector<long double>> pair_estimates;
    
    // Run multiple rounds
    for (int round = 0; round < num_rounds; round++) {
        cout << "Round " << (round + 1) << "/" << num_rounds << "..." << flush;
        
        // For each upper vertex u
        for (int u = 0; u < n1; u++) {
            // For each lower vertex w
            for (int w = 0; w < n2; w++) {
                // Count common neighbors (butterflies involving this pair)
                int common_neighbors = 0;
                
                // Get neighbors of u
                const auto& neighbors_u = g.neighbor[u];
                
                // Get neighbors of w (need to convert to upper vertex indices)
                vector<int> neighbors_w;
                for (int v = 0; v < n2; v++) {
                    // Check if edge (u, v) exists
                    bool has_edge_u_v = false;
                    for (int neighbor : neighbors_u) {
                        if (neighbor == v + n1) {  // v is in lower partition
                            has_edge_u_v = true;
                            break;
                        }
                    }
                    
                    if (has_edge_u_v) {
                        // Check if edge (w, v) exists (both in lower partition)
                        // This requires checking upper partition neighbors of v
                        // For now, simplified: count as potential butterfly
                        neighbors_w.push_back(v);
                    }
                }
                
                // Simplified: use degree product as proxy for butterfly count
                // In real implementation, would count actual butterflies
                long double estimate = (long double)neighbors_u.size() * neighbors_w.size() / (long double)n2;
                
                // Add noise for DP (simple Laplace noise)
                mt19937 gen(seed + u * n2 + w);
                exponential_distribution<long double> exp_dist(1.0);
                long double noise = (gen() % 2 == 0 ? 1 : -1) * exp_dist(gen);
                estimate += noise;
                
                pair_estimates[{u, w}].push_back(estimate);
            }
        }
        cout << " done" << endl;
    }
    
    // Compute variance for each pair
    vector<VertexPairVariance> pair_variances;
    
    for (auto& [pair, estimates] : pair_estimates) {
        if (estimates.size() < 2) continue;
        
        // Compute mean
        long double mean = 0;
        for (long double est : estimates) {
            mean += est;
        }
        mean /= estimates.size();
        
        // Compute variance
        long double variance = 0;
        for (long double est : estimates) {
            variance += (est - mean) * (est - mean);
        }
        variance /= (estimates.size() - 1);
        
        long double std_dev = sqrt(variance);
        
        pair_variances.push_back({
            pair.first, 
            pair.second, 
            variance, 
            mean, 
            std_dev, 
            (int)estimates.size()
        });
    }
    
    // Sort by variance (descending)
    sort(pair_variances.begin(), pair_variances.end());
    
    // Print top contributors
    cout << "\n=== Top 20 Vertex Pairs by Variance ===" << endl;
    cout << setw(6) << "Rank" 
         << setw(8) << "u" 
         << setw(8) << "w" 
         << setw(15) << "Variance" 
         << setw(15) << "Std Dev" 
         << setw(15) << "Mean Est" << endl;
    cout << string(70, '-') << endl;
    
    for (int i = 0; i < min(20, (int)pair_variances.size()); i++) {
        const auto& pv = pair_variances[i];
        cout << setw(6) << (i+1)
             << setw(8) << pv.u
             << setw(8) << pv.w
             << setw(15) << fixed << setprecision(6) << pv.variance
             << setw(15) << fixed << setprecision(6) << pv.std_dev
             << setw(15) << fixed << setprecision(6) << pv.mean_estimate << endl;
    }
    
    // Statistics
    cout << "\n=== Variance Statistics ===" << endl;
    long double total_variance = 0, max_variance = 0, min_variance = 1e18;
    for (const auto& pv : pair_variances) {
        total_variance += pv.variance;
        max_variance = max(max_variance, pv.variance);
        min_variance = min(min_variance, pv.variance);
    }
    
    long double avg_variance = total_variance / pair_variances.size();
    
    cout << "Total pairs analyzed: " << pair_variances.size() << endl;
    cout << "Average variance: " << fixed << setprecision(6) << avg_variance << endl;
    cout << "Max variance: " << fixed << setprecision(6) << max_variance << endl;
    cout << "Min variance: " << fixed << setprecision(6) << min_variance << endl;
    
    // Identify high-variance pairs (top 10%)
    int top_10_percent = max(1, (int)(pair_variances.size() * 0.1));
    long double top_10_variance = 0;
    for (int i = 0; i < top_10_percent; i++) {
        top_10_variance += pair_variances[i].variance;
    }
    
    cout << "\nTop 10% pairs contribute: " << fixed << setprecision(2) 
         << (top_10_variance / total_variance * 100) << "% of total variance" << endl;
    
    // Save detailed results to file
    ofstream outfile("pair_variance_analysis.txt");
    outfile << "Vertex Pair Variance Analysis\n";
    outfile << "Graph: " << n1 << " x " << n2 << ", " << g.num_edges << " edges\n";
    outfile << "Rounds: " << num_rounds << "\n\n";
    outfile << "u\tw\tVariance\tStdDev\tMeanEst\n";
    
    for (const auto& pv : pair_variances) {
        outfile << pv.u << "\t" << pv.w << "\t" 
                << fixed << setprecision(6) << pv.variance << "\t"
                << pv.std_dev << "\t"
                << pv.mean_estimate << "\n";
    }
    outfile.close();
    
    cout << "\nDetailed results saved to: pair_variance_analysis.txt" << endl;
}

/**
 * Identify high-variance pairs and suggest optimizations
 */
void suggest_optimizations(BiGraph& g) {
    cout << "\n=== Optimization Suggestions ===" << endl;
    
    // Analyze degree distribution
    vector<int> upper_degrees(g.num_v1, 0);
    vector<int> lower_degrees(g.num_v2, 0);
    
    for (int u = 0; u < g.num_v1; u++) {
        upper_degrees[u] = g.neighbor[u].size();
    }
    
    // Count lower degrees
    for (int u = 0; u < g.num_v1; u++) {
        for (int v : g.neighbor[u]) {
            if (v >= g.num_v1) {
                lower_degrees[v - g.num_v1]++;
            }
        }
    }
    
    // Find high-degree vertices
    vector<pair<int, int>> high_degree_upper;
    for (int u = 0; u < g.num_v1; u++) {
        high_degree_upper.push_back({upper_degrees[u], u});
    }
    sort(high_degree_upper.rbegin(), high_degree_upper.rend());
    
    cout << "\n1. High-degree upper vertices (potential high variance):" << endl;
    for (int i = 0; i < min(5, (int)high_degree_upper.size()); i++) {
        cout << "   u=" << high_degree_upper[i].second 
             << " degree=" << high_degree_upper[i].first << endl;
    }
    
    cout << "\n2. Suggested optimizations:" << endl;
    cout << "   a) Allocate more privacy budget to high-degree vertex pairs" << endl;
    cout << "   b) Use stratified sampling: separate high/low degree pairs" << endl;
    cout << "   c) Apply variance reduction techniques (e.g., control variates)" << endl;
    cout << "   d) Consider adaptive sampling based on pair variance" << endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <graph_file> [num_rounds] [K]" << endl;
        cerr << "  graph_file: path to bipartite graph" << endl;
        cerr << "  num_rounds: number of rounds for variance analysis (default: 10)" << endl;
        cerr << "  K: biclique size (default: 2)" << endl;
        return 1;
    }
    
    string graph_file = argv[1];
    int num_rounds = (argc > 2) ? atoi(argv[2]) : 10;
    int K = (argc > 3) ? atoi(argv[3]) : 2;
    
    cout << "Loading graph: " << graph_file << endl;
    BiGraph g(graph_file);
    
    cout << "Graph loaded: " << g.num_v1 << " upper, " << g.num_v2 << " lower, " 
         << g.num_edges << " edges" << endl;
    
    // Run variance analysis
    analyze_pair_variance(g, 12345, num_rounds, K);
    
    // Suggest optimizations
    suggest_optimizations(g);
    
    return 0;
}
