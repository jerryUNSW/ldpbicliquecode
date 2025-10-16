/*
 * test_batch_with_ground_truth.cpp
 * 
 * Test the batch function with proper ground truth comparison
 * Runs 10 rounds and calculates relative errors for Q = [4,5,6,7,8,9,10]
 */

#include "biclique.h"
#include "bigraph.h"
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <chrono>
#include <cstdlib>
#include <cmath>

using namespace std;

// Global variables (same as main.cpp) - only define the ones not already in biclique.cpp
long double Eps, Eps0, Eps1, Eps2, p;
bool use_probability_filtering;
vector<int> priv_deg;
int priv_dmax_1, priv_dmax_2, iteration, num_rounds;
vector<long double> estis, naive_estis;
long double RR_time, server_side_time, naive_server_side, local_count_time, deg_esti_time, communication_cost;
long double real = 0.0;  // Use long double consistently for ground truth
long double real_ld = 0.0;  // High precision ground truth for large numbers (kept for compatibility)
int P___, K___;
vector<vector<int>> up_options, lo_options;
long double avg_estimated_variance;

// Variables defined in biclique.cpp - declare as extern
extern bool multi_estimator_switch;
extern bool two_noisy_graph_switch;

// Function to get ground truth for a specific (P,Q) combination
long double get_ground_truth(string dataset, BiGraph& g, int P, int Q) {
    // Set global variables for the ground truth computation
    P___ = P;
    K___ = Q;
    
    // Call the function that computes or fetches ground truth
    fetch_or_compute_biclique_count(P___, K___, dataset, g);
    
    // Return the high precision value if available, otherwise the regular value
    return (real_ld > 0) ? real_ld : static_cast<long double>(real);
}

// Function to get ground truth for all Q values [4,5,6,7,8,9,10]
vector<long double> get_all_ground_truth(string dataset, BiGraph& g, int P) {
    vector<long double> ground_truth(7);
    
    for (int i = 0; i < 7; i++) {
        int Q = i + 4;  // Q = 4, 5, 6, 7, 8, 9, 10
        ground_truth[i] = get_ground_truth(dataset, g, P, Q);
    }
    
    return ground_truth;
}

// Function to run a specific algorithm for all Q values
vector<vector<long double>> run_algorithm(BiGraph& g, int algorithm, int num_rounds, double epsilon) {
    vector<vector<long double>> results(num_rounds, vector<long double>(7));
    
    // Set algorithm-specific parameters
    if (algorithm == 0) {
        // Naive algorithm
        cout << "Running Naive algorithm..." << endl;
    } else if (algorithm == 1) {
        // oneR
        cout << "Running oneR algorithm..." << endl;
    } else if (algorithm == 2) {
        // ADV
        cout << "Running ADV algorithm..." << endl;
        multi_estimator_switch = false;
        two_noisy_graph_switch = false;
    } else if (algorithm == 3) {
        // ADV+
        cout << "Running ADV+ algorithm..." << endl;
        multi_estimator_switch = true;
        two_noisy_graph_switch = false;
    } else if (algorithm == 4) {
        // ADV++
        cout << "Running ADV++ algorithm..." << endl;
        multi_estimator_switch = true;
        two_noisy_graph_switch = true;
    }
    
    for (int round = 0; round < num_rounds; round++) {
        unsigned long seed = 12345 + round;
        
    if (algorithm == 0) {
        // Run naive batch algorithm (builds noisy graph once, handles all Q values)
        vector<long double> batch_results = naive_biclique_batch(g, seed, 2);
        for (int i = 0; i < 7; i++) {
            results[round][i] = batch_results[i];
        }
    } else if (algorithm == 1) {
        // Run one_round batch algorithm (builds noisy graph once, handles all Q values)
        vector<long double> batch_results = one_round_biclique_2_K_batch(g, seed);
        for (int i = 0; i < 7; i++) {
            results[round][i] = batch_results[i];
        }
    } else {
        // Run batch function for ADV algorithms
        vector<long double> batch_results = wedge_based_two_round_2_K_biclique_batch(g, seed);
        for (int i = 0; i < 7; i++) {
            results[round][i] = batch_results[i];
        }
    }
    }
    
    return results;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cout << "Usage: " << argv[0] << " <dataset_path> <epsilon> <num_rounds>" << endl;
        cout << "Example: " << argv[0] << " ../bidata/co 1.0 10" << endl;
        cout << "Algorithms tested: 0=Naive, 1=oneR, 2=ADV, 3=ADV+, 4=ADV++" << endl;
        return 1;
    }
    
    string dataset_path = argv[1];
    double epsilon = atof(argv[2]);
    int num_rounds = atoi(argv[3]);
    
    cout << "=== Algorithm Comparison: (2,Q)-Biclique Counting ===" << endl;
    cout << "Dataset: " << dataset_path << endl;
    cout << "Epsilon: " << epsilon << endl;
    cout << "Rounds: " << num_rounds << endl;
    cout << "Q values: [2, 3, 4]" << endl;
    cout << "Algorithms: 0=Naive, 1=oneR, 2=ADV, 3=ADV+, 4=ADV++" << endl;
    cout << "=====================================================" << endl;
    
    // Load the graph
    BiGraph g(dataset_path);
    
    cout << "Graph loaded successfully:" << endl;
    cout << "  Vertices: " << g.num_nodes() << " (V1: " << g.num_v1 << ", V2: " << g.num_v2 << ")" << endl;
    cout << "  Edges: " << g.num_edges << endl;
    cout << endl;
    
    // Set global epsilon
    Eps = epsilon;
    use_probability_filtering = false; // No filtering for batch version
    
    // Get ground truth for all Q values (2,3,4)
    cout << "Computing ground truth for (2,Q)-bicliques..." << endl;
    vector<long double> ground_truth = get_all_ground_truth(dataset_path, g, 2);
    for (int i = 0; i < 3; i++) {
        int Q = i + 2;
        cout << "  (2," << Q << ")-bicliques: " << ground_truth[i] << endl;
    }
    cout << endl;
    
    // Test all algorithms
    vector<int> algorithms = {0, 1, 2, 3, 4}; // Naive, oneR, ADV, ADV+, ADV++
    vector<string> algorithm_names = {"Naive", "oneR", "ADV", "ADV+", "ADV++"};
    
    // Store results for all algorithms
    vector<vector<vector<long double>>> all_algorithm_results(5);
    vector<vector<vector<long double>>> all_algorithm_errors(5);
    
    auto total_start_time = chrono::high_resolution_clock::now();
    
    for (int alg_idx = 0; alg_idx < algorithms.size(); alg_idx++) {
        int algorithm = algorithms[alg_idx];
        cout << "\n=== Testing " << algorithm_names[alg_idx] << " (Algorithm " << algorithm << ") ===" << endl;
        
        // Run algorithm
        vector<vector<long double>> results = run_algorithm(g, algorithm, num_rounds, epsilon);
        all_algorithm_results[alg_idx] = results;
        
        // Calculate relative errors using high precision ground truth
        vector<vector<long double>> relative_errors(num_rounds, vector<long double>(7));
        for (int round = 0; round < num_rounds; round++) {
            for (int i = 0; i < 7; i++) {
                if (ground_truth[i] > 0) {
                    relative_errors[round][i] = abs(results[round][i] - ground_truth[i]) / ground_truth[i];
                } else {
                    relative_errors[round][i] = (results[round][i] > 0) ? 1.0 : 0.0;
                }
            }
        }
        all_algorithm_errors[alg_idx] = relative_errors;
    }
    
    auto total_end_time = chrono::high_resolution_clock::now();
    auto total_duration = chrono::duration_cast<chrono::milliseconds>(total_end_time - total_start_time);
    
    cout << endl;
    cout << "=== Algorithm Comparison Results ===" << endl;
    cout << "Total execution time: " << total_duration.count() << " ms" << endl;
    cout << "Average time per algorithm: " << (total_duration.count() / algorithms.size()) << " ms" << endl;
    cout << endl;
    
    // Calculate and display statistics for each algorithm
    for (int alg_idx = 0; alg_idx < algorithms.size(); alg_idx++) {
        cout << "\n=== " << algorithm_names[alg_idx] << " Results ===" << endl;
        cout << "Q\tGround Truth\tMean Estimate\tMean Rel Error\tStd Rel Error" << endl;
        cout << "-\t------------\t-------------\t-------------\t-------------" << endl;
        
        for (int i = 0; i < 3; i++) {
            int Q = i + 2;
            
            // Calculate mean estimate
            long double mean_estimate = 0;
            for (int round = 0; round < num_rounds; round++) {
                mean_estimate += all_algorithm_results[alg_idx][round][i];
            }
            mean_estimate /= num_rounds;
            
            // Calculate mean relative error
            long double mean_rel_error = 0;
            for (int round = 0; round < num_rounds; round++) {
                mean_rel_error += all_algorithm_errors[alg_idx][round][i];
            }
            mean_rel_error /= num_rounds;
            
            // Calculate standard deviation of relative error
            long double std_rel_error = 0;
            for (int round = 0; round < num_rounds; round++) {
                long double diff = all_algorithm_errors[alg_idx][round][i] - mean_rel_error;
                std_rel_error += diff * diff;
            }
            std_rel_error = sqrt(std_rel_error / num_rounds);
            
            cout << Q << "\t" << ground_truth[i] << "\t" 
                 << fixed << setprecision(2) << mean_estimate << "\t"
                 << scientific << setprecision(3) << mean_rel_error << "\t"
                 << scientific << setprecision(3) << std_rel_error << endl;
        }
    }
    
    // Algorithm comparison summary
    cout << "\n=== Algorithm Comparison Summary ===" << endl;
    cout << "Q\tNaive\t\toneR\t\tADV\t\tADV+\t\tADV++" << endl;
    cout << "-\t-----\t\t----\t\t---\t\t----\t\t-----" << endl;
    
    for (int i = 0; i < 7; i++) {
        int Q = i + 4;
        cout << Q;
        
        for (int alg_idx = 0; alg_idx < algorithms.size(); alg_idx++) {
            // Calculate mean relative error for this algorithm and Q
            long double mean_rel_error = 0;
            for (int round = 0; round < num_rounds; round++) {
                mean_rel_error += all_algorithm_errors[alg_idx][round][i];
            }
            mean_rel_error /= num_rounds;
            
            cout << "\t" << scientific << setprecision(2) << mean_rel_error;
        }
        cout << endl;
    }
    
    return 0;
}
