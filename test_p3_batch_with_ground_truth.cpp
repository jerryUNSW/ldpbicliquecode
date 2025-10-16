/*
 * test_p3_batch_with_ground_truth.cpp
 * 
 * Test P=3 algorithms with proper ground truth comparison
 * Runs 10 rounds and calculates relative errors for Q = [4,5,6,7,8,9,10]
 * Only tests Q values where (3,Q)-biclique count > 0
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
long double real = 0.0;
long double real_ld = 0.0;  // High precision ground truth for large numbers
int P___, K___;
vector<vector<int>> up_options, lo_options;
long double avg_estimated_variance;

// Variables defined in biclique.cpp - declare as extern
extern bool multi_estimator_switch;
extern bool two_noisy_graph_switch;
extern vector<long double> deg_estis;

// Function to get ground truth for a specific (P,Q) combination
long double get_ground_truth(string dataset, BiGraph& g, int P, int Q) {
    // Set global variables for the ground truth computation
    P___ = P;
    K___ = Q;
    
    // Call the function that computes or fetches ground truth
    fetch_or_compute_biclique_count(P___, K___, dataset, g);
    
    // Return the high precision value if available, otherwise the regular value
    return real;  // Now both real and real_ld are long double
}

// Function to get ground truth for all Q values [2,3,4] and filter out zero counts
vector<long double> get_all_ground_truth_p3(string dataset, BiGraph& g, int P) {
    vector<long double> ground_truth(3);
    vector<bool> valid_q(3, false);
    
    cout << "Checking ground truth for (3,Q)-bicliques..." << endl;
    
    for (int i = 0; i < 3; i++) {
        int Q = i + 2;  // Q = 2, 3, 4
        ground_truth[i] = get_ground_truth(dataset, g, P, Q);
        
        if (ground_truth[i] > 0) {
            valid_q[i] = true;
            cout << "  (3," << Q << ")-bicliques: " << ground_truth[i] << " [VALID]" << endl;
        } else {
            cout << "  (3," << Q << ")-bicliques: " << ground_truth[i] << " [SKIP - zero count]" << endl;
        }
    }
    
    // Filter out zero counts
    vector<long double> filtered_ground_truth;
    vector<int> valid_q_values;
    
    for (int i = 0; i < 3; i++) {
        if (valid_q[i]) {
            filtered_ground_truth.push_back(ground_truth[i]);
            valid_q_values.push_back(i + 2);  // Q values
        }
    }
    
    cout << "Valid Q values to test: ";
    for (int q_val : valid_q_values) {
        cout << q_val << " ";
    }
    cout << endl;
    
    return filtered_ground_truth;
}

// Function to run a specific algorithm for valid Q values only
vector<vector<long double>> run_algorithm_p3(BiGraph& g, int algorithm, int num_rounds, double epsilon, 
                                           const vector<int>& valid_q_values) {
    int num_valid_q = valid_q_values.size();
    vector<vector<long double>> results(num_rounds, vector<long double>(num_valid_q));
    
    // Set algorithm-specific parameters
    if (algorithm == 0) {
        // Naive algorithm
        cout << "Running Naive algorithm..." << endl;
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
            vector<long double> batch_results = naive_biclique_with_vertex_sampling_batch(g, seed, 3, 0.1, 20, round);
            for (int q_idx = 0; q_idx < num_valid_q; q_idx++) {
                int Q = valid_q_values[q_idx];
                int q_array_idx = Q - 2;  // Convert Q to array index (Q=2 -> idx=0, Q=3 -> idx=1, Q=4 -> idx=2)
                results[round][q_idx] = batch_results[q_array_idx];
            }
        } else {
            // Run ADV batch algorithm (builds noisy graph once, handles all Q values)
            vector<long double> batch_results = wedge_based_two_round_3_K_biclique_rejection_sampling_batch(g, seed);
            for (int q_idx = 0; q_idx < num_valid_q; q_idx++) {
                int Q = valid_q_values[q_idx];
                int q_array_idx = Q - 2;  // Convert Q to array index (Q=2 -> idx=0, Q=3 -> idx=1, Q=4 -> idx=2)
                results[round][q_idx] = batch_results[q_array_idx];
            }
        }
    }
    
    return results;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cout << "Usage: " << argv[0] << " <dataset_path> <epsilon> <num_rounds>" << endl;
        cout << "Example: " << argv[0] << " ../bidata/librec-filmtrust-ratings 1.0 10" << endl;
        cout << "Algorithms tested: 0=Naive, 2=ADV, 3=ADV+, 4=ADV++" << endl;
        return 1;
    }
    
    string dataset_path = argv[1];
    double epsilon = atof(argv[2]);
    int num_rounds = atoi(argv[3]);
    
    cout << "=== Algorithm Comparison: (3,Q)-Biclique Counting ===" << endl;
    cout << "Dataset: " << dataset_path << endl;
    cout << "Epsilon: " << epsilon << endl;
    cout << "Rounds: " << num_rounds << endl;
    cout << "Q values: [2, 3, 4] (filtered by non-zero ground truth)" << endl;
    cout << "Algorithms: 0=Naive, 2=ADV, 3=ADV+, 4=ADV++" << endl;
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
    
    // Get ground truth for all Q values and filter out zero counts
    cout << "Computing ground truth for (3,Q)-bicliques..." << endl;
    vector<long double> ground_truth = get_all_ground_truth_p3(dataset_path, g, 3);
    
    if (ground_truth.empty()) {
        cout << "ERROR: No valid (3,Q)-bicliques found for any Q value!" << endl;
        return 1;
    }
    
    // Determine valid Q values (Q = 2,3,4)
    vector<int> valid_q_values;
    for (int i = 0; i < 3; i++) {
        int Q = i + 2;
        long double gt = get_ground_truth(dataset_path, g, 3, Q);
        if (gt > 0) {
            valid_q_values.push_back(Q);
        }
    }
    
    cout << endl;
    
    // Test all algorithms
    vector<int> algorithms = {0, 2, 3, 4}; // Naive, ADV, ADV+, ADV++
    vector<string> algorithm_names = {"Naive", "ADV", "ADV+", "ADV++"};
    
    // Store results for all algorithms
    vector<vector<vector<long double>>> all_algorithm_results(4);
    vector<vector<vector<long double>>> all_algorithm_errors(4);
    
    auto total_start_time = chrono::high_resolution_clock::now();
    
    for (int alg_idx = 0; alg_idx < algorithms.size(); alg_idx++) {
        int algorithm = algorithms[alg_idx];
        cout << "\n=== Testing " << algorithm_names[alg_idx] << " (Algorithm " << algorithm << ") ===" << endl;
        
        // Run algorithm
        vector<vector<long double>> results = run_algorithm_p3(g, algorithm, num_rounds, epsilon, valid_q_values);
        all_algorithm_results[alg_idx] = results;
        
        // Calculate relative errors using high precision ground truth
        vector<vector<long double>> relative_errors(num_rounds, vector<long double>(valid_q_values.size()));
        for (int round = 0; round < num_rounds; round++) {
            for (int q_idx = 0; q_idx < valid_q_values.size(); q_idx++) {
                if (ground_truth[q_idx] > 0) {
                    relative_errors[round][q_idx] = abs(results[round][q_idx] - ground_truth[q_idx]) / ground_truth[q_idx];
                } else {
                    relative_errors[round][q_idx] = (results[round][q_idx] > 0) ? 1.0 : 0.0;
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
        
        for (int q_idx = 0; q_idx < valid_q_values.size(); q_idx++) {
            int Q = valid_q_values[q_idx];
            
            // Calculate mean estimate
            long double mean_estimate = 0;
            for (int round = 0; round < num_rounds; round++) {
                mean_estimate += all_algorithm_results[alg_idx][round][q_idx];
            }
            mean_estimate /= num_rounds;
            
            // Calculate mean relative error
            long double mean_rel_error = 0;
            for (int round = 0; round < num_rounds; round++) {
                mean_rel_error += all_algorithm_errors[alg_idx][round][q_idx];
            }
            mean_rel_error /= num_rounds;
            
            // Calculate standard deviation of relative error
            long double std_rel_error = 0;
            for (int round = 0; round < num_rounds; round++) {
                long double diff = all_algorithm_errors[alg_idx][round][q_idx] - mean_rel_error;
                std_rel_error += diff * diff;
            }
            std_rel_error = sqrt(std_rel_error / num_rounds);
            
            cout << Q << "\t" << ground_truth[q_idx] << "\t" 
                 << fixed << setprecision(2) << mean_estimate << "\t"
                 << scientific << setprecision(3) << mean_rel_error << "\t"
                 << scientific << setprecision(3) << std_rel_error << endl;
        }
    }
    
    // Algorithm comparison summary
    cout << "\n=== Algorithm Comparison Summary ===" << endl;
    cout << "Q\tNaive\t\tADV\t\tADV+\t\tADV++" << endl;
    cout << "-\t-----\t\t---\t\t----\t\t-----" << endl;
    
    for (int q_idx = 0; q_idx < valid_q_values.size(); q_idx++) {
        int Q = valid_q_values[q_idx];
        cout << Q;
        
        for (int alg_idx = 0; alg_idx < algorithms.size(); alg_idx++) {
            // Calculate mean relative error for this algorithm and Q
            long double mean_rel_error = 0;
            for (int round = 0; round < num_rounds; round++) {
                mean_rel_error += all_algorithm_errors[alg_idx][round][q_idx];
            }
            mean_rel_error /= num_rounds;
            
            cout << "\t" << scientific << setprecision(2) << mean_rel_error;
        }
        cout << endl;
    }
    
    return 0;
}
