#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include "biclique.h"
#include "fixed_ground_truth.cpp"

// Global variables (same as in test_batch_with_ground_truth.cpp)
double Eps = 1.0;
bool use_probability_filtering = false;
bool multi_estimator_switch = false;
bool two_noisy_graph_switch = false;
int P___ = 2;
int K___ = 4;
unsigned long long real = 0;

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " <dataset_path>" << endl;
        cout << "Example: " << argv[0] << " ../bidata/lrcwiki" << endl;
        return 1;
    }
    
    string dataset_path = argv[1];
    
    cout << "=== Testing Fixed Ground Truth Computation ===" << endl;
    cout << "Dataset: " << dataset_path << endl;
    cout << "=============================================" << endl;
    
    // Load the graph
    BiGraph g(dataset_path);
    
    cout << "Graph loaded successfully:" << endl;
    cout << "  Vertices: " << g.num_nodes() << " (V1: " << g.num_v1 << ", V2: " << g.num_v2 << ")" << endl;
    cout << "  Edges: " << g.num_edges << endl;
    cout << endl;
    
    // Get ground truth for all Q values with fixed precision
    cout << "Computing ground truth for (2,Q)-bicliques with fixed precision..." << endl;
    vector<long double> ground_truth = get_all_ground_truth_fixed(dataset_path, g, 2);
    
    cout << "\n=== Fixed Ground Truth Results ===" << endl;
    for (int i = 0; i < 7; i++) {
        int Q = i + 4;
        cout << "  (2," << Q << ")-bicliques: " << scientific << setprecision(6) << ground_truth[i] << endl;
    }
    cout << endl;
    
    // Compare with old ground truth
    cout << "=== Comparison with Old Ground Truth ===" << endl;
    cout << "Q\tOld Ground Truth\t\tFixed Ground Truth\t\tDifference" << endl;
    cout << "-\t---------------\t\t------------------\t\t----------" << endl;
    
    // Get old ground truth for comparison
    vector<unsigned long long> old_ground_truth(7);
    for (int i = 0; i < 7; i++) {
        int Q = i + 4;
        P___ = 2;
        K___ = Q;
        fetch_or_compute_biclique_count(P___, K___, dataset_path, g);
        old_ground_truth[i] = real;
    }
    
    for (int i = 0; i < 7; i++) {
        int Q = i + 4;
        long double difference = abs(ground_truth[i] - static_cast<long double>(old_ground_truth[i]));
        cout << Q << "\t" << scientific << setprecision(6) << old_ground_truth[i] 
             << "\t\t" << ground_truth[i] 
             << "\t\t" << difference << endl;
    }
    
    return 0;
}
