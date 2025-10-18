#include "biclique.h"
#include "bigraph.h"
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace std;

// Global variables needed for compilation
long double Eps = 1.0, Eps0 = 0.05, Eps1 = 0.3, Eps2 = 0.65, p = 0.1;
long double m3__ = 0, m2__ = 0, m1__ = 0, m0__ = 0, communication_cost = 0;
long double alpha0 = 0.05, alpha1 = 0.3, alpha2 = 0.65;
vector<int> priv_deg;
int priv_dmax_1 = 0, priv_dmax_2 = 0;
vector<vector<int>> up_options, lo_options;
long double real_ld = 0;
bool use_probability_filtering = false;
bool two_noisy_graph_switch = false;

// Simple implementation of missing functions
long double power(long double base, int exponent) {
    long double result = 1.0;
    for (int i = 0; i < exponent; i++) {
        result *= base;
    }
    return result;
}

void print_dash(int n) {
    for (int i = 0; i < n; i++) {
        cout << "-";
    }
    cout << endl;
}

void test_f_distribution_real(string dataset, int num_samples) {
    cout << "=== Testing f_u distribution with REAL implementation ===" << endl;
    
    // Load dataset
    BiGraph g;
    string filename = "data/" + dataset;
    g.loadGraph(filename);
    cout << "Loaded dataset: " << dataset << " with " << g.num_nodes() << " nodes" << endl;
    
    // Set up parameters
    p = 0.1;  // RR parameter
    Eps = 1.0;  // Privacy parameter
    Eps0 = 0.05;  // Degree estimation
    Eps1 = 0.3;   // First round
    Eps2 = 0.65;  // Second round
    
    // Collect f values for analysis
    vector<double> f_values_q2, f_values_q3, f_values_q4;
    
    cout << "Sampling " << num_samples << " f values using REAL implementation..." << endl;
    
    // Initialize random seed
    srand(42);
    
    for (int sample = 0; sample < num_samples; sample++) {
        // Create noisy graph using REAL implementation
        BiGraph g2(g);
        construct_noisy_graph(g, g2, sample);
        
        // Sample random vertex pairs
        int u = rand() % g.num_v1;
        int w = g.num_v1 + (rand() % g.num_v2);
        
        // Compute f using REAL locally_compute_f_given_q_and_x function
        try {
            double f_q2 = locally_compute_f_given_q_and_x(2, u, g, g2);
            double f_q3 = locally_compute_f_given_q_and_x(3, u, g, g2);
            double f_q4 = locally_compute_f_given_q_and_x(4, u, g, g2);
            
            f_values_q2.push_back(f_q2);
            f_values_q3.push_back(f_q3);
            f_values_q4.push_back(f_q4);
        } catch (const exception& e) {
            cout << "Error computing f for sample " << sample << ": " << e.what() << endl;
            continue;
        }
    }
    
    // Compute statistics for q=2
    if (!f_values_q2.empty()) {
        double mean_f2 = accumulate(f_values_q2.begin(), f_values_q2.end(), 0.0) / f_values_q2.size();
        double var_f2 = 0.0;
        for (double f : f_values_q2) {
            var_f2 += (f - mean_f2) * (f - mean_f2);
        }
        var_f2 /= f_values_q2.size();
        
        // Compute skewness and kurtosis
        double skew_f2 = 0.0, kurt_f2 = 0.0;
        for (double f : f_values_q2) {
            double normalized = (f - mean_f2) / sqrt(var_f2);
            skew_f2 += normalized * normalized * normalized;
            kurt_f2 += normalized * normalized * normalized * normalized;
        }
        skew_f2 /= f_values_q2.size();
        kurt_f2 /= f_values_q2.size();
        kurt_f2 -= 3.0;  // Excess kurtosis
        
        cout << "\n=== REAL f_u Distribution Statistics ===" << endl;
        cout << "q=2: mean=" << mean_f2 << ", var=" << var_f2 << ", std=" << sqrt(var_f2) << endl;
        cout << "q=2: skewness=" << skew_f2 << ", excess_kurtosis=" << kurt_f2 << endl;
    }
    
    // Compute statistics for q=3
    if (!f_values_q3.empty()) {
        double mean_f3 = accumulate(f_values_q3.begin(), f_values_q3.end(), 0.0) / f_values_q3.size();
        double var_f3 = 0.0;
        for (double f : f_values_q3) {
            var_f3 += (f - mean_f3) * (f - mean_f3);
        }
        var_f3 /= f_values_q3.size();
        
        double skew_f3 = 0.0, kurt_f3 = 0.0;
        for (double f : f_values_q3) {
            double normalized = (f - mean_f3) / sqrt(var_f3);
            skew_f3 += normalized * normalized * normalized;
            kurt_f3 += normalized * normalized * normalized * normalized;
        }
        skew_f3 /= f_values_q3.size();
        kurt_f3 /= f_values_q3.size();
        kurt_f3 -= 3.0;
        
        cout << "q=3: mean=" << mean_f3 << ", var=" << var_f3 << ", std=" << sqrt(var_f3) << endl;
        cout << "q=3: skewness=" << skew_f3 << ", excess_kurtosis=" << kurt_f3 << endl;
    }
    
    // Compute statistics for q=4
    if (!f_values_q4.empty()) {
        double mean_f4 = accumulate(f_values_q4.begin(), f_values_q4.end(), 0.0) / f_values_q4.size();
        double var_f4 = 0.0;
        for (double f : f_values_q4) {
            var_f4 += (f - mean_f4) * (f - mean_f4);
        }
        var_f4 /= f_values_q4.size();
        
        double skew_f4 = 0.0, kurt_f4 = 0.0;
        for (double f : f_values_q4) {
            double normalized = (f - mean_f4) / sqrt(var_f4);
            skew_f4 += normalized * normalized * normalized;
            kurt_f4 += normalized * normalized * normalized * normalized;
        }
        skew_f4 /= f_values_q4.size();
        kurt_f4 /= f_values_q4.size();
        kurt_f4 -= 3.0;
        
        cout << "q=4: mean=" << mean_f4 << ", var=" << var_f4 << ", std=" << sqrt(var_f4) << endl;
        cout << "q=4: skewness=" << skew_f4 << ", excess_kurtosis=" << kurt_f4 << endl;
    }
    
    cout << "\n=== REAL Implementation Analysis ===" << endl;
    cout << "Using actual locally_compute_f_given_q_and_x function" << endl;
    cout << "Higher q values should show more Gaussian behavior due to:" << endl;
    cout << "1. More terms in the sum (CLT applies)" << endl;
    cout << "2. Independence of noise sources" << endl;
    cout << "3. Bounded influence of individual terms" << endl;
    
    // Output data for Python visualization
    cout << "\n=== Data for Python Visualization ===" << endl;
    cout << "f_values_q2: ";
    for (int i = 0; i < min(20, (int)f_values_q2.size()); i++) {
        cout << f_values_q2[i] << " ";
    }
    cout << endl;
    
    cout << "f_values_q3: ";
    for (int i = 0; i < min(20, (int)f_values_q3.size()); i++) {
        cout << f_values_q3[i] << " ";
    }
    cout << endl;
    
    cout << "f_values_q4: ";
    for (int i = 0; i < min(20, (int)f_values_q4.size()); i++) {
        cout << f_values_q4[i] << " ";
    }
    cout << endl;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <dataset> [num_samples]" << endl;
        cout << "Example: " << argv[0] << " unicode 1000" << endl;
        return 1;
    }
    
    string dataset = argv[1];
    int num_samples = 1000;
    
    if (argc >= 3) {
        num_samples = atoi(argv[2]);
    }
    
    cout << "Running REAL f distribution test..." << endl;
    cout << "Dataset: " << dataset << endl;
    cout << "Number of samples: " << num_samples << endl;
    
    // Initialize random seed
    srand(42);  // Fixed seed for reproducibility
    
    try {
        test_f_distribution_real(dataset, num_samples);
        cout << "\nTest completed successfully!" << endl;
    } catch (const exception& e) {
        cout << "Error during test: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}
