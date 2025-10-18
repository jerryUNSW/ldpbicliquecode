#include "biclique.h"
#include <iostream>
#include <string>

using namespace std;

// Global variables needed for compilation (copied from main.cpp)
long double p, Eps, Eps0, Eps1, Eps2;
long double alpha0, alpha1, alpha2;
vector<int> priv_deg;
int priv_dmax_1, priv_dmax_2, iteration, num_rounds;
long double m3__ = 0, m2__ = 0, m1__ = 0, m0__ = 0, real_stars = 0,
            RR_time, server_side_time, naive_server_side, local_count_time = 0, deg_esti_time = 0,
            communication_cost = 0;
long double real = 0.0;
long double real_ld = 0.0;
int P___, K___;
vector<vector<int>> up_options, lo_options;
// bool two_noisy_graph_switch = false;  // Already defined in biclique.cpp
// bool multi_estimator_switch = false;  // Already defined in biclique.cpp
bool use_probability_filtering = false;
vector<long double> estis, naive_estis;

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <dataset> [num_samples]" << endl;
        cout << "Example: " << argv[0] << " librec-filmtrust-ratings 5000" << endl;
        return 1;
    }
    
    string dataset = argv[1];
    int num_samples = 5000;
    
    if (argc >= 3) {
        num_samples = atoi(argv[2]);
    }
    
    cout << "Running f_avg distribution test..." << endl;
    cout << "Dataset: " << dataset << endl;
    cout << "Number of samples: " << num_samples << endl;
    
    // Initialize random seed
    srand(42);  // Fixed seed for reproducibility
    
    try {
        test_f_distribution_p2(dataset, num_samples);
        cout << "\nTest completed successfully!" << endl;
    } catch (const exception& e) {
        cout << "Error during test: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}
