#include "biclique.h"
#include <iostream>
#include <string>

using namespace std;

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
    
    cout << "Running f distribution test..." << endl;
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
