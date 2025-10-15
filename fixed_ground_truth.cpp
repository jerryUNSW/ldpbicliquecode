#include <iostream>
#include <string>
#include <vector>
#include <sqlite3.h>
#include <sstream>
#include <iomanip>
#include "biclique.h"

// Global variables
extern unsigned long long real;
extern int P___;
extern int K___;

// Function to parse scientific notation string to long double
long double parseScientificNotation(const std::string& str) {
    if (str.empty()) return 0.0;
    
    // Handle regular numbers
    if (str.find('e') == std::string::npos && str.find('E') == std::string::npos) {
        return std::stold(str);
    }
    
    // Handle scientific notation
    std::istringstream iss(str);
    long double value;
    iss >> value;
    return value;
}

// Fixed version that handles large numbers properly
void fetch_or_compute_biclique_count_fixed(int P___, int K___, 
    string dataset, BiGraph& g) {

    sqlite3* db;
    if (sqlite3_open("../biclq_counts.db", &db) != SQLITE_OK) {
        std::cerr << "Error opening database: " << sqlite3_errmsg(db) << std::endl;
        exit(1);
    }
    
    // Dataset, p, and q values to filter
    size_t found = dataset.find_last_of("/\\");
    std::string dataset_to_find = dataset.substr(found + 1);

    int p_to_find = P___;
    int q_to_find = K___;

    // Query to retrieve from the fixed table
    const char* sql = "SELECT dataset, count FROM pqbiclique_counts_fixed WHERE dataset = ? AND p = ? AND q = ?;";
    sqlite3_stmt* stmt;

    // Prepare statement
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        std::cerr << "Error preparing statement: " << sqlite3_errmsg(db) << std::endl;
        sqlite3_close(db);
        exit(1);
    }

    // Bind parameters
    sqlite3_bind_text(stmt, 1, dataset_to_find.c_str(), -1, SQLITE_STATIC);
    sqlite3_bind_int(stmt, 2, p_to_find);
    sqlite3_bind_int(stmt, 3, q_to_find);

    // Execute the query
    if (sqlite3_step(stmt) == SQLITE_ROW) {
        std::string dataset_name = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
        std::string count_str = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1));
        
        long double count_value = parseScientificNotation(count_str);
        
        std::cout << "Dataset: " << dataset_name << ", p = " << P___ << ", q = " << K___ 
                  << ", biclique count = " << std::scientific << std::setprecision(6) << count_value << std::endl;
        
        // Store as long double to preserve precision
        real = static_cast<unsigned long long>(count_value);
        
    } else {
        std::cout << "No matching data found in fixed table. Computing ground truth..." << std::endl;

        // Compute ground truth using exact counting
        biGraph convertedGraph = convertBiGraphTobiGraph(g);
        std::cout << "Converted graph: n1=" << convertedGraph.n1 
                  << ", n2=" << convertedGraph.n2 
                  << ", m=" << convertedGraph.m << std::endl;
        
        BCListPlusPlus* counter = new BCListPlusPlus(&convertedGraph, P___, K___);
        long double exact_count = counter->exactCount();
        delete counter;
        
        std::cout << "Exact count = " << std::scientific << std::setprecision(6) << exact_count << std::endl;

        // Insert into both tables
        const char* insert_sql = "INSERT INTO pqbiclique_counts_fixed (dataset, p, q, count) VALUES (?, ?, ?, ?);";
        sqlite3_stmt* insert_stmt;

        if (sqlite3_prepare_v2(db, insert_sql, -1, &insert_stmt, nullptr) != SQLITE_OK) {
            std::cerr << "Error preparing insert statement: " << sqlite3_errmsg(db) << std::endl;
            sqlite3_close(db);
            exit(1);
        }

        // Convert to scientific notation string for storage
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(6) << exact_count;
        std::string count_str = oss.str();

        sqlite3_bind_text(insert_stmt, 1, dataset_to_find.c_str(), -1, SQLITE_STATIC);
        sqlite3_bind_int(insert_stmt, 2, P___);
        sqlite3_bind_int(insert_stmt, 3, K___);
        sqlite3_bind_text(insert_stmt, 4, count_str.c_str(), -1, SQLITE_STATIC);

        if (sqlite3_step(insert_stmt) != SQLITE_DONE) {
            std::cerr << "Error inserting data: " << sqlite3_errmsg(db) << std::endl;
        } else {
            std::cout << "Inserted new biclique count into fixed database." << std::endl;
        }
        
        sqlite3_finalize(insert_stmt);
        real = static_cast<unsigned long long>(exact_count);
    }
    
    sqlite3_finalize(stmt);
    sqlite3_close(db);
}

// Function to get ground truth for all Q values with fixed precision
std::vector<long double> get_all_ground_truth_fixed(string dataset, BiGraph& g, int P) {
    std::vector<long double> ground_truth(7);
    
    for (int i = 0; i < 7; i++) {
        int Q = i + 4;  // Q = 4, 5, 6, 7, 8, 9, 10
        
        // Set global variables
        P___ = P;
        K___ = Q;
        
        // Get ground truth with fixed precision
        fetch_or_compute_biclique_count_fixed(P___, K___, dataset, g);
        ground_truth[i] = static_cast<long double>(real);
    }
    
    return ground_truth;
}
