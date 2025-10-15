#include <iostream>
#include <sqlite3.h>
#include <iomanip>

int main() {
    sqlite3* db;
    if (sqlite3_open("../biclq_counts.db", &db) != SQLITE_OK) {
        std::cerr << "Error opening database: " << sqlite3_errmsg(db) << std::endl;
        return 1;
    }
    
    std::cout << "=== Testing Ground Truth Retrieval Fix ===" << std::endl;
    
    // Test query for lrcwiki Q=6-9
    const char* sql = "SELECT dataset, count FROM pqbiclique_counts WHERE dataset = 'lrcwiki' AND p = 2 AND q >= 6 ORDER BY q;";
    sqlite3_stmt* stmt;
    
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        std::cerr << "Error preparing statement: " << sqlite3_errmsg(db) << std::endl;
        sqlite3_close(db);
        return 1;
    }
    
    std::cout << "Q\tOld Method (int64)\tNew Method (double)\tDifference" << std::endl;
    std::cout << "-\t------------------\t------------------\t----------" << std::endl;
    
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        std::string dataset = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
        
        // Old method (causes overflow)
        long long old_count = sqlite3_column_int64(stmt, 1);
        
        // New method (handles REAL values)
        int column_type = sqlite3_column_type(stmt, 1);
        long double new_count;
        
        if (column_type == SQLITE_INTEGER) {
            new_count = static_cast<long double>(sqlite3_column_int64(stmt, 1));
        } else if (column_type == SQLITE_FLOAT) {
            new_count = static_cast<long double>(sqlite3_column_double(stmt, 1));
        } else {
            new_count = 0;
        }
        
        // Get Q value from the query (we need to modify the query to include it)
        // For now, let's just show the results
        std::cout << "?\t" << old_count << "\t\t\t" 
                  << std::scientific << std::setprecision(6) << new_count << "\t\t"
                  << (new_count - old_count) << std::endl;
    }
    
    sqlite3_finalize(stmt);
    sqlite3_close(db);
    
    return 0;
}
