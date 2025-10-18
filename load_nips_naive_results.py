#!/usr/bin/env python3
"""
Load nips dataset Naive results into the algorithm_comparison_p3 table.
"""

import sqlite3

def load_nips_naive_results():
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # Naive results for nips (extracted from p3_nips_new.txt)
    naive_results = [
        (4, 1.729841e+14, 0.0, 1.557307e+03, 0.0),  # Ground truth, estimate (unknown), relative error, std error
        (5, 2.026990e+15, 0.0, 7.404205e+03, 0.0),
        (6, 2.054109e+16, 0.0, 3.392706e+04, 0.0),
        (7, 1.848309e+17, 0.0, 1.983840e+05, 0.0),
        (8, 1.505491e+18, 0.0, 6.806663e+05, 0.0),
        (9, 1.126738e+19, 0.0, 4.067774e+06, 0.0),
        (10, 7.844682e+19, 0.0, 1.627487e+07, 0.0)
    ]
    
    dataset = "nips"
    epsilon = 1.0
    p = 3
    T_samples = -1  # All triplets used
    
    # Insert Naive results
    print("Inserting Naive results for nips...")
    for q, ground_truth, estimate, relative_error, std_relative_error in naive_results:
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'Naive', q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p))
            print(f"  Q={q}: Ground Truth={ground_truth:.2e}, Rel Error={relative_error:.3e}")
        except Exception as e:
            print(f"Error inserting Naive Q={q}: {e}")
    
    # Commit changes
    conn.commit()
    
    # Verify data
    print(f"\nVerifying data for {dataset}...")
    cursor.execute('''
        SELECT algorithm, q, relative_error 
        FROM algorithm_comparison_p3 
        WHERE dataset = ? AND p = ? AND epsilon = ?
        ORDER BY algorithm, q
    ''', (dataset, p, epsilon))
    
    results = cursor.fetchall()
    for algorithm, q, rel_error in results:
        print(f"  {algorithm} Q={q}: Rel Error={rel_error:.3e}")
    
    conn.close()
    print(f"\nSuccessfully loaded {dataset} Naive results into database!")

if __name__ == "__main__":
    load_nips_naive_results()

