#!/usr/bin/env python3
"""
Add librec-filmtrust-ratings Naive results to the algorithm_comparison_p3 table.
"""

import sqlite3

def add_filmtrust_naive_results():
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # Ground truth values from the output
    ground_truth_values = {
        4: 1.110059e+12,
        5: 8.618833e+12,
        6: 5.800279e+13,
        7: 3.381429e+14,
        8: 1.719092e+15,
        9: 7.683863e+15,
        10: 3.043229e+16
    }
    
    # Naive results from the experiments
    naive_results = {
        4: 1.126060e+14,
        5: 8.357680e+14,
        6: 4.348975e+15,
        7: 4.716876e+16,
        8: 2.166225e+17,
        9: 1.407367e+18,
        10: 4.391287e+18
    }
    
    dataset = "librec-filmtrust-ratings"
    epsilon = 1.0
    p = 3
    T_samples = 10  # 10 samples per iteration (vertex sampling)
    
    # Insert Naive results
    print(f"Adding {dataset} Naive results...")
    for q in range(4, 11):
        ground_truth = ground_truth_values[q]
        estimate = naive_results[q]
        
        # Calculate relative error
        if ground_truth > 0:
            rel_error = abs(estimate - ground_truth) / ground_truth
        else:
            rel_error = float('inf')
        
        # Insert into database
        cursor.execute('''
            INSERT OR REPLACE INTO algorithm_comparison_p3 
            (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (dataset, 'Naive', q, ground_truth, estimate, rel_error, 0.0, T_samples, epsilon, p))
        
        print(f"  Naive Q={q}: GT={ground_truth:.2e}, Est={estimate:.2e}, Rel Error={rel_error:.3f}")
    
    conn.commit()
    conn.close()
    print(f"\nSuccessfully added {dataset} Naive results to database!")

if __name__ == "__main__":
    add_filmtrust_naive_results()
