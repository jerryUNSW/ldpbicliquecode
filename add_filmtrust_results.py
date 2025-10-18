#!/usr/bin/env python3
"""
Add librec-filmtrust-ratings ADV and ADV+ results to the algorithm_comparison_p3 table.
"""

import sqlite3

def add_filmtrust_results():
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
    
    # ADV results
    adv_results = {
        4: 1.716344e+13,
        5: 5.464494e+13,
        6: 2.351488e+15,
        7: 9.168473e+15,
        8: 5.454362e+17,
        9: 1.141966e+18,
        10: 9.990508e+19
    }
    
    # ADV+ results
    adv_plus_results = {
        4: 1.047784e+13,
        5: 1.858006e+13,
        6: -8.420763e+14,
        7: -7.151986e+14,
        8: 5.570240e+16,
        9: 1.353813e+15,
        10: -3.151237e+18
    }
    
    dataset = "librec-filmtrust-ratings"
    epsilon = 1.0
    p = 3
    T_samples = -1  # No sampling (all triplets used)
    
    # Insert ADV results
    print(f"Adding {dataset} ADV results...")
    for q in range(4, 11):
        ground_truth = ground_truth_values[q]
        estimate = adv_results[q]
        
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
        ''', (dataset, 'ADV', q, ground_truth, estimate, rel_error, 0.0, T_samples, epsilon, p))
        
        print(f"  ADV Q={q}: GT={ground_truth:.2e}, Est={estimate:.2e}, Rel Error={rel_error:.3f}")
    
    # Insert ADV+ results
    print(f"\nAdding {dataset} ADV+ results...")
    for q in range(4, 11):
        ground_truth = ground_truth_values[q]
        estimate = adv_plus_results[q]
        
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
        ''', (dataset, 'ADV+', q, ground_truth, estimate, rel_error, 0.0, T_samples, epsilon, p))
        
        print(f"  ADV+ Q={q}: GT={ground_truth:.2e}, Est={estimate:.2e}, Rel Error={rel_error:.3f}")
    
    conn.commit()
    conn.close()
    print(f"\nSuccessfully added {dataset} results to database!")

if __name__ == "__main__":
    add_filmtrust_results()
