#!/usr/bin/env python3
"""
Add librec-filmtrust-ratings ADV++ results for p=3, q=4-10 to the database.
"""

import sqlite3

def add_filmtrust_adv_plus_plus_results():
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # Ground truth values
    ground_truth_values = {
        4: 1.110059e+12,
        5: 8.618833e+12,
        6: 5.800279e+13,
        7: 3.381429e+14,
        8: 1.719092e+15,
        9: 7.683863e+15,
        10: 3.043229e+16
    }
    
    # ADV++ results from the comparison table
    adv_plus_plus_results = {
        4: 5.741829e+12,   # Mean estimate
        5: 2.942153e+12,   # Mean estimate
        6: -3.200351e+14,  # Mean estimate (negative)
        7: 1.058325e+14,   # Mean estimate
        8: 1.355410e+16,   # Mean estimate
        9: -1.484289e+16,  # Mean estimate (negative)
        10: -4.466144e+17  # Mean estimate (negative)
    }
    
    # Relative errors from the comparison table
    relative_errors = {
        4: 4.173e+00,      # 4.17
        5: 6.586e-01,      # 0.66
        6: 6.518e+00,      # 6.52
        7: 6.870e-01,      # 0.69
        8: 6.884e+00,      # 6.88
        9: 2.932e+00,      # 2.93
        10: 1.568e+01      # 15.68
    }
    
    dataset = "librec-filmtrust-ratings"
    epsilon = 1.0
    p = 3
    T_samples = -1  # No sampling for ADV++
    
    # Insert ADV++ results for q=4-10
    print(f"Adding {dataset} ADV++ results for p=3, q=4-10...")
    for q in [4, 5, 6, 7, 8, 9, 10]:
        ground_truth = ground_truth_values[q]
        estimate = adv_plus_plus_results[q]
        rel_error = relative_errors[q]
        
        # Insert into database
        cursor.execute('''
            INSERT OR REPLACE INTO algorithm_comparison_p3 
            (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (dataset, 'ADV++', q, ground_truth, estimate, rel_error, 0.0, T_samples, epsilon, p))
        
        print(f"  ADV++ Q={q}: GT={ground_truth:.2e}, Est={estimate:.2e}, Rel Error={rel_error:.3e}")
    
    conn.commit()
    conn.close()
    print(f"\nSuccessfully added {dataset} ADV++ results for p=3, q=4-10!")

if __name__ == "__main__":
    add_filmtrust_adv_plus_plus_results()
