#!/usr/bin/env python3
"""
Update librec-filmtrust-ratings Naive results for p=3, q=4-10 with new values from testing.
"""

import sqlite3

def update_filmtrust_naive_p3_all_q():
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
    
    # New Naive results from re-computation (10-round averages)
    new_naive_results = {
        4: 1.092668e+17,  # Average from 10 rounds
        5: 5.345897e+18,  # Average from 10 rounds
        6: 2.386802e+20,  # Average from 10 rounds
        7: 8.817101e+21,  # Average from 10 rounds
        8: 2.977100e+23,  # Average from 10 rounds
        9: 8.715427e+24,  # Average from 10 rounds
        10: 2.219288e+26  # Average from 10 rounds
    }
    
    dataset = "librec-filmtrust-ratings"
    epsilon = 1.0
    p = 3
    T_samples = 10  # 10 samples per iteration (vertex sampling)
    
    # Update Naive results for q=4-10
    print(f"Updating {dataset} Naive results for p=3, q=4-10...")
    for q in [4, 5, 6, 7, 8, 9, 10]:
        ground_truth = ground_truth_values[q]
        estimate = new_naive_results[q]
        
        # Calculate relative error
        if ground_truth > 0:
            rel_error = abs(estimate - ground_truth) / ground_truth
        else:
            rel_error = float('inf')
        
        # Update in database
        cursor.execute('''
            UPDATE algorithm_comparison_p3 
            SET estimate = ?, relative_error = ?
            WHERE dataset = ? AND algorithm = ? AND q = ? AND epsilon = ? AND p = ?
        ''', (estimate, rel_error, dataset, 'Naive', q, epsilon, p))
        
        print(f"  Naive Q={q}: GT={ground_truth:.2e}, Est={estimate:.2e}, Rel Error={rel_error:.3e}")
    
    conn.commit()
    conn.close()
    print(f"\nSuccessfully updated {dataset} Naive results for p=3, q=4-10!")

if __name__ == "__main__":
    update_filmtrust_naive_p3_all_q()
