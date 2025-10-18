#!/usr/bin/env python3
"""
Update librec-filmtrust-ratings Naive results for q=8,9,10 with new values.
"""

import sqlite3

def update_filmtrust_naive_q8_10():
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # Ground truth values
    ground_truth_values = {
        8: 1.719092e+15,
        9: 7.683863e+15,
        10: 3.043229e+16
    }
    
    # New Naive results from re-computation
    new_naive_results = {
        8: 3.535875e+17,
        9: 1.217735e+18,
        10: 7.390881e+18
    }
    
    dataset = "librec-filmtrust-ratings"
    epsilon = 1.0
    p = 3
    T_samples = 10  # 10 samples per iteration (vertex sampling)
    
    # Update Naive results for q=8,9,10
    print(f"Updating {dataset} Naive results for q=8,9,10...")
    for q in [8, 9, 10]:
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
        
        print(f"  Naive Q={q}: GT={ground_truth:.2e}, Est={estimate:.2e}, Rel Error={rel_error:.3f}")
    
    conn.commit()
    conn.close()
    print(f"\nSuccessfully updated {dataset} Naive results for q=8,9,10!")

if __name__ == "__main__":
    update_filmtrust_naive_q8_10()
