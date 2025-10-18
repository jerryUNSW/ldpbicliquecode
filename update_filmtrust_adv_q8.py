#!/usr/bin/env python3
"""
Update librec-filmtrust-ratings ADV result for q=8 with new value.
"""

import sqlite3

def update_filmtrust_adv_q8():
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # Ground truth value
    ground_truth = 1.719092e+15
    
    # New ADV result from re-computation
    new_adv_estimate = 4.476957e+17
    new_adv_rel_error = 259.4257382893
    
    dataset = "librec-filmtrust-ratings"
    epsilon = 1.0
    p = 3
    q = 8
    T_samples = 2000000  # Rejection sampling with T=2M
    
    # Update ADV result for q=8
    print(f"Updating {dataset} ADV result for q=8...")
    
    cursor.execute('''
        UPDATE algorithm_comparison_p3 
        SET estimate = ?, relative_error = ?, T_samples = ?
        WHERE dataset = ? AND algorithm = ? AND q = ? AND epsilon = ? AND p = ?
    ''', (new_adv_estimate, new_adv_rel_error, T_samples, dataset, 'ADV', q, epsilon, p))
    
    print(f"  ADV Q={q}: GT={ground_truth:.2e}, Est={new_adv_estimate:.2e}, Rel Error={new_adv_rel_error:.3f}")
    
    conn.commit()
    conn.close()
    print(f"\nSuccessfully updated {dataset} ADV result for q=8!")

if __name__ == "__main__":
    update_filmtrust_adv_q8()
