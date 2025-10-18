#!/usr/bin/env python3
"""
Add csbwiki Naive results to the algorithm_comparison_p3 table.
"""

import sqlite3

def add_csbwiki_naive_results():
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # Ground truth values from the file
    ground_truth_values = {
        4: 3.520414e+13,
        5: 7.532651e+15,
        6: 1.487379e+18,
        7: 2.677936e+20,
        8: 4.398261e+22,
        9: 6.621362e+24,
        10: 9.191803e+26
    }
    
    # Naive algorithm estimates and relative errors from p3_csbwiki_2M.txt
    naive_results = [
        (4, 1.305684e+16, 369.89),  # (q, estimate, relative_error)
        (5, 4.531345e+17, 59.16),
        (6, 1.389676e+19, 8.34),
        (7, 3.796647e+20, 0.42),
        (8, 1.062635e+22, 0.76),
        (9, 3.077620e+23, 0.95),
        (10, 1.077190e+25, 0.99)
    ]
    
    dataset = "csbwiki"
    epsilon = 1.0
    p = 3
    T_samples = -1  # All triplets used
    std_relative_error = 0.0  # Single run, no standard deviation
    
    # Insert Naive results
    print("Inserting Naive results for csbwiki...")
    for q, estimate, relative_error in naive_results:
        ground_truth = ground_truth_values[q]
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'Naive', q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p))
            print(f"  Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.2f}")
        except Exception as e:
            print(f"Error inserting Naive Q={q}: {e}")
    
    # Commit changes
    conn.commit()
    
    # Verify data
    print(f"\nVerifying all csbwiki data...")
    cursor.execute('''
        SELECT algorithm, q, ground_truth, estimate, relative_error 
        FROM algorithm_comparison_p3 
        WHERE dataset = ? AND p = ? AND epsilon = ?
        ORDER BY algorithm, q
    ''', (dataset, p, epsilon))
    
    results = cursor.fetchall()
    for algorithm, q, ground_truth, estimate, rel_error in results:
        print(f"  {algorithm} Q={q}: GT={ground_truth:.2e}, Est={estimate:.2e}, Rel Error={rel_error:.2f}")
    
    conn.close()
    print(f"\nSuccessfully added csbwiki Naive results to database!")

if __name__ == "__main__":
    add_csbwiki_naive_results()
