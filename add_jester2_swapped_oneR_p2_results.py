#!/usr/bin/env python3
"""
Script to add jester2-swapped OneR results to the p=2 database.
"""

import sqlite3

def add_oneR_results():
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # OneR results for jester2-swapped (p=2, q=4-10)
    oneR_results = [
        (4, 8.869e+18, 8.824e+18, 0.0050),
        (5, 8.248e+22, 8.197e+22, 0.0061),
        (6, 6.744e+26, 6.822e+26, 0.0115),
        (7, 4.821e+30, 4.878e+30, 0.0118),
        (8, 3.039e+34, 2.980e+34, 0.0194),
        (9, 1.707e+38, 1.739e+38, 0.0184),
        (10, 8.645e+41, 8.367e+41, 0.0322)
    ]
    
    # Insert OneR results
    print("Inserting OneR results for jester2-swapped p=2...")
    for q, ground_truth, estimate, relative_error in oneR_results:
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                'jester2-swapped',
                'OneR',
                q,
                ground_truth,
                estimate,
                relative_error,
                0.0,   # std_relative_error (not available from single run)
                -1,    # T_samples (no sampling)
                1.0,   # epsilon
                2      # p
            ))
            print(f"Inserted jester2-swapped OneR q={q}: rel_err={relative_error:.4f}")
        except Exception as e:
            print(f"Error inserting OneR q={q}: {e}")
    
    # Commit changes
    conn.commit()
    conn.close()
    print("\nSuccessfully added jester2-swapped OneR results to p=3 database")

if __name__ == "__main__":
    add_oneR_results()
