#!/usr/bin/env python3
"""
Script to add jester2-swapped naive results to the database.
"""

import sqlite3

def add_naive_results():
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # Naive results for jester2-swapped (p=3, q=4-10)
    results = [
        (4, 2.132051e+19, 1.578799e+18, 0.9259493),  # q=4
        (5, 1.756765e+23, 3.622652e+21, 0.9793789),  # q=5
        (6, 1.379549e+27, 8.239859e+24, 0.9940271), # q=6
        (7, 9.725885e+30, 1.871544e+28, 0.9980757), # q=7
        (8, 6.098389e+34, 3.991373e+31, 0.9993455), # q=8
        (9, 3.419788e+38, 8.078718e+34, 0.9997638), # q=9
        (10, 1.730112e+42, 1.538356e+38, 0.9999111) # q=10
    ]
    
    for q, ground_truth, estimate, relative_error in results:
        try:
            # Insert naive result
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                'jester2-swapped',
                'Naive',
                q,
                ground_truth,
                estimate,
                relative_error,
                0.0,  # std_relative_error (single run)
                -1,   # T_samples (no sampling)
                1.0,  # epsilon
                3     # p
            ))
            print(f"Inserted jester2-swapped Naive q={q}: ground_truth={ground_truth:.2e}, estimate={estimate:.2e}, rel_err={relative_error:.6f}")
            
        except Exception as e:
            print(f"Error inserting q={q}: {e}")
    
    # Commit changes
    conn.commit()
    conn.close()
    print("Successfully added jester2-swapped naive results to database")

if __name__ == "__main__":
    add_naive_results()
