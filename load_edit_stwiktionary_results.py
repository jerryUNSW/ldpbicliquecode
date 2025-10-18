#!/usr/bin/env python3
"""
Load ADV, ADV+, ADV++ results for edit-stwiktionary dataset into the algorithm_comparison_p3 table.
"""

import sqlite3
import re

def load_edit_stwiktionary_results():
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # Create table if it doesn't exist
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS algorithm_comparison_p3 (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            dataset TEXT NOT NULL,
            algorithm TEXT NOT NULL,
            q INTEGER NOT NULL,
            ground_truth REAL NOT NULL,
            estimate REAL NOT NULL,
            relative_error REAL NOT NULL,
            std_relative_error REAL NOT NULL,
            T_samples INTEGER DEFAULT -1,
            epsilon REAL DEFAULT 1.0,
            p INTEGER DEFAULT 3,
            UNIQUE(dataset, algorithm, q, T_samples, epsilon, p)
        )
    ''')
    
    # ADV results
    adv_results = [
        (4, 2.399290e+09, 53988615005.89, 2.150e+01, 0.000e+00),
        (5, 1.072e+11, -129811073435.20, 2.211e+00, 0.000e+00),
        (6, 4.214e+12, -15668343642074.99, 4.718e+00, 0.000e+00),
        (7, 1.464e+14, 36110876344621.77, 7.533e-01, 0.000e+00),
        (8, 4.539e+15, 4342437993995599.38, 4.331e-02, 0.000e+00),
        (9, 1.270e+17, -7965261049618416.74, 1.063e+00, 0.000e+00),
        (10, 3.237e+18, -854541237284586556.25, 1.264e+00, 0.000e+00)
    ]
    
    # ADV+ results
    adv_plus_results = [
        (4, 2.399e+09, 52349309064.78, 2.082e+01, 0.000e+00),
        (5, 1.072e+11, -108904519488.07, 2.016e+00, 0.000e+00),
        (6, 4.214e+12, -14614388644791.15, 4.468e+00, 0.000e+00),
        (7, 1.464e+14, 42567104070908.53, 7.092e-01, 0.000e+00),
        (8, 4.539e+15, 4091340478569225.31, 9.863e-02, 0.000e+00),
        (9, 1.270e+17, -17834753843471278.54, 1.140e+00, 0.000e+00),
        (10, 3.237e+18, -1001301793099988431.25, 1.309e+00, 0.000e+00)
    ]
    
    # ADV++ results
    adv_plus_plus_results = [
        (4, 2.399e+09, 16489831131.21, 5.873e+00, 0.000e+00),
        (5, 1.072e+11, 23870470113.80, 7.773e-01, 0.000e+00),
        (6, 4.214e+12, 826820887468.73, 8.038e-01, 0.000e+00),
        (7, 1.464e+14, 64723948822609.38, 5.579e-01, 0.000e+00),
        (8, 4.539e+15, 1814331233632382.61, 6.003e-01, 0.000e+00),
        (9, 1.270e+17, 42795593377455772.52, 6.630e-01, 0.000e+00),
        (10, 3.237e+18, 923789891217098126.06, 7.146e-01, 0.000e+00)
    ]
    
    dataset = "edit-stwiktionary"
    epsilon = 1.0
    p = 3
    T_samples = -1  # All triplets used
    
    # Insert ADV results
    print("Inserting ADV results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in adv_results:
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'ADV', q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p))
            print(f"  Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3e}")
        except Exception as e:
            print(f"Error inserting ADV Q={q}: {e}")
    
    # Insert ADV+ results
    print("\nInserting ADV+ results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in adv_plus_results:
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'ADV+', q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p))
            print(f"  Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3e}")
        except Exception as e:
            print(f"Error inserting ADV+ Q={q}: {e}")
    
    # Insert ADV++ results
    print("\nInserting ADV++ results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in adv_plus_plus_results:
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'ADV++', q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p))
            print(f"  Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3e}")
        except Exception as e:
            print(f"Error inserting ADV++ Q={q}: {e}")
    
    # Commit changes
    conn.commit()
    
    # Verify insertion
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
    print(f"\nSuccessfully loaded {dataset} results into database!")

if __name__ == "__main__":
    load_edit_stwiktionary_results()

