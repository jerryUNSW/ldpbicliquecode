#!/usr/bin/env python3
"""
Load actual NIPS results from testing-nips-no-sampling.txt into the database.
"""

import sqlite3
import re

def parse_test_results():
    """Parse the actual results from testing-nips-no-sampling.txt"""
    
    # Ground truth values
    ground_truth = {
        4: 1.729841e+14,
        5: 2.026990e+15,
        6: 2.054109e+16,
        7: 1.848309e+17,
        8: 1.505491e+18,
        9: 1.126738e+19,
        10: 7.844682e+19
    }
    
    # ADV results from line 51
    adv_results = {
        4: 4.442019e+14,
        5: 8.447451e+15,
        6: 1.490508e+17,
        7: 2.456005e+18,
        8: 3.857563e+19,
        9: 5.732922e+20,
        10: 8.076765e+21
    }
    
    # ADV+ results from line 62
    adv_plus_results = {
        4: -2.413750e+14,  # Negative result!
        5: 1.220725e+16,
        6: 2.141426e+17,
        7: -3.416424e+18,  # Negative result!
        8: -9.742705e+19,  # Negative result!
        9: 5.527845e+20,
        10: 3.078710e+22
    }
    
    return ground_truth, adv_results, adv_plus_results

def compute_relative_error(estimated, true_value):
    """Compute relative error as |estimated - true_value| / |true_value| * 100"""
    return abs(estimated - true_value) / abs(true_value) * 100

def load_actual_results():
    """Load actual results into database"""
    
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # Ensure the table exists
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
    conn.commit()
    
    # Parse actual results
    ground_truth, adv_results, adv_plus_results = parse_test_results()
    
    dataset = "nips"
    epsilon = 1.0
    p = 3
    T_samples = -1  # All triplets used
    
    # Insert ADV results
    print("Inserting ACTUAL ADV results...")
    for q in sorted(ground_truth.keys()):
        gt = ground_truth[q]
        estimate = adv_results[q]
        rel_error = compute_relative_error(estimate, gt)
        
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'ADV', q, gt, estimate, rel_error, 0.0, T_samples, epsilon, p))
            print(f"  Q={q}: GT={gt:.2e}, Est={estimate:.2e}, Rel Error={rel_error:.2f}%")
        except Exception as e:
            print(f"Error inserting ADV Q={q}: {e}")
    
    # Insert ADV+ results
    print("\nInserting ACTUAL ADV+ results...")
    for q in sorted(ground_truth.keys()):
        gt = ground_truth[q]
        estimate = adv_plus_results[q]
        rel_error = compute_relative_error(estimate, gt)
        
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'ADV+', q, gt, estimate, rel_error, 0.0, T_samples, epsilon, p))
            print(f"  Q={q}: GT={gt:.2e}, Est={estimate:.2e}, Rel Error={rel_error:.2f}%")
        except Exception as e:
            print(f"Error inserting ADV+ Q={q}: {e}")
    
    # Commit changes
    conn.commit()
    
    # Verify data
    print(f"\nVerifying ACTUAL data for {dataset}...")
    cursor.execute('''
        SELECT algorithm, q, ground_truth, estimate, relative_error 
        FROM algorithm_comparison_p3 
        WHERE dataset = ? AND p = ? AND epsilon = ? AND algorithm IN ('ADV', 'ADV+')
        ORDER BY algorithm, q
    ''', (dataset, p, epsilon))
    
    results = cursor.fetchall()
    for algorithm, q, gt, estimate, rel_error in results:
        print(f"  {algorithm} Q={q}: GT={gt:.2e}, Est={estimate:.2e}, Rel Error={rel_error:.2f}%")
    
    conn.close()
    print(f"\nSuccessfully loaded ACTUAL {dataset} results into database!")

if __name__ == "__main__":
    load_actual_results()
