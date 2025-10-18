#!/usr/bin/env python3
"""
Load all p=3 results (lrcwiki and edit-stwiktionary) into the algorithm_comparison_p3 table.
"""

import sqlite3
import re

def load_all_p3_results():
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
    
    # Load lrcwiki results
    print("Loading lrcwiki results...")
    load_lrcwiki_results(cursor)
    
    # Load edit-stwiktionary results
    print("\nLoading edit-stwiktionary results...")
    load_edit_stwiktionary_results(cursor)
    
    # Commit changes
    conn.commit()
    
    # Verify insertion
    print(f"\nVerifying all data...")
    cursor.execute('''
        SELECT dataset, algorithm, COUNT(*) as count 
        FROM algorithm_comparison_p3 
        WHERE p = 3 AND epsilon = 1.0
        GROUP BY dataset, algorithm 
        ORDER BY dataset, algorithm
    ''')
    
    results = cursor.fetchall()
    for dataset, algorithm, count in results:
        print(f"  {dataset} {algorithm}: {count} entries")
    
    conn.close()
    print(f"\nSuccessfully loaded all p=3 results into database!")

def load_lrcwiki_results(cursor):
    """Load lrcwiki results from the provided data"""
    
    # Naive results (from the file)
    naive_results = [
        (4, 5.949427e+10, 1.267050e+14, 2.128701e+03, 0.0),
        (5, 1.261650e+13, 5.647905e+15, 4.466601e+02, 0.0),
        (6, 2.263535e+15, 4.768782e+17, 2.096785e+02, 0.0),
        (7, 3.491078e+17, 4.861602e+19, 1.382579e+02, 0.0),
        (8, 4.711424e+19, 1.817259e+20, 2.857134e+00, 0.0),
        (9, 5.648001e+21, 1.326744e+22, 1.349051e+00, 0.0),
        (10, 6.088413e+23, 1.273033e+25, 1.990910e+01, 0.0)
    ]
    
    # ADV results (from the file)
    adv_results = [
        (4, 5.949427e+10, 5.414138e+10, 6.924796e-01, 0.0),
        (5, 1.261650e+13, 1.339127e+13, 6.569245e-01, 0.0),
        (6, 2.263535e+15, 2.421166e+15, 6.685328e-01, 0.0),
        (7, 3.491078e+17, 4.795089e+17, 5.734291e-01, 0.0),
        (8, 4.711424e+19, 3.609108e+19, 7.659310e-01, 0.0),
        (9, 5.648001e+21, 5.198696e+21, 7.204184e-01, 0.0),
        (10, 6.088413e+23, 8.317748e+23, 5.661590e-01, 0.0)
    ]
    
    # ADV+ results (from the file)
    adv_plus_results = [
        (4, 5.949427e+10, 8.044460e+10, 5.508349e-01, 0.0),
        (5, 1.261650e+13, 1.142170e+13, 6.997400e-01, 0.0),
        (6, 2.263535e+15, 1.726937e+15, 7.616923e-01, 0.0),
        (7, 3.491078e+17, 3.243931e+17, 7.288061e-01, 0.0),
        (8, 4.711424e+19, 5.037519e+19, 6.691541e-01, 0.0),
        (9, 5.648001e+21, 6.878506e+21, 6.178577e-01, 0.0),
        (10, 6.088413e+23, 5.532886e+23, 7.087502e-01, 0.0)
    ]
    
    # ADV++ results (from the file)
    adv_plus_plus_results = [
        (4, 5.949427e+10, 6.316425e+10, 6.456784e-01, 0.0),
        (5, 1.261650e+13, 1.148664e+13, 7.039371e-01, 0.0),
        (6, 2.263535e+15, 2.050630e+15, 7.043396e-01, 0.0),
        (7, 3.491078e+17, 3.201092e+17, 7.166136e-01, 0.0),
        (8, 4.711424e+19, 2.884199e+19, 8.120369e-01, 0.0),
        (9, 5.648001e+21, 6.002237e+21, 6.627011e-01, 0.0),
        (10, 6.088413e+23, 7.427150e+23, 6.198808e-01, 0.0)
    ]
    
    dataset = "lrcwiki"
    epsilon = 1.0
    p = 3
    T_samples = -1  # All triplets used
    
    # Insert Naive results
    print("  Inserting Naive results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in naive_results:
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'Naive', q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p))
            print(f"    Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3e}")
        except Exception as e:
            print(f"Error inserting Naive Q={q}: {e}")
    
    # Insert ADV results
    print("  Inserting ADV results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in adv_results:
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'ADV', q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p))
            print(f"    Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3e}")
        except Exception as e:
            print(f"Error inserting ADV Q={q}: {e}")
    
    # Insert ADV+ results
    print("  Inserting ADV+ results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in adv_plus_results:
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'ADV+', q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p))
            print(f"    Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3e}")
        except Exception as e:
            print(f"Error inserting ADV+ Q={q}: {e}")
    
    # Insert ADV++ results
    print("  Inserting ADV++ results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in adv_plus_plus_results:
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'ADV++', q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p))
            print(f"    Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3e}")
        except Exception as e:
            print(f"Error inserting ADV++ Q={q}: {e}")

def load_edit_stwiktionary_results(cursor):
    """Load edit-stwiktionary results"""
    
    # Naive results (from the user's output)
    naive_results = [
        (4, 2.399290e+09, 2.637855e+11, 1.089432e+02, 0.0),
        (5, 1.071954e+11, 3.526413e+12, 3.189705e+01, 0.0),
        (6, 4.214460e+12, 5.439186e+13, 1.190601e+01, 0.0),
        (7, 1.463854e+14, 8.564192e+14, 4.850442e+00, 0.0),
        (8, 4.539023e+15, 1.060288e+16, 1.335939e+00, 0.0),
        (9, 1.269991e+17, 2.239479e+17, 7.633822e-01, 0.0),
        (10, 3.236539e+18, 3.775316e+18, 1.664670e-01, 0.0)
    ]
    
    # ADV results
    adv_results = [
        (4, 2.399290e+09, 5.398862e+10, 2.150e+01, 0.000e+00),
        (5, 1.072e+11, -1.298111e+11, 2.211e+00, 0.000e+00),
        (6, 4.214e+12, -1.566834e+13, 4.718e+00, 0.000e+00),
        (7, 1.464e+14, 3.611088e+13, 7.533e-01, 0.000e+00),
        (8, 4.539e+15, 4.342438e+15, 4.331e-02, 0.000e+00),
        (9, 1.270e+17, -7.965261e+15, 1.063e+00, 0.000e+00),
        (10, 3.237e+18, -8.545412e+17, 1.264e+00, 0.000e+00)
    ]
    
    # ADV+ results
    adv_plus_results = [
        (4, 2.399e+09, 5.234931e+10, 2.082e+01, 0.000e+00),
        (5, 1.072e+11, -1.089045e+11, 2.016e+00, 0.000e+00),
        (6, 4.214e+12, -1.461439e+13, 4.468e+00, 0.000e+00),
        (7, 1.464e+14, 4.256710e+13, 7.092e-01, 0.000e+00),
        (8, 4.539e+15, 4.091340e+15, 9.863e-02, 0.000e+00),
        (9, 1.270e+17, -1.783475e+16, 1.140e+00, 0.000e+00),
        (10, 3.237e+18, -1.001302e+18, 1.309e+00, 0.000e+00)
    ]
    
    # ADV++ results
    adv_plus_plus_results = [
        (4, 2.399e+09, 1.648983e+10, 5.873e+00, 0.000e+00),
        (5, 1.072e+11, 2.387047e+10, 7.773e-01, 0.000e+00),
        (6, 4.214e+12, 8.268209e+11, 8.038e-01, 0.000e+00),
        (7, 1.464e+14, 6.472395e+13, 5.579e-01, 0.000e+00),
        (8, 4.539e+15, 1.814331e+15, 6.003e-01, 0.000e+00),
        (9, 1.270e+17, 4.279559e+16, 6.630e-01, 0.000e+00),
        (10, 3.237e+18, 9.237899e+17, 7.146e-01, 0.000e+00)
    ]
    
    dataset = "edit-stwiktionary"
    epsilon = 1.0
    p = 3
    T_samples = -1  # All triplets used
    
    # Insert Naive results
    print("  Inserting Naive results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in naive_results:
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'Naive', q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p))
            print(f"    Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3e}")
        except Exception as e:
            print(f"Error inserting Naive Q={q}: {e}")
    
    # Insert ADV results
    print("  Inserting ADV results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in adv_results:
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'ADV', q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p))
            print(f"    Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3e}")
        except Exception as e:
            print(f"Error inserting ADV Q={q}: {e}")
    
    # Insert ADV+ results
    print("  Inserting ADV+ results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in adv_plus_results:
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'ADV+', q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p))
            print(f"    Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3e}")
        except Exception as e:
            print(f"Error inserting ADV+ Q={q}: {e}")
    
    # Insert ADV++ results
    print("  Inserting ADV++ results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in adv_plus_plus_results:
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'ADV++', q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p))
            print(f"    Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3e}")
        except Exception as e:
            print(f"Error inserting ADV++ Q={q}: {e}")

if __name__ == "__main__":
    load_all_p3_results()
