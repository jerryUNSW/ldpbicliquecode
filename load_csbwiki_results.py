#!/usr/bin/env python3
"""
Load csbwiki dataset ADV and ADV+ results into the algorithm_comparison_p3 table.
"""

import sqlite3

def load_csbwiki_results():
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
    
    # ADV algorithm estimates and relative errors
    adv_results = [
        (4, 6.788692e+13, 0.928),  # (q, estimate, relative_error)
        (5, 7.661287e+15, 0.017),
        (6, 1.518377e+18, 0.021),
        (7, 2.737018e+20, 0.022),
        (8, 4.395603e+22, 0.0006),
        (9, 6.384268e+24, 0.036),
        (10, 8.436928e+26, 0.082)
    ]
    
    # ADV+ algorithm estimates and relative errors
    adv_plus_results = [
        (4, 6.491104e+13, 0.844),  # (q, estimate, relative_error)
        (5, 6.980601e+15, 0.073),
        (6, 1.346859e+18, 0.094),
        (7, 2.362330e+20, 0.118),
        (8, 3.684946e+22, 0.162),
        (9, 5.203504e+24, 0.214),
        (10, 6.701875e+26, 0.271)
    ]
    
    dataset = "csbwiki"
    epsilon = 1.0
    p = 3
    T_samples = -1  # All triplets used
    std_relative_error = 0.0  # Single run, no standard deviation
    
    # Insert ADV results
    print("Inserting ADV results...")
    for q, estimate, relative_error in adv_results:
        ground_truth = ground_truth_values[q]
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'ADV', q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p))
            print(f"  Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3f}")
        except Exception as e:
            print(f"Error inserting ADV Q={q}: {e}")
    
    # Insert ADV+ results
    print("\nInserting ADV+ results...")
    for q, estimate, relative_error in adv_plus_results:
        ground_truth = ground_truth_values[q]
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'ADV+', q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p))
            print(f"  Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3f}")
        except Exception as e:
            print(f"Error inserting ADV+ Q={q}: {e}")
    
    # Commit changes
    conn.commit()
    
    # Verify data
    print(f"\nVerifying data for {dataset}...")
    cursor.execute('''
        SELECT algorithm, q, ground_truth, estimate, relative_error 
        FROM algorithm_comparison_p3 
        WHERE dataset = ? AND p = ? AND epsilon = ?
        ORDER BY algorithm, q
    ''', (dataset, p, epsilon))
    
    results = cursor.fetchall()
    for algorithm, q, ground_truth, estimate, rel_error in results:
        print(f"  {algorithm} Q={q}: GT={ground_truth:.2e}, Est={estimate:.2e}, Rel Error={rel_error:.3f}")
    
    conn.close()
    print(f"\nSuccessfully loaded {dataset} results into database!")

if __name__ == "__main__":
    load_csbwiki_results()
