#!/usr/bin/env python3
"""
Update csbwiki dataset results in the algorithm_comparison_p3 table with new data from testing-csbwiki-no-sampling.txt
This includes ADV, ADV+, and ADV++ results for p=3 with T=-1 (all triplets used).
"""

import sqlite3

def update_csbwiki_p3_results():
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
    
    # Ground truth values from the testing file
    ground_truth_values = {
        4: 3.520414e+13,
        5: 7.532651e+15,
        6: 1.487379e+18,
        7: 2.677936e+20,
        8: 4.398261e+22,
        9: 6.621362e+24,
        10: 9.191803e+26
    }
    
    # ADV algorithm estimates and relative errors (from testing file)
    adv_results = [
        (4, 6.788692e+13, 0.928),  # (q, estimate, relative_error)
        (5, 7.661287e+15, 0.017),
        (6, 1.518377e+18, 0.021),
        (7, 2.737018e+20, 0.022),
        (8, 4.395603e+22, 0.0006),
        (9, 6.384268e+24, 0.036),
        (10, 8.436928e+26, 0.082)
    ]
    
    # ADV+ algorithm estimates and relative errors (from testing file)
    adv_plus_results = [
        (4, 6.491104e+13, 0.844),  # (q, estimate, relative_error)
        (5, 6.980601e+15, 0.073),
        (6, 1.346859e+18, 0.094),
        (7, 2.362330e+20, 0.118),
        (8, 3.684946e+22, 0.162),
        (9, 5.203504e+24, 0.214),
        (10, 6.701875e+26, 0.271)
    ]
    
    # ADV++ algorithm estimates and relative errors (from testing file)
    adv_plus_plus_results = [
        (4, 3.764505e+13, 0.069),  # (q, estimate, relative_error)
        (5, 6.296393e+15, 0.164),
        (6, 1.177493e+18, 0.208),
        (7, 2.007289e+20, 0.250),
        (8, 3.126076e+22, 0.289),
        (9, 4.473512e+24, 0.324),
        (10, 5.917465e+26, 0.356)
    ]
    
    dataset = "csbwiki"
    epsilon = 1.0
    p = 3
    T_samples = -1  # All triplets used (no sampling)
    std_relative_error = 0.0  # Single run, no standard deviation
    
    # Insert ADV results
    print("Updating ADV results...")
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
            print(f"Error updating ADV Q={q}: {e}")
    
    # Insert ADV+ results
    print("\nUpdating ADV+ results...")
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
            print(f"Error updating ADV+ Q={q}: {e}")
    
    # Insert ADV++ results
    print("\nUpdating ADV++ results...")
    for q, estimate, relative_error in adv_plus_plus_results:
        ground_truth = ground_truth_values[q]
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset, 'ADV++', q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p))
            print(f"  Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3f}")
        except Exception as e:
            print(f"Error updating ADV++ Q={q}: {e}")
    
    # Commit changes
    conn.commit()
    
    # Verify data
    print(f"\nVerifying updated data for {dataset}...")
    cursor.execute('''
        SELECT algorithm, q, ground_truth, estimate, relative_error 
        FROM algorithm_comparison_p3 
        WHERE dataset = ? AND p = ? AND epsilon = ? AND T_samples = ?
        ORDER BY algorithm, q
    ''', (dataset, p, epsilon, T_samples))
    
    results = cursor.fetchall()
    for algorithm, q, ground_truth, estimate, rel_error in results:
        print(f"  {algorithm} Q={q}: GT={ground_truth:.2e}, Est={estimate:.2e}, Rel Error={rel_error:.3f}")
    
    # Show summary statistics
    print(f"\nSummary for {dataset} (p={p}, ε={epsilon}, T={T_samples}):")
    cursor.execute('''
        SELECT algorithm, 
               COUNT(*) as num_q_values,
               AVG(relative_error) as avg_rel_error,
               MIN(relative_error) as min_rel_error,
               MAX(relative_error) as max_rel_error
        FROM algorithm_comparison_p3 
        WHERE dataset = ? AND p = ? AND epsilon = ? AND T_samples = ?
        GROUP BY algorithm
        ORDER BY algorithm
    ''', (dataset, p, epsilon, T_samples))
    
    summary_results = cursor.fetchall()
    for algorithm, num_q, avg_error, min_error, max_error in summary_results:
        print(f"  {algorithm}: {num_q} Q-values, Avg Rel Error: {avg_error:.3f}, Range: [{min_error:.3f}, {max_error:.3f}]")
    
    conn.close()
    print(f"\nSuccessfully updated {dataset} results in database!")

if __name__ == "__main__":
    update_csbwiki_p3_results()
