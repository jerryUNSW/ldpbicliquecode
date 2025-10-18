#!/usr/bin/env python3
"""
Load nips dataset ADV, ADV+, ADV++ results into the algorithm_comparison_p3 table.
"""

import sqlite3

def load_nips_results():
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
    
    # ADV results
    adv_results = [
        (4, 1.729841e+14, 6647543525249683.81, 3.743e+01, 0.000e+00),
        (5, 2.027e+15, 995038755842054561.88, 4.899e+02, 0.000e+00),
        (6, 2.054e+16, 147525980191338103680.00, 7.181e+03, 0.000e+00),
        (7, 1.848e+17, 25656960504668958218240.00, 1.388e+05, 0.000e+00),
        (8, 1.505e+18, 5024786113263302170116096.00, 3.338e+06, 0.000e+00),
        (9, 1.127e+19, 1002442740784710342036946944.00, 8.897e+07, 0.000e+00),
        (10, 7.845e+19, 191017389698057575864232574976.00, 2.435e+09, 0.000e+00)
    ]
    
    # ADV+ results
    adv_plus_results = [
        (4, 1.730e+14, 6951139104037008.83, 3.918e+01, 0.000e+00),
        (5, 2.027e+15, 1242035828274722241.50, 6.117e+02, 0.000e+00),
        (6, 2.054e+16, 260104870456352081440.00, 1.266e+04, 0.000e+00),
        (7, 1.848e+17, 66163295374015044558848.00, 3.580e+05, 0.000e+00),
        (8, 1.505e+18, 17407010157280706664857600.00, 1.156e+07, 0.000e+00),
        (9, 1.127e+19, 4324353206448457028066082816.00, 3.838e+08, 0.000e+00),
        (10, 7.845e+19, 988330874395095321549618544640.00, 1.260e+10, 0.000e+00)
    ]
    
    # ADV++ results
    adv_plus_plus_results = [
        (4, 1.730e+14, 753055543727766.94, 3.353e+00, 0.000e+00),
        (5, 2.027e+15, 70420084664150883.19, 3.374e+01, 0.000e+00),
        (6, 2.054e+16, 6280456114820272327.50, 3.048e+02, 0.000e+00),
        (7, 1.848e+17, 649065167740753498240.00, 3.511e+03, 0.000e+00),
        (8, 1.505e+18, 76276083922232982618112.00, 5.066e+04, 0.000e+00),
        (9, 1.127e+19, 9159495316735466761355264.00, 8.129e+05, 0.000e+00),
        (10, 7.845e+19, 1048029047164569731125477376.00, 1.336e+07, 0.000e+00)
    ]
    
    dataset = "nips"
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
    
    # Verify data
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
    load_nips_results()

