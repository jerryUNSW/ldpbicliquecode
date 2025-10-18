#!/usr/bin/env python3
"""
Update edit-stwiktionary ADV, ADV+, ADV++ results in the algorithm_comparison_p3 table.
"""

import sqlite3

def update_edit_stwiktionary_results():
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # New ADV results
    adv_results = [
        (4, 2.399290e+09, 51969522804.36, 2.066e+01, 6.596e-01),
        (5, 1.072e+11, -55027378173.21, 1.513e+00, 4.139e-01),
        (6, 4.214e+12, -13468068839458.76, 4.196e+00, 4.340e-01),
        (7, 1.464e+14, 45755145317133.32, 6.874e-01, 1.662e-01),
        (8, 4.539e+15, 3992561521027777.45, 1.384e-01, 5.981e-02),
        (9, 1.270e+17, -25123766835771459.03, 1.198e+00, 1.259e-01),
        (10, 3.237e+18, -1396058120188330220.00, 1.431e+00, 1.256e-01)
    ]
    
    # New ADV+ results
    adv_plus_results = [
        (4, 2.399e+09, 51768423716.31, 2.058e+01, 5.019e-01),
        (5, 1.072e+11, -62905886032.74, 1.587e+00, 3.599e-01),
        (6, 4.214e+12, -13514483546913.48, 4.207e+00, 3.419e-01),
        (7, 1.464e+14, 41883355234178.43, 7.285e-01, 2.331e-01),
        (8, 4.539e+15, 3761660786913025.54, 2.583e-01, 1.591e-01),
        (9, 1.270e+17, -26718160220866925.00, 1.210e+00, 2.829e-01),
        (10, 3.237e+18, -1121426102690657786.00, 1.346e+00, 1.767e-01)
    ]
    
    # New ADV++ results
    adv_plus_plus_results = [
        (4, 2.399e+09, 16896567708.97, 6.042e+00, 3.851e-01),
        (5, 1.072e+11, 47578785115.21, 6.129e-01, 2.309e-01),
        (6, 4.214e+12, 1509620845493.20, 7.027e-01, 2.363e-01),
        (7, 1.464e+14, 75711671080999.92, 5.825e-01, 1.853e-01),
        (8, 4.539e+15, 1821529543369854.17, 6.707e-01, 1.905e-01),
        (9, 1.270e+17, 35955597670649350.64, 7.505e-01, 2.349e-01),
        (10, 3.237e+18, 646686753194772483.56, 8.002e-01, 2.794e-01)
    ]
    
    dataset = "edit-stwiktionary"
    epsilon = 1.0
    p = 3
    T_samples = -1  # All triplets used
    
    # Update ADV results
    print("Updating ADV results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in adv_results:
        try:
            cursor.execute('''
                UPDATE algorithm_comparison_p3 
                SET ground_truth = ?, estimate = ?, relative_error = ?, std_relative_error = ?
                WHERE dataset = ? AND algorithm = ? AND q = ? AND p = ? AND epsilon = ?
            ''', (ground_truth, estimate, relative_error, std_relative_error, dataset, 'ADV', q, p, epsilon))
            print(f"  Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3e}")
        except Exception as e:
            print(f"Error updating ADV Q={q}: {e}")
    
    # Update ADV+ results
    print("\nUpdating ADV+ results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in adv_plus_results:
        try:
            cursor.execute('''
                UPDATE algorithm_comparison_p3 
                SET ground_truth = ?, estimate = ?, relative_error = ?, std_relative_error = ?
                WHERE dataset = ? AND algorithm = ? AND q = ? AND p = ? AND epsilon = ?
            ''', (ground_truth, estimate, relative_error, std_relative_error, dataset, 'ADV+', q, p, epsilon))
            print(f"  Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3e}")
        except Exception as e:
            print(f"Error updating ADV+ Q={q}: {e}")
    
    # Update ADV++ results
    print("\nUpdating ADV++ results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in adv_plus_plus_results:
        try:
            cursor.execute('''
                UPDATE algorithm_comparison_p3 
                SET ground_truth = ?, estimate = ?, relative_error = ?, std_relative_error = ?
                WHERE dataset = ? AND algorithm = ? AND q = ? AND p = ? AND epsilon = ?
            ''', (ground_truth, estimate, relative_error, std_relative_error, dataset, 'ADV++', q, p, epsilon))
            print(f"  Q={q}: Ground Truth={ground_truth:.2e}, Estimate={estimate:.2e}, Rel Error={relative_error:.3e}")
        except Exception as e:
            print(f"Error updating ADV++ Q={q}: {e}")
    
    # Commit changes
    conn.commit()
    
    # Verify updates
    print(f"\nVerifying updated data for {dataset}...")
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
    print(f"\nSuccessfully updated {dataset} results in database!")

if __name__ == "__main__":
    update_edit_stwiktionary_results()

