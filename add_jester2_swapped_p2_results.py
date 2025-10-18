#!/usr/bin/env python3
"""
Script to add jester2-swapped p=2 results to the database.
"""

import sqlite3

def add_p2_results():
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # Ground truth values for jester2-swapped p=2, q=4-10
    ground_truth_values = [
        8.868889e+18,  # q=4
        8.248e+22,     # q=5
        6.744e+26,     # q=6
        4.821e+30,     # q=7
        3.039e+34,     # q=8
        1.707e+38,     # q=9
        8.645e+41      # q=10
    ]
    
    # Naive results (mean estimates from the experiment)
    naive_estimates = [
        1.461967e+18,  # q=4
        5.745043e+21,  # q=5
        2.114416e+25,  # q=6
        7.220588e+28,  # q=7
        2.268215e+32,  # q=8
        6.528326e+35,  # q=9
        1.721855e+39   # q=10
    ]
    
    # ADV results (mean estimates from the experiment)
    adv_estimates = [
        8.888541e+18,  # q=4
        8.268102e+22,  # q=5
        6.763138e+26,  # q=6
        4.836707e+30,  # q=7
        3.049681e+34,  # q=8
        1.714453e+38,  # q=9
        8.685348e+41   # q=10
    ]
    
    # ADV+ results (mean estimates from the experiment)
    adv_plus_estimates = [
        8.904776e+18,  # q=4
        8.293455e+22,  # q=5
        6.790651e+26,  # q=6
        4.860529e+30,  # q=7
        3.067074e+34,  # q=8
        1.725489e+38,  # q=9
        8.747418e+41   # q=10
    ]
    
    # ADV++ results (mean estimates from the experiment)
    adv_plus_plus_estimates = [
        8.895011e+18,  # q=4
        8.280800e+22,  # q=5
        6.777448e+26,  # q=6
        4.849124e+30,  # q=7
        3.058682e+34,  # q=8
        1.720106e+38,  # q=9
        8.716795e+41   # q=10
    ]
    
    # Relative errors from the experiment
    naive_errors = [0.8352, 0.9303, 0.9686, 0.9850, 0.9925, 0.9962, 0.9980]
    adv_errors = [0.00979, 0.01268, 0.01534, 0.01790, 0.02043, 0.02295, 0.02545]
    adv_plus_errors = [0.00957, 0.01320, 0.01647, 0.01951, 0.02244, 0.02531, 0.02816]
    adv_plus_plus_errors = [0.00652, 0.00848, 0.01037, 0.01220, 0.01399, 0.01577, 0.01754]
    
    # Insert results for each algorithm
    algorithms = [
        ("Naive", naive_estimates, naive_errors),
        ("ADV", adv_estimates, adv_errors),
        ("ADV+", adv_plus_estimates, adv_plus_errors),
        ("ADV++", adv_plus_plus_estimates, adv_plus_plus_errors)
    ]
    
    for algorithm_name, estimates, errors in algorithms:
        for i, q in enumerate(range(4, 11)):
            try:
                cursor.execute('''
                    INSERT OR REPLACE INTO algorithm_comparison_p3 
                    (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    'jester2-swapped',
                    algorithm_name,
                    q,
                    ground_truth_values[i],
                    estimates[i],
                    errors[i],
                    0.0,  # std_relative_error (single run)
                    -1,   # T_samples (no sampling)
                    1.0,  # epsilon
                    2     # p=2
                ))
                print(f"Inserted jester2-swapped {algorithm_name} p=2 q={q}: rel_err={errors[i]:.4f}")
            except Exception as e:
                print(f"Error inserting {algorithm_name} q={q}: {e}")
    
    # Commit changes
    conn.commit()
    conn.close()
    print("Successfully added jester2-swapped p=2 results to database")

if __name__ == "__main__":
    add_p2_results()
