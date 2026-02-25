#!/usr/bin/env python3
"""
Script to add jester2-swapped p=3 eps2 results to the database.
Loads data from jester-p-3-exp-eps2 directory and computes statistics.
"""

import sqlite3
import re
import os
from pathlib import Path
import numpy as np

def extract_data_from_file(filepath):
    """Extract key metrics from a single result file."""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        
        # Extract Q value from filename
        q_match = re.search(r'q(\d+)_', filepath)
        q = int(q_match.group(1)) if q_match else None
        
        # Extract algorithm from filename
        if 'naive' in filepath:
            algorithm = 'Naive'
        elif 'advplusplus' in filepath:
            algorithm = 'ADV++'
        elif 'advplus' in filepath:
            algorithm = 'ADV+'
        elif 'adv' in filepath:
            algorithm = 'ADV'
        else:
            algorithm = 'Unknown'
        
        # Extract ground truth (appears once per file)
        ground_truth_match = re.search(r'real count = ([\d\.e\+\-]+)', content)
        if not ground_truth_match:
            # Try alternative pattern
            ground_truth_match = re.search(r'biclique count = ([\d\.e\+\-]+)', content)
        ground_truth = float(ground_truth_match.group(1)) if ground_truth_match else None
        
        # Extract all relative errors (may be multiple runs)
        rel_errors = []
        estimates = []
        
        # Pattern 1: "adv rel err = ..." (for naive, appears once)
        adv_rel_err_match = re.search(r'adv rel err = ([\d\.e\-\+]+)', content)
        if adv_rel_err_match:
            rel_errors.append(float(adv_rel_err_match.group(1)))
        
        # Pattern 2: "relative error = ..." (for ADV/ADV+/ADV++, may appear multiple times)
        rel_err_matches = re.findall(r'relative error = ([\d\.e\-\+]+)', content)
        for match in rel_err_matches:
            rel_errors.append(float(match))
        
        # Extract estimates
        # Pattern 1: "estimate = ..." (may appear multiple times)
        estimate_matches = re.findall(r'estimate = ([\d\.e\+\-]+)', content)
        for match in estimate_matches:
            estimates.append(float(match))
        
        # Pattern 2: "# Mean = ..." (for naive)
        mean_match = re.search(r'# Mean = ([\d\.e\+\-]+)', content)
        if mean_match:
            estimates.append(float(mean_match.group(1)))
        
        # Compute statistics
        if rel_errors:
            mean_rel_error = np.mean(rel_errors)
            std_rel_error = np.std(rel_errors) if len(rel_errors) > 1 else 0.0
        else:
            mean_rel_error = None
            std_rel_error = None
        
        if estimates:
            mean_estimate = np.mean(estimates)
        else:
            mean_estimate = None
        
        return {
            'q': q,
            'algorithm': algorithm,
            'ground_truth': ground_truth,
            'estimate': mean_estimate,
            'relative_error': mean_rel_error,
            'std_relative_error': std_rel_error,
            'num_runs': len(rel_errors) if rel_errors else 0
        }
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return None

def load_all_data(data_dir_str="jester-p-3-exp-eps2"):
    """Load data from all result files."""
    data_dir = Path(data_dir_str)
    all_data = []
    
    if not data_dir.exists():
        print(f"Directory {data_dir} not found!")
        return []
    
    for filepath in data_dir.glob("*.txt"):
        data = extract_data_from_file(str(filepath))
        if data:
            all_data.append(data)
    
    return all_data

def add_p3_eps2_results():
    """Add jester2-swapped p=3 eps2 results to the database."""
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # Load data from directory
    print("Loading data from jester-p-3-exp-eps2 directory...")
    all_data = load_all_data("jester-p-3-exp-eps2")
    
    if not all_data:
        print("No data found!")
        conn.close()
        return
    
    print(f"Loaded {len(all_data)} data points")
    
    # Group by algorithm and q
    results_by_alg_q = {}
    for data in all_data:
        key = (data['algorithm'], data['q'])
        if key not in results_by_alg_q:
            results_by_alg_q[key] = []
        results_by_alg_q[key].append(data)
    
    # Insert results into database
    print("\nInserting results into database...")
    inserted_count = 0
    
    for (algorithm, q), data_list in sorted(results_by_alg_q.items()):
        # Use the first entry's ground truth (should be the same for all)
        ground_truth = data_list[0]['ground_truth']
        
        # Compute mean estimate and relative error across all entries
        estimates = [d['estimate'] for d in data_list if d['estimate'] is not None]
        rel_errors = [d['relative_error'] for d in data_list if d['relative_error'] is not None]
        
        if not rel_errors:
            print(f"Warning: No relative errors found for {algorithm} q={q}")
            continue
        
        mean_estimate = np.mean(estimates) if estimates else None
        mean_rel_error = np.mean(rel_errors)
        std_rel_error = np.std(rel_errors) if len(rel_errors) > 1 else 0.0
        
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                'jester2-swapped',
                algorithm,
                q,
                ground_truth,
                mean_estimate,
                mean_rel_error,
                std_rel_error,
                -1,   # T_samples (no sampling)
                2.0,  # epsilon = 2.0
                3     # p = 3
            ))
            print(f"Inserted jester2-swapped {algorithm} p=3 eps=2.0 q={q}: "
                  f"rel_err={mean_rel_error:.6f}, std={std_rel_error:.6f}, runs={len(rel_errors)}")
            inserted_count += 1
        except Exception as e:
            print(f"Error inserting {algorithm} q={q}: {e}")
    
    # Commit changes
    conn.commit()
    conn.close()
    print(f"\nSuccessfully added {inserted_count} jester2-swapped p=3 eps2 results to database")

if __name__ == "__main__":
    add_p3_eps2_results()



