#!/usr/bin/env python3

import sqlite3
import re
import os
import glob

def parse_naive_results(file_path):
    """Parse naive algorithm results from a _new.txt file"""
    results = []
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Extract dataset name from file path
    dataset = os.path.basename(file_path).replace('p3_', '').replace('_new.txt', '')
    
    # Find all ALG=0 sections
    sections = re.findall(r'DATASET: (\w+)\s+ALG=0\s+Q=(\d+)(.*?)(?=DATASET:|$)', content, re.DOTALL)
    
    for section in sections:
        dataset_name, q_value, section_content = section
        
        # Extract ground truth
        real_match = re.search(r'real count = ([\d\.e\+\-]+)', section_content)
        if not real_match:
            continue
        ground_truth = float(real_match.group(1))
        
        # Extract estimate (mean)
        mean_match = re.search(r'# Mean = ([\d\.e\+\-]+)', section_content)
        if not mean_match:
            continue
        estimate = float(mean_match.group(1))
        
        # Extract relative error
        rel_err_match = re.search(r'adv rel err = ([\d\.e\+\-]+)', section_content)
        if not rel_err_match:
            continue
        relative_error = float(rel_err_match.group(1))
        
        results.append({
            'dataset': dataset,
            'epsilon': 1.0,
            'algorithm_name': 'Naive',
            'q_value': int(q_value),
            'ground_truth': ground_truth,
            'estimate': estimate,
            'relative_error': relative_error,
            'num_rounds': 20,  # Based on the vertex sampling approach mentioned
            'T_samples': None  # Naive algorithm doesn't use T_samples
        })
    
    return results

def load_to_database(results, db_path):
    """Load results into the algorithm_comparison_p3 table"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    for result in results:
        try:
            cursor.execute("""
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, epsilon, algorithm_name, q_value, ground_truth, estimate, relative_error, num_rounds, T_samples)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                result['dataset'],
                result['epsilon'],
                result['algorithm_name'],
                result['q_value'],
                result['ground_truth'],
                result['estimate'],
                result['relative_error'],
                result['num_rounds'],
                result['T_samples']
            ))
        except Exception as e:
            print(f"Error inserting {result}: {e}")
    
    conn.commit()
    conn.close()

def main():
    # Find all _new.txt files
    new_files = glob.glob('/data/yizhangh/ldp-pq/larger_datasets_batch_results_FIXED/p3_*_new.txt')
    
    print(f"Found {len(new_files)} _new.txt files:")
    for file in new_files:
        print(f"  - {os.path.basename(file)}")
    
    all_results = []
    
    for file_path in new_files:
        print(f"\nProcessing {os.path.basename(file_path)}...")
        results = parse_naive_results(file_path)
        print(f"  Found {len(results)} naive algorithm results")
        
        for result in results:
            print(f"    Q={result['q_value']}: GT={result['ground_truth']:.2e}, Est={result['estimate']:.2e}, RelErr={result['relative_error']:.2f}%")
        
        all_results.extend(results)
    
    print(f"\nTotal results to load: {len(all_results)}")
    
    # Load to database
    db_path = '/data/yizhangh/ldp-pq/../biclq_counts.db'
    load_to_database(all_results, db_path)
    
    print(f"\nSuccessfully loaded {len(all_results)} naive algorithm results to algorithm_comparison_p3 table")

if __name__ == "__main__":
    main()
