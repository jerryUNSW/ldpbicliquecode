#!/usr/bin/env python3
"""
Correctly update lrcwiki results from p3_lrcwiki_2M.txt file.
Extract the final "adv rel err" values for each algorithm and Q.
"""

import sqlite3
import re

def parse_lrcwiki_2M_file():
    """Parse the p3_lrcwiki_2M.txt file and extract the final relative error values"""
    
    file_path = "larger_datasets_batch_results_FIXED/p3_lrcwiki_2M.txt"
    
    # Dictionary to store results by algorithm and Q
    results = {
        'Naive': {},  # ALG=0
        'ADV': {},    # ALG=2
        'ADV+': {},   # ALG=3
        'ADV++': {}   # ALG=4
    }
    
    # Algorithm mapping
    alg_mapping = {0: 'Naive', 2: 'ADV', 3: 'ADV+', 4: 'ADV++'}
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Split by algorithm sections
    sections = re.split(r'DATASET: lrcwiki  ALG=(\d+)  Q=(\d+)', content)
    
    for i in range(1, len(sections), 3):
        alg_num = int(sections[i])
        q_value = int(sections[i+1])
        section_content = sections[i+2]
        
        if alg_num not in alg_mapping:
            continue
            
        algorithm = alg_mapping[alg_num]
        
        # Extract ground truth
        gt_match = re.search(r'biclique count = ([\d.e+-]+)', section_content)
        if not gt_match:
            continue
        ground_truth = float(gt_match.group(1))
        
        # Extract the final "adv rel err" value (this is the key difference)
        rel_err_match = re.search(r'adv rel err = ([\d.e+-]+)', section_content)
        if not rel_err_match:
            continue
        relative_error = float(rel_err_match.group(1))
        
        # Extract estimate from the final mean
        mean_match = re.search(r'# Mean = ([\d.e+-]+)', section_content)
        estimate = float(mean_match.group(1)) if mean_match else 0.0
        
        results[algorithm][q_value] = {
            'ground_truth': ground_truth,
            'estimate': estimate,
            'relative_error': relative_error,
            'std_relative_error': 0.0  # Not available in this file
        }
    
    return results

def update_lrcwiki_database():
    """Update the database with the corrected lrcwiki results"""
    
    # Parse the file
    results = parse_lrcwiki_2M_file()
    
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    dataset = "lrcwiki"
    epsilon = 1.0
    p = 3
    T_samples = -1  # All triplets used
    
    # Update results for each algorithm
    for algorithm, q_data in results.items():
        print(f"Updating {algorithm} results...")
        for q, data in q_data.items():
            try:
                cursor.execute('''
                    UPDATE algorithm_comparison_p3 
                    SET ground_truth = ?, estimate = ?, relative_error = ?, std_relative_error = ?
                    WHERE dataset = ? AND algorithm = ? AND q = ? AND p = ? AND epsilon = ?
                ''', (data['ground_truth'], data['estimate'], data['relative_error'], 
                      data['std_relative_error'], dataset, algorithm, q, p, epsilon))
                print(f"  Q={q}: Ground Truth={data['ground_truth']:.2e}, Estimate={data['estimate']:.2e}, Rel Error={data['relative_error']:.3e}")
            except Exception as e:
                print(f"Error updating {algorithm} Q={q}: {e}")
    
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
    update_lrcwiki_database()

