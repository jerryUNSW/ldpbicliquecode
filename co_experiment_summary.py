#!/usr/bin/env python3

import os
import glob
import re
import pandas as pd

def parse_log_file(log_file):
    """Parse a single log file and extract key metrics"""
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        
        # Extract parameters from filename
        filename = os.path.basename(log_file)
        match = re.match(r'co_eps(\d+)_p(\d+)_q(\d+)_alg(\d+)\.log', filename)
        if not match:
            return None
            
        epsilon = int(match.group(1))
        p = int(match.group(2))
        q = int(match.group(3))
        algorithm = int(match.group(4))
        
        # Map algorithm numbers to names
        alg_names = {0: 'Naive', 1: 'ADV', 2: 'ADV+', 3: 'ADV++', 4: 'ADV+++'}
        algorithm_name = alg_names.get(algorithm, f'Unknown-{algorithm}')
        
        # Extract metrics
        estimate_match = re.search(r'estimate = ([\d.e+-]+)', content)
        real_match = re.search(r'real count = ([\d.]+)', content)
        rel_error_match = re.search(r'adv rel err = ([\d.e+-]+|nan|-nan)', content)
        time_match = re.search(r'time:([\d.]+)', content)
        
        if estimate_match and real_match and rel_error_match:
            estimate = float(estimate_match.group(1))
            real_count = float(real_match.group(1))
            
            # Handle -nan values
            rel_error_str = rel_error_match.group(1)
            if rel_error_str in ['-nan', 'nan']:
                relative_error = float('inf')  # Use infinity for nan values
            else:
                relative_error = float(rel_error_str)
                
            time_taken = float(time_match.group(1)) if time_match else None
            
            return {
                'epsilon': epsilon,
                'p': p,
                'q': q,
                'algorithm': algorithm,
                'algorithm_name': algorithm_name,
                'estimate': estimate,
                'real_count': real_count,
                'relative_error': relative_error,
                'time': time_taken,
                'log_file': log_file
            }
        else:
            return None
            
    except Exception as e:
        return None

def main():
    logs_dir = 'larger_p_q_experiment/logs'
    
    # Parse all log files
    log_files = glob.glob(os.path.join(logs_dir, "*.log"))
    results = []
    
    for log_file in log_files:
        result = parse_log_file(log_file)
        if result:
            results.append(result)
    
    if not results:
        print("No valid results found!")
        return
    
    df = pd.DataFrame(results)
    
    # Create summary tables
    print("=== MULTI-DATASET EXPERIMENT SUMMARY ===\n")
    
    # 1. Overall statistics
    print("1. OVERALL STATISTICS")
    print(f"Total experiments: {len(df)}")
    print(f"Successful experiments: {len(df[df['relative_error'] != float('inf')])}")
    print(f"Failed experiments: {len(df[df['relative_error'] == float('inf')])}")
    print(f"Datasets: {sorted(df['dataset'].unique())}")
    print()
    
    # 2. Algorithm performance summary
    print("2. ALGORITHM PERFORMANCE (Mean Relative Error %)")
    df_finite = df[df['relative_error'] != float('inf')].copy()
    alg_summary = df_finite.groupby('algorithm_name')['relative_error'].agg(['mean', 'std', 'count']).round(2)
    print(alg_summary)
    print()
    
    # 3. Dataset performance summary
    print("3. DATASET PERFORMANCE (Mean Relative Error %)")
    dataset_summary = df_finite.groupby('dataset')['relative_error'].agg(['mean', 'std', 'count']).round(2)
    print(dataset_summary)
    print()
    
    # 4. Epsilon performance summary
    print("4. EPSILON PERFORMANCE (Mean Relative Error %)")
    eps_summary = df_finite.groupby('epsilon')['relative_error'].agg(['mean', 'std', 'count']).round(2)
    print(eps_summary)
    print()
    
    # 5. Best performing combinations
    print("5. TOP 10 BEST PERFORMING COMBINATIONS")
    best_results = df_finite.nsmallest(10, 'relative_error')[['dataset', 'epsilon', 'p', 'q', 'algorithm_name', 'relative_error', 'real_count', 'estimate']]
    print(best_results.to_string(index=False))
    print()
    
    # 6. Worst performing combinations
    print("6. TOP 10 WORST PERFORMING COMBINATIONS")
    worst_results = df_finite.nlargest(10, 'relative_error')[['dataset', 'epsilon', 'p', 'q', 'algorithm_name', 'relative_error', 'real_count', 'estimate']]
    print(worst_results.to_string(index=False))
    print()
    
    # 7. Performance by (P,Q) combinations
    print("7. PERFORMANCE BY (P,Q) COMBINATIONS")
    df_finite['pq_label'] = df_finite['p'].astype(str) + ',' + df_finite['q'].astype(str)
    pq_summary = df_finite.groupby('pq_label')['relative_error'].agg(['mean', 'std', 'count']).round(2)
    print(pq_summary)
    print()
    
    # 8. Algorithm comparison by epsilon
    print("8. ALGORITHM COMPARISON BY EPSILON")
    comparison = df_finite.groupby(['algorithm_name', 'epsilon'])['relative_error'].mean().unstack().round(2)
    print(comparison)
    print()
    
    # 9. Failed experiments analysis
    failed_df = df[df['relative_error'] == float('inf')]
    if len(failed_df) > 0:
        print("9. FAILED EXPERIMENTS ANALYSIS")
        print("Failed experiments by dataset:")
        print(failed_df['dataset'].value_counts())
        print("\nFailed experiments by algorithm:")
        print(failed_df['algorithm_name'].value_counts())
        print("\nFailed experiments by epsilon:")
        print(failed_df['epsilon'].value_counts())
        print("\nFailed experiments by (P,Q):")
        failed_df['pq_label'] = failed_df['p'].astype(str) + ',' + failed_df['q'].astype(str)
        print(failed_df['pq_label'].value_counts())
    else:
        print("9. FAILED EXPERIMENTS ANALYSIS")
        print("No failed experiments!")
    
    print("\n=== END OF SUMMARY ===")

if __name__ == "__main__":
    main()
