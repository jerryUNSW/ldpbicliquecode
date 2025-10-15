#!/usr/bin/env python3

import os
import re
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def parse_log_file(log_file):
    """Parse a single log file and extract results"""
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        
        filename = os.path.basename(log_file)
        parts = filename.replace('.log', '').split('_')
        
        dataset = parts[0]
        epsilon = int(parts[1].replace('eps', ''))
        p = int(parts[2].replace('p', ''))
        q = int(parts[3].replace('q', ''))
        algorithm = int(parts[4].replace('alg', ''))
        
        # Map algorithm numbers to names
        alg_names = {0: 'Naive', 1: 'ADV', 2: 'ADV+', 3: 'ADV++', 4: 'ADV+++'}
        algorithm_name = alg_names.get(algorithm, f'Unknown-{algorithm}')
        
        # Extract relative errors
        rel_error_matches = re.findall(r'adv rel err = ([\d.-]+)', content)
        
        if rel_error_matches:
            rel_errors = []
            for x in rel_error_matches:
                if x == '-nan' or x == 'nan':
                    rel_errors.append(float('inf'))
                else:
                    rel_errors.append(float(x))
            
            # Calculate statistics
            finite_errors = [e for e in rel_errors if e != float('inf')]
            if finite_errors:
                mean_rel_error = np.mean(finite_errors)
                std_rel_error = np.std(finite_errors)
            else:
                mean_rel_error = float('inf')
                std_rel_error = 0.0
            
            # Extract other metrics
            estimate_match = re.search(r'# Mean = ([\d.]+)', content)
            real_match = re.search(r'real count = ([\d.]+)', content)
            time_match = re.search(r'time:([\d.]+)', content)
            
            estimate = float(estimate_match.group(1)) if estimate_match else None
            real_count = float(real_match.group(1)) if real_match else None
            time_taken = float(time_match.group(1)) if time_match else None
            
            return {
                'dataset': dataset,
                'epsilon': epsilon,
                'p': p,
                'q': q,
                'algorithm': algorithm,
                'algorithm_name': algorithm_name,
                'estimate': estimate,
                'real_count': real_count,
                'relative_error': mean_rel_error,
                'relative_error_std': std_rel_error,
                'num_rounds': len(rel_errors),
                'time': time_taken,
                'log_file': log_file
            }
        else:
            print(f"Warning: No relative errors found in {log_file}")
            return None
            
    except Exception as e:
        print(f"Error parsing {log_file}: {e}")
        return None

def main():
    # Analyze results
    logs_dir = 'larger_p_q_experiment/logs'
    results_dir = 'larger_p_q_experiment'
    
    # Get all log files
    log_files = glob.glob(os.path.join(logs_dir, "*.log"))
    results = []
    
    print(f"Found {len(log_files)} log files to parse...")
    
    for log_file in log_files:
        result = parse_log_file(log_file)
        if result:
            results.append(result)
    
    if not results:
        print("No valid results found!")
        return
    
    df = pd.DataFrame(results)
    
    print(f"Successfully parsed {len(results)} experiments")
    print(f"Datasets: {sorted(df['dataset'].unique())}")
    print(f"Algorithms: {sorted(df['algorithm_name'].unique())}")
    print(f"Epsilons: {sorted(df['epsilon'].unique())}")
    print(f"P values: {sorted(df['p'].unique())}")
    print(f"Q values: {sorted(df['q'].unique())}")
    
    # Filter out infinite relative errors for analysis
    df_finite = df[df['relative_error'] != float('inf')].copy()
    print(f"Experiments with finite relative error: {len(df_finite)}")
    
    if df_finite.empty:
        print("No experiments with finite relative error found!")
        return
    
    # Set up plotting
    plt.style.use('default')
    algorithm_order = ['Naive', 'ADV', 'ADV+', 'ADV++', 'ADV+++']
    
    # Generate separate plots for each (dataset, p, epsilon) combination
    datasets = sorted(df['dataset'].unique())
    epsilons = sorted(df['epsilon'].unique())
    p_values = sorted(df['p'].unique())
    
    for dataset in datasets:
        dataset_df = df_finite[df_finite['dataset'] == dataset]
        if dataset_df.empty:
            continue
            
        for p in p_values:
            for eps in epsilons:
                fig, ax = plt.subplots(1, 1, figsize=(10, 6))
                
                # Filter for this specific (dataset, p, epsilon) combination
                plot_df = dataset_df[(dataset_df['p'] == p) & (dataset_df['epsilon'] == eps)]
                if plot_df.empty:
                    plt.close()
                    continue
                
                # Group by algorithm and Q, calculate mean relative error
                grouped = plot_df.groupby(['algorithm_name', 'q'])['relative_error'].mean().reset_index()
                pivot_df = grouped.pivot(index='q', columns='algorithm_name', values='relative_error')
                
                # Ensure consistent algorithm ordering
                available_algs = [alg for alg in algorithm_order if alg in pivot_df.columns]
                pivot_df = pivot_df[available_algs]
                
                if pivot_df.empty:
                    plt.close()
                    continue
                
                # Create the plot
                pivot_df.plot(kind='bar', ax=ax, width=0.8)
                ax.set_title(f'Algorithm Comparison by Q Value - {dataset.upper()} Dataset\nP = {p}, Epsilon = {eps}', 
                            fontsize=14, fontweight='bold')
                ax.set_xlabel('Q Value', fontsize=12)
                ax.set_ylabel('Mean Relative Error (%)', fontsize=12)
                ax.set_yscale('log')
                ax.legend(title='Algorithm', fontsize=10, loc='upper right')
                ax.grid(True, alpha=0.3)
                ax.tick_params(axis='x', rotation=0, labelsize=11)
                ax.tick_params(axis='y', labelsize=11)
                
                plt.tight_layout()
                
                # Save with descriptive filename
                filename = f'algorithm_comparison_{dataset}_p{p}_eps{eps}.pdf'
                plt.savefig(os.path.join(results_dir, filename), 
                           bbox_inches='tight', dpi=300)
                plt.close()
                
                print(f"Generated plot: {filename}")
    
    print(f"\nAll plots saved to: {results_dir}/")
    print("Plots generated for each (dataset, P, epsilon) combination showing Q values")

if __name__ == "__main__":
    main()
