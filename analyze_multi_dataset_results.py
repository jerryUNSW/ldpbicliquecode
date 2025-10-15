#!/usr/bin/env python3

import os
import glob
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def parse_log_file(log_file):
    """Parse a single log file and extract key metrics"""
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        
        # Extract parameters from filename
        filename = os.path.basename(log_file)
        match = re.match(r'(\w+)_eps(\d+)_p(\d+)_q(\d+)_alg(\d+)\.log', filename)
        if not match:
            return None
            
        dataset = match.group(1)
        epsilon = int(match.group(2))
        p = int(match.group(3))
        q = int(match.group(4))
        algorithm = int(match.group(5))
        
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
                'dataset': dataset,
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
            print(f"Warning: Could not parse {log_file}")
            return None
            
    except Exception as e:
        print(f"Error parsing {log_file}: {e}")
        return None

def main():
    logs_dir = 'larger_p_q_experiment/logs'
    results_dir = 'larger_p_q_experiment'
    
    # Parse all log files
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
    
    # Filter out infinite relative errors for some plots
    df_finite = df[df['relative_error'] != float('inf')].copy()
    print(f"Experiments with finite relative error: {len(df_finite)}")
    
    # Create plots
    plt.style.use('default')
    
    # Algorithm Comparison by Q value - Separate plots for each dataset and epsilon
    # Define consistent algorithm order
    algorithm_order = ['Naive', 'ADV', 'ADV+', 'ADV++', 'ADV+++']
    
    datasets = sorted(df['dataset'].unique())
    epsilons = sorted(df['epsilon'].unique())
    
    for dataset in datasets:
        dataset_df = df_finite[df_finite['dataset'] == dataset]
        if dataset_df.empty:
            continue
            
        for eps in epsilons:
            # Create individual plot for each dataset and epsilon
            fig, ax = plt.subplots(1, 1, figsize=(10, 6))
            
            eps_df = dataset_df[dataset_df['epsilon'] == eps]
            if eps_df.empty:
                plt.close()
                continue
            
            # Group by algorithm and q, compute mean relative error
            grouped = eps_df.groupby(['algorithm_name', 'q'])['relative_error'].mean().reset_index()
            
            # Pivot for plotting
            pivot_df = grouped.pivot(index='q', columns='algorithm_name', values='relative_error')
            
            # Reorder columns to match desired order
            available_algs = [alg for alg in algorithm_order if alg in pivot_df.columns]
            pivot_df = pivot_df[available_algs]
            
            # Plot
            pivot_df.plot(kind='bar', ax=ax, width=0.8)
            ax.set_title(f'Algorithm Comparison by Q Value - {dataset.upper()} Dataset (Epsilon = {eps})', fontsize=14, fontweight='bold')
            ax.set_xlabel('Q Value', fontsize=12)
            ax.set_ylabel('Mean Relative Error (%)', fontsize=12)
            ax.set_yscale('log')
            ax.legend(title='Algorithm', fontsize=10)
            ax.grid(True, alpha=0.3)
            ax.tick_params(axis='x', rotation=0)
            
            plt.tight_layout()
            plt.savefig(os.path.join(results_dir, f'algorithm_comparison_by_q_{dataset}_eps{eps}.pdf'), 
                        bbox_inches='tight')
            plt.close()
    
    # Summary statistics
    print("\n=== SUMMARY STATISTICS ===")
    print("\nRelative Error by Algorithm (mean ± std):")
    for alg in sorted(df_finite['algorithm_name'].unique()):
        alg_data = df_finite[df_finite['algorithm_name'] == alg]['relative_error']
        print(f"{alg}: {alg_data.mean():.2f} ± {alg_data.std():.2f}")
    
    print("\nRelative Error by Dataset (mean ± std):")
    for dataset in sorted(df_finite['dataset'].unique()):
        dataset_data = df_finite[df_finite['dataset'] == dataset]['relative_error']
        print(f"{dataset}: {dataset_data.mean():.2f} ± {dataset_data.std():.2f}")
    
    print("\nRelative Error by Epsilon (mean ± std):")
    for eps in sorted(df_finite['epsilon'].unique()):
        eps_data = df_finite[df_finite['epsilon'] == eps]['relative_error']
        print(f"Epsilon {eps}: {eps_data.mean():.2f} ± {eps_data.std():.2f}")
    
    print("\nExperiments with infinite relative error (failed):")
    failed_df = df[df['relative_error'] == float('inf')]
    if len(failed_df) > 0:
        print(failed_df[['dataset', 'epsilon', 'p', 'q', 'algorithm_name']].value_counts().sort_index())
    else:
        print("None")
    
    print(f"\nPlots saved to: {results_dir}/")
    for dataset in sorted(df['dataset'].unique()):
        for eps in sorted(df['epsilon'].unique()):
            print(f"- algorithm_comparison_by_q_{dataset}_eps{eps}.pdf")

if __name__ == "__main__":
    main()
