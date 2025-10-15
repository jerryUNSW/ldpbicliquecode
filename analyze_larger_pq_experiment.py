#!/usr/bin/env python3

import os
import re
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

def parse_log_file(log_file):
    """Parse a single log file and extract key metrics"""
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        
        # Extract key information from filename
        filename = os.path.basename(log_file)
        parts = filename.replace('.log', '').split('_')
        
        dataset = parts[0]
        epsilon = int(parts[1].replace('eps', ''))
        p = int(parts[2].replace('p', ''))
        q = int(parts[3].replace('q', ''))
        algorithm = int(parts[4].replace('alg', ''))
        
        # Algorithm names
        alg_names = {0: 'Naive', 1: 'ADV', 2: 'ADV+', 3: 'ADV++'}
        algorithm_name = alg_names.get(algorithm, f'Unknown-{algorithm}')
        
        # Extract metrics from log content
        estimate_match = re.search(r'# Mean = ([\d.]+)', content)
        real_match = re.search(r'real count = ([\d.]+)', content)
        rel_error_match = re.search(r'adv rel err = ([\d.]+)', content)
        time_match = re.search(r'time:([\d.]+)', content)
        
        if estimate_match and real_match and rel_error_match:
            estimate = float(estimate_match.group(1))
            real_count = float(real_match.group(1))
            # Handle -nan values
            rel_error_str = rel_error_match.group(1)
            if rel_error_str == '-nan' or rel_error_str == 'nan':
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

def parse_all_logs(logs_dir):
    """Parse all log files in the directory"""
    log_files = glob.glob(os.path.join(logs_dir, "*.log"))
    results = []
    
    print(f"Found {len(log_files)} log files to parse...")
    
    for log_file in log_files:
        result = parse_log_file(log_file)
        if result:
            results.append(result)
    
    return pd.DataFrame(results)

def create_relative_error_plots(df, output_dir):
    """Create plots showing relative error by algorithm, epsilon, p, q"""
    
    # Set up the plotting style
    plt.style.use('default')
    fig_size = (15, 10)
    
    # 1. Relative error by algorithm and epsilon (averaged over p,q)
    fig, axes = plt.subplots(2, 2, figsize=fig_size)
    fig.suptitle('Relative Error Analysis - Larger PQ Experiment', fontsize=16)
    
    datasets = df['dataset'].unique()
    
    for i, dataset in enumerate(datasets):
        ax = axes[i//2, i%2] if len(datasets) > 1 else axes[0]
        
        dataset_df = df[df['dataset'] == dataset]
        
        # Group by algorithm and epsilon, compute mean relative error
        grouped = dataset_df.groupby(['algorithm_name', 'epsilon'])['relative_error'].mean().reset_index()
        
        # Pivot for plotting
        pivot_df = grouped.pivot(index='epsilon', columns='algorithm_name', values='relative_error')
        
        # Plot
        pivot_df.plot(kind='bar', ax=ax, width=0.8)
        ax.set_title(f'Dataset: {dataset}')
        ax.set_xlabel('Epsilon')
        ax.set_ylabel('Mean Relative Error (%)')
        ax.set_yscale('log')
        ax.legend(title='Algorithm')
        ax.grid(True, alpha=0.3)
        
        # Rotate x-axis labels
        ax.tick_params(axis='x', rotation=0)
    
    # Hide unused subplots if only one dataset
    if len(datasets) == 1:
        axes[0, 1].set_visible(False)
        axes[1, 0].set_visible(False)
        axes[1, 1].set_visible(False)
    
    plt.tight_layout()
    # Skip saving this plot as requested
    plt.close()
    
    # 2. Relative error by (p,q) combinations
    fig, axes = plt.subplots(2, 2, figsize=fig_size)
    fig.suptitle('Relative Error by (P,Q) Combinations', fontsize=16)
    
    for i, dataset in enumerate(datasets):
        ax = axes[i//2, i%2] if len(datasets) > 1 else axes[0]
        
        dataset_df = df[df['dataset'] == dataset]
        
        # Create (p,q) labels
        dataset_df['pq_label'] = dataset_df['p'].astype(str) + ',' + dataset_df['q'].astype(str)
        
        # Group by algorithm and (p,q), compute mean relative error
        grouped = dataset_df.groupby(['algorithm_name', 'pq_label'])['relative_error'].mean().reset_index()
        
        # Pivot for plotting
        pivot_df = grouped.pivot(index='pq_label', columns='algorithm_name', values='relative_error')
        
        # Plot
        pivot_df.plot(kind='bar', ax=ax, width=0.8)
        ax.set_title(f'Dataset: {dataset}')
        ax.set_xlabel('(P,Q)')
        ax.set_ylabel('Mean Relative Error (%)')
        ax.set_yscale('log')
        ax.legend(title='Algorithm')
        ax.grid(True, alpha=0.3)
        
        # Rotate x-axis labels
        ax.tick_params(axis='x', rotation=45)
    
    # Hide unused subplots if only one dataset
    if len(datasets) == 1:
        axes[0, 1].set_visible(False)
        axes[1, 0].set_visible(False)
        axes[1, 1].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'relative_error_by_pq_combinations.png'), 
                dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'relative_error_by_pq_combinations.pdf'), 
                bbox_inches='tight')
    plt.close()
    
    # 3. Heatmap of relative error by algorithm and (p,q) for each dataset
    for dataset in datasets:
        dataset_df = df[df['dataset'] == dataset]
        
        fig, axes = plt.subplots(1, len(df['epsilon'].unique()), figsize=(5*len(df['epsilon'].unique()), 8))
        if len(df['epsilon'].unique()) == 1:
            axes = [axes]
        
        fig.suptitle(f'Relative Error Heatmap - Dataset: {dataset}', fontsize=16)
        
        for i, epsilon in enumerate(sorted(df['epsilon'].unique())):
            ax = axes[i]
            
            eps_df = dataset_df[dataset_df['epsilon'] == epsilon]
            
            # Create pivot table
            pivot_df = eps_df.pivot_table(
                values='relative_error', 
                index='q', 
                columns=['p', 'algorithm_name'], 
                aggfunc='mean'
            )
            
            # Flatten column names for better display
            pivot_df.columns = [f'P{p}-{alg}' for p, alg in pivot_df.columns]
            
            # Create heatmap
            im = ax.imshow(pivot_df.values, cmap='YlOrRd', aspect='auto')
            
            # Set labels
            ax.set_xticks(range(len(pivot_df.columns)))
            ax.set_xticklabels(pivot_df.columns, rotation=45, ha='right')
            ax.set_yticks(range(len(pivot_df.index)))
            ax.set_yticklabels(pivot_df.index)
            
            ax.set_xlabel('Algorithm and P')
            ax.set_ylabel('Q')
            ax.set_title(f'Epsilon = {epsilon}')
            
            # Add colorbar
            plt.colorbar(im, ax=ax, label='Relative Error (%)')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'heatmap_{dataset}.png'), 
                    dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, f'heatmap_{dataset}.pdf'), 
                    bbox_inches='tight')
        plt.close()

def create_summary_table(df, output_dir):
    """Create summary tables of results"""
    
    # Overall summary by algorithm
    summary = df.groupby(['algorithm_name', 'dataset']).agg({
        'relative_error': ['mean', 'std', 'count'],
        'time': ['mean', 'std']
    }).round(2)
    
    # Save to CSV
    summary.to_csv(os.path.join(output_dir, 'summary_by_algorithm.csv'))
    
    # Detailed results
    detailed = df[['dataset', 'epsilon', 'p', 'q', 'algorithm_name', 
                   'estimate', 'real_count', 'relative_error', 'time']].copy()
    detailed = detailed.sort_values(['dataset', 'epsilon', 'p', 'q', 'algorithm_name'])
    detailed.to_csv(os.path.join(output_dir, 'detailed_results.csv'), index=False)
    
    # Print summary
    print("\n=== SUMMARY BY ALGORITHM ===")
    print(summary)
    
    print(f"\nDetailed results saved to: {os.path.join(output_dir, 'detailed_results.csv')}")
    print(f"Summary saved to: {os.path.join(output_dir, 'summary_by_algorithm.csv')}")

def main():
    # Configuration
    logs_dir = "larger_p_q_experiment/logs"
    results_dir = "larger_p_q_experiment"
    
    # Create results directory
    os.makedirs(results_dir, exist_ok=True)
    
    # Parse all log files
    print("Parsing log files...")
    df = parse_all_logs(logs_dir)
    
    if df.empty:
        print("No valid results found!")
        return
    
    print(f"Successfully parsed {len(df)} experiments")
    print(f"Datasets: {df['dataset'].unique()}")
    print(f"Algorithms: {df['algorithm_name'].unique()}")
    print(f"Epsilons: {sorted(df['epsilon'].unique())}")
    print(f"P values: {sorted(df['p'].unique())}")
    print(f"Q values: {sorted(df['q'].unique())}")
    
    # Create plots
    print("\nCreating plots...")
    create_relative_error_plots(df, results_dir)
    
    # Create summary tables
    print("\nCreating summary tables...")
    create_summary_table(df, results_dir)
    
    print(f"\nAnalysis complete! Results saved to: {results_dir}")

if __name__ == "__main__":
    main()
