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
    # Analyze unicode results
    logs_dir = 'larger_p_q_experiment/logs'
    results_dir = 'larger_p_q_experiment'
    
    # Get all unicode log files
    unicode_log_files = glob.glob(os.path.join(logs_dir, "unicode_*.log"))
    results = []
    
    print(f"Found {len(unicode_log_files)} unicode log files to parse...")
    
    for log_file in unicode_log_files:
        result = parse_log_file(log_file)
        if result:
            results.append(result)
    
    if not results:
        print("No valid results found!")
        return
    
    df = pd.DataFrame(results)
    
    print(f"Successfully parsed {len(results)} unicode experiments")
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
    
    # Generate algorithm comparison by Q plots for each epsilon
    epsilons = sorted(df['epsilon'].unique())
    
    for eps in epsilons:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        eps_df = df_finite[df_finite['epsilon'] == eps]
        if eps_df.empty:
            plt.close()
            continue
        
        # Group by algorithm and Q, calculate mean relative error
        grouped = eps_df.groupby(['algorithm_name', 'q'])['relative_error'].mean().reset_index()
        pivot_df = grouped.pivot(index='q', columns='algorithm_name', values='relative_error')
        
        # Ensure consistent algorithm ordering
        available_algs = [alg for alg in algorithm_order if alg in pivot_df.columns]
        pivot_df = pivot_df[available_algs]
        
        # Create the plot
        pivot_df.plot(kind='bar', ax=ax, width=0.8)
        ax.set_title(f'Algorithm Comparison by Q Value - UNICODE Dataset (Epsilon = {eps})', 
                    fontsize=16, fontweight='bold')
        ax.set_xlabel('Q Value', fontsize=14)
        ax.set_ylabel('Mean Relative Error (%)', fontsize=14)
        ax.set_yscale('log')
        ax.legend(title='Algorithm', fontsize=12, loc='upper right')
        ax.grid(True, alpha=0.3)
        ax.tick_params(axis='x', rotation=0, labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
        
        plt.tight_layout()
        plt.savefig(os.path.join(results_dir, f'unicode_algorithm_comparison_by_q_eps{eps}.pdf'), 
                   bbox_inches='tight', dpi=300)
        plt.close()
        
        print(f"Generated plot: unicode_algorithm_comparison_by_q_eps{eps}.pdf")
    
    # Generate summary statistics
    print("\n" + "="*60)
    print("UNICODE DATASET ANALYSIS SUMMARY")
    print("="*60)
    
    print(f"\nTotal experiments: {len(df)}")
    print(f"Successful experiments: {len(df_finite)}")
    print(f"Failed experiments: {len(df) - len(df_finite)}")
    
    print("\nRelative Error by Algorithm (mean ± std):")
    for alg in sorted(df_finite['algorithm_name'].unique()):
        alg_data = df_finite[df_finite['algorithm_name'] == alg]['relative_error']
        print(f"  {alg}: {alg_data.mean():.2f} ± {alg_data.std():.2f}%")
    
    print("\nRelative Error by Epsilon (mean ± std):")
    for eps in sorted(df_finite['epsilon'].unique()):
        eps_data = df_finite[df_finite['epsilon'] == eps]['relative_error']
        print(f"  Epsilon {eps}: {eps_data.mean():.2f} ± {eps_data.std():.2f}%")
    
    print("\nRelative Error by P Value (mean ± std):")
    for p in sorted(df_finite['p'].unique()):
        p_data = df_finite[df_finite['p'] == p]['relative_error']
        print(f"  P = {p}: {p_data.mean():.2f} ± {p_data.std():.2f}%")
    
    print("\nRelative Error by Q Value (mean ± std):")
    for q in sorted(df_finite['q'].unique()):
        q_data = df_finite[df_finite['q'] == q]['relative_error']
        print(f"  Q = {q}: {q_data.mean():.2f} ± {q_data.std():.2f}%")
    
    # Find best and worst performing combinations
    print("\nTop 5 Best Performing Combinations:")
    best_combinations = df_finite.nsmallest(5, 'relative_error')[['epsilon', 'p', 'q', 'algorithm_name', 'relative_error']]
    for idx, row in best_combinations.iterrows():
        print(f"  ε={row['epsilon']}, p={row['p']}, q={row['q']}, {row['algorithm_name']}: {row['relative_error']:.2f}%")
    
    print("\nTop 5 Worst Performing Combinations:")
    worst_combinations = df_finite.nlargest(5, 'relative_error')[['epsilon', 'p', 'q', 'algorithm_name', 'relative_error']]
    for idx, row in worst_combinations.iterrows():
        print(f"  ε={row['epsilon']}, p={row['p']}, q={row['q']}, {row['algorithm_name']}: {row['relative_error']:.2f}%")
    
    # Save detailed results to CSV
    df_finite.to_csv(os.path.join(results_dir, 'unicode_detailed_results.csv'), index=False)
    print(f"\nDetailed results saved to: unicode_detailed_results.csv")
    
    print(f"\nPlots saved to: {results_dir}/")
    for eps in sorted(df['epsilon'].unique()):
        print(f"  - unicode_algorithm_comparison_by_q_eps{eps}.pdf")

if __name__ == "__main__":
    main()

