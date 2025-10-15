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
            
            mean_rel_error = np.mean([e for e in rel_errors if e != float('inf')]) if any(e != float('inf') for e in rel_errors) else float('inf')
            std_rel_error = np.std([e for e in rel_errors if e != float('inf')]) if any(e != float('inf') for e in rel_errors) else 0.0
            
            return {
                'dataset': dataset,
                'epsilon': epsilon,
                'p': p,
                'q': q,
                'algorithm': algorithm,
                'algorithm_name': algorithm_name,
                'relative_error': mean_rel_error,
                'relative_error_std': std_rel_error,
                'num_rounds': len(rel_errors),
                'log_file': log_file
            }
        else:
            return None
            
    except Exception as e:
        print(f"Error parsing {log_file}: {e}")
        return None

def main():
    logs_dir = 'larger_p_q_experiment/logs'
    results_dir = 'larger_p_q_experiment'
    
    # Focus only on CO dataset with P=3 and P=4
    co_log_files = glob.glob(os.path.join(logs_dir, "co_*p3*.log")) + glob.glob(os.path.join(logs_dir, "co_*p4*.log"))
    
    print(f"Found {len(co_log_files)} CO log files for P=3 and P=4...")
    
    results = []
    for log_file in co_log_files:
        result = parse_log_file(log_file)
        if result:
            results.append(result)
    
    if not results:
        print("No valid results found!")
        return
    
    df = pd.DataFrame(results)
    df_finite = df[df['relative_error'] != float('inf')].copy()
    
    print(f"Successfully parsed {len(results)} experiments")
    print(f"Experiments with finite relative error: {len(df_finite)}")
    print(f"P values: {sorted(df['p'].unique())}")
    print(f"Epsilons: {sorted(df['epsilon'].unique())}")
    print(f"Q values: {sorted(df['q'].unique())}")
    
    plt.style.use('default')
    algorithm_order = ['Naive', 'ADV', 'ADV+', 'ADV++', 'ADV+++']
    
    # Generate plots for CO dataset, P=3 and P=4
    for p_val in [3, 4]:
        p_df = df_finite[df_finite['p'] == p_val]
        if p_df.empty:
            print(f"No data for CO P={p_val}")
            continue
            
        for eps in sorted(p_df['epsilon'].unique()):
            fig, ax = plt.subplots(1, 1, figsize=(10, 6))
            
            eps_df = p_df[p_df['epsilon'] == eps]
            if eps_df.empty:
                plt.close()
                continue
            
            # Group by algorithm and Q value, calculate mean relative error
            grouped = eps_df.groupby(['algorithm_name', 'q'])['relative_error'].mean().reset_index()
            pivot_df = grouped.pivot(index='q', columns='algorithm_name', values='relative_error')
            
            # Ensure consistent algorithm ordering
            available_algs = [alg for alg in algorithm_order if alg in pivot_df.columns]
            if not available_algs:
                plt.close()
                continue
                
            pivot_df = pivot_df[available_algs]
            
            # Create the plot
            pivot_df.plot(kind='bar', ax=ax, width=0.8)
            ax.set_title(f'Algorithm Comparison by Q Value - CO Dataset (P={p_val}, Epsilon = {eps})', 
                        fontsize=14, fontweight='bold')
            ax.set_xlabel('Q Value', fontsize=12)
            ax.set_ylabel('Mean Relative Error (%)', fontsize=12)
            ax.set_yscale('log')
            ax.legend(title='Algorithm', fontsize=10)
            ax.grid(True, alpha=0.3)
            ax.tick_params(axis='x', rotation=0)
            
            plt.tight_layout()
            output_file = os.path.join(results_dir, f'algorithm_comparison_co_p{p_val}_eps{eps}.pdf')
            plt.savefig(output_file, bbox_inches='tight')
            plt.close()
            
            print(f"Generated plot: algorithm_comparison_co_p{p_val}_eps{eps}.pdf")
    
    print(f"\nAll plots saved to: {results_dir}/")

if __name__ == "__main__":
    main()
