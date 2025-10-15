#!/usr/bin/env python3

import os
import re
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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
        alg_names = {0: 'Naive', 2: 'ADV+', 3: 'ADV++', 4: 'ADV+++'}
        algorithm_name = alg_names.get(algorithm, f'Unknown-{algorithm}')
        
        # Extract all relative errors from the log (for multi-round experiments)
        rel_error_matches = re.findall(r'adv rel err = ([\d.]+)', content)
        
        if rel_error_matches:
            # Convert to float and compute statistics
            rel_errors = [float(x) for x in rel_error_matches]
            mean_rel_error = np.mean(rel_errors)
            std_rel_error = np.std(rel_errors)
            
            # Get other metrics from the last run
            estimate_match = re.search(r'# Mean = ([\d.]+)', content)
            real_match = re.search(r'real count = ([\d.]+)', content)
            
            estimate = float(estimate_match.group(1)) if estimate_match else None
            real_count = float(real_match.group(1)) if real_match else None
            
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
                'num_rounds': len(rel_errors)
            }
        else:
            return None
            
    except Exception as e:
        print(f"Error parsing {log_file}: {e}")
        return None

def parse_all_logs(logs_dir):
    """Parse all log files in the directory"""
    log_files = glob.glob(os.path.join(logs_dir, "*.log"))
    results = []
    
    for log_file in log_files:
        result = parse_log_file(log_file)
        if result:
            results.append(result)
    
    return pd.DataFrame(results)

def create_algorithm_comparison_plots(df, output_dir):
    """Create focused plots comparing algorithms as Q increases"""
    
    if df.empty:
        print("No data to plot!")
        return
    
    # Filter for algorithms 0, 2, 3, 4 and focus on P=2, epsilon=1
    focus_df = df[(df['algorithm'].isin([0, 2, 3, 4])) & 
                  (df['p'] == 2) & 
                  (df['epsilon'] == 1)].copy()
    
    if focus_df.empty:
        print("No data for P=2, epsilon=1 with algorithms 0,2,3,4")
        return
    
    # Set up the plotting style
    plt.style.use('default')
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Algorithm Comparison: Relative Error vs Q (P=2, Epsilon=1)', fontsize=16)
    
    # Plot 1: Line plot showing relative error vs Q for each algorithm
    ax1 = axes[0, 0]
    algorithms = [0, 2, 3, 4]
    alg_names = {0: 'Naive', 2: 'ADV+', 3: 'ADV++', 4: 'ADV+++'}
    colors = {0: 'red', 2: 'blue', 3: 'green', 4: 'orange'}
    markers = {0: 'o', 2: 's', 3: '^', 4: 'D'}
    
    for alg in algorithms:
        alg_data = focus_df[focus_df['algorithm'] == alg].sort_values('q')
        if not alg_data.empty:
            ax1.plot(alg_data['q'], alg_data['relative_error'], 
                    marker=markers[alg], label=alg_names[alg], 
                    color=colors[alg], linewidth=2, markersize=8)
    
    ax1.set_xlabel('Q')
    ax1.set_ylabel('Relative Error (%)')
    ax1.set_yscale('log')
    ax1.set_title('Relative Error vs Q')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(sorted(focus_df['q'].unique()))
    
    # Plot 2: Bar plot comparing algorithms for each Q
    ax2 = axes[0, 1]
    q_values = sorted(focus_df['q'].unique())
    x = np.arange(len(q_values))
    width = 0.2
    
    for i, alg in enumerate(algorithms):
        alg_data = focus_df[focus_df['algorithm'] == alg]
        if not alg_data.empty:
            rel_errors = []
            for q in q_values:
                q_data = alg_data[alg_data['q'] == q]
                if not q_data.empty:
                    rel_errors.append(q_data['relative_error'].iloc[0])
                else:
                    rel_errors.append(0)  # No data for this Q
            
            ax2.bar(x + i*width, rel_errors, width, 
                   label=alg_names[alg], color=colors[alg], alpha=0.8)
    
    ax2.set_xlabel('Q')
    ax2.set_ylabel('Relative Error (%)')
    ax2.set_yscale('log')
    ax2.set_title('Relative Error by Algorithm and Q')
    ax2.set_xticks(x + width * 1.5)
    ax2.set_xticklabels(q_values)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Heatmap of relative errors
    ax3 = axes[1, 0]
    pivot_data = focus_df.pivot_table(
        values='relative_error', 
        index='q', 
        columns='algorithm_name', 
        aggfunc='mean'
    )
    
    # Reorder columns to match our desired order
    column_order = ['Naive', 'ADV+', 'ADV++', 'ADV+++']
    pivot_data = pivot_data.reindex(columns=column_order)
    
    im = ax3.imshow(pivot_data.values, cmap='YlOrRd', aspect='auto')
    ax3.set_xticks(range(len(pivot_data.columns)))
    ax3.set_xticklabels(pivot_data.columns)
    ax3.set_yticks(range(len(pivot_data.index)))
    ax3.set_yticklabels(pivot_data.index)
    ax3.set_xlabel('Algorithm')
    ax3.set_ylabel('Q')
    ax3.set_title('Relative Error Heatmap')
    
    # Add text annotations
    for i in range(len(pivot_data.index)):
        for j in range(len(pivot_data.columns)):
            if not np.isnan(pivot_data.iloc[i, j]):
                ax3.text(j, i, f'{pivot_data.iloc[i, j]:.1f}', 
                        ha='center', va='center', fontweight='bold')
    
    plt.colorbar(im, ax=ax3, label='Relative Error (%)')
    
    # Plot 4: Algorithm ranking by Q
    ax4 = axes[1, 1]
    ranking_data = []
    for q in q_values:
        q_data = focus_df[focus_df['q'] == q]
        if not q_data.empty:
            # Sort by relative error (lower is better)
            ranked = q_data.sort_values('relative_error')
            for rank, (_, row) in enumerate(ranked.iterrows(), 1):
                ranking_data.append({
                    'Q': q,
                    'Algorithm': row['algorithm_name'],
                    'Rank': rank,
                    'Relative_Error': row['relative_error']
                })
    
    ranking_df = pd.DataFrame(ranking_data)
    
    # Create ranking plot
    for alg in ['Naive', 'ADV+', 'ADV++', 'ADV+++']:
        alg_ranks = ranking_df[ranking_df['Algorithm'] == alg]
        if not alg_ranks.empty:
            ax4.plot(alg_ranks['Q'], alg_ranks['Rank'], 
                    marker=markers[list(alg_names.keys())[list(alg_names.values()).index(alg)]], 
                    label=alg, linewidth=2, markersize=8)
    
    ax4.set_xlabel('Q')
    ax4.set_ylabel('Rank (1=Best)')
    ax4.set_title('Algorithm Ranking by Q')
    ax4.set_xticks(q_values)
    ax4.set_yticks([1, 2, 3, 4])
    ax4.invert_yaxis()  # Lower rank (1) at top
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'algorithm_comparison_by_q.png'), 
                dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'algorithm_comparison_by_q.pdf'), 
                bbox_inches='tight')
    plt.close()
    
    # Print summary table
    print("\n=== ALGORITHM COMPARISON BY Q ===")
    print("P=2, Epsilon=1")
    print("\nRelative Error (%) by Algorithm and Q:")
    
    summary_table = focus_df.pivot_table(
        values='relative_error', 
        index='algorithm_name', 
        columns='q', 
        aggfunc='mean'
    ).round(2)
    
    # Reorder rows and columns
    summary_table = summary_table.reindex(['Naive', 'ADV+', 'ADV++', 'ADV+++'])
    summary_table = summary_table.reindex(columns=sorted(summary_table.columns))
    
    print(summary_table)
    
    # Save detailed results
    focus_df.to_csv(os.path.join(output_dir, 'algorithm_comparison_detailed.csv'), index=False)
    summary_table.to_csv(os.path.join(output_dir, 'algorithm_comparison_summary.csv'))
    
    print(f"\nDetailed results saved to: {os.path.join(output_dir, 'algorithm_comparison_detailed.csv')}")
    print(f"Summary table saved to: {os.path.join(output_dir, 'algorithm_comparison_summary.csv')}")
    print(f"Plots saved to: {os.path.join(output_dir, 'algorithm_comparison_by_q.png')}")

def main():
    # Configuration
    logs_dir = "larger_p_q_experiment/logs"
    results_dir = "larger_p_q_experiment"
    
    # Create results directory
    os.makedirs(results_dir, exist_ok=True)
    
    # Parse all log files
    print("Parsing log files for algorithm comparison...")
    df = parse_all_logs(logs_dir)
    
    if df.empty:
        print("No valid results found!")
        return
    
    print(f"Successfully parsed {len(df)} experiments")
    
    # Create focused comparison plots
    print("\nCreating algorithm comparison plots...")
    create_algorithm_comparison_plots(df, results_dir)
    
    print(f"\nAlgorithm comparison analysis complete!")

if __name__ == "__main__":
    main()

