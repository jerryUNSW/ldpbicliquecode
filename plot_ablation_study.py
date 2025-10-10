#!/usr/bin/env python3
"""
Ablation Study Visualization
Compares ADV (base), ADV+ (+ multi-estimator), and ADV++ (+ two noisy graphs)
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pathlib import Path

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (15, 5)
plt.rcParams['font.size'] = 11

def load_ablation_data(csv_path):
    """Load ablation study results from CSV"""
    df = pd.read_csv(csv_path, names=['dataset', 'alg_id', 'alg_name', 'mean', 'rel_err', 'time'])
    return df

def plot_ablation_comparison(df, output_dir='ablation_study_results'):
    """Create comprehensive ablation study plots"""
    
    datasets = df['dataset'].unique()
    alg_names = ['ADV (base MRCN)', 'ADV+ (+ multi-estimator)', 'ADV++ (+ two noisy graphs)']
    colors = ['#e74c3c', '#3498db', '#2ecc71']  # Red, Blue, Green
    
    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # 1. Relative Error Comparison
    ax = axes[0]
    x = np.arange(len(datasets))
    width = 0.25
    
    for i, (alg, color) in enumerate(zip(['alg2', 'alg3', 'alg4'], colors)):
        alg_data = df[df['alg_id'] == alg]
        errors = [alg_data[alg_data['dataset'] == d]['rel_err'].values[0] for d in datasets]
        ax.bar(x + i*width, errors, width, label=alg_names[i], color=color, alpha=0.8)
    
    ax.set_xlabel('Dataset', fontweight='bold')
    ax.set_ylabel('Relative Error', fontweight='bold')
    ax.set_title('(a) Relative Error Comparison', fontweight='bold')
    ax.set_xticks(x + width)
    ax.set_xticklabels(datasets, rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # 2. Runtime Comparison
    ax = axes[1]
    for i, (alg, color) in enumerate(zip(['alg2', 'alg3', 'alg4'], colors)):
        alg_data = df[df['alg_id'] == alg]
        times = [alg_data[alg_data['dataset'] == d]['time'].values[0] for d in datasets]
        ax.bar(x + i*width, times, width, label=alg_names[i], color=color, alpha=0.8)
    
    ax.set_xlabel('Dataset', fontweight='bold')
    ax.set_ylabel('Runtime (seconds)', fontweight='bold')
    ax.set_title('(b) Runtime Comparison', fontweight='bold')
    ax.set_xticks(x + width)
    ax.set_xticklabels(datasets, rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    # 3. Relative Improvement (normalized to ADV baseline)
    ax = axes[2]
    base_errors = df[df['alg_id'] == 'alg2'].set_index('dataset')['rel_err']
    
    improvements_alg3 = []
    improvements_alg4 = []
    
    for dataset in datasets:
        base_err = base_errors[dataset]
        alg3_err = df[(df['dataset'] == dataset) & (df['alg_id'] == 'alg3')]['rel_err'].values[0]
        alg4_err = df[(df['dataset'] == dataset) & (df['alg_id'] == 'alg4')]['rel_err'].values[0]
        
        # Improvement = (baseline - new) / baseline * 100
        improvements_alg3.append((base_err - alg3_err) / base_err * 100)
        improvements_alg4.append((base_err - alg4_err) / base_err * 100)
    
    ax.bar(x, improvements_alg3, width*1.5, label='ADV+ improvement', color=colors[1], alpha=0.8)
    ax.bar(x + width*1.5, improvements_alg4, width*1.5, label='ADV++ improvement', color=colors[2], alpha=0.8)
    
    ax.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xlabel('Dataset', fontweight='bold')
    ax.set_ylabel('Error Reduction (%)', fontweight='bold')
    ax.set_title('(c) Relative Improvement over Base MRCN', fontweight='bold')
    ax.set_xticks(x + width*0.75)
    ax.set_xticklabels(datasets, rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/ablation_comparison.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_dir}/ablation_comparison.png', dpi=300, bbox_inches='tight')
    print(f"Saved ablation comparison plot to {output_dir}/ablation_comparison.pdf")
    
    plt.show()

def create_ablation_table(df, output_dir='ablation_study_results'):
    """Create a LaTeX table for the paper"""
    
    datasets = df['dataset'].unique()
    
    latex_table = r"""\begin{table}[t]
\centering
\caption{Ablation Study: Component Contributions to MRCN Performance}
\label{tab:ablation}
\begin{tabular}{l|ccc|ccc}
\hline
\multirow{2}{*}{\textbf{Dataset}} & \multicolumn{3}{c|}{\textbf{Relative Error}} & \multicolumn{3}{c}{\textbf{Runtime (s)}} \\
& ADV & ADV+ & ADV++ & ADV & ADV+ & ADV++ \\
\hline
"""
    
    for dataset in datasets:
        row = f"{dataset}"
        
        for alg in ['alg2', 'alg3', 'alg4']:
            data = df[(df['dataset'] == dataset) & (df['alg_id'] == alg)]
            rel_err = data['rel_err'].values[0]
            row += f" & {rel_err:.4f}"
        
        for alg in ['alg2', 'alg3', 'alg4']:
            data = df[(df['dataset'] == dataset) & (df['alg_id'] == alg)]
            time_val = data['time'].values[0]
            row += f" & {time_val:.2f}"
        
        row += r" \\"
        latex_table += row + "\n"
    
    latex_table += r"""\hline
\end{tabular}
\end{table}
"""
    
    with open(f'{output_dir}/ablation_table.tex', 'w') as f:
        f.write(latex_table)
    
    print(f"Saved LaTeX table to {output_dir}/ablation_table.tex")
    
    # Also create a summary text
    print("\n" + "="*60)
    print("ABLATION STUDY SUMMARY")
    print("="*60)
    
    for dataset in datasets:
        print(f"\n{dataset}:")
        base_err = df[(df['dataset'] == dataset) & (df['alg_id'] == 'alg2')]['rel_err'].values[0]
        alg3_err = df[(df['dataset'] == dataset) & (df['alg_id'] == 'alg3')]['rel_err'].values[0]
        alg4_err = df[(df['dataset'] == dataset) & (df['alg_id'] == 'alg4')]['rel_err'].values[0]
        
        improvement_3 = (base_err - alg3_err) / base_err * 100
        improvement_4 = (base_err - alg4_err) / base_err * 100
        
        print(f"  Base (ADV) error: {base_err:.4f}")
        print(f"  ADV+ error: {alg3_err:.4f} ({improvement_3:+.1f}%)")
        print(f"  ADV++ error: {alg4_err:.4f} ({improvement_4:+.1f}%)")

def main():
    csv_path = 'ablation_study_results/ablation_summary.csv'
    
    if not Path(csv_path).exists():
        print(f"Error: {csv_path} not found!")
        print("Please run the ablation study first: ./run_ablation_study.sh")
        return
    
    # Load data
    df = load_ablation_data(csv_path)
    
    # Create visualizations
    plot_ablation_comparison(df)
    
    # Create LaTeX table
    create_ablation_table(df)

if __name__ == '__main__':
    main()

