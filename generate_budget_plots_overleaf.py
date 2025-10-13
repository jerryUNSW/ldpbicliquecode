#!/usr/bin/env python3
"""
Generate budget allocation plots for Overleaf in the style of Figure 9.
Creates plots for ε=1.0, p=2, q=3 and ε=2.0, p=2, q=2 configurations.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import shutil

# Set publication-quality plot parameters (matching legacy scripts)
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 300
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams["legend.framealpha"] = 0
plt.rcParams["legend.handletextpad"] = 0.1
plt.rcParams["legend.columnspacing"] = 0.2
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 8

def plot_dataset_multi_alg(csv_files_dict, dataset_name, ax, epsilon=1.0, p=2, q=2):
    """
    Plot multiple algorithms for a single dataset on the given axis.
    Matches the style of Figure 9 (budget allocation analysis).
    """
    # Algorithm colors and markers (matching Figure 9 style)
    alg_styles = {
        'MRCN': {'color': 'blue', 'marker': 'o', 'linestyle': '-'},
        'MRCN+': {'color': 'green', 'marker': 's', 'linestyle': '-'},
        'MRCN++': {'color': 'red', 'marker': '^', 'linestyle': '-'}
    }
    
    has_data = False
    alpha1_values = None
    
    # Plot each algorithm
    for alg_name, csv_file in sorted(csv_files_dict.items()):
        df = pd.read_csv(csv_file)
        
        if len(df) == 0:
            print(f"Warning: No data for {dataset_name} - {alg_name}")
            continue
        
        has_data = True
        alpha1_values = df['alpha1'].values
        rel_errors = df['relative_error'].values
        
        style = alg_styles.get(alg_name, {'color': 'black', 'marker': 'o', 'linestyle': '-'})
        
        # Plot with algorithm-specific style
        ax.plot(alpha1_values, rel_errors, 
                marker=style['marker'], 
                markerfacecolor='white', 
                markeredgecolor=style['color'],
                markeredgewidth=1.5, 
                color=style['color'], 
                linewidth=2,
                linestyle=style['linestyle'],
                label=alg_name)
    
    if not has_data:
        return False
    
    # Set log scale on y-axis
    ax.set_yscale('log')
    
    # Labels (matching Figure 9 style)
    ax.set_xlabel(r'$\alpha_1$', fontsize=20)
    ax.set_ylabel('relative error', fontsize=20)
    
    # Title format matching Figure 9
    if epsilon == 1.0:
        title = f'\\texttt{{{dataset_name}}}'
    else:
        title = f'\\texttt{{{dataset_name}}}'
    
    ax.set_title(title, fontsize=16)
    
    # Legend
    ax.legend(loc='best', fontsize=14, ncol=1)
    
    # Grid
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    
    # Format x-axis to show alpha1 values nicely
    if alpha1_values is not None:
        ax.set_xticks(alpha1_values)
        ax.set_xticklabels([f'{v:.1f}' for v in alpha1_values], rotation=0, fontsize=20)
        ax.tick_params(axis='y', labelsize=20)
    
    return True

def generate_plots_for_config(results_dir, epsilon, p, q, overleaf_dir):
    """
    Generate plots for a specific configuration and copy to Overleaf directory.
    """
    print(f"\n{'='*60}")
    print(f"Processing: ε={epsilon}, p={p}, q={q}")
    print(f"{'='*60}")
    
    if not os.path.exists(results_dir):
        print(f"Error: Directory {results_dir} not found")
        return False
    
    # Find all CSV files
    csv_files = [f for f in os.listdir(results_dir) if f.endswith('_budget_allocation.csv')]
    
    if not csv_files:
        print("No CSV files found")
        return False
    
    # Dataset name mapping (matching Figure 9 datasets)
    dataset_names = {
        'bag-kos': 'bag-kos',
        'unicode': 'unicode',
        'lrcwiki': 'lrcwiki',
        'csbwiki': 'csbwiki',
        'librec-filmtrust-ratings': 'librec-filmtrust-ratings',
        'rmwiki': 'rmwiki',
        'bpywiki': 'bpywiki',
        'co': 'co',
        'to': 'to'
    }
    
    # Algorithm name mapping
    algorithm_names = {
        'alg2': 'MRCN',
        'alg3': 'MRCN+',
        'alg4': 'MRCN++'
    }
    
    # Group CSV files by dataset
    datasets_dict = {}
    
    for csv_file in sorted(csv_files):
        # Parse filename: dataset_alg#_budget_allocation.csv
        if '_alg' in csv_file:
            parts = csv_file.replace('_budget_allocation.csv', '').split('_alg')
            dataset_key = parts[0]
            alg_key = 'alg' + parts[1]
        else:
            continue
        
        csv_path = os.path.join(results_dir, csv_file)
        
        # Check if dataset has complete data
        df = pd.read_csv(csv_path)
        if len(df) == 0:
            print(f"Skipping {dataset_key} - {alg_key}: no data")
            continue
        
        alg_name = algorithm_names.get(alg_key, alg_key.upper())
        
        # Add to datasets dictionary
        if dataset_key not in datasets_dict:
            datasets_dict[dataset_key] = {}
        datasets_dict[dataset_key][alg_name] = csv_path
    
    # Filter to only include datasets that have all three algorithms
    complete_datasets = {}
    expected_algorithms = {'MRCN', 'MRCN+', 'MRCN++'}
    
    for dataset_key, alg_dict in datasets_dict.items():
        available_algs = set(alg_dict.keys())
        if available_algs == expected_algorithms:
            complete_datasets[dataset_key] = alg_dict
        else:
            missing = expected_algorithms - available_algs
            print(f"Skipping {dataset_key}: missing algorithms {missing}")
    
    if not complete_datasets:
        print("No complete datasets to plot (need all 3 algorithms per dataset)")
        return False
    
    # Focus on the datasets shown in Figure 9
    figure9_datasets = ['bag-kos', 'unicode', 'lrcwiki', 'csbwiki', 'librec-filmtrust-ratings', 'rmwiki']
    
    # Create individual plots for Figure 9 datasets
    print(f"\nGenerating plots for Figure 9 datasets...")
    generated_files = []
    
    for dataset_key in figure9_datasets:
        if dataset_key not in complete_datasets:
            print(f"Skipping {dataset_key}: not available in this configuration")
            continue
            
        alg_csv_dict = complete_datasets[dataset_key]
        dataset_name = dataset_names.get(dataset_key, dataset_key)
        
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        
        success = plot_dataset_multi_alg(alg_csv_dict, dataset_name, ax, epsilon=epsilon, p=p, q=q)
        
        if success:
            plt.tight_layout()
            
            # Generate filename matching Figure 9 format
            eps_str = f"{int(epsilon)}" if epsilon == int(epsilon) else f"{epsilon:.1f}"
            output_filename = f'{dataset_key}-budget-allocation-eps{eps_str}-p{p}q{q}.pdf'
            output_pdf = os.path.join(overleaf_dir, output_filename)
            
            plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
            plt.close()
            
            generated_files.append(output_filename)
            
            # Calculate statistics for each algorithm
            print(f"\n{dataset_name}:")
            for alg_name, csv_path in sorted(alg_csv_dict.items()):
                df = pd.read_csv(csv_path)
                min_error = df['relative_error'].min()
                best_alpha1 = df.loc[df['relative_error'].idxmin(), 'alpha1']
                print(f"  {alg_name}: Best α₁ = {best_alpha1:.2f}, RelErr = {min_error:.4f}")
            print(f"  Generated: {output_filename}")
    
    print(f"\n✓ Generated {len(generated_files)} plots for Overleaf")
    return True

def main():
    # Overleaf directory for budget allocation plots
    overleaf_dir = "overleaf-paper/revision-plots"
    
    # Create overleaf directory if it doesn't exist
    os.makedirs(overleaf_dir, exist_ok=True)
    
    # Configuration 2: ε=1.0, p=2, q=3
    config2_dir = "budget_allocation_results_eps1.0_p2q3"
    success2 = generate_plots_for_config(config2_dir, 1.0, 2, 3, overleaf_dir)
    
    # Configuration 3: ε=2.0, p=2, q=2  
    config3_dir = "budget_allocation_results_eps2.0_p2q2"
    success3 = generate_plots_for_config(config3_dir, 2.0, 2, 2, overleaf_dir)
    
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    
    if success2:
        print("✅ ε=1.0, p=2, q=3: Plots generated for Overleaf")
    else:
        print("❌ ε=1.0, p=2, q=3: Failed to generate plots")
        
    if success3:
        print("✅ ε=2.0, p=2, q=2: Plots generated for Overleaf")
    else:
        print("❌ ε=2.0, p=2, q=2: Failed to generate plots")
    
    if success2 or success3:
        print(f"\n📁 Plots saved to: {overleaf_dir}")
        print("📝 Ready for inclusion in Overleaf Figure 9 style")

if __name__ == '__main__':
    main()
