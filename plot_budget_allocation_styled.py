#!/usr/bin/env python3
"""
Plot privacy budget allocation results in the style of Figure 6.
Creates publication-quality plots with log scale and markers.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

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
    
    Args:
        csv_files_dict: Dictionary mapping algorithm names to CSV file paths
        dataset_name: Name of the dataset for title
        ax: Matplotlib axis object
        epsilon: Total privacy budget (for title)
        p: First parameter of (p,q)-biclique
        q: Second parameter of (p,q)-biclique
    """
    # Algorithm colors and markers
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
    
    # Labels (matching legacy script font sizes)
    ax.set_xlabel(r'$\alpha_1$', fontsize=20)
    ax.set_ylabel('relative error', fontsize=20)
    ax.set_title(f'{dataset_name}, $\\epsilon = {epsilon}$, $p = {p}$, $q = {q}$', fontsize=16)
    
    # Legend
    ax.legend(loc='best', fontsize=14, ncol=1)
    
    # Grid
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    
    # Format x-axis to show alpha1 values nicely (matching legacy script)
    if alpha1_values is not None:
        ax.set_xticks(alpha1_values)
        ax.set_xticklabels([f'{v:.1f}' for v in alpha1_values], rotation=0, fontsize=20)
        ax.tick_params(axis='y', labelsize=20)
    
    return True

def main():
    # Get directory from command line or use default
    if len(sys.argv) > 1:
        results_dir = sys.argv[1]
    else:
        results_dir = 'budget_allocation_results'
    
    # Get epsilon, P and Q values from command line (optional, defaults: 1.0, 2, 2)
    epsilon = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    p = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    q = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    
    if not os.path.exists(results_dir):
        print(f"Error: Directory {results_dir} not found")
        return
    
    # Find all CSV files
    csv_files = [f for f in os.listdir(results_dir) if f.endswith('_budget_allocation.csv')]
    
    if not csv_files:
        print("No CSV files found")
        return
    
    # Dataset name mapping for better titles
    dataset_names = {
        'co': 'CO',
        'to': 'TO', 
        'unicode': 'UN',
        'lrcwiki': 'LR',
        'csbwiki': 'CS',
        'librec-filmtrust-ratings': 'LIBREC-FILMTRUST-RATINGS',
        'rmwiki': 'RM',
        'bag-kos': 'BAG-KOS',
        'bpywiki': 'BP'
    }
    
    # Algorithm name mapping
    algorithm_names = {
        'alg2': 'MRCN',
        'alg3': 'MRCN+',
        'alg4': 'MRCN++'
    }
    
    # Group CSV files by dataset
    datasets_dict = {}  # dataset_key -> {alg_name: csv_path}
    
    for csv_file in sorted(csv_files):
        # Parse filename: dataset_alg#_budget_allocation.csv or dataset_budget_allocation.csv
        if '_alg' in csv_file:
            parts = csv_file.replace('_budget_allocation.csv', '').split('_alg')
            dataset_key = parts[0]
            alg_key = 'alg' + parts[1]
        else:
            dataset_key = csv_file.replace('_budget_allocation.csv', '')
            alg_key = 'alg4'  # default to MRCN++ for old format
        
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
    # This ensures we only plot datasets that were completely processed by the test script
    complete_datasets = {}
    expected_algorithms = {'MRCN', 'MRCN+', 'MRCN++'}
    
    # Exclude datasets that are in progress (not fully verified)
    exclude_datasets = set()  # All datasets are now complete
    
    for dataset_key, alg_dict in datasets_dict.items():
        if dataset_key in exclude_datasets:
            print(f"Skipping {dataset_key}: excluded (in progress or not verified)")
            continue
            
        available_algs = set(alg_dict.keys())
        if available_algs == expected_algorithms:
            complete_datasets[dataset_key] = alg_dict
        else:
            missing = expected_algorithms - available_algs
            print(f"Skipping {dataset_key}: missing algorithms {missing}")
    
    if not complete_datasets:
        print("No complete datasets to plot (need all 3 algorithms per dataset)")
        return
    
    datasets_dict = complete_datasets
    
    # Create individual plots for each dataset (with all algorithms on one plot)
    print(f"\nGenerating individual plots for {len(datasets_dict)} datasets...")
    for dataset_key, alg_csv_dict in sorted(datasets_dict.items()):
        dataset_name = dataset_names.get(dataset_key, dataset_key.upper())
        
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        
        success = plot_dataset_multi_alg(alg_csv_dict, dataset_name, ax, epsilon=epsilon, p=p, q=q)
        
        if success:
            plt.tight_layout()
            
            # Save as PDF only with budget allocation naming including epsilon
            # Format epsilon as integer if it's a whole number, otherwise as decimal
            eps_str = f"{int(epsilon)}" if epsilon == int(epsilon) else f"{epsilon:.1f}".replace('.', 'p')
            output_pdf = os.path.join(results_dir, f'{dataset_key}-budget-allocation-eps{eps_str}-p{p}q{q}.pdf')
            
            plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
            plt.close()
            
            # Calculate statistics for each algorithm
            print(f"\n{dataset_name}:")
            for alg_name, csv_path in sorted(alg_csv_dict.items()):
                df = pd.read_csv(csv_path)
                min_error = df['relative_error'].min()
                best_alpha1 = df.loc[df['relative_error'].idxmin(), 'alpha1']
                print(f"  {alg_name}: Best α₁ = {best_alpha1:.2f}, RelErr = {min_error:.4f}")
            print(f"  Saved: {output_pdf}")
    
    # Create combined plot (grid size based on number of datasets)
    n_datasets = len(datasets_dict)
    if n_datasets == 1:
        n_rows, n_cols = 1, 1
        figsize = (6, 5)
    elif n_datasets == 2:
        n_rows, n_cols = 1, 2
        figsize = (12, 5)
    elif n_datasets <= 4:
        n_rows, n_cols = 2, 2
        figsize = (12, 10)
    elif n_datasets <= 6:
        n_rows, n_cols = 2, 3
        figsize = (18, 10)
    else:
        n_rows, n_cols = 3, 3
        figsize = (18, 15)
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_datasets == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    print(f"\nGenerating combined comparison plot...")
    for idx, (dataset_key, alg_csv_dict) in enumerate(sorted(datasets_dict.items())):
        if idx < len(axes):
            dataset_name = dataset_names.get(dataset_key, dataset_key.upper())
            plot_dataset_multi_alg(alg_csv_dict, dataset_name, axes[idx], epsilon=epsilon, p=p, q=q)
    
    # Hide unused subplots
    for idx in range(len(datasets_dict), len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    
    # Save combined plot as PDF only with epsilon in filename
    eps_str = f"{int(epsilon)}" if epsilon == int(epsilon) else f"{epsilon:.1f}".replace('.', 'p')
    output_combined_pdf = os.path.join(results_dir, f'budget-allocation-comparison-eps{eps_str}-p{p}q{q}.pdf')
    
    plt.savefig(output_combined_pdf, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nCombined plot saved: {output_combined_pdf}")
    print("\n✓ All plots generated successfully!")

if __name__ == '__main__':
    main()

