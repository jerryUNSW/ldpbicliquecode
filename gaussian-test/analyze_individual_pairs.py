#!/usr/bin/env python3
"""
Analyze f_u and f_w distributions for each individual vertex pair.
This provides better insight into the Gaussian assumption per pair.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import jarque_bera
import re
import os

def parse_cpp_output(filename):
    """Parse the C++ output file to extract f_u and f_w values."""
    f_u_data = {}  # pair_idx -> list of f_u values
    f_w_data = {}  # pair_idx -> list of f_w values
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Pair_') and '_p2_f_u:' in line:
                # Parse: Pair_0_p2_f_u: val1 val2 val3 ...
                match = re.match(r'Pair_(\d+)_p2_f_u:\s*(.+)', line)
                if match:
                    pair_idx = int(match.group(1))
                    values = [float(x) for x in match.group(2).split()]
                    f_u_data[pair_idx] = values
                    
            elif line.startswith('Pair_') and '_p2_f_w:' in line:
                # Parse: Pair_0_p2_f_w: val1 val2 val3 ...
                match = re.match(r'Pair_(\d+)_p2_f_w:\s*(.+)', line)
                if match:
                    pair_idx = int(match.group(1))
                    values = [float(x) for x in match.group(2).split()]
                    f_w_data[pair_idx] = values
    
    return f_u_data, f_w_data

def analyze_individual_pair(pair_idx, f_u_vals, f_w_vals):
    """Analyze a single vertex pair."""
    f_u_vals = np.array(f_u_vals)
    f_w_vals = np.array(f_w_vals)
    f_avg_vals = (f_u_vals + f_w_vals) / 2
    
    # Compute statistics
    f_u_stats = {
        'mean': np.mean(f_u_vals),
        'std': np.std(f_u_vals),
        'skewness': stats.skew(f_u_vals),
        'kurtosis': stats.kurtosis(f_u_vals),
        'jb_pvalue': jarque_bera(f_u_vals)[1]
    }
    
    f_w_stats = {
        'mean': np.mean(f_w_vals),
        'std': np.std(f_w_vals),
        'skewness': stats.skew(f_w_vals),
        'kurtosis': stats.kurtosis(f_w_vals),
        'jb_pvalue': jarque_bera(f_w_vals)[1]
    }
    
    f_avg_stats = {
        'mean': np.mean(f_avg_vals),
        'std': np.std(f_avg_vals),
        'skewness': stats.skew(f_avg_vals),
        'kurtosis': stats.kurtosis(f_avg_vals),
        'jb_pvalue': jarque_bera(f_avg_vals)[1]
    }
    
    return f_u_stats, f_w_stats, f_avg_stats, f_avg_vals

def plot_individual_pairs(f_u_data, f_w_data, output_dir):
    """Create plots for each individual vertex pair."""
    
    # Create a large subplot grid for all pairs
    n_pairs = len(f_u_data)
    fig, axes = plt.subplots(n_pairs, 3, figsize=(18, 4*n_pairs))
    if n_pairs == 1:
        axes = axes.reshape(1, -1)
    
    fig.suptitle('Individual Vertex Pair Analysis for Gaussian Assumption', fontsize=16, fontweight='bold')
    
    all_f_avg_stats = []
    
    for pair_idx in sorted(f_u_data.keys()):
        if pair_idx in f_w_data:
            f_u_vals = f_u_data[pair_idx]
            f_w_vals = f_w_data[pair_idx]
            f_avg_vals = [(f_u + f_w) / 2 for f_u, f_w in zip(f_u_vals, f_w_vals)]
            
            # Analyze this pair
            f_u_stats, f_w_stats, f_avg_stats, f_avg_vals = analyze_individual_pair(
                pair_idx, f_u_vals, f_w_vals
            )
            all_f_avg_stats.append(f_avg_stats)
            
            # Plot for this pair
            row = pair_idx
            
            # f_u histogram
            axes[row, 0].hist(f_u_vals, bins=30, alpha=0.7, color='blue', density=True)
            axes[row, 0].set_title(f'Pair {pair_idx}: f_u (JB p={f_u_stats["jb_pvalue"]:.3f})')
            axes[row, 0].set_xlabel('f_u values')
            axes[row, 0].set_ylabel('Density')
            axes[row, 0].grid(True, alpha=0.3)
            
            # f_w histogram
            axes[row, 1].hist(f_w_vals, bins=30, alpha=0.7, color='red', density=True)
            axes[row, 1].set_title(f'Pair {pair_idx}: f_w (JB p={f_w_stats["jb_pvalue"]:.3f})')
            axes[row, 1].set_xlabel('f_w values')
            axes[row, 1].set_ylabel('Density')
            axes[row, 1].grid(True, alpha=0.3)
            
            # f_avg histogram
            axes[row, 2].hist(f_avg_vals, bins=30, alpha=0.7, color='green', density=True)
            axes[row, 2].set_title(f'Pair {pair_idx}: f_avg (JB p={f_avg_stats["jb_pvalue"]:.3f})')
            axes[row, 2].set_xlabel('f_avg values')
            axes[row, 2].set_ylabel('Density')
            axes[row, 2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'individual_pairs_analysis.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    return all_f_avg_stats

def plot_qq_plots(f_u_data, f_w_data, output_dir):
    """Create Q-Q plots for each pair."""
    n_pairs = len(f_u_data)
    fig, axes = plt.subplots(n_pairs, 3, figsize=(18, 4*n_pairs))
    if n_pairs == 1:
        axes = axes.reshape(1, -1)
    
    fig.suptitle('Q-Q Plots for Individual Vertex Pairs', fontsize=16, fontweight='bold')
    
    for pair_idx in sorted(f_u_data.keys()):
        if pair_idx in f_w_data:
            f_u_vals = f_u_data[pair_idx]
            f_w_vals = f_w_data[pair_idx]
            f_avg_vals = [(f_u + f_w) / 2 for f_u, f_w in zip(f_u_vals, f_w_vals)]
            
            row = pair_idx
            
            # f_u Q-Q plot
            stats.probplot(f_u_vals, dist="norm", plot=axes[row, 0])
            axes[row, 0].set_title(f'Pair {pair_idx}: f_u Q-Q Plot')
            axes[row, 0].grid(True, alpha=0.3)
            
            # f_w Q-Q plot
            stats.probplot(f_w_vals, dist="norm", plot=axes[row, 1])
            axes[row, 1].set_title(f'Pair {pair_idx}: f_w Q-Q Plot')
            axes[row, 1].grid(True, alpha=0.3)
            
            # f_avg Q-Q plot
            stats.probplot(f_avg_vals, dist="norm", plot=axes[row, 2])
            axes[row, 2].set_title(f'Pair {pair_idx}: f_avg Q-Q Plot')
            axes[row, 2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'individual_pairs_qq.png'), dpi=300, bbox_inches='tight')
    plt.close()

def print_individual_statistics(all_f_avg_stats):
    """Print statistics for each pair."""
    print("\n" + "="*80)
    print("INDIVIDUAL VERTEX PAIR STATISTICS")
    print("="*80)
    
    for i, stats in enumerate(all_f_avg_stats):
        print(f"\nPair {i}:")
        print(f"  Mean: {stats['mean']:.4f}")
        print(f"  Std: {stats['std']:.4f}")
        print(f"  Skewness: {stats['skewness']:.4f}")
        print(f"  Kurtosis: {stats['kurtosis']:.4f}")
        print(f"  Jarque-Bera p-value: {stats['jb_pvalue']:.6f}")
        
        # Interpretation
        if stats['jb_pvalue'] > 0.05:
            print(f"  → Gaussian assumption ACCEPTED (p > 0.05)")
        else:
            print(f"  → Gaussian assumption REJECTED (p ≤ 0.05)")
            
        if abs(stats['skewness']) < 0.5:
            print(f"  → Skewness is small (|{stats['skewness']:.3f}| < 0.5)")
        else:
            print(f"  → Skewness is moderate/large (|{stats['skewness']:.3f}| ≥ 0.5)")

def main():
    """Main analysis function."""
    print("=== Individual Vertex Pair Analysis for Gaussian Assumption ===")
    
    # Parse the C++ output
    data_file = '/data/yizhangh/ldp-pq/gaussian-test/f_distribution_unicode_500.txt'
    f_u_data, f_w_data = parse_cpp_output(data_file)
    
    print(f"Loaded data for {len(f_u_data)} pairs")
    
    # Create output directory
    output_dir = '/data/yizhangh/ldp-pq/gaussian-test'
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate individual pair plots
    all_f_avg_stats = plot_individual_pairs(f_u_data, f_w_data, output_dir)
    
    # Generate Q-Q plots
    plot_qq_plots(f_u_data, f_w_data, output_dir)
    
    # Print individual statistics
    print_individual_statistics(all_f_avg_stats)
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    accepted_pairs = sum(1 for stats in all_f_avg_stats if stats['jb_pvalue'] > 0.05)
    total_pairs = len(all_f_avg_stats)
    
    print(f"Pairs with Gaussian assumption accepted: {accepted_pairs}/{total_pairs}")
    print(f"Pairs with Gaussian assumption rejected: {total_pairs - accepted_pairs}/{total_pairs}")
    
    # Average statistics
    avg_skewness = np.mean([stats['skewness'] for stats in all_f_avg_stats])
    avg_kurtosis = np.mean([stats['kurtosis'] for stats in all_f_avg_stats])
    
    print(f"\nAverage across all pairs:")
    print(f"  Mean skewness: {avg_skewness:.4f}")
    print(f"  Mean kurtosis: {avg_kurtosis:.4f}")
    
    print(f"\nPlots saved to:")
    print(f"  {output_dir}/individual_pairs_analysis.png")
    print(f"  {output_dir}/individual_pairs_qq.png")

if __name__ == "__main__":
    main()
