#!/usr/bin/env python3
"""
Analyze f_u and f_w distributions from C++ output to test Gaussian assumption.
This addresses Reviewer Q3 about why f_avg = (f_u + f_w)/2 can be assumed Gaussian.
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

def compute_statistics(values, name):
    """Compute statistical measures for a distribution."""
    values = np.array(values)
    
    stats_dict = {
        'mean': np.mean(values),
        'std': np.std(values),
        'skewness': stats.skew(values),
        'kurtosis': stats.kurtosis(values),
        'jarque_bera_stat': jarque_bera(values)[0],
        'jarque_bera_pvalue': jarque_bera(values)[1],
        'count': len(values)
    }
    
    print(f"\n{name} Statistics:")
    print(f"  Count: {stats_dict['count']}")
    print(f"  Mean: {stats_dict['mean']:.4f}")
    print(f"  Std: {stats_dict['std']:.4f}")
    print(f"  Skewness: {stats_dict['skewness']:.4f}")
    print(f"  Kurtosis: {stats_dict['kurtosis']:.4f}")
    print(f"  Jarque-Bera p-value: {stats_dict['jarque_bera_pvalue']:.6f}")
    
    return stats_dict

def plot_distributions(f_u_data, f_w_data, output_dir):
    """Create comprehensive distribution plots."""
    
    # Collect all f_u and f_w values
    all_f_u = []
    all_f_w = []
    all_f_avg = []
    
    for pair_idx in f_u_data:
        if pair_idx in f_w_data:
            f_u_vals = f_u_data[pair_idx]
            f_w_vals = f_w_data[pair_idx]
            f_avg_vals = [(f_u + f_w) / 2 for f_u, f_w in zip(f_u_vals, f_w_vals)]
            
            all_f_u.extend(f_u_vals)
            all_f_w.extend(f_w_vals)
            all_f_avg.extend(f_avg_vals)
    
    all_f_u = np.array(all_f_u)
    all_f_w = np.array(all_f_w)
    all_f_avg = np.array(all_f_avg)
    
    # Create 6-panel plot
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Distribution Analysis for Gaussian Assumption (Reviewer Q3)', fontsize=16, fontweight='bold')
    
    # Panel 1: f_u histogram
    axes[0,0].hist(all_f_u, bins=50, alpha=0.7, color='blue', density=True, label='f_u')
    axes[0,0].set_title('f_u Distribution')
    axes[0,0].set_xlabel('f_u values')
    axes[0,0].set_ylabel('Density')
    axes[0,0].grid(True, alpha=0.3)
    
    # Panel 2: f_w histogram  
    axes[0,1].hist(all_f_w, bins=50, alpha=0.7, color='red', density=True, label='f_w')
    axes[0,1].set_title('f_w Distribution')
    axes[0,1].set_xlabel('f_w values')
    axes[0,1].set_ylabel('Density')
    axes[0,1].grid(True, alpha=0.3)
    
    # Panel 3: f_avg histogram
    axes[0,2].hist(all_f_avg, bins=50, alpha=0.7, color='green', density=True, label='f_avg')
    axes[0,2].set_title('f_avg = (f_u + f_w)/2 Distribution')
    axes[0,2].set_xlabel('f_avg values')
    axes[0,2].set_ylabel('Density')
    axes[0,2].grid(True, alpha=0.3)
    
    # Panel 4: f_u Q-Q plot
    stats.probplot(all_f_u, dist="norm", plot=axes[1,0])
    axes[1,0].set_title('f_u Q-Q Plot vs Normal')
    axes[1,0].grid(True, alpha=0.3)
    
    # Panel 5: f_w Q-Q plot
    stats.probplot(all_f_w, dist="norm", plot=axes[1,1])
    axes[1,1].set_title('f_w Q-Q Plot vs Normal')
    axes[1,1].grid(True, alpha=0.3)
    
    # Panel 6: f_avg Q-Q plot
    stats.probplot(all_f_avg, dist="norm", plot=axes[1,2])
    axes[1,2].set_title('f_avg Q-Q Plot vs Normal')
    axes[1,2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'f_distribution_analysis.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    return all_f_u, all_f_w, all_f_avg

def main():
    """Main analysis function."""
    print("=== Gaussian Assumption Analysis for Reviewer Q3 ===")
    print("Testing why f_avg = (f_u + f_w)/2 can be assumed Gaussian")
    
    # Parse the C++ output
    data_file = '/data/yizhangh/ldp-pq/gaussian-test/f_distribution_unicode_500.txt'
    f_u_data, f_w_data = parse_cpp_output(data_file)
    
    print(f"Loaded data for {len(f_u_data)} pairs")
    
    # Create output directory
    output_dir = '/data/yizhangh/ldp-pq/gaussian-test'
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate plots and get combined data
    all_f_u, all_f_w, all_f_avg = plot_distributions(f_u_data, f_w_data, output_dir)
    
    # Compute statistics
    f_u_stats = compute_statistics(all_f_u, "f_u")
    f_w_stats = compute_statistics(all_f_w, "f_w") 
    f_avg_stats = compute_statistics(all_f_avg, "f_avg (Gaussian assumption)")
    
    # Summary for Reviewer Q3
    print("\n" + "="*60)
    print("SUMMARY FOR REVIEWER Q3")
    print("="*60)
    print("Question: Why can we assume f_avg = (f_u + f_w)/2 is approximately Gaussian?")
    print()
    print("Key Findings:")
    print(f"1. f_avg has Jarque-Bera p-value: {f_avg_stats['jarque_bera_pvalue']:.6f}")
    if f_avg_stats['jarque_bera_pvalue'] > 0.05:
        print("   → p > 0.05: Cannot reject Gaussian hypothesis")
    else:
        print("   → p ≤ 0.05: Evidence against Gaussian hypothesis")
    
    print(f"2. f_avg skewness: {f_avg_stats['skewness']:.4f} (0 = symmetric)")
    print(f"3. f_avg kurtosis: {f_avg_stats['kurtosis']:.4f} (0 = normal kurtosis)")
    print(f"4. Sample size: {f_avg_stats['count']} observations")
    
    print("\nTheoretical Justification:")
    print("- f_u and f_w are both sums of many independent random variables (edge noise)")
    print("- By Central Limit Theorem, both f_u and f_w should be approximately normal")
    print("- f_avg = (f_u + f_w)/2 is a linear combination of two approximately normal variables")
    print("- Linear combinations of normal variables are also normal")
    
    print(f"\nPlots saved to: {output_dir}/f_distribution_analysis.png")

if __name__ == "__main__":
    main()
