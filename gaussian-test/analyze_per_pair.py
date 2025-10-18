#!/usr/bin/env python3
"""
Analyze f_u, f_w, and f_avg distributions for each individual vertex pair.
For each pair (u,w), we have 500 samples of each estimator.
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
                match = re.match(r'Pair_(\d+)_p2_f_u:\s*(.+)', line)
                if match:
                    pair_idx = int(match.group(1))
                    values = [float(x) for x in match.group(2).split()]
                    f_u_data[pair_idx] = values
                    
            elif line.startswith('Pair_') and '_p2_f_w:' in line:
                match = re.match(r'Pair_(\d+)_p2_f_w:\s*(.+)', line)
                if match:
                    pair_idx = int(match.group(1))
                    values = [float(x) for x in match.group(2).split()]
                    f_w_data[pair_idx] = values
    
    return f_u_data, f_w_data

def analyze_single_pair(pair_idx, f_u_vals, f_w_vals):
    """Analyze a single vertex pair (u,w) with 500 samples each."""
    f_u_vals = np.array(f_u_vals)
    f_w_vals = np.array(f_w_vals)
    f_avg_vals = (f_u_vals + f_w_vals) / 2
    
    print(f"\n=== Pair {pair_idx} Analysis ===")
    print(f"Sample size: {len(f_u_vals)} for each estimator")
    
    # Statistics for f_u
    f_u_mean = np.mean(f_u_vals)
    f_u_std = np.std(f_u_vals)
    f_u_skew = stats.skew(f_u_vals)
    f_u_kurt = stats.kurtosis(f_u_vals)
    f_u_jb_p = jarque_bera(f_u_vals)[1]
    
    print(f"f_u: mean={f_u_mean:.3f}, std={f_u_std:.3f}, skew={f_u_skew:.3f}, kurt={f_u_kurt:.3f}, JB_p={f_u_jb_p:.6f}")
    
    # Statistics for f_w
    f_w_mean = np.mean(f_w_vals)
    f_w_std = np.std(f_w_vals)
    f_w_skew = stats.skew(f_w_vals)
    f_w_kurt = stats.kurtosis(f_w_vals)
    f_w_jb_p = jarque_bera(f_w_vals)[1]
    
    print(f"f_w: mean={f_w_mean:.3f}, std={f_w_std:.3f}, skew={f_w_skew:.3f}, kurt={f_w_kurt:.3f}, JB_p={f_w_jb_p:.6f}")
    
    # Statistics for f_avg
    f_avg_mean = np.mean(f_avg_vals)
    f_avg_std = np.std(f_avg_vals)
    f_avg_skew = stats.skew(f_avg_vals)
    f_avg_kurt = stats.kurtosis(f_avg_vals)
    f_avg_jb_p = jarque_bera(f_avg_vals)[1]
    
    print(f"f_avg: mean={f_avg_mean:.3f}, std={f_avg_std:.3f}, skew={f_avg_skew:.3f}, kurt={f_avg_kurt:.3f}, JB_p={f_avg_jb_p:.6f}")
    
    # Interpretation
    print(f"Gaussian assumption for f_avg: {'ACCEPTED' if f_avg_jb_p > 0.05 else 'REJECTED'} (p={f_avg_jb_p:.6f})")
    
    return {
        'f_u': {'mean': f_u_mean, 'std': f_u_std, 'skew': f_u_skew, 'kurt': f_u_kurt, 'jb_p': f_u_jb_p},
        'f_w': {'mean': f_w_mean, 'std': f_w_std, 'skew': f_w_skew, 'kurt': f_w_kurt, 'jb_p': f_w_jb_p},
        'f_avg': {'mean': f_avg_mean, 'std': f_avg_std, 'skew': f_avg_skew, 'kurt': f_avg_kurt, 'jb_p': f_avg_jb_p}
    }

def plot_single_pair(pair_idx, f_u_vals, f_w_vals, f_avg_vals, output_dir):
    """Create plots for a single vertex pair."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f'Vertex Pair {pair_idx} Analysis (500 samples each)', fontsize=16, fontweight='bold')
    
    # Top row: Histograms
    axes[0,0].hist(f_u_vals, bins=30, alpha=0.7, color='blue', density=True)
    axes[0,0].set_title(f'f_u Distribution (n={len(f_u_vals)})')
    axes[0,0].set_xlabel('f_u values')
    axes[0,0].set_ylabel('Density')
    axes[0,0].grid(True, alpha=0.3)
    
    axes[0,1].hist(f_w_vals, bins=30, alpha=0.7, color='red', density=True)
    axes[0,1].set_title(f'f_w Distribution (n={len(f_w_vals)})')
    axes[0,1].set_xlabel('f_w values')
    axes[0,1].set_ylabel('Density')
    axes[0,1].grid(True, alpha=0.3)
    
    axes[0,2].hist(f_avg_vals, bins=30, alpha=0.7, color='green', density=True)
    axes[0,2].set_title(f'f_avg Distribution (n={len(f_avg_vals)})')
    axes[0,2].set_xlabel('f_avg values')
    axes[0,2].set_ylabel('Density')
    axes[0,2].grid(True, alpha=0.3)
    
    # Bottom row: Q-Q plots
    stats.probplot(f_u_vals, dist="norm", plot=axes[1,0])
    axes[1,0].set_title('f_u Q-Q Plot vs Normal')
    axes[1,0].grid(True, alpha=0.3)
    
    stats.probplot(f_w_vals, dist="norm", plot=axes[1,1])
    axes[1,1].set_title('f_w Q-Q Plot vs Normal')
    axes[1,1].grid(True, alpha=0.3)
    
    stats.probplot(f_avg_vals, dist="norm", plot=axes[1,2])
    axes[1,2].set_title('f_avg Q-Q Plot vs Normal')
    axes[1,2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'pair_{pair_idx}_analysis.png'), dpi=300, bbox_inches='tight')
    plt.close()

def main():
    """Main analysis function."""
    print("=== Per-Pair Analysis for Gaussian Assumption ===")
    print("For each vertex pair (u,w): 500 f_u, 500 f_w, 500 f_avg samples")
    
    # Parse the C++ output
    data_file = '/data/yizhangh/ldp-pq/gaussian-test/f_distribution_unicode_500.txt'
    f_u_data, f_w_data = parse_cpp_output(data_file)
    
    print(f"Loaded data for {len(f_u_data)} pairs")
    
    # Create output directory
    output_dir = '/data/yizhangh/ldp-pq/gaussian-test'
    os.makedirs(output_dir, exist_ok=True)
    
    # Analyze each pair individually
    all_results = []
    
    for pair_idx in sorted(f_u_data.keys()):
        if pair_idx in f_w_data:
            f_u_vals = f_u_data[pair_idx]
            f_w_vals = f_w_data[pair_idx]
            f_avg_vals = [(f_u + f_w) / 2 for f_u, f_w in zip(f_u_vals, f_w_vals)]
            
            # Analyze this pair
            results = analyze_single_pair(pair_idx, f_u_vals, f_w_vals)
            all_results.append((pair_idx, results))
            
            # Create individual plot
            plot_single_pair(pair_idx, f_u_vals, f_w_vals, f_avg_vals, output_dir)
    
    # Summary across all pairs
    print("\n" + "="*80)
    print("SUMMARY ACROSS ALL PAIRS")
    print("="*80)
    
    f_avg_accepted = sum(1 for _, results in all_results if results['f_avg']['jb_p'] > 0.05)
    f_avg_rejected = len(all_results) - f_avg_accepted
    
    print(f"f_avg Gaussian assumption:")
    print(f"  Accepted: {f_avg_accepted}/{len(all_results)} pairs")
    print(f"  Rejected: {f_avg_rejected}/{len(all_results)} pairs")
    
    # Average statistics
    avg_f_avg_skew = np.mean([results['f_avg']['skew'] for _, results in all_results])
    avg_f_avg_kurt = np.mean([results['f_avg']['kurt'] for _, results in all_results])
    
    print(f"\nAverage f_avg statistics across all pairs:")
    print(f"  Mean skewness: {avg_f_avg_skew:.4f}")
    print(f"  Mean kurtosis: {avg_f_avg_kurt:.4f}")
    
    print(f"\nIndividual pair plots saved to: {output_dir}/pair_X_analysis.png")

if __name__ == "__main__":
    main()
