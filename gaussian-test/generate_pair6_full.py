#!/usr/bin/env python3
"""
Generate full 6-subplot analysis figure for pair 6 (for response letter).
Based on analyze_with_common_neighbors.py but outputs PDF format.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import jarque_bera
import re
import os

def parse_cpp_output(filename):
    """Parse the C++ output file to extract f_u, f_w values, common neighbor counts, and degrees."""
    f_u_data = {}
    f_w_data = {}
    common_neighbors = {}
    degree_u = {}
    degree_w = {}
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Pair_') and '_real_common_neighbors:' in line:
                match = re.match(r'Pair_(\d+)_real_common_neighbors:\s*(\d+)', line)
                if match:
                    pair_idx = int(match.group(1))
                    count = int(match.group(2))
                    common_neighbors[pair_idx] = count
                    
            elif line.startswith('Pair_') and '_degree_u:' in line:
                match = re.match(r'Pair_(\d+)_degree_u:\s*(\d+)', line)
                if match:
                    pair_idx = int(match.group(1))
                    degree = int(match.group(2))
                    degree_u[pair_idx] = degree
                    
            elif line.startswith('Pair_') and '_degree_w:' in line:
                match = re.match(r'Pair_(\d+)_degree_w:\s*(\d+)', line)
                if match:
                    pair_idx = int(match.group(1))
                    degree = int(match.group(2))
                    degree_w[pair_idx] = degree
                    
            elif line.startswith('Pair_') and '_p2_f_u:' in line:
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
    
    return f_u_data, f_w_data, common_neighbors, degree_u, degree_w

def plot_pair6_full(f_u_vals, f_w_vals, f_avg_vals, real_common_neighbors, degree_u, degree_w, output_file):
    """Create full 6-subplot figure for pair 6."""
    # Convert to numpy arrays
    f_u_vals = np.array(f_u_vals)
    f_w_vals = np.array(f_w_vals)
    f_avg_vals = np.array(f_avg_vals)
    
    # Compute normality scores
    f_u_jb_p = jarque_bera(f_u_vals)[1]
    f_w_jb_p = jarque_bera(f_w_vals)[1]
    f_avg_jb_p = jarque_bera(f_avg_vals)[1]
    
    fig, axes = plt.subplots(1, 6, figsize=(18, 4))
    
    # All 6 plots in one row: 3 histograms + 3 Q-Q plots
    # f_u histogram
    axes[0].hist(f_u_vals, bins=20, alpha=0.7, color='blue', density=True, label='$f_u$ data')
    axes[0].axvline(x=real_common_neighbors, color='red', linestyle='--', linewidth=2, label=f'Real={real_common_neighbors}')
    
    # Overlay normal distribution
    f_u_mean = np.mean(f_u_vals)
    f_u_std = np.std(f_u_vals)
    x_range = np.linspace(f_u_vals.min(), f_u_vals.max(), 100)
    normal_curve = stats.norm.pdf(x_range, f_u_mean, f_u_std)
    axes[0].plot(x_range, normal_curve, 'k-', linewidth=2, label=f'N({f_u_mean:.1f}, {f_u_std:.1f})')
    
    axes[0].set_title(f'$f_u$ (deg={degree_u}, Jarque-Bera score={f_u_jb_p:.3f})', fontsize=10)
    axes[0].set_xlabel('$f_u$ values', fontsize=10)
    axes[0].set_ylabel('Density', fontsize=8)
    axes[0].legend(loc='upper left', frameon=False, fontsize=7)
    axes[0].grid(True, alpha=0.3)
    # Set y-axis limit to avoid legend interference
    axes[0].set_ylim(0, max(normal_curve) * 1.3)
    
    # f_w histogram
    axes[1].hist(f_w_vals, bins=20, alpha=0.7, color='red', density=True, label='$f_w$ data')
    axes[1].axvline(x=real_common_neighbors, color='blue', linestyle='--', linewidth=2, label=f'Real={real_common_neighbors}')
    
    # Overlay normal distribution
    f_w_mean = np.mean(f_w_vals)
    f_w_std = np.std(f_w_vals)
    x_range = np.linspace(f_w_vals.min(), f_w_vals.max(), 100)
    normal_curve = stats.norm.pdf(x_range, f_w_mean, f_w_std)
    axes[1].plot(x_range, normal_curve, 'k-', linewidth=2, label=f'N({f_w_mean:.1f}, {f_w_std:.1f})')
    
    axes[1].set_title(f'$f_w$ (deg={degree_w}, Jarque-Bera score={f_w_jb_p:.3f})', fontsize=10)
    axes[1].set_xlabel('$f_w$ values', fontsize=10)
    axes[1].set_ylabel('Density', fontsize=8)
    axes[1].legend(loc='upper left', frameon=False, fontsize=7)
    axes[1].grid(True, alpha=0.3)
    # Set y-axis limit to avoid legend interference
    axes[1].set_ylim(0, max(normal_curve) * 1.3)
    
    # f_avg histogram
    axes[2].hist(f_avg_vals, bins=20, alpha=0.7, color='green', density=True, label='$f_{avg}$ data')
    axes[2].axvline(x=real_common_neighbors, color='red', linestyle='--', linewidth=2, label=f'Real={real_common_neighbors}')
    
    # Overlay normal distribution
    f_avg_mean = np.mean(f_avg_vals)
    f_avg_std = np.std(f_avg_vals)
    x_range = np.linspace(f_avg_vals.min(), f_avg_vals.max(), 100)
    normal_curve = stats.norm.pdf(x_range, f_avg_mean, f_avg_std)
    axes[2].plot(x_range, normal_curve, 'k-', linewidth=2, label=f'N({f_avg_mean:.1f}, {f_avg_std:.1f})')
    
    axes[2].set_title(f'$f_{{avg}}$ (Jarque-Bera score={f_avg_jb_p:.3f})', fontsize=10)
    axes[2].set_xlabel('$f_{avg}$ values', fontsize=10)
    axes[2].set_ylabel('Density', fontsize=8)
    axes[2].legend(loc='upper left', frameon=False, fontsize=7)
    axes[2].grid(True, alpha=0.3)
    # Set y-axis limit to avoid legend interference
    axes[2].set_ylim(0, max(normal_curve) * 1.3)
    
    # Q-Q plots
    stats.probplot(f_u_vals, dist="norm", plot=axes[3])
    axes[3].set_title('$f_u$ Q-Q Plot', fontsize=10)
    axes[3].set_xlabel('Theoretical Quantiles', fontsize=10)
    axes[3].set_ylabel('Sample Quantiles', fontsize=8)
    axes[3].grid(True, alpha=0.3)
    
    stats.probplot(f_w_vals, dist="norm", plot=axes[4])
    axes[4].set_title('$f_w$ Q-Q Plot', fontsize=10)
    axes[4].set_xlabel('Theoretical Quantiles', fontsize=10)
    axes[4].set_ylabel('Sample Quantiles', fontsize=8)
    axes[4].grid(True, alpha=0.3)
    
    stats.probplot(f_avg_vals, dist="norm", plot=axes[5])
    axes[5].set_title('$f_{avg}$ Q-Q Plot', fontsize=10)
    axes[5].set_xlabel('Theoretical Quantiles', fontsize=10)
    axes[5].set_ylabel('Sample Quantiles', fontsize=8)
    axes[5].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, format='pdf', bbox_inches='tight')
    plt.close()
    
    print(f"Full 6-subplot figure saved to: {output_file}")
    print(f"Statistics: deg_u={degree_u}, deg_w={degree_w}, Real CN={real_common_neighbors}")
    print(f"Jarque-Bera p-values: f_u={f_u_jb_p:.4f}, f_w={f_w_jb_p:.4f}, f_avg={f_avg_jb_p:.4f}")

def main():
    """Generate full 6-subplot figure for pair 6."""
    print("=== Generating Full 6-Subplot Figure for Pair 6 ===")
    
    # Parse the C++ output
    data_file = '/data/yizhangh/ldp-pq/gaussian-test/f_distribution_lrcwiki_high_cn_with_degrees_1000.txt'
    f_u_data, f_w_data, common_neighbors, degree_u, degree_w = parse_cpp_output(data_file)
    
    # Extract pair 6 data
    pair_idx = 6
    if pair_idx not in f_u_data or pair_idx not in f_w_data:
        print(f"ERROR: Pair {pair_idx} not found in data!")
        return
    
    f_u_vals = f_u_data[pair_idx]
    f_w_vals = f_w_data[pair_idx]
    f_avg_vals = [(f_u + f_w) / 2 for f_u, f_w in zip(f_u_vals, f_w_vals)]
    real_cn = common_neighbors[pair_idx]
    deg_u = degree_u[pair_idx]
    deg_w = degree_w[pair_idx]
    
    # Generate PDF figure
    output_file = '/data/yizhangh/ldp-pq/gaussian-test/pair_6_full_analysis.pdf'
    plot_pair6_full(f_u_vals, f_w_vals, f_avg_vals, real_cn, deg_u, deg_w, output_file)

if __name__ == "__main__":
    main()

