#!/usr/bin/env python3
"""
Generate 2-subplot figure for pair 6 (f_avg only) for the revised paper.
Based on analyze_with_common_neighbors.py but outputs PDF format with only f_avg plots.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import jarque_bera
import re
import os

def plt_settings():
    """Apply styling from biclique plots."""
    plt.style.use('default')
    plt.rcParams['savefig.dpi'] = 300 
    plt.rcParams['figure.dpi'] = 300 
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams["legend.framealpha"] = 0
    plt.rcParams["legend.handletextpad"] = 0.1
    plt.rcParams["legend.columnspacing"] = 0.2
    plt.rcParams['pdf.fonttype'] = 42

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

def plot_pair6_favg_only(f_avg_vals, real_common_neighbors, degree_u, degree_w, output_file):
    """Create 2-subplot figure for pair 6 (f_avg only)."""
    plt_settings()
    
    # Convert to numpy array
    f_avg_vals = np.array(f_avg_vals)
    
    # Compute normality score
    f_avg_jb_p = jarque_bera(f_avg_vals)[1]
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left plot: f_avg histogram with normal overlay
    axes[0].hist(f_avg_vals, bins=25, alpha=0.7, color='green', density=True, label='$f_{avg}$ data', edgecolor='black', linewidth=0.5)
    axes[0].axvline(x=real_common_neighbors, color='red', linestyle='--', linewidth=2.5, label=f'Real={real_common_neighbors}')
    
    # Overlay normal distribution
    f_avg_mean = np.mean(f_avg_vals)
    f_avg_std = np.std(f_avg_vals)
    x_range = np.linspace(f_avg_vals.min(), f_avg_vals.max(), 100)
    normal_curve = stats.norm.pdf(x_range, f_avg_mean, f_avg_std)
    axes[0].plot(x_range, normal_curve, 'k-', linewidth=2.5, label=f'N({f_avg_mean:.1f}, {f_avg_std:.1f})')
    
    axes[0].set_title(f'(deg$_u$={degree_u}, deg$_w$={degree_w}, Jarque-Bera score={f_avg_jb_p:.3f})', fontsize=16)
    axes[0].set_xlabel('$f_{avg}$ values', fontsize=16)
    axes[0].set_ylabel('Density', fontsize=16)
    axes[0].legend(loc='upper left', frameon=False, fontsize=14)
    axes[0].grid(True, alpha=0.3)
    axes[0].set_ylim(0, max(normal_curve) * 1.3)
    axes[0].tick_params(axis='both', which='major', labelsize=14)
    
    # Right plot: f_avg Q-Q plot
    stats.probplot(f_avg_vals, dist="norm", plot=axes[1])
    axes[1].set_title(f'$f_{{avg}}$ Q-Q Plot vs Normal', fontsize=16)
    axes[1].set_xlabel('Theoretical Quantiles', fontsize=16)
    axes[1].set_ylabel('Sample Quantiles', fontsize=16)
    axes[1].grid(True, alpha=0.3)
    axes[1].tick_params(axis='both', which='major', labelsize=14)
    
    plt.tight_layout()
    plt.savefig(output_file, format='pdf', bbox_inches='tight')
    plt.close()
    
    print(f"2-subplot figure saved to: {output_file}")
    print(f"Statistics: deg_u={degree_u}, deg_w={degree_w}, Real CN={real_common_neighbors}")
    print(f"Jarque-Bera p-value for f_avg: {f_avg_jb_p:.4f}")

def main():
    """Generate 2-subplot figure for pair 6 (f_avg only)."""
    print("=== Generating 2-Subplot Figure for Pair 6 (f_avg only) ===")
    
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
    output_file = '/data/yizhangh/ldp-pq/gaussian-test/pair_6_favg_only.pdf'
    plot_pair6_favg_only(f_avg_vals, real_cn, deg_u, deg_w, output_file)

if __name__ == "__main__":
    main()
