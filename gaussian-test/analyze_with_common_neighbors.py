#!/usr/bin/env python3
"""
Analysis script for f_u, f_w, f_avg distributions with degree annotations.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import jarque_bera
import re
import os

def parse_cpp_output(filename):
    """Parse the C++ output file to extract f_u, f_w values, common neighbor counts, and degrees."""
    f_u_data = {}  # pair_idx -> list of f_u values
    f_w_data = {}  # pair_idx -> list of f_w values
    common_neighbors = {}  # pair_idx -> real common neighbor count
    degree_u = {}  # pair_idx -> degree of vertex u
    degree_w = {}  # pair_idx -> degree of vertex w
    
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

def analyze_single_pair(pair_idx, f_u_vals, f_w_vals, real_common_neighbors):
    """Analyze a single vertex pair (u,w) with 500 samples each."""
    f_u_vals = np.array(f_u_vals)
    f_w_vals = np.array(f_w_vals)
    f_avg_vals = (f_u_vals + f_w_vals) / 2
    
    # Basic statistics
    f_u_mean, f_u_std = np.mean(f_u_vals), np.std(f_u_vals)
    f_w_mean, f_w_std = np.mean(f_w_vals), np.std(f_w_vals)
    f_avg_mean, f_avg_std = np.mean(f_avg_vals), np.std(f_avg_vals)
    
    # Normality tests
    f_u_jb_p = jarque_bera(f_u_vals)[1]
    f_w_jb_p = jarque_bera(f_w_vals)[1]
    f_avg_jb_p = jarque_bera(f_avg_vals)[1]
    
    # Skewness and kurtosis
    f_u_skew, f_u_kurt = stats.skew(f_u_vals), stats.kurtosis(f_u_vals)
    f_w_skew, f_w_kurt = stats.skew(f_w_vals), stats.kurtosis(f_w_vals)
    f_avg_skew, f_avg_kurt = stats.skew(f_avg_vals), stats.kurtosis(f_avg_vals)
    
    return {
        'f_u': {'mean': f_u_mean, 'std': f_u_std, 'skew': f_u_skew, 'kurt': f_u_kurt, 'jb_p': f_u_jb_p},
        'f_w': {'mean': f_w_mean, 'std': f_w_std, 'skew': f_w_skew, 'kurt': f_w_kurt, 'jb_p': f_w_jb_p},
        'f_avg': {'mean': f_avg_mean, 'std': f_avg_std, 'skew': f_avg_skew, 'kurt': f_avg_kurt, 'jb_p': f_avg_jb_p}
    }

def plot_single_pair(pair_idx, f_u_vals, f_w_vals, f_avg_vals, real_common_neighbors, degree_u, degree_w, output_dir):
    """Create plots for a single vertex pair."""
    # Convert to numpy arrays
    f_u_vals = np.array(f_u_vals)
    f_w_vals = np.array(f_w_vals)
    f_avg_vals = np.array(f_avg_vals)
    
    # Compute normality scores
    f_u_jb_p = jarque_bera(f_u_vals)[1]
    f_w_jb_p = jarque_bera(f_w_vals)[1]
    f_avg_jb_p = jarque_bera(f_avg_vals)[1]
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f'Vertex Pair {pair_idx} Analysis (Real CN={real_common_neighbors}, 1000 samples each)', fontsize=16, fontweight='bold')
    
    # Top row: Histograms with normal overlays
    # f_u histogram
    axes[0,0].hist(f_u_vals, bins=30, alpha=0.7, color='blue', density=True, label='fu data')
    axes[0,0].axvline(x=real_common_neighbors, color='red', linestyle='--', linewidth=2, label=f'Real CN={real_common_neighbors}')
    
    # Overlay normal distribution
    f_u_mean = np.mean(f_u_vals)
    f_u_std = np.std(f_u_vals)
    x_range = np.linspace(f_u_vals.min(), f_u_vals.max(), 100)
    normal_curve = stats.norm.pdf(x_range, f_u_mean, f_u_std)
    axes[0,0].plot(x_range, normal_curve, 'k-', linewidth=2, label=f'Normal(μ={f_u_mean:.1f}, σ={f_u_std:.1f})')
    
    axes[0,0].set_title(f'fu Distribution (n={len(f_u_vals)}, deg_u={degree_u})\nJarque-Bera p-value={f_u_jb_p:.4f}', fontsize=10)
    axes[0,0].set_xlabel('fu values', fontsize=8)
    axes[0,0].set_ylabel('Density', fontsize=8)
    axes[0,0].legend(loc='upper left', frameon=False, fontsize=8)
    axes[0,0].grid(True, alpha=0.3)
    
    # f_w histogram
    axes[0,1].hist(f_w_vals, bins=30, alpha=0.7, color='red', density=True, label='fw data')
    axes[0,1].axvline(x=real_common_neighbors, color='blue', linestyle='--', linewidth=2, label=f'Real CN={real_common_neighbors}')
    
    # Overlay normal distribution
    f_w_mean = np.mean(f_w_vals)
    f_w_std = np.std(f_w_vals)
    x_range = np.linspace(f_w_vals.min(), f_w_vals.max(), 100)
    normal_curve = stats.norm.pdf(x_range, f_w_mean, f_w_std)
    axes[0,1].plot(x_range, normal_curve, 'k-', linewidth=2, label=f'Normal(μ={f_w_mean:.1f}, σ={f_w_std:.1f})')
    
    axes[0,1].set_title(f'fw Distribution (n={len(f_w_vals)}, deg_w={degree_w})\nJarque-Bera p-value={f_w_jb_p:.4f}', fontsize=10)
    axes[0,1].set_xlabel('fw values', fontsize=8)
    axes[0,1].set_ylabel('Density', fontsize=8)
    axes[0,1].legend(loc='upper left', frameon=False, fontsize=8)
    axes[0,1].grid(True, alpha=0.3)
    
    # f_avg histogram
    axes[0,2].hist(f_avg_vals, bins=30, alpha=0.7, color='green', density=True, label='favg data')
    axes[0,2].axvline(x=real_common_neighbors, color='red', linestyle='--', linewidth=2, label=f'Real CN={real_common_neighbors}')
    
    # Overlay normal distribution
    f_avg_mean = np.mean(f_avg_vals)
    f_avg_std = np.std(f_avg_vals)
    x_range = np.linspace(f_avg_vals.min(), f_avg_vals.max(), 100)
    normal_curve = stats.norm.pdf(x_range, f_avg_mean, f_avg_std)
    axes[0,2].plot(x_range, normal_curve, 'k-', linewidth=2, label=f'Normal(μ={f_avg_mean:.1f}, σ={f_avg_std:.1f})')
    
    axes[0,2].set_title(f'favg Distribution (n={len(f_avg_vals)})\nJarque-Bera p-value={f_avg_jb_p:.4f}', fontsize=10)
    axes[0,2].set_xlabel('favg values', fontsize=8)
    axes[0,2].set_ylabel('Density', fontsize=8)
    axes[0,2].legend(loc='upper left', frameon=False, fontsize=8)
    axes[0,2].grid(True, alpha=0.3)
    
    # Bottom row: Q-Q plots
    stats.probplot(f_u_vals, dist="norm", plot=axes[1,0])
    axes[1,0].set_title('fu Q-Q Plot vs Normal', fontsize=10)
    axes[1,0].grid(True, alpha=0.3)
    
    stats.probplot(f_w_vals, dist="norm", plot=axes[1,1])
    axes[1,1].set_title('fw Q-Q Plot vs Normal', fontsize=10)
    axes[1,1].grid(True, alpha=0.3)
    
    stats.probplot(f_avg_vals, dist="norm", plot=axes[1,2])
    axes[1,2].set_title('favg Q-Q Plot vs Normal', fontsize=10)
    axes[1,2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'pair_{pair_idx}_analysis_with_cn.png'), dpi=300, bbox_inches='tight')
    plt.close()

def main():
    """Main analysis function."""
    print("=== Per-Pair Analysis with Real Common Neighbor Counts and Degrees ===")
    print("For each vertex pair (u,w): 1000 f_u, 1000 f_w, 1000 f_avg samples + real CN count + degrees")
    
    # Parse the C++ output
    data_file = '/data/yizhangh/ldp-pq/gaussian-test/f_distribution_lrcwiki_high_cn_with_degrees_1000.txt'
    f_u_data, f_w_data, common_neighbors, degree_u, degree_w = parse_cpp_output(data_file)
    
    print(f"Loaded data for {len(f_u_data)} pairs")
    print(f"Common neighbor counts: {dict(common_neighbors)}")
    print(f"Degrees u: {dict(degree_u)}")
    print(f"Degrees w: {dict(degree_w)}")
    
    # Create output directory
    output_dir = '/data/yizhangh/ldp-pq/gaussian-test'
    os.makedirs(output_dir, exist_ok=True)
    
    # Analyze each pair individually
    all_results = []
    for pair_idx in sorted(f_u_data.keys()):
        if pair_idx in f_w_data and pair_idx in common_neighbors:
            f_u_vals = f_u_data[pair_idx]
            f_w_vals = f_w_data[pair_idx]
            real_cn = common_neighbors[pair_idx]
            
            print(f"\n--- Pair {pair_idx} ---")
            print(f"Real common neighbors: {real_cn}")
            print(f"Degree u: {degree_u.get(pair_idx, 'N/A')}")
            print(f"Degree w: {degree_w.get(pair_idx, 'N/A')}")
            print(f"f_u samples: {len(f_u_vals)}")
            print(f"f_w samples: {len(f_w_vals)}")
            
            # Analyze this pair
            results = analyze_single_pair(pair_idx, f_u_vals, f_w_vals, real_cn)
            all_results.append((pair_idx, results))
            
            # Create individual plot
            f_avg_vals = [(f_u + f_w) / 2 for f_u, f_w in zip(f_u_vals, f_w_vals)]
            deg_u = degree_u.get(pair_idx, 0)
            deg_w = degree_w.get(pair_idx, 0)
            plot_single_pair(pair_idx, f_u_vals, f_w_vals, f_avg_vals, real_cn, deg_u, deg_w, output_dir)
    
    # Summary across all pairs
    print(f"\n=== Summary across {len(all_results)} pairs ===")
    print("Pair\tReal CN\tDeg_u\tDeg_w\tfu_JB_p\tfw_JB_p\tfavg_JB_p")
    print("-" * 60)
    
    for pair_idx, results in all_results:
        real_cn = common_neighbors.get(pair_idx, 0)
        deg_u = degree_u.get(pair_idx, 0)
        deg_w = degree_w.get(pair_idx, 0)
        fu_jb = results['f_u']['jb_p']
        fw_jb = results['f_w']['jb_p']
        favg_jb = results['f_avg']['jb_p']
        
        print(f"{pair_idx}\t{real_cn}\t{deg_u}\t{deg_w}\t{fu_jb:.4f}\t{fw_jb:.4f}\t{favg_jb:.4f}")
    
    print(f"\nPlots saved to: {output_dir}")
    print("Each plot shows: fu, fw, favg distributions + Q-Q plots + real CN + degrees")

if __name__ == "__main__":
    main()
