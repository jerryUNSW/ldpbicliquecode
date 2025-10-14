#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

# Read the data
df = pd.read_csv('/data/yizhangh/ldp-pq/f_distribution_data.csv')

print("=== Multi-Source f Distribution Analysis ===")
print(f"Total data points: {len(df)}")
print(f"Number of vertex pairs: {df['pair_idx'].nunique()}")
print(f"Runs per pair: {len(df[df['pair_idx'] == 0])}")

# Analysis for each vertex pair
for pair_idx in sorted(df['pair_idx'].unique()):
    pair_data = df[df['pair_idx'] == pair_idx]
    u, w = pair_data['u'].iloc[0], pair_data['w'].iloc[0]
    true_common = pair_data['true_common'].iloc[0]
    
    print(f"\n--- Pair {pair_idx} (u={u}, w={w}, true={true_common}) ---")
    
    # Individual estimates
    f_u_w = pair_data['f_u_w'].values
    f_w_u = pair_data['f_w_u'].values
    f_avg = pair_data['f_avg'].values
    
    # Statistics for individual estimates
    print(f"f_u_w: mean={f_u_w.mean():.3f}, std={f_u_w.std():.3f}, skew={stats.skew(f_u_w):.3f}, kurt={stats.kurtosis(f_u_w):.3f}")
    print(f"f_w_u: mean={f_w_u.mean():.3f}, std={f_w_u.std():.3f}, skew={stats.skew(f_w_u):.3f}, kurt={stats.kurtosis(f_w_u):.3f}")
    print(f"f_avg: mean={f_avg.mean():.3f}, std={f_avg.std():.3f}, skew={stats.skew(f_avg):.3f}, kurt={stats.kurtosis(f_avg):.3f}")
    
    # Comprehensive normality assessment
    _, p_u_w = stats.normaltest(f_u_w)
    _, p_w_u = stats.normaltest(f_w_u)
    _, p_avg = stats.normaltest(f_avg)
    
    # Additional normality tests
    shapiro_u_w, shapiro_p_u_w = stats.shapiro(f_u_w) if len(f_u_w) <= 5000 else (np.nan, np.nan)
    shapiro_w_u, shapiro_p_w_u = stats.shapiro(f_w_u) if len(f_w_u) <= 5000 else (np.nan, np.nan)
    shapiro_avg, shapiro_p_avg = stats.shapiro(f_avg) if len(f_avg) <= 5000 else (np.nan, np.nan)
    
    # Kolmogorov-Smirnov test
    ks_stat_u_w, ks_p_u_w = stats.kstest(f_u_w, 'norm', args=(f_u_w.mean(), f_u_w.std()))
    ks_stat_w_u, ks_p_w_u = stats.kstest(f_w_u, 'norm', args=(f_w_u.mean(), f_w_u.std()))
    ks_stat_avg, ks_p_avg = stats.kstest(f_avg, 'norm', args=(f_avg.mean(), f_avg.std()))
    
    # Calculate normality scores (0-100, higher is more normal)
    def normality_score(p_values):
        """Convert p-values to normality score (0-100)"""
        valid_p = [p for p in p_values if not np.isnan(p)]
        if not valid_p:
            return np.nan
        # Use the minimum p-value among tests, convert to score
        min_p = min(valid_p)
        # p > 0.2 = excellent (100), p > 0.1 = very good (90-100), p > 0.05 = good (70-90)
        if min_p > 0.2:
            return 100
        elif min_p > 0.1:
            return 90 + (min_p - 0.1) * 100  # 90-100
        elif min_p > 0.05:
            return 70 + (min_p - 0.05) * 400  # 70-90
        elif min_p > 0.01:
            return 40 + (min_p - 0.01) * 750  # 40-70
        else:
            return min_p * 4000  # 0-40
    
    score_u_w = normality_score([p_u_w, shapiro_p_u_w, ks_p_u_w])
    score_w_u = normality_score([p_w_u, shapiro_p_w_u, ks_p_w_u])
    score_avg = normality_score([p_avg, shapiro_p_avg, ks_p_avg])
    
    print(f"Normality tests (p-values): f_u_w={p_u_w:.3f}, f_w_u={p_w_u:.3f}, f_avg={p_avg:.3f}")
    print(f"Shapiro-Wilk (p-values): f_u_w={shapiro_p_u_w:.3f}, f_w_u={shapiro_p_w_u:.3f}, f_avg={shapiro_p_avg:.3f}")
    print(f"KS test (p-values): f_u_w={ks_p_u_w:.3f}, f_w_u={ks_p_w_u:.3f}, f_avg={ks_p_avg:.3f}")
    print(f"Normality Scores (0-100): f_u_w={score_u_w:.1f}, f_w_u={score_w_u:.1f}, f_avg={score_avg:.1f}")

# Create comprehensive distribution visualizations for all pairs
num_pairs = df['pair_idx'].nunique()

# Create a large figure with subplots for each vertex pair
fig, axes = plt.subplots(num_pairs, 3, figsize=(18, 4*num_pairs))
if num_pairs == 1:
    axes = axes.reshape(1, -1)

for pair_idx in sorted(df['pair_idx'].unique()):
    pair_data = df[df['pair_idx'] == pair_idx]
    u, w = pair_data['u'].iloc[0], pair_data['w'].iloc[0]
    true_common = pair_data['true_common'].iloc[0]
    
    f_u_w = pair_data['f_u_w'].values
    f_w_u = pair_data['f_w_u'].values
    f_avg = pair_data['f_avg'].values
    
    # Histograms with normal overlay
    ax1 = axes[pair_idx, 0]
    ax1.hist(f_u_w, bins=20, alpha=0.7, density=True, color='blue', 
             label=f'f_u_w (mean={f_u_w.mean():.2f}, std={f_u_w.std():.2f})')
    
    # Overlay normal distribution
    x_range = np.linspace(f_u_w.min(), f_u_w.max(), 100)
    normal_pdf = stats.norm.pdf(x_range, f_u_w.mean(), f_u_w.std())
    ax1.plot(x_range, normal_pdf, 'r-', linewidth=2, label='Normal fit')
    ax1.set_title(f'f_u_w Distribution (Pair {pair_idx}: u={u}, w={w}, true={true_common})')
    ax1.set_xlabel('f_u_w value')
    ax1.set_ylabel('Density')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    ax2 = axes[pair_idx, 1]
    ax2.hist(f_w_u, bins=20, alpha=0.7, density=True, color='red',
             label=f'f_w_u (mean={f_w_u.mean():.2f}, std={f_w_u.std():.2f})')
    
    # Overlay normal distribution
    x_range = np.linspace(f_w_u.min(), f_w_u.max(), 100)
    normal_pdf = stats.norm.pdf(x_range, f_w_u.mean(), f_w_u.std())
    ax2.plot(x_range, normal_pdf, 'r-', linewidth=2, label='Normal fit')
    ax2.set_title(f'f_w_u Distribution (Pair {pair_idx}: u={u}, w={w}, true={true_common})')
    ax2.set_xlabel('f_w_u value')
    ax2.set_ylabel('Density')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    ax3 = axes[pair_idx, 2]
    ax3.hist(f_avg, bins=20, alpha=0.7, density=True, color='green',
             label=f'f_avg (mean={f_avg.mean():.2f}, std={f_avg.std():.2f})')
    
    # Overlay normal distribution
    x_range = np.linspace(f_avg.min(), f_avg.max(), 100)
    normal_pdf = stats.norm.pdf(x_range, f_avg.mean(), f_avg.std())
    ax3.plot(x_range, normal_pdf, 'r-', linewidth=2, label='Normal fit')
    ax3.set_title(f'f_avg Distribution (Pair {pair_idx}: u={u}, w={w}, true={true_common})')
    ax3.set_xlabel('f_avg value')
    ax3.set_ylabel('Density')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/data/yizhangh/ldp-pq/f_distributions_all_pairs.png', dpi=300, bbox_inches='tight')
print(f"\nDistribution plots saved to f_distributions_all_pairs.png")

# Create Q-Q plots for normality assessment
fig, axes = plt.subplots(num_pairs, 3, figsize=(18, 4*num_pairs))
if num_pairs == 1:
    axes = axes.reshape(1, -1)

for pair_idx in sorted(df['pair_idx'].unique()):
    pair_data = df[df['pair_idx'] == pair_idx]
    u, w = pair_data['u'].iloc[0], pair_data['w'].iloc[0]
    true_common = pair_data['true_common'].iloc[0]
    
    f_u_w = pair_data['f_u_w'].values
    f_w_u = pair_data['f_w_u'].values
    f_avg = pair_data['f_avg'].values
    
    # Q-Q plots for normality testing
    stats.probplot(f_u_w, dist="norm", plot=axes[pair_idx, 0])
    axes[pair_idx, 0].set_title(f'Q-Q Plot: f_u_w (Pair {pair_idx}, true={true_common})')
    axes[pair_idx, 0].grid(True, alpha=0.3)
    
    stats.probplot(f_w_u, dist="norm", plot=axes[pair_idx, 1])
    axes[pair_idx, 1].set_title(f'Q-Q Plot: f_w_u (Pair {pair_idx}, true={true_common})')
    axes[pair_idx, 1].grid(True, alpha=0.3)
    
    stats.probplot(f_avg, dist="norm", plot=axes[pair_idx, 2])
    axes[pair_idx, 2].set_title(f'Q-Q Plot: f_avg (Pair {pair_idx}, true={true_common})')
    axes[pair_idx, 2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/data/yizhangh/ldp-pq/qq_plots_all_pairs.png', dpi=300, bbox_inches='tight')
print(f"Q-Q plots saved to qq_plots_all_pairs.png")

# Correlation analysis
print(f"\n=== Correlation Analysis ===")
print(f"Correlation between f_u_w and f_w_u: {np.corrcoef(f_u_w, f_w_u)[0,1]:.3f}")

# Check if averaging improves normality
print(f"\n=== Normality Improvement Analysis ===")
print("For each vertex pair, compare individual vs averaged estimates:")

for pair_idx in sorted(df['pair_idx'].unique()):
    pair_data = df[df['pair_idx'] == pair_idx]
    f_u_w = pair_data['f_u_w'].values
    f_w_u = pair_data['f_w_u'].values
    f_avg = pair_data['f_avg'].values
    
    # Normality tests
    _, p_u_w = stats.normaltest(f_u_w)
    _, p_w_u = stats.normaltest(f_w_u)
    _, p_avg = stats.normaltest(f_avg)
    
    # Skewness and kurtosis
    skew_u_w = abs(stats.skew(f_u_w))
    skew_w_u = abs(stats.skew(f_w_u))
    skew_avg = abs(stats.skew(f_avg))
    
    kurt_u_w = abs(stats.kurtosis(f_u_w))
    kurt_w_u = abs(stats.kurtosis(f_w_u))
    kurt_avg = abs(stats.kurtosis(f_avg))
    
    print(f"Pair {pair_idx}:")
    print(f"  Individual: f_u_w p={p_u_w:.3f}, skew={skew_u_w:.3f}, kurt={kurt_u_w:.3f}")
    print(f"  Individual: f_w_u p={p_w_u:.3f}, skew={skew_w_u:.3f}, kurt={kurt_w_u:.3f}")
    print(f"  Averaged:   f_avg p={p_avg:.3f}, skew={skew_avg:.3f}, kurt={kurt_avg:.3f}")
    
    # Check if averaging improves normality
    avg_improves = (skew_avg < min(skew_u_w, skew_w_u)) and (kurt_avg < min(kurt_u_w, kurt_w_u))
    print(f"  Averaging improves normality: {'✓' if avg_improves else '✗'}")
    print()

print("\n=== Summary ===")
print("The multi-source approach computes f estimates from both perspectives (u→w and w→u)")
print("and averages them. This should theoretically improve the distribution properties")
print("by reducing bias and variance through averaging independent estimates.")

print("\n=== Normality Scoring Guide ===")
print("Normality Score (0-100 scale):")
print("- 90-100: Excellent normality (p > 0.1)")
print("- 70-90:  Good normality (p > 0.05)")  
print("- 40-70:  Moderate normality (p > 0.01)")
print("- 0-40:   Poor normality (p < 0.01)")
print()
print("Higher scores indicate better adherence to normal distribution assumptions")
print("required by the moment correction technique.")
