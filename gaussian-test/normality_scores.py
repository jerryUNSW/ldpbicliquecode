#!/usr/bin/env python3
"""
Comprehensive normality scoring for each vertex pair.
Uses multiple metrics to gauge how close each distribution is to normal.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import jarque_bera, normaltest, shapiro, kstest
import re
import os

def parse_cpp_output(filename):
    """Parse the C++ output file to extract f_u and f_w values."""
    f_u_data = {}
    f_w_data = {}
    
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

def compute_normality_scores(values, name):
    """Compute multiple normality scores for a distribution."""
    values = np.array(values)
    
    scores = {}
    
    # 1. Jarque-Bera test (p-value)
    jb_stat, jb_p = jarque_bera(values)
    scores['jarque_bera_p'] = jb_p
    scores['jarque_bera_stat'] = jb_stat
    
    # 2. D'Agostino-Pearson test (p-value)
    try:
        dp_stat, dp_p = normaltest(values)
        scores['dagostino_pearson_p'] = dp_p
        scores['dagostino_pearson_stat'] = dp_stat
    except:
        scores['dagostino_pearson_p'] = np.nan
        scores['dagostino_pearson_stat'] = np.nan
    
    # 3. Shapiro-Wilk test (p-value) - only for n <= 5000
    if len(values) <= 5000:
        try:
            sw_stat, sw_p = shapiro(values)
            scores['shapiro_wilk_p'] = sw_p
            scores['shapiro_wilk_stat'] = sw_stat
        except:
            scores['shapiro_wilk_p'] = np.nan
            scores['shapiro_wilk_stat'] = np.nan
    else:
        scores['shapiro_wilk_p'] = np.nan
        scores['shapiro_wilk_stat'] = np.nan
    
    # 4. Kolmogorov-Smirnov test vs normal
    ks_stat, ks_p = kstest(values, 'norm', args=(np.mean(values), np.std(values)))
    scores['kolmogorov_smirnov_p'] = ks_p
    scores['kolmogorov_smirnov_stat'] = ks_stat
    
    # 5. Skewness (closer to 0 = more normal)
    skewness = stats.skew(values)
    scores['skewness'] = skewness
    scores['skewness_abs'] = abs(skewness)
    
    # 6. Kurtosis (closer to 0 = more normal)
    kurtosis = stats.kurtosis(values)
    scores['kurtosis'] = kurtosis
    scores['kurtosis_abs'] = abs(kurtosis)
    
    # 7. Combined normality score (0-1, higher = more normal)
    # Based on p-values and shape measures
    p_score = np.mean([jb_p, dp_p if not np.isnan(dp_p) else 1.0, 
                       sw_p if not np.isnan(sw_p) else 1.0, ks_p])
    shape_score = 1.0 - min(1.0, (abs(skewness) + abs(kurtosis)) / 4.0)  # Normalize to 0-1
    scores['combined_score'] = (p_score + shape_score) / 2.0
    
    # 8. Effect size measures
    scores['skewness_effect'] = 'small' if abs(skewness) < 0.5 else 'moderate' if abs(skewness) < 1.0 else 'large'
    scores['kurtosis_effect'] = 'small' if abs(kurtosis) < 0.5 else 'moderate' if abs(kurtosis) < 1.0 else 'large'
    
    return scores

def print_normality_scores(pair_idx, f_u_scores, f_w_scores, f_avg_scores):
    """Print comprehensive normality scores for a pair."""
    print(f"\n=== Pair {pair_idx} Normality Scores ===")
    
    print(f"\nf_u scores:")
    print(f"  Jarque-Bera p-value: {f_u_scores['jarque_bera_p']:.6f}")
    print(f"  D'Agostino-Pearson p-value: {f_u_scores['dagostino_pearson_p']:.6f}")
    print(f"  Shapiro-Wilk p-value: {f_u_scores['shapiro_wilk_p']:.6f}")
    print(f"  Kolmogorov-Smirnov p-value: {f_u_scores['kolmogorov_smirnov_p']:.6f}")
    print(f"  Skewness: {f_u_scores['skewness']:.4f} (effect: {f_u_scores['skewness_effect']})")
    print(f"  Kurtosis: {f_u_scores['kurtosis']:.4f} (effect: {f_u_scores['kurtosis_effect']})")
    print(f"  Combined normality score: {f_u_scores['combined_score']:.4f}")
    
    print(f"\nf_w scores:")
    print(f"  Jarque-Bera p-value: {f_w_scores['jarque_bera_p']:.6f}")
    print(f"  D'Agostino-Pearson p-value: {f_w_scores['dagostino_pearson_p']:.6f}")
    print(f"  Shapiro-Wilk p-value: {f_w_scores['shapiro_wilk_p']:.6f}")
    print(f"  Kolmogorov-Smirnov p-value: {f_w_scores['kolmogorov_smirnov_p']:.6f}")
    print(f"  Skewness: {f_w_scores['skewness']:.4f} (effect: {f_w_scores['skewness_effect']})")
    print(f"  Kurtosis: {f_w_scores['kurtosis']:.4f} (effect: {f_w_scores['kurtosis_effect']})")
    print(f"  Combined normality score: {f_w_scores['combined_score']:.4f}")
    
    print(f"\nf_avg scores (Gaussian assumption):")
    print(f"  Jarque-Bera p-value: {f_avg_scores['jarque_bera_p']:.6f}")
    print(f"  D'Agostino-Pearson p-value: {f_avg_scores['dagostino_pearson_p']:.6f}")
    print(f"  Shapiro-Wilk p-value: {f_avg_scores['shapiro_wilk_p']:.6f}")
    print(f"  Kolmogorov-Smirnov p-value: {f_avg_scores['kolmogorov_smirnov_p']:.6f}")
    print(f"  Skewness: {f_avg_scores['skewness']:.4f} (effect: {f_avg_scores['skewness_effect']})")
    print(f"  Kurtosis: {f_avg_scores['kurtosis']:.4f} (effect: {f_avg_scores['kurtosis_effect']})")
    print(f"  Combined normality score: {f_avg_scores['combined_score']:.4f}")
    
    # Interpretation
    if f_avg_scores['combined_score'] > 0.7:
        interpretation = "STRONGLY NORMAL"
    elif f_avg_scores['combined_score'] > 0.5:
        interpretation = "MODERATELY NORMAL"
    elif f_avg_scores['combined_score'] > 0.3:
        interpretation = "WEAKLY NORMAL"
    else:
        interpretation = "NOT NORMAL"
    
    print(f"\n  → {interpretation} (score: {f_avg_scores['combined_score']:.3f})")

def create_normality_heatmap(all_scores, output_dir):
    """Create a heatmap of normality scores across all pairs."""
    n_pairs = len(all_scores)
    
    # Prepare data for heatmap
    estimators = ['f_u', 'f_w', 'f_avg']
    metrics = ['jarque_bera_p', 'dagostino_pearson_p', 'shapiro_wilk_p', 
               'kolmogorov_smirnov_p', 'skewness_abs', 'kurtosis_abs', 'combined_score']
    
    # Create matrix
    data_matrix = np.zeros((n_pairs * len(estimators), len(metrics)))
    row_labels = []
    
    for i, (pair_idx, scores) in enumerate(all_scores):
        for j, estimator in enumerate(estimators):
            row_idx = i * len(estimators) + j
            row_labels.append(f"Pair{pair_idx}_{estimator}")
            
            for k, metric in enumerate(metrics):
                if metric in ['skewness_abs', 'kurtosis_abs']:
                    # For these, lower is better (closer to 0)
                    data_matrix[row_idx, k] = 1.0 - min(1.0, scores[estimator][metric])
                else:
                    # For p-values and combined_score, higher is better
                    data_matrix[row_idx, k] = scores[estimator][metric]
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(12, 8))
    im = ax.imshow(data_matrix, cmap='RdYlGn', aspect='auto')
    
    # Set labels
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(metrics, rotation=45, ha='right')
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Normality Score (0=not normal, 1=perfectly normal)')
    
    plt.title('Normality Scores Heatmap Across All Pairs', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'normality_scores_heatmap.png'), dpi=300, bbox_inches='tight')
    plt.close()

def main():
    """Main analysis function."""
    print("=== Comprehensive Normality Scoring Analysis ===")
    print("Multiple metrics to gauge how close each distribution is to normal")
    
    # Parse the C++ output
    data_file = '/data/yizhangh/ldp-pq/gaussian-test/f_distribution_unicode_high_cn_500.txt'
    f_u_data, f_w_data = parse_cpp_output(data_file)
    
    print(f"Loaded data for {len(f_u_data)} pairs")
    
    # Create output directory
    output_dir = '/data/yizhangh/ldp-pq/gaussian-test'
    os.makedirs(output_dir, exist_ok=True)
    
    # Analyze each pair
    all_scores = []
    
    for pair_idx in sorted(f_u_data.keys()):
        if pair_idx in f_w_data:
            f_u_vals = f_u_data[pair_idx]
            f_w_vals = f_w_data[pair_idx]
            f_avg_vals = [(f_u + f_w) / 2 for f_u, f_w in zip(f_u_vals, f_w_vals)]
            
            # Compute scores for each estimator
            f_u_scores = compute_normality_scores(f_u_vals, 'f_u')
            f_w_scores = compute_normality_scores(f_w_vals, 'f_w')
            f_avg_scores = compute_normality_scores(f_avg_vals, 'f_avg')
            
            # Print detailed scores
            print_normality_scores(pair_idx, f_u_scores, f_w_scores, f_avg_scores)
            
            # Store for summary
            all_scores.append((pair_idx, {
                'f_u': f_u_scores,
                'f_w': f_w_scores,
                'f_avg': f_avg_scores
            }))
    
    # Create heatmap
    create_normality_heatmap(all_scores, output_dir)
    
    # Summary statistics
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    
    f_avg_scores = [scores['f_avg']['combined_score'] for _, scores in all_scores]
    
    print(f"f_avg combined normality scores:")
    print(f"  Mean: {np.mean(f_avg_scores):.4f}")
    print(f"  Std: {np.std(f_avg_scores):.4f}")
    print(f"  Min: {np.min(f_avg_scores):.4f}")
    print(f"  Max: {np.max(f_avg_scores):.4f}")
    
    # Categorize pairs
    strongly_normal = sum(1 for score in f_avg_scores if score > 0.7)
    moderately_normal = sum(1 for score in f_avg_scores if 0.5 < score <= 0.7)
    weakly_normal = sum(1 for score in f_avg_scores if 0.3 < score <= 0.5)
    not_normal = sum(1 for score in f_avg_scores if score <= 0.3)
    
    print(f"\nCategorization of f_avg distributions:")
    print(f"  Strongly normal (>0.7): {strongly_normal}/{len(f_avg_scores)} pairs")
    print(f"  Moderately normal (0.5-0.7): {moderately_normal}/{len(f_avg_scores)} pairs")
    print(f"  Weakly normal (0.3-0.5): {weakly_normal}/{len(f_avg_scores)} pairs")
    print(f"  Not normal (≤0.3): {not_normal}/{len(f_avg_scores)} pairs")
    
    print(f"\nHeatmap saved to: {output_dir}/normality_scores_heatmap.png")

if __name__ == "__main__":
    main()
