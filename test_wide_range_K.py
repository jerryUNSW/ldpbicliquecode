#!/usr/bin/env python3

import numpy as np
import math
from scipy import stats
import matplotlib.pyplot as plt

def factorial(n):
    """Compute factorial"""
    if n < 0:
        return 0
    return math.factorial(n)

def binomial(n, k):
    """Compute binomial coefficient C(n, k)"""
    if k > n or k < 0:
        return 0
    return factorial(n) // (factorial(k) * factorial(n - k))

def stirling_first(n, k):
    """Compute Stirling numbers of the first kind"""
    if n == 0 and k == 0:
        return 1
    if n == 0 or k == 0:
        return 0
    if k > n:
        return 0
    if n == k:
        return 1
    
    # Use recurrence relation: s(n,k) = (n-1)*s(n-1,k) + s(n-1,k-1)
    dp = np.zeros((n+1, k+1))
    dp[0][0] = 1
    
    for i in range(1, n+1):
        for j in range(1, min(i+1, k+1)):
            dp[i][j] = (i-1) * dp[i-1][j] + dp[i-1][j-1]
    
    return dp[n][k]

def estimate_choose_K_moment_correction(K, f_prime, variance):
    """
    Estimate C(f, K) using moment correction assuming f' ~ N(f, variance)
    """
    if K == 0:
        return 1.0
    if K == 1:
        return f_prime
    
    # Construct corrected moments up to K
    moments = np.zeros(K + 1)
    moments[1] = f_prime
    
    # Compute higher moments using normal distribution properties
    for k in range(2, K + 1):
        moments[k] = f_prime ** k
        
        # Subtract bias terms due to variance
        for i in range(1, k//2 + 1):
            if k - 2*i >= 0:
                # Double factorial: (2i-1)!! = (2i-1)(2i-3)...1
                double_fact = 1
                for j in range(2*i-1, 0, -2):
                    double_fact *= j
                
                coeff = binomial(k, 2*i) * double_fact * (variance ** i)
                moments[k] -= coeff * moments[k - 2*i]
    
    # Use inclusion-exclusion to get C(f, K)
    result = 0
    for i in range(K + 1):
        sign = 1 if i % 2 == 0 else -1
        coeff = stirling_first(K, K - i)
        result += sign * coeff * moments[K - i]
    
    return result / factorial(K)

def estimate_choose_K_naive(K, f_prime):
    """
    Naive approach: directly compute C(f', K) without correction
    """
    if f_prime < K:
        return 0
    
    result = 1.0
    for i in range(K):
        result *= (f_prime - i)
    return result / factorial(K)

def test_wide_range_K():
    """
    Test moment correction for a wide range of K values
    """
    print("=== Testing Moment Correction for Wide Range of K Values ===\n")
    
    # Test parameters - 10 diverse f values, K from 2 to f/2
    true_f_values = [4, 6, 8, 10, 15, 20, 25, 30, 40, 50]  # Small to large f values
    variance_values = [1.0, 2.0, 4.0]  # Different noise levels
    
    results = []
    n_samples = 5000  # Reduced for faster computation
    
    print(f"Testing {len(true_f_values)} f values with K from 2 to f/2 × {len(variance_values)} variances")
    print(f"Each case uses {n_samples} Monte Carlo samples")
    print("=" * 100)
    
    for true_f in true_f_values:
        print(f"\nTesting with true f = {true_f}")
        print("-" * 60)
        
        # Generate K values from 2 to f/2
        max_K = true_f // 2
        K_values = list(range(2, max_K + 1))
        
        for variance in variance_values:
            print(f"\n  Variance σ² = {variance}")
            
            # Generate samples once for this f and variance
            np.random.seed(42)
            f_prime_samples = np.random.normal(true_f, np.sqrt(variance), n_samples)
            
            # Test all K values for this f and variance
            for K in K_values:
                    
                true_choose_k = binomial(true_f, K)
                
                # Estimate C(f, K) using both methods
                moment_correction_estimates = []
                naive_estimates = []
                
                for f_prime in f_prime_samples:
                    # Moment correction method
                    mc_est = estimate_choose_K_moment_correction(K, f_prime, variance)
                    moment_correction_estimates.append(mc_est)
                    
                    # Naive method
                    naive_est = estimate_choose_K_naive(K, f_prime)
                    naive_estimates.append(naive_est)
                
                # Calculate statistics
                mc_mean = np.mean(moment_correction_estimates)
                mc_bias = mc_mean - true_choose_k
                mc_rmse = np.sqrt(np.mean((np.array(moment_correction_estimates) - true_choose_k)**2))
                
                naive_mean = np.mean(naive_estimates)
                naive_bias = naive_mean - true_choose_k
                naive_rmse = np.sqrt(np.mean((np.array(naive_estimates) - true_choose_k)**2))
                
                # Calculate ratios
                bias_ratio = abs(naive_bias)/abs(mc_bias) if mc_bias != 0 else float('inf')
                rmse_ratio = naive_rmse/mc_rmse if mc_rmse != 0 else float('inf')
                
                # Print results for this K
                print(f"    K={K:2d}: True={true_choose_k:8.0f} | MC Bias={mc_bias:8.2f}, Naive Bias={naive_bias:8.2f} | Bias Ratio={bias_ratio:6.1f}x | RMSE Ratio={rmse_ratio:5.2f}x")
                
                results.append({
                    'true_f': true_f,
                    'K': K,
                    'variance': variance,
                    'true_choose_k': true_choose_k,
                    'mc_bias': mc_bias,
                    'mc_rmse': mc_rmse,
                    'naive_bias': naive_bias,
                    'naive_rmse': naive_rmse,
                    'bias_ratio': bias_ratio,
                    'rmse_ratio': rmse_ratio
                })
    
    return results

def plot_wide_range_results(results):
    """
    Plot results for wide range of K values
    """
    # Extract data for plotting
    true_f_vals = sorted(list(set([r['true_f'] for r in results])))
    variance_vals = sorted(list(set([r['variance'] for r in results])))
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Moment Correction Performance Across Wide Range of K Values', fontsize=16)
    
    # Plot 1: Bias improvement vs K
    ax1 = axes[0, 0]
    for f in true_f_vals:
        for var in variance_vals:
            data = [r for r in results if r['true_f'] == f and r['variance'] == var]
            if data:
                K_vals = [r['K'] for r in data]
                bias_ratios = [r['bias_ratio'] for r in data if r['bias_ratio'] != float('inf')]
                K_vals_filtered = [K_vals[i] for i, r in enumerate(data) if r['bias_ratio'] != float('inf')]
                
                ax1.plot(K_vals_filtered, bias_ratios, 'o-', label=f'f={f}, σ²={var}', alpha=0.7, markersize=4)
    
    ax1.set_xlabel('K')
    ax1.set_ylabel('Bias Improvement Ratio (Naive/MC)')
    ax1.set_title('Bias Improvement vs K')
    ax1.set_yscale('log')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: RMSE improvement vs K
    ax2 = axes[0, 1]
    for f in true_f_vals:
        for var in variance_vals:
            data = [r for r in results if r['true_f'] == f and r['variance'] == var]
            if data:
                K_vals = [r['K'] for r in data]
                rmse_ratios = [r['rmse_ratio'] for r in data]
                
                ax2.plot(K_vals, rmse_ratios, 'o-', label=f'f={f}, σ²={var}', alpha=0.7, markersize=4)
    
    ax2.set_xlabel('K')
    ax2.set_ylabel('RMSE Improvement Ratio (Naive/MC)')
    ax2.set_title('RMSE Improvement vs K')
    ax2.axhline(y=1, color='r', linestyle='--', alpha=0.5, label='No improvement')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Absolute bias comparison
    ax3 = axes[1, 0]
    for f in true_f_vals:
        for var in variance_vals:
            data = [r for r in results if r['true_f'] == f and r['variance'] == var]
            if data:
                K_vals = [r['K'] for r in data]
                mc_bias_abs = [abs(r['mc_bias']) for r in data]
                naive_bias_abs = [abs(r['naive_bias']) for r in data]
                
                ax3.plot(K_vals, mc_bias_abs, 'o-', label=f'MC f={f}, σ²={var}', alpha=0.7, markersize=4)
                ax3.plot(K_vals, naive_bias_abs, 's--', label=f'Naive f={f}, σ²={var}', alpha=0.7, markersize=4)
    
    ax3.set_xlabel('K')
    ax3.set_ylabel('Absolute Bias')
    ax3.set_title('Absolute Bias Comparison')
    ax3.set_yscale('log')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: RMSE comparison
    ax4 = axes[1, 1]
    for f in true_f_vals:
        for var in variance_vals:
            data = [r for r in results if r['true_f'] == f and r['variance'] == var]
            if data:
                K_vals = [r['K'] for r in data]
                mc_rmse = [r['mc_rmse'] for r in data]
                naive_rmse = [r['naive_rmse'] for r in data]
                
                ax4.plot(K_vals, mc_rmse, 'o-', label=f'MC f={f}, σ²={var}', alpha=0.7, markersize=4)
                ax4.plot(K_vals, naive_rmse, 's--', label=f'Naive f={f}, σ²={var}', alpha=0.7, markersize=4)
    
    ax4.set_xlabel('K')
    ax4.set_ylabel('RMSE')
    ax4.set_title('RMSE Comparison')
    ax4.set_yscale('log')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('/data/yizhangh/ldp-pq/wide_range_K_results.png', dpi=300, bbox_inches='tight')
    plt.show()

def analyze_wide_range_results(results):
    """
    Analyze results across wide range of K values
    """
    print(f"\n=== Analysis Across Wide Range of K Values ===")
    
    # Group by K ranges
    small_K = [r for r in results if 2 <= r['K'] <= 5]
    medium_K = [r for r in results if 6 <= r['K'] <= 10]
    large_K = [r for r in results if 11 <= r['K'] <= 20]
    
    def analyze_group(group, name):
        if not group:
            return
        bias_ratios = [r['bias_ratio'] for r in group if r['bias_ratio'] != float('inf')]
        rmse_ratios = [r['rmse_ratio'] for r in group]
        
        mc_better_bias = sum(1 for r in group if abs(r['mc_bias']) < abs(r['naive_bias']))
        mc_better_rmse = sum(1 for r in group if r['mc_rmse'] < r['naive_rmse'])
        
        print(f"\n{name} (K range):")
        print(f"  Cases: {len(group)}")
        print(f"  Avg bias improvement: {np.mean(bias_ratios):.1f}x" if bias_ratios else "  Avg bias improvement: N/A")
        print(f"  Avg RMSE improvement: {np.mean(rmse_ratios):.2f}x")
        print(f"  MC better bias: {mc_better_bias}/{len(group)} ({mc_better_bias/len(group)*100:.1f}%)")
        print(f"  MC better RMSE: {mc_better_rmse}/{len(group)} ({mc_better_rmse/len(group)*100:.1f}%)")
    
    analyze_group(small_K, "Small K")
    analyze_group(medium_K, "Medium K")
    analyze_group(large_K, "Large K")
    
    # Overall statistics
    all_bias_ratios = [r['bias_ratio'] for r in results if r['bias_ratio'] != float('inf')]
    all_rmse_ratios = [r['rmse_ratio'] for r in results]
    
    mc_better_bias_all = sum(1 for r in results if abs(r['mc_bias']) < abs(r['naive_bias']))
    mc_better_rmse_all = sum(1 for r in results if r['mc_rmse'] < r['naive_rmse'])
    
    print(f"\nOverall (All K values):")
    print(f"  Total cases: {len(results)}")
    print(f"  Avg bias improvement: {np.mean(all_bias_ratios):.1f}x")
    print(f"  Avg RMSE improvement: {np.mean(all_rmse_ratios):.2f}x")
    print(f"  MC better bias: {mc_better_bias_all}/{len(results)} ({mc_better_bias_all/len(results)*100:.1f}%)")
    print(f"  MC better RMSE: {mc_better_rmse_all}/{len(results)} ({mc_better_rmse_all/len(results)*100:.1f}%)")

if __name__ == "__main__":
    # Run the wide range test
    results = test_wide_range_K()
    
    # Analyze results
    analyze_wide_range_results(results)
    
    # Plot results
    plot_wide_range_results(results)
    
    print(f"\nResults saved to: wide_range_K_results.png")
