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
        # For normal distribution: E[X^k] = sum_{i=0}^{k/2} C(k,2i) * (2i-1)!! * μ^(k-2i) * σ^(2i)
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

def test_moment_correction():
    """
    Test the accuracy of moment correction vs naive approach
    """
    print("=== Testing Moment Correction for C(f, K) Estimation ===\n")
    
    # Test parameters
    true_f_values = [5, 10, 15, 20, 25, 30]
    K_values = [2, 3, 4, 5, 6]
    variance_values = [0.5, 1.0, 2.0, 4.0]
    
    results = []
    
    for true_f in true_f_values:
        print(f"Testing with true f = {true_f}")
        print("-" * 50)
        
        for K in K_values:
            if K > true_f:  # Skip if K > f (would be 0)
                continue
                
            true_choose_k = binomial(true_f, K)
            print(f"  K={K}: True C({true_f},{K}) = {true_choose_k}")
            
            for variance in variance_values:
                # Generate many samples of f' ~ N(f, variance)
                n_samples = 10000
                f_prime_samples = np.random.normal(true_f, np.sqrt(variance), n_samples)
                
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
                
                print(f"    σ²={variance}:")
                print(f"      Moment Correction: bias={mc_bias:.3f}, RMSE={mc_rmse:.3f}")
                print(f"      Naive:           bias={naive_bias:.3f}, RMSE={naive_rmse:.3f}")
                print(f"      Improvement:     bias_ratio={abs(naive_bias)/abs(mc_bias):.2f}x, RMSE_ratio={naive_rmse/mc_rmse:.2f}x")
                
                results.append({
                    'true_f': true_f,
                    'K': K,
                    'variance': variance,
                    'true_choose_k': true_choose_k,
                    'mc_bias': mc_bias,
                    'mc_rmse': mc_rmse,
                    'naive_bias': naive_bias,
                    'naive_rmse': naive_rmse
                })
        
        print()
    
    return results

def plot_results(results):
    """
    Plot the results comparing moment correction vs naive approach
    """
    # Extract data for plotting
    variances = sorted(list(set([r['variance'] for r in results])))
    K_vals = sorted(list(set([r['K'] for r in results])))
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Moment Correction vs Naive C(f,K) Estimation', fontsize=16)
    
    # Plot 1: Bias comparison
    ax1 = axes[0, 0]
    for K in K_vals:
        mc_biases = [r['mc_bias'] for r in results if r['K'] == K]
        naive_biases = [r['naive_bias'] for r in results if r['K'] == K]
        variances_plot = [r['variance'] for r in results if r['K'] == K]
        
        ax1.plot(variances_plot, np.abs(mc_biases), 'o-', label=f'MC K={K}', alpha=0.7)
        ax1.plot(variances_plot, np.abs(naive_biases), 's--', label=f'Naive K={K}', alpha=0.7)
    
    ax1.set_xlabel('Variance (σ²)')
    ax1.set_ylabel('Absolute Bias')
    ax1.set_title('Bias Comparison')
    ax1.set_yscale('log')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: RMSE comparison
    ax2 = axes[0, 1]
    for K in K_vals:
        mc_rmse = [r['mc_rmse'] for r in results if r['K'] == K]
        naive_rmse = [r['naive_rmse'] for r in results if r['K'] == K]
        variances_plot = [r['variance'] for r in results if r['K'] == K]
        
        ax2.plot(variances_plot, mc_rmse, 'o-', label=f'MC K={K}', alpha=0.7)
        ax2.plot(variances_plot, naive_rmse, 's--', label=f'Naive K={K}', alpha=0.7)
    
    ax2.set_xlabel('Variance (σ²)')
    ax2.set_ylabel('RMSE')
    ax2.set_title('RMSE Comparison')
    ax2.set_yscale('log')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Improvement ratio (bias)
    ax3 = axes[1, 0]
    for K in K_vals:
        bias_ratios = [abs(r['naive_bias'])/abs(r['mc_bias']) if r['mc_bias'] != 0 else np.nan 
                      for r in results if r['K'] == K]
        variances_plot = [r['variance'] for r in results if r['K'] == K]
        
        ax3.plot(variances_plot, bias_ratios, 'o-', label=f'K={K}', alpha=0.7)
    
    ax3.set_xlabel('Variance (σ²)')
    ax3.set_ylabel('Bias Improvement Ratio (Naive/MC)')
    ax3.set_title('Bias Improvement from Moment Correction')
    ax3.axhline(y=1, color='r', linestyle='--', alpha=0.5, label='No improvement')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Improvement ratio (RMSE)
    ax4 = axes[1, 1]
    for K in K_vals:
        rmse_ratios = [r['naive_rmse']/r['mc_rmse'] for r in results if r['K'] == K]
        variances_plot = [r['variance'] for r in results if r['K'] == K]
        
        ax4.plot(variances_plot, rmse_ratios, 'o-', label=f'K={K}', alpha=0.7)
    
    ax4.set_xlabel('Variance (σ²)')
    ax4.set_ylabel('RMSE Improvement Ratio (Naive/MC)')
    ax4.set_title('RMSE Improvement from Moment Correction')
    ax4.axhline(y=1, color='r', linestyle='--', alpha=0.5, label='No improvement')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('/data/yizhangh/ldp-pq/moment_correction_accuracy.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Run the test
    results = test_moment_correction()
    
    # Plot results
    plot_results(results)
    
    # Summary statistics
    print("\n=== Summary ===")
    mc_bias_improvements = [abs(r['naive_bias'])/abs(r['mc_bias']) for r in results if r['mc_bias'] != 0]
    mc_rmse_improvements = [r['naive_rmse']/r['mc_rmse'] for r in results]
    
    print(f"Moment Correction Bias Improvement: {np.mean(mc_bias_improvements):.2f}x (avg)")
    print(f"Moment Correction RMSE Improvement: {np.mean(mc_rmse_improvements):.2f}x (avg)")
    print(f"Cases where MC has lower bias: {sum(1 for r in results if abs(r['mc_bias']) < abs(r['naive_bias']))}/{len(results)}")
    print(f"Cases where MC has lower RMSE: {sum(1 for r in results if r['mc_rmse'] < r['naive_rmse'])}/{len(results)}")
