#!/usr/bin/env python3

import numpy as np
import math
from sympy import symbols, expand, simplify

# Define symbols: fp = f' (noisy estimate), s = variance
fp, s = symbols('fp s')

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

def estimate_choose_K_moment_correction_symbolic(K, f_prime, variance):
    """
    Estimate C(f, K) using moment correction - symbolic version
    This follows the exact same logic as the numeric version in test_wide_range_K.py
    """
    if K == 0:
        return 1.0
    if K == 1:
        return f_prime
    
    # Construct corrected moments up to K
    moments = [0] * (K + 1)
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
    
    return simplify(result / factorial(K))

if __name__ == "__main__":
    print("=== Moment Correction Formulas for K=1 to 10 ===")
    print("Using the recursive implementation from test_wide_range_K.py\n")
    
    for K in range(1, 11):
        formula = estimate_choose_K_moment_correction_symbolic(K, fp, s)
        expanded = expand(formula)
        print(f"K={K}: {expanded}")
        print()
