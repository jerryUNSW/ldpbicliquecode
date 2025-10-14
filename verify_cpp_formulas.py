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

def compute_local_res_cpp_style(K, f_prime, variance):
    """
    Exact implementation of the C++ compute_local_res function
    """
    moment = [0] * (max(K, 10) + 1)  # Ensure enough space for all moments
    moment[1] = f_prime
    moment[2] = f_prime**2 - variance

    if K == 2:
        return (moment[2] - moment[1]) / 2

    if K >= 3:
        moment[3] = f_prime**3 - 3 * f_prime * variance
        if K == 3:
            return (moment[3] - 3 * moment[2] + 2 * moment[1]) / 6

    if K >= 4:
        moment[4] = f_prime**4 - 6 * moment[2] * variance - 3 * variance**2
        if K == 4:
            return (moment[4] - 6 * moment[3] + 11 * moment[2] - 6 * moment[1]) / 24

    if K >= 5:
        moment[5] = f_prime**5 - 10 * moment[3] * variance - 15 * f_prime * variance**2
        if K == 5:
            return (moment[5] - 10 * moment[4] + 35 * moment[3] - 50 * moment[2] + 24 * moment[1]) / 120

    if K >= 6:
        moment[6] = f_prime**6 - 15 * moment[4] * variance - 45 * moment[2] * variance**2 - 15 * variance**3
        if K == 6:
            return (moment[6] - 15 * moment[5] + 85 * moment[4] - 225 * moment[3] + 274 * moment[2] - 120 * moment[1]) / 720

    if K >= 7:
        moment[7] = f_prime**7 + 21 * moment[5] * variance + 105 * moment[3] * variance**2 + 105 * f_prime * variance**3
        if K == 7:
            return (moment[7] - 21 * moment[6] + 105 * moment[5] - 210 * moment[4] + 252 * moment[3] - 140 * moment[2] + 24 * moment[1]) / 5040

    if K >= 8:
        moment[8] = f_prime**8 - 28 * moment[6] * variance - 140 * moment[4] * variance**2 - 210 * moment[2] * variance**3 - 105 * variance**4
        if K == 8:
            return (moment[8] - 28 * moment[7] + 140 * moment[6] - 364 * moment[5] + 560 * moment[4] - 560 * moment[3] + 336 * moment[2] - 70 * moment[1]) / 40320

    if K >= 9:
        moment[9] = f_prime**9 - 36 * moment[7] * variance - 210 * moment[5] * variance**2 - 420 * moment[3] * variance**3 - 315 * moment[1] * variance**4
        if K == 9:
            return (moment[9] - 36 * moment[8] + 168 * moment[7] - 504 * moment[6] + 1260 * moment[5] - 2520 * moment[4] + 3024 * moment[3] - 2016 * moment[2] + 504 * moment[1]) / 362880

    if K == 10:
        moment[10] = f_prime**10
        return (moment[10] - 45 * moment[9] + 210 * moment[8] - 630 * moment[7] + 1260 * moment[6] - 2520 * moment[5] + 3024 * moment[4] - 2520 * moment[3] + 1260 * moment[2] - 210 * moment[1]) / 3628800

    return 0

if __name__ == "__main__":
    print("=== C++ Implementation Formulas for K=1 to 10 ===\n")
    
    for K in range(1, 11):
        formula = compute_local_res_cpp_style(K, fp, s)
        expanded = expand(formula)
        print(f"K={K}: {expanded}")
        print()
