#!/usr/bin/env python3
"""
Compute relative errors for ADV and ADV+ algorithms from NIPS testing results.
"""

import numpy as np

# Ground truth values from the testing results
ground_truth = {
    4: 1.729841e+14,
    5: 2.026990e+15,
    6: 2.054109e+16,
    7: 1.848309e+17,
    8: 1.505491e+18,
    9: 1.126738e+19,
    10: 7.844682e+19
}

# ADV algorithm results
adv_results = {
    4: 4.442019e+14,
    5: 8.447451e+15,
    6: 1.490508e+17,
    7: 2.456005e+18,
    8: 3.857563e+19,
    9: 5.732922e+20,
    10: 8.076765e+21
}

# ADV+ algorithm results
adv_plus_results = {
    4: -2.413750e+14,
    5: 1.220725e+16,
    6: 2.141426e+17,
    7: -3.416424e+18,
    8: -9.742705e+19,
    9: 5.527845e+20,
    10: 3.078710e+22
}

# Naive algorithm results (not available in current test - would need to run separately)
# Note: The test file shows "Naive skipped" - these would need to be obtained from a separate run
naive_results = {
    4: None,  # Would need actual results
    5: None,
    6: None,
    7: None,
    8: None,
    9: None,
    10: None
}

def compute_relative_error(estimated, true_value):
    """Compute relative error as |estimated - true_value| / |true_value| * 100"""
    return abs(estimated - true_value) / abs(true_value) * 100

print("=== Relative Error Analysis for NIPS Dataset (3,Q)-Biclique Counting ===")
print("Epsilon: 1.0")
print("=" * 100)
print(f"{'Q':<3} {'Ground Truth':<15} {'ADV':<15} {'ADV+':<15} {'Naive':<15} {'ADV Rel Err (%)':<15} {'ADV+ Rel Err (%)':<15} {'Naive Rel Err (%)':<15}")
print("-" * 100)

for q in sorted(ground_truth.keys()):
    gt = ground_truth[q]
    adv = adv_results[q]
    adv_plus = adv_plus_results[q]
    naive = naive_results[q]
    
    adv_error = compute_relative_error(adv, gt)
    adv_plus_error = compute_relative_error(adv_plus, gt)
    
    if naive is not None:
        naive_error = compute_relative_error(naive, gt)
        naive_str = f"{naive:<15.2e}"
        naive_error_str = f"{naive_error:<15.2f}"
    else:
        naive_str = "N/A"
        naive_error_str = "N/A"
    
    print(f"{q:<3} {gt:<15.2e} {adv:<15.2e} {adv_plus:<15.2e} {naive_str:<15} {adv_error:<15.2f} {adv_plus_error:<15.2f} {naive_error_str:<15}")

print("\n=== Summary Statistics ===")
adv_errors = [compute_relative_error(adv_results[q], ground_truth[q]) for q in ground_truth.keys()]
adv_plus_errors = [compute_relative_error(adv_plus_results[q], ground_truth[q]) for q in ground_truth.keys()]

print(f"ADV Algorithm:")
print(f"  Mean Relative Error: {np.mean(adv_errors):.2f}%")
print(f"  Median Relative Error: {np.median(adv_errors):.2f}%")
print(f"  Min Relative Error: {np.min(adv_errors):.2f}%")
print(f"  Max Relative Error: {np.max(adv_errors):.2f}%")
print(f"  Std Dev: {np.std(adv_errors):.2f}%")

print(f"\nADV+ Algorithm:")
print(f"  Mean Relative Error: {np.mean(adv_plus_errors):.2f}%")
print(f"  Median Relative Error: {np.median(adv_plus_errors):.2f}%")
print(f"  Min Relative Error: {np.min(adv_plus_errors):.2f}%")
print(f"  Max Relative Error: {np.max(adv_plus_errors):.2f}%")
print(f"  Std Dev: {np.std(adv_plus_errors):.2f}%")

print(f"\n=== Algorithm Comparison ===")
better_count = sum(1 for i in range(len(adv_errors)) if adv_plus_errors[i] < adv_errors[i])
total_count = len(adv_errors)
print(f"ADV+ performs better than ADV in {better_count}/{total_count} cases ({better_count/total_count*100:.1f}%)")

# Check for negative estimates (which indicate algorithm instability)
adv_negative = sum(1 for q in adv_results.keys() if adv_results[q] < 0)
adv_plus_negative = sum(1 for q in adv_plus_results.keys() if adv_plus_results[q] < 0)

print(f"\nNegative Estimates:")
print(f"ADV: {adv_negative}/{len(adv_results)} cases")
print(f"ADV+: {adv_plus_negative}/{len(adv_plus_results)} cases")

# Performance analysis by Q value
print(f"\n=== Performance by Q Value ===")
for q in sorted(ground_truth.keys()):
    adv_err = compute_relative_error(adv_results[q], ground_truth[q])
    adv_plus_err = compute_relative_error(adv_plus_results[q], ground_truth[q])
    winner = "ADV+" if adv_plus_err < adv_err else "ADV"
    print(f"Q={q}: ADV={adv_err:.1f}%, ADV+={adv_plus_err:.1f}% (Winner: {winner})")

print(f"\n=== Note on Naive Algorithm ===")
print("The Naive algorithm was skipped in this test run.")
print("To include Naive in the comparison, you would need to:")
print("1. Run the test with Naive algorithm enabled")
print("2. Update the naive_results dictionary with the actual results")
print("3. Re-run this analysis script")
