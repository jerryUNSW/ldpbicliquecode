# Bug Report: Algorithm Switch Logic Not Implemented in Batch Functions

## Summary
The wedge-based 3K biclique counting algorithms have a critical bug where the `multi_estimator_switch` logic is not properly implemented in batch and rejection sampling versions. This causes Algorithm 2 (ADV) to behave identically to Algorithm 3 (ADV+), making the algorithm comparison results invalid.

## Affected Functions
The following functions have the bug:
1. `wedge_based_two_round_3_K_biclique_no_sampling_batch` (line 1471)
2. `wedge_based_two_round_3_K_biclique_rejection_sampling` (line 1704)  
3. `wedge_based_two_round_3_K_biclique_rejection_sampling_batch` (line 1918)

## Root Cause
All affected functions are missing the `if(multi_estimator_switch)` conditional structure entirely. They always compute f1, f2, f3 and use the multi-estimator approach regardless of the switch setting:

```cpp
// Always computes f1, f2, f3 (missing multi_estimator_switch check)
// Always averages them
fuvw = (f1 + f2 + f3)/3;

// Always uses multi-estimator variance calculation
bool improvement = true;
if(improvement){
    // Multi-estimator variance calculation
}
```

The `bool improvement = true;` flag itself is not problematic - it's just an internal optimization within the multi-estimator approach. The real issue is the complete absence of the `multi_estimator_switch` conditional logic.

## Expected Behavior
- **Algorithm 2 (ADV)**: `multi_estimator_switch = false` → Should use single estimator (f1 only)
- **Algorithm 3 (ADV+)**: `multi_estimator_switch = true` → Should use multi estimator (f1, f2, f3 averaged)
- **Algorithm 4 (ADV++)**: `multi_estimator_switch = true` + `two_noisy_graph_switch = true` → Should use multi estimator + two noisy graphs

## Actual Behavior
All algorithms run the same multi-estimator code regardless of the `multi_estimator_switch` setting, causing:
- Algorithm 2 to perform similarly to Algorithm 3
- Algorithm 3 and 4 to work correctly (but Algorithm 2 is wrong)
- Invalid algorithm comparison results

## Correct Implementation
Only `wedge_based_two_round_3_K_biclique` (the main non-batch function) correctly implements the switch logic:

```cpp
if(multi_estimator_switch){
    // Multi-estimator: compute f1, f2, f3 and average
    fuvw = (f1 + f2 + f3)/3;
    // Complex variance calculation with cross-terms
} else {
    // Single-estimator: only compute f1
    fuvw = f1;  // No averaging!
    // Simple variance calculation
}
```

## Impact
- **Algorithm comparison results are invalid** for p=3 experiments
- **Algorithm 2 (ADV) performance is artificially inflated** due to running multi-estimator code
- **Research conclusions may be incorrect** due to unfair algorithm comparison

## Files Affected
- `biclique.cpp` (lines 1471, 1704, 1918)
- `test_p3_batch_with_ground_truth.cpp` (uses the broken batch function)
- All p=3 algorithm comparison results and plots

## Fix Required
Add the missing `if(multi_estimator_switch)` conditional structure to all three affected functions, following this pattern:

```cpp
if(multi_estimator_switch){
    // Multi-source estimator logic
    // Compute f1, f2, f3
    fuvw = (f1 + f2 + f3)/3;
    
    bool improvement = true;  // This stays inside the multi-estimator branch
    if(improvement){
        // Multi-estimator variance calculation
    }
} else {
    // Single-source estimator logic
    // Only compute f1
    fuvw = f1;
    
    // Simple variance calculation (no improvement flag needed)
}
```

## Test Results Affected
- `testing-csbwiki-no-sampling.txt` - Results are invalid due to this bug
- `algorithm_comparison_csbwiki_eps1_FIXED_p3.pdf` - Plot shows incorrect algorithm comparison
- Database entries in `algorithm_comparison_p3.db` - Need to be regenerated after fix

## Date Discovered
October 18, 2024

## Severity
**CRITICAL** - Invalidates all p=3 algorithm comparison results and research conclusions.
