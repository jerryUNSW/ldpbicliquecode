# Sampling Strategy Analysis: Degree-Proportional vs Low-Degree Prioritization

## The Problem You've Identified

You're right to be confused - there's a fundamental tension:

**Option A: Sample high-degree pairs more often**
- Intuition: High-degree pairs have high variance, so sample them more
- Problem: We still need to compute the estimator for each sample
- Cost: Same computational cost per pair regardless of degree
- Benefit: More samples of high-variance pairs → better estimates for them

**Option B: Sample low-degree pairs more often**
- Intuition: Low-degree pairs already have low variance, don't need many samples
- Problem: Wastes samples on pairs that don't need refinement
- Cost: Same computational cost per pair
- Benefit: Allocate limited samples to where they help most

## The Key Insight: What Are We Optimizing?

The answer depends on what you want to minimize:

### Scenario 1: Minimize Total Variance (All Pairs)
```
Total Variance = Σ variance(u,w) for all pairs

If we have T total samples to allocate:
- High-degree pairs: high variance, benefit most from extra samples
- Low-degree pairs: low variance, already accurate

STRATEGY: Sample high-degree pairs more often
REASON: Reduces the largest contributors to total variance
```

### Scenario 2: Minimize Maximum Error (Worst-Case)
```
Max Error = max(error(u,w)) for all pairs

If we have T total samples:
- High-degree pairs: already high variance, hard to fix
- Low-degree pairs: already low error, don't need help

STRATEGY: Sample low-degree pairs more often (or not at all)
REASON: Ensures all pairs meet minimum accuracy threshold
```

### Scenario 3: Minimize Computational Cost for Target Accuracy
```
Goal: All pairs have error < threshold

If we have T total samples:
- High-degree pairs: need many samples to reach threshold
- Low-degree pairs: need few samples to reach threshold

STRATEGY: Adaptive sampling based on current error
REASON: Allocate samples to pairs that haven't reached target yet
```

## The Real Issue: Computational Cost

You're right that computing the estimator is expensive. But here's the key:

**The cost is NOT in computing the estimator itself** - it's in:
1. Sampling the pair (u,w)
2. Running the algorithm on that pair
3. Computing the estimate

If you sample a pair multiple times, you're running the algorithm multiple times on that pair, which is expensive.

## Recommended Strategy: Adaptive Sampling

Instead of uniform or degree-proportional sampling, use **adaptive sampling**:

```
Algorithm:
1. Initial round: Sample all pairs once (or random subset)
2. Compute variance for each pair
3. Identify pairs with high variance
4. Allocate remaining samples to high-variance pairs
5. Repeat until budget exhausted

Pseudocode:
  variance_estimate = {}
  samples_per_pair = {}
  
  # Phase 1: Initial sampling
  for each pair (u,w):
    estimate = run_algorithm(u, w)
    variance_estimate[u,w] = estimate_variance(estimate)
    samples_per_pair[u,w] = 1
  
  # Phase 2: Adaptive refinement
  remaining_budget = T - num_pairs
  while remaining_budget > 0:
    # Find pair with highest variance
    worst_pair = argmax(variance_estimate)
    
    # Sample it again
    estimate = run_algorithm(worst_pair)
    variance_estimate[worst_pair] = update_variance(estimate)
    samples_per_pair[worst_pair] += 1
    remaining_budget -= 1
```

## Why This Works Better

1. **Targets high-variance pairs**: Focuses computational budget where it helps most
2. **Adaptive**: Responds to actual variance, not just degree
3. **Efficient**: Doesn't waste samples on already-accurate pairs
4. **Practical**: Can stop early if variance is acceptable

## Connection to Your Degree-Product Finding

Your degree product finding is still valuable:

```
Degree product tells us:
- Which pairs WILL LIKELY have high variance (prediction)
- Where to focus initial sampling effort
- Which pairs to prioritize in adaptive refinement

But it's not perfect because:
- Actual variance depends on graph structure, not just degrees
- Some high-degree pairs might have low variance (sparse neighborhoods)
- Some low-degree pairs might have high variance (dense neighborhoods)
```

## Practical Implementation

### Option 1: Degree-Proportional Sampling (Simple)
```cpp
// Sample pairs with probability proportional to degree product
for (int sample = 0; sample < T; sample++) {
    // Sample pair (u,w) with probability ∝ degree[u] * degree[w]
    pair<int,int> (u,w) = sample_pair_by_degree_product();
    
    // Run algorithm and accumulate estimate
    estimate += run_algorithm(u, w);
}

// Average over all samples
final_estimate = estimate / T;
```

**Pros**: Simple, focuses on high-variance pairs
**Cons**: Might oversample some pairs, undersample others

### Option 2: Adaptive Sampling (Better)
```cpp
// Phase 1: Sample all pairs once
map<pair<int,int>, vector<double>> estimates;
for (auto& pair : all_pairs) {
    estimates[pair].push_back(run_algorithm(pair));
}

// Phase 2: Adaptively refine high-variance pairs
for (int sample = 0; sample < remaining_budget; sample++) {
    // Find pair with highest variance
    auto worst_pair = find_highest_variance_pair(estimates);
    
    // Sample it again
    estimates[worst_pair].push_back(run_algorithm(worst_pair));
}

// Compute final estimates with confidence intervals
for (auto& pair : all_pairs) {
    final_estimate[pair] = mean(estimates[pair]);
    confidence[pair] = std_dev(estimates[pair]);
}
```

**Pros**: Targets actual variance, efficient, adaptive
**Cons**: More complex, requires multiple rounds

## Answer to Your Confusion

**You should NOT prioritize low-degree pairs.** Here's why:

1. **Low-degree pairs already have low variance** - they don't need refinement
2. **High-degree pairs have high variance** - they benefit most from extra samples
3. **Computational budget is limited** - allocate to where it helps most

**BUT** - don't use pure degree-proportional sampling either. Use **adaptive sampling** that:
- Starts with degree product as a hint
- Measures actual variance
- Allocates remaining budget to highest-variance pairs

## Summary

| Strategy | Pros | Cons | When to Use |
|----------|------|------|------------|
| Uniform sampling | Simple | Wastes budget on low-variance pairs | Baseline only |
| Degree-proportional | Targets high-variance | Doesn't adapt to actual variance | Quick approximation |
| Adaptive sampling | Optimal budget allocation | More complex | Production code |

**Recommendation**: Start with degree-proportional (simple), then upgrade to adaptive sampling (better results).
