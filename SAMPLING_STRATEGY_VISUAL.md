# Sampling Strategy Comparison: Visual Guide

## The Core Question

You have a budget of T samples to allocate across all (u,w) pairs. How should you allocate them?

## Visual Comparison

### Strategy 1: Uniform Sampling
```
Pair variance:  [LOW]  [LOW]  [HIGH] [HIGH] [HIGH]
Samples:        [1]    [1]    [1]    [1]    [1]
Result:         ✓      ✓      ✗      ✗      ✗
                Good   Good   Bad    Bad    Bad
```
**Problem**: Wastes samples on already-accurate pairs

### Strategy 2: Degree-Proportional Sampling
```
Pair variance:  [LOW]  [LOW]  [HIGH] [HIGH] [HIGH]
Degree product: [10]   [15]   [100]  [120]  [150]
Samples:        [1]    [1]    [3]    [4]    [5]
Result:         ✓      ✓      ✓      ✓      ✓
                Good   Good   Good   Good   Good
```
**Benefit**: Allocates more samples to high-variance pairs

### Strategy 3: Adaptive Sampling (Best)
```
Round 1:
Pair variance:  [LOW]  [LOW]  [HIGH] [HIGH] [HIGH]
Samples:        [1]    [1]    [1]    [1]    [1]
Measured var:   [0.1]  [0.2]  [2.5]  [2.8]  [3.2]

Round 2 (allocate remaining budget to highest variance):
Pair variance:  [LOW]  [LOW]  [HIGH] [HIGH] [HIGH]
Additional:     [0]    [0]    [2]    [2]    [3]
Total samples:  [1]    [1]    [3]    [3]    [4]
Result:         ✓      ✓      ✓      ✓      ✓
                Good   Good   Better Better Better
```
**Benefit**: Responds to actual measured variance, not just degree

## Why NOT Low-Degree Prioritization?

```
❌ WRONG: Prioritize low-degree pairs
Pair variance:  [LOW]  [LOW]  [HIGH] [HIGH] [HIGH]
Samples:        [5]    [4]    [1]    [1]    [1]
Result:         ✓✓✓    ✓✓✓    ✗      ✗      ✗
                Overkill Overkill Bad   Bad   Bad

Problem: Wastes budget on pairs that don't need it
```

## The Math Behind It

### Total Variance Minimization

If you have T samples to allocate:

```
Total Variance = Σ variance(u,w)

For each pair, variance decreases with more samples:
variance(u,w) with n samples ≈ base_variance(u,w) / √n

To minimize total variance:
- Allocate more samples to pairs with high base_variance
- Allocate fewer samples to pairs with low base_variance

This is exactly what degree-proportional sampling does!
```

### Optimal Allocation (Neyman Allocation)

```
Optimal samples for pair (u,w):
n(u,w) = T × √(variance(u,w)) / Σ√(variance(u,w))

This means:
- Pairs with high variance get more samples
- Pairs with low variance get fewer samples
- Budget is allocated optimally

Degree product is a good proxy for variance!
```

## Practical Example

Suppose you have 100 samples to allocate across 5 pairs:

```
Pair    Degree  Degree_Product  Estimated_Variance  Optimal_Samples
(u1,w1) 5       50              0.5                 5
(u2,w2) 8       80              0.8                 8
(u3,w3) 20      200             2.0                 20
(u4,w4) 25      250             2.5                 25
(u5,w5) 30      300             3.0                 30
        ---     ---             ---                 ---
Total                           8.8                 88 (+ 12 for rounding)

Allocation:
- Low-degree pairs (u1,w1), (u2,w2): 5-8 samples each
- High-degree pairs (u3,w3), (u4,w4), (u5,w5): 20-30 samples each
```

## Implementation Recommendation

### Phase 1: Use Degree-Proportional (Simple)
```python
# Sample pairs with probability ∝ degree product
total_samples = 100
for _ in range(total_samples):
    # Sample pair (u,w) with probability ∝ degree[u] * degree[w]
    u, w = sample_pair_by_degree_product()
    estimate = run_algorithm(u, w)
    accumulate_estimate(u, w, estimate)

# Average estimates
final_estimate = compute_average()
```

**Pros**: Simple, effective, uses your degree product finding
**Cons**: Might not be perfectly optimal

### Phase 2: Upgrade to Adaptive (Better)
```python
# Phase 1: Initial sampling
estimates = {}
for u, w in all_pairs:
    estimates[u,w] = [run_algorithm(u, w)]

# Phase 2: Adaptive refinement
remaining_budget = 50  # Allocate 50% of budget to refinement
for _ in range(remaining_budget):
    # Find pair with highest variance
    worst_pair = max(estimates.keys(), 
                     key=lambda p: variance(estimates[p]))
    
    # Sample it again
    estimates[worst_pair].append(run_algorithm(worst_pair))

# Compute final estimates
final_estimate = {p: mean(estimates[p]) for p in all_pairs}
```

**Pros**: Optimal allocation, responds to actual variance
**Cons**: More complex, requires two phases

## Answer to Your Confusion

**Q: Should we prioritize low-degree pairs?**
**A: No. Prioritize high-degree pairs.**

**Q: Why?**
**A: Because high-degree pairs have high variance and benefit most from extra samples.**

**Q: But we compute the same estimator for each pair...**
**A: Exactly! So allocate your computational budget to pairs that need it most (high-variance ones).**

**Q: How do we know which pairs have high variance?**
**A: Use degree product as a proxy (your finding!), or measure it adaptively.**

## Key Takeaway

Your degree product finding is valuable because it tells you:
- **Which pairs to sample more often** (high degree product = high variance)
- **Which pairs to sample less often** (low degree product = low variance)

This is exactly the right intuition. Don't second-guess it!
