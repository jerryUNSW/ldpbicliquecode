# Butterfly (2,2)-Biclique Counting: Lean Implementation

## Overview

This document describes the Lean 4 implementation for formally verifying the butterfly counting formula in the non-interactive DP setting. The butterfly is a special case of (2,2)-biclique.

## The Formula

For a bipartite graph G with privacy budget ε, let G' be the noisy graph after randomized response. The unbiased estimator for butterfly count is:

$$f_{BTF}(G') = \frac{1}{(\mu-1)^4} \left( \mu^4 n'_{\mathcal{M}_6} - \mu^3 n'_{\mathcal{M}_5} + \mu^2 (n'_{\mathcal{M}_4} + n'_{\mathcal{M}_3}) - \mu n'_{\mathcal{M}_2} + n'_{\mathcal{M}_1} \right)$$

where:
- μ = e^ε
- n'_{M_i} = count of motif type M_i in the noisy graph G'
- M₆ is the butterfly (complete 2×2 biclique)

**Theorem**: $\mathbb{E}[f_{BTF}(G')] = BTF_G$ (unbiasedness)

## Implementation Structure

### File: `ButterflyCounting_2_2.lean`

The Lean file is organized into 10 parts:

#### Part 1: Basic Definitions
- `BipartiteGraph`: Structure for bipartite graphs
- `MotifType22`: The 6 distinct motif types (M₁ through M₆)
- `count_motif_type`: Count motifs of each type in a graph
- `butterfly_count`: Count butterflies (M₆) in a graph

#### Part 2: Randomized Response
- `μ(ε)`: Privacy parameter μ = e^ε
- `flip_prob(ε)`: Flipping probability p = 1/(1+μ)
- `keep_prob(ε)`: Keeping probability 1-p = μ/(1+μ)

#### Part 3: Transformation Probabilities
- `transformation_prob_22`: Compute probability that motif M_i becomes M_j
- Algorithm: For each motif, consider all ways to flip 0-4 edges

#### Part 4: Transformation Matrix
- `transformation_matrix_explicit`: The 6×6 matrix T from the paper
- Matrix structure (normalized by 1/(1+μ)⁴):
  ```
  T = 1/(1+μ)⁴ × [
    [μ⁴,     4μ³,     2μ²,     4μ²,     4μ,     1    ],
    [μ³, μ⁴+3μ²,  μ³+μ,  2μ³+2μ, 3μ²+1,   μ    ],
    [μ², 2μ³+2μ,  μ⁴+1,     4μ², 2μ³+2μ,  μ²   ],
    [μ², 2μ³+2μ,   2μ², μ⁴+2μ²+1, 2μ³+2μ,  μ²   ],
    [μ,  3μ²+1,   μ³+μ,  2μ³+2μ, μ⁴+3μ²,  μ³   ],
    [1,     4μ,     2μ²,     4μ²,   4μ³,  μ⁴   ]
  ]
  ```

#### Part 5: Expected Relationship
- `expected_motif_counts`: Theorem stating E[n'_j] = Σ_i T[i][j] · n_i
- This is the key relationship that links true and noisy counts

#### Part 6: Inverse Matrix and Estimator
- `inverse_entry_butterfly`: The last column of T⁻¹
  - T⁻¹[i][5] = (-1)^i · μ^i / (μ-1)⁴
- `butterfly_estimator`: The estimator formula

#### Part 7: Main Unbiasedness Theorem
- `unbiased_butterfly_estimator`: Proves E[estimator] = true count
- This is the main result we want to verify

#### Part 8: Verification of Inverse Structure
- `inverse_structure_correct`: Verifies the closed-form inverse formula
- Proves that T⁻¹ has the claimed structure

#### Part 9: Example Verification
- `transformation_prob_M4_to_M4`: Verifies T[4,4] computation from paper
- Shows how transformation probabilities are computed

#### Part 10: Helper Functions
- `apply_randomized_response`: Placeholder for randomized response mechanism

## The 6 Motif Types

For (2,2) case with 2 upper and 2 lower vertices:

1. **M₁**: 0 edges (empty)
2. **M₂**: 1 edge
3. **M₃**: 2 edges, non-adjacent (share no vertex)
4. **M₄**: 2 edges, adjacent (share one vertex)
5. **M₅**: 3 edges
6. **M₆**: 4 edges (butterfly/complete biclique)

## Proof Strategy

### Step 1: Transformation Matrix Construction
Prove that `transformation_matrix_explicit` correctly computes the transformation probabilities by:
- Enumerating all ways to flip 0-4 edges
- Summing probabilities for each transformation pattern
- Matching the matrix structure from the paper

### Step 2: Expected Relationship
Prove `expected_motif_counts` using:
- Law of total expectation
- Linearity of expectation
- Independence of edge flips in randomized response

### Step 3: Matrix Inversion
Prove `inverse_structure_correct` by:
- Computing T⁻¹ symbolically
- Verifying the last column matches the formula
- Using properties of the transformation matrix structure

### Step 4: Unbiasedness
Prove `unbiased_butterfly_estimator` by:
1. Using `expected_motif_counts` to express E[n'_i] in terms of n_j
2. Substituting into the estimator formula
3. Using T · T⁻¹ = I to show the formula extracts n_{M₆} correctly
4. Simplifying to show E[estimator] = n_{M₆} = BTF_G

## Current Status

### ✅ Completed
- [x] Basic structure and definitions
- [x] Transformation matrix explicit formula
- [x] Estimator formula
- [x] Theorem statements

### 🔄 In Progress
- [ ] Proof of transformation matrix correctness
- [ ] Proof of expected relationship
- [ ] Proof of inverse structure
- [ ] Proof of unbiasedness

### 📋 Next Steps

1. **Complete Transformation Probability Computation**
   - Implement `can_transform_by_flips` to check if M_i can become M_j by flipping x edges
   - Prove that `transformation_prob_22` matches `transformation_matrix_explicit`

2. **Prove Expected Relationship**
   - Formalize randomized response as a probability measure
   - Prove the linear relationship between true and noisy counts

3. **Prove Inverse Structure**
   - Compute T⁻¹ symbolically for the 6×6 case
   - Verify the last column formula

4. **Complete Unbiasedness Proof**
   - Chain together all the lemmas
   - Final verification that the estimator is unbiased

## Example: T[4,4] Computation

From the paper's example, T[4,4] (M₄ to M₄) is computed as:

- **x = 0**: No edges flip → (1-p)⁴ = (μ/(1+μ))⁴
- **x = 2**: Two non-adjacent edges flip → 2 · p²(1-p)² = 2 · (1/(1+μ))²(μ/(1+μ))²
- **x = 4**: All edges flip → p⁴ = (1/(1+μ))⁴

Sum = (μ²+1)²/(1+μ)⁴

This matches the matrix entry: T[4,4] = (μ⁴+2μ²+1)/(1+μ)⁴ = (μ²+1)²/(1+μ)⁴

## Testing Strategy

1. **Small Graph Examples**
   - Test on graphs with known butterfly counts
   - Verify the estimator formula matches expected values

2. **Symbolic Verification**
   - Verify transformation matrix entries symbolically
   - Check inverse matrix computation

3. **Numerical Verification**
   - For specific ε values, compute the matrix and verify properties
   - Check that T · T⁻¹ = I

## References

- **Paper**: `overleaf-paper/non-interactive.tex` - Butterfly counting section
- **Theorem**: Theorem 1 (thm:btf) - Unbiased estimator formula
- **Matrix**: Transformation matrix T (6×6) and its structure
- **Code**: `biclique.cpp` - Implementation reference

## Benefits of This Formalization

1. **Mathematical Rigor**: Prove correctness, not just empirical validation
2. **Reproducibility**: Complete derivation that can be verified by anyone
3. **Foundation**: Base for extending to (2,q) and general (p,q) cases
4. **Research Value**: Formal proofs strengthen the theoretical contribution
5. **Error Detection**: Catch mistakes in combinatorial calculations early

## Notes

- The current implementation uses `sorry` for proofs that need to be completed
- The randomized response mechanism is a placeholder and needs proper probability measure formalization
- Some helper functions (like `can_transform_by_flips`) need implementation
- The file compiles but proofs need to be filled in
