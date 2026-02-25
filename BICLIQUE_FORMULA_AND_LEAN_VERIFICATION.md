# Mathematical Formula for (p,q)-Biclique Estimation in Non-Interactive DP

## Overview

This document describes the mathematical formula for estimating (p,q)-biclique counts in a non-interactive differential privacy setting, and discusses how Lean 4 can be used to formally verify the derivation.

## The Core Formula

### General (p,q)-Biclique Formula

For a bipartite graph $G$ with privacy budget $\varepsilon$, let $G'_{\varepsilon}$ be the noisy graph after applying randomized response. The unbiased estimator for the number of $(p,q)$-bicliques is given by:

**Theorem (from `3oneroundapproach.tex`):**

$$\hat{n}_{\mathcal{B}_{pq}} = \frac{ \sum_{i=0}^{p \times q} (-1)^{i} \mu^i n_{\mathcal{B}_i}' } { (1-\mu)^{p \times q}}$$

where:
- $\mu = e^{\varepsilon}$ (the privacy parameter)
- $n_{\mathcal{B}_i}'$ is the count of motifs with exactly $i$ edges in the noisy graph $G'_{\varepsilon}$
- $\mathcal{B}_{pq}$ represents the complete $(p,q)$-biclique (with all $pq$ edges present)
- The formula satisfies: $\mathbb{E}[\hat{n}_{\mathcal{B}_{pq}}] = C_{p,q}(G)$ (unbiasedness)

### Special Case: (2,q)-Biclique

For the special case where $p=2$, the formula simplifies to:

$$\hat{n}_{\mathcal{B}_{2q}} = \frac{1}{(1-\mu)^{2q}} \sum_{i=0}^{2q} (-\mu)^i \cdot m_i$$

where $m_i$ counts motifs with $i$ edges between two upper vertices and $q$ lower vertices in $G'_{\varepsilon}$.

This is implemented in the code at `biclique.cpp:3254-3258`:
```cpp
long double res = 0, mu = std::exp(Eps);
for (size_t i = 0; i <= 2 * K; ++i) {
    res += std::pow(-mu, i) * m__[i];
}
res /= std::pow(1 - mu, 2 * K);
```

## Mathematical Derivation

### Step 1: Transformation Probability Matrix

The derivation is based on the **transformation probability matrix** $\mathcal{T}$, where:
- $\mathcal{T}[i,j] = TP(\mathcal{B}_i, \mathcal{B}_j)$ = probability that a motif with $i$ edges in $G$ becomes a motif with $j$ edges in $G'$
- For randomized response with flipping probability $p = \frac{1}{1+e^{\varepsilon}}$:
  - An edge flips with probability $p$
  - An edge stays the same with probability $1-p = \frac{\mu}{1+\mu}$ where $\mu = e^{\varepsilon}$

### Step 2: Linear System

The relationship between true and noisy motif counts is:

$$\mathbb{E}(n_{\mathcal{B}_0}', n_{\mathcal{B}_1}', \cdots , n_{\mathcal{B}_{pq}}') = (n_{\mathcal{B}_0}, n_{\mathcal{B}_1}, \cdots , n_{\mathcal{B}_{pq}}) \times \mathcal{T}$$

### Step 3: Inverse Transformation

To recover the true counts, we invert the transformation matrix:

$$(\hat{n}_{\mathcal{B}_0}, \hat{n}_{\mathcal{B}_{1}},  \cdots , \hat{n}_{\mathcal{B}_{pq}}) = (n_{\mathcal{B}_0}', n_{\mathcal{B}_1}', \cdots , n_{\mathcal{B}_{pq}}') \times \mathcal{T}^{-1}$$

### Step 4: Closed-Form Solution

For the special structure of randomized response, the inverse matrix $\mathcal{T}^{-1}$ has a closed form that leads to the formula above. The key insight is that the transformation probabilities follow a binomial pattern based on the number of edges flipped.

## Why This Formula Works

1. **Inclusion-Exclusion Principle**: The formula accounts for all possible ways a $(p,q)$-biclique can appear in the noisy graph, including cases where some edges are flipped.

2. **Unbiasedness**: The expectation of the estimator equals the true count because the transformation probabilities are correctly accounted for in the inverse matrix.

3. **Combinatorial Structure**: The formula leverages the fact that randomized response creates a linear relationship between true and noisy motif counts.

## Using Lean 4 to Verify the Derivation

### Why Lean 4?

1. **Formal Verification**: Prove mathematically (not just empirically) that the estimator is unbiased
2. **Combinatorial Correctness**: Verify that the inclusion-exclusion terms are exhaustive
3. **Automated Proof Search**: Use LLM-assisted proof generation with Lean as the verifier

### Key Components to Verify in Lean

#### 1. Transformation Probability Matrix

```lean
-- Define the transformation probability
def transformation_prob (p q : ℕ) (i j : ℕ) (μ : ℝ) : ℝ :=
  -- Probability that a motif with i edges becomes one with j edges
  -- This involves binomial coefficients and powers of μ/(1+μ) and 1/(1+μ)
```

#### 2. Linear System Relationship

```lean
-- Prove that expected noisy counts = true counts × transformation matrix
theorem expected_noisy_counts (G : BipartiteGraph) (ε : ℝ) :
  let μ := Real.exp ε
  let G' := apply_randomized_response G ε
  let T := transformation_matrix p q μ
  𝔼[count_motifs G' i] = ∑ j, T[i][j] * count_motifs G j
```

#### 3. Unbiasedness Theorem

```lean
-- The core theorem: the estimator is unbiased
theorem unbiased_biclique_estimator (G : BipartiteGraph) (p q : ℕ) (ε : ℝ) :
  let μ := Real.exp ε
  let G' := apply_randomized_response G ε
  let estimator := (∑ i : Fin (p * q + 1), (-1)^i * μ^i * count_motifs G' i) / (1 - μ)^(p * q)
  𝔼[estimator] = count_bicliques G p q
```

#### 4. Inverse Matrix Properties

```lean
-- Prove that T^{-1} exists and has the claimed structure
theorem inverse_transformation_exists (p q : ℕ) (μ : ℝ) (hμ : μ > 0) (hμ_ne_one : μ ≠ 1) :
  ∃ T_inv, T_inv * transformation_matrix p q μ = 1 ∧
  T_inv[i][pq] = (-1)^i * μ^i / (1 - μ)^(p * q)
```

### Implementation Strategy

#### Phase 1: Basic Definitions
- Define bipartite graphs in Lean
- Define randomized response mechanism
- Define motif counting functions

#### Phase 2: Transformation Probabilities
- Prove the transformation probability formula
- Verify the matrix structure

#### Phase 3: Linear System
- Prove the relationship between true and noisy counts
- Verify the matrix multiplication

#### Phase 4: Unbiasedness
- Prove the main unbiasedness theorem
- Verify the closed-form formula

#### Phase 5: Special Cases
- Verify the (2,q) formula specifically
- Prove correctness for small values (p=2, q=2,3,4)

### Example Lean Proof Sketch

```lean
import Mathlib.Probability.ProbabilityMassFunction
import Mathlib.Data.Matrix.Basic

-- Define a bipartite graph
structure BipartiteGraph where
  upper_vertices : Finset ℕ
  lower_vertices : Finset ℕ
  edges : Finset (ℕ × ℕ)

-- Randomized response
def randomized_response (edge : Bool) (ε : ℝ) : PMF Bool :=
  let p := 1 / (1 + Real.exp ε)
  if edge then PMF.ofFintype (fun b => if b then 1 - p else p)
  else PMF.ofFintype (fun b => if b then p else 1 - p)

-- Count motifs with i edges
def count_motifs (G : BipartiteGraph) (p q i : ℕ) : ℕ :=
  -- Count all (p,q) vertex sets with exactly i edges

-- Transformation probability
def transformation_prob (p q i j : ℕ) (μ : ℝ) : ℝ :=
  -- Binomial probability: choose (j-i) edges to flip (if j >= i)
  -- This is simplified; actual formula is more complex

-- Main theorem
theorem unbiased_estimator (G : BipartiteGraph) (p q : ℕ) (ε : ℝ) :
  let μ := Real.exp ε
  let G' := apply_rr G ε
  let est := (∑ i, (-1)^i * μ^i * count_motifs G' p q i) / (1 - μ)^(p * q)
  𝔼[est] = count_bicliques G p q := by
  -- Proof using transformation matrix inversion
  sorry
```

## Connection to Existing Code

The formula is implemented in several places:

1. **`biclique.cpp:3254-3258`**: One-round (2,q)-biclique estimation
   ```cpp
   long double res = 0, mu = std::exp(Eps);
   for (size_t i = 0; i <= 2 * K; ++i) {
       res += std::pow(-mu, i) * m__[i];
   }
   res /= std::pow(1 - mu, 2 * K);
   ```

2. **`biclique.cpp:2397-2500`**: `compute_local_res` function for (2,K) estimation using moment-based correction (different approach for multi-round)

3. **`overleaf-paper/3oneroundapproach.tex`**: Theorem statement and proof sketch

## Challenges for Lean Verification

1. **Combinatorial Complexity**: The transformation matrix for large (p,q) has $2^{pq}$ possible motif types
2. **Matrix Inversion**: Proving the closed-form inverse exists and has the claimed structure
3. **Expectation Calculations**: Handling probability distributions and expectations in Lean
4. **Binomial Identities**: Many combinatorial identities needed for the proof

## Next Steps

1. **Start with Small Cases**: Verify (2,2) and (2,3) first
2. **Build Transformation Matrix**: Prove the structure for small cases
3. **Generalize**: Use induction or pattern matching to extend to arbitrary (p,q)
4. **Integration**: Connect Lean proofs with the C++ implementation for certified correctness

## References

- Paper: `overleaf-paper/3oneroundapproach.tex` - Theorem for (p,q)-biclique
- Paper: `overleaf-paper/non-interactive.tex` - Butterfly (2,2)-biclique case
- Code: `biclique.cpp` - Implementation of the formula
- Documentation: `LEAN4_DP_MOTIF_VERIFICATION.md` - General framework for Lean verification
- **Guide: `LEAN_REPRODUCTION_GUIDE.md`** - Step-by-step guide for reproducing the formula in Lean
- **Lean File: `BicliqueTransformationMatrix.lean`** - Complete Lean formalization skeleton
