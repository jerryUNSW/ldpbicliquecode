# Using Lean 4 to Reproduce the (p,q)-Biclique Estimation Formula

## Overview

This guide explains how to use Lean 4 to formally derive and verify the formula for estimating (p,q)-biclique counts in the non-interactive DP setting. The derivation follows these steps:

1. **Build noisy graph G'** using randomized response
2. **Consider all motifs** with p upper and q lower vertices
3. **Compute transformation probabilities** between motifs
4. **Build transformation matrix T** and derive its inverse
5. **Prove the unbiased estimator formula**

## Step-by-Step Derivation in Lean

### Step 1: Define Motifs and Edge Counts

In the non-interactive approach, we consider all possible motifs with p upper vertices and q lower vertices. A motif is defined by:
- A set of p upper vertices
- A set of q lower vertices  
- The edges present between them (can be 0 to p×q edges)

```lean
structure Motif (p q : ℕ) where
  upper_set : Finset ℕ
  lower_set : Finset ℕ
  upper_card : upper_set.card = p
  lower_card : lower_set.card = q

def motif_edge_count (G : BipartiteGraph) (M : Motif p q) : ℕ :=
  (M.upper_set ×ˢ M.lower_set).filter (fun (u, v) => (u, v) ∈ G.edges) |>.card
```

We group motifs by their edge count: **B_i** = set of motifs with exactly i edges (0 ≤ i ≤ p×q).

### Step 2: Transformation Probabilities

The key insight: when we apply randomized response to graph G to get G', each edge flips independently with probability p = 1/(1+e^ε).

**Transformation Probability**: Given a motif M with i edges in G, what is the probability it becomes a motif with j edges in G'?

This depends on:
- How many of the i present edges flip (become absent)
- How many of the (p×q - i) absent edges flip (become present)

**Mathematical Formula**:

For a motif with i edges in G:
- It has i edges present, (p×q - i) edges absent
- To get j edges in G': we need exactly j edges present after flipping

Let:
- k = number of present edges that flip (0 ≤ k ≤ i)
- ℓ = number of absent edges that flip (0 ≤ ℓ ≤ p×q - i)

Constraint: j = i - k + ℓ, so ℓ = j - i + k

The transformation probability is:

$$TP(\mathcal{B}_i, \mathcal{B}_j) = \sum_{k=0}^{i} \binom{i}{k} p^k (1-p)^{i-k} \binom{pq-i}{\ell} p^\ell (1-p)^{pq-i-\ell}$$

where ℓ = j - i + k (if this is valid, i.e., 0 ≤ ℓ ≤ pq - i), otherwise 0.

**In Lean**:

```lean
def transformation_prob (p q i j : ℕ) (ε : ℝ) : ℝ :=
  let total_edges := p * q
  let p_flip := 1 / (1 + Real.exp ε)  -- flipping probability
  let p_keep := 1 - p_flip            -- keeping probability
  
  ∑ k : Fin (i + 1),
    let ℓ := (j : ℤ) - (i : ℤ) + (k : ℤ)
    if 0 ≤ ℓ ∧ ℓ ≤ total_edges - i then
      let ℓ_nat := ℓ.natAbs
      (Nat.choose i k : ℝ) * (p_flip ^ k) * (p_keep ^ (i - k)) *
      (Nat.choose (total_edges - i) ℓ_nat : ℝ) * 
      (p_flip ^ ℓ_nat) * (p_keep ^ (total_edges - i - ℓ_nat))
    else 0
```

### Step 3: Build Transformation Matrix

The transformation matrix **T** is a (p×q+1) × (p×q+1) matrix where:
- T[i][j] = TP(B_i, B_j) = probability that a motif with i edges becomes one with j edges

```lean
def transformation_matrix (p q : ℕ) (ε : ℝ) : 
    Matrix (Fin (p * q + 1)) (Fin (p * q + 1)) ℝ :=
  fun i j => transformation_prob p q i.val j.val ε
```

**Key Property**: The expected relationship is:

$$\mathbb{E}[n'_{\mathcal{B}_j}] = \sum_{i=0}^{pq} T[i][j] \cdot n_{\mathcal{B}_i}$$

In matrix form:
$$\mathbb{E}[\mathbf{n}'] = \mathbf{n} \times T$$

where:
- **n** = vector of true motif counts [n_B₀, n_B₁, ..., n_B_{pq}]
- **n'** = vector of noisy motif counts [n'_{B₀}, n'_{B₁}, ..., n'_{B_{pq}}]

### Step 4: Derive the Inverse Matrix

To recover true counts from noisy counts, we need:
$$\hat{\mathbf{n}} = \mathbf{n}' \times T^{-1}$$

The key insight: for randomized response, the inverse matrix **T⁻¹** has a closed-form structure.

**Theorem**: The last column of T⁻¹ (which gives the estimate for n_{B_{pq}}, i.e., the biclique count) is:

$$T^{-1}[i][pq] = \frac{(-1)^i \mu^i}{(1-\mu)^{pq}}$$

where μ = e^ε.

**In Lean**:

```lean
def inverse_entry (p q i : ℕ) (ε : ℝ) : ℝ :=
  let μ_val := Real.exp ε
  let denom := (1 - μ_val) ^ (p * q)
  (-1 : ℝ) ^ i * μ_val ^ i / denom
```

### Step 5: Prove the Estimator Formula

Using the inverse matrix structure, the estimator for (p,q)-biclique count is:

$$\hat{n}_{\mathcal{B}_{pq}} = \sum_{i=0}^{pq} T^{-1}[i][pq] \cdot n'_{\mathcal{B}_i}$$

Substituting the formula for T⁻¹[i][pq]:

$$\hat{n}_{\mathcal{B}_{pq}} = \sum_{i=0}^{pq} \frac{(-1)^i \mu^i}{(1-\mu)^{pq}} \cdot n'_{\mathcal{B}_i}$$

$$= \frac{1}{(1-\mu)^{pq}} \sum_{i=0}^{pq} (-1)^i \mu^i \cdot n'_{\mathcal{B}_i}$$

**In Lean**:

```lean
def biclique_estimator (G' : BipartiteGraph) (p q : ℕ) (ε : ℝ) : ℝ :=
  let μ_val := Real.exp ε
  let denom := (1 - μ_val) ^ (p * q)
  let numerator := ∑ i : Fin (p * q + 1),
    (-1 : ℝ) ^ (i : ℕ) * μ_val ^ (i : ℕ) * (count_motifs G' p q i : ℝ)
  numerator / denom
```

### Step 6: Prove Unbiasedness

**Main Theorem**: The estimator is unbiased:

$$\mathbb{E}[\hat{n}_{\mathcal{B}_{pq}}] = n_{\mathcal{B}_{pq}} = C_{p,q}(G)$$

**Proof Strategy in Lean**:

1. Use the expected relationship: E[n'_j] = Σ_i T[i][j] · n_i
2. Substitute into the estimator formula
3. Use matrix multiplication: T · T⁻¹ = I
4. Show that the formula extracts the (pq)-th component correctly

```lean
theorem unbiased_biclique_estimator (G : BipartiteGraph) (p q : ℕ) (ε : ℝ)
    (hε_pos : ε > 0) (hμ_ne_one : μ ε ≠ 1) :
  let G' := apply_randomized_response G ε
  𝔼[biclique_estimator G' p q ε] = count_motifs G p q (p * q) := by
  -- Proof using transformation matrix inversion
  sorry
```

## Verification for Specific Cases

### Case 1: (2,2) - Butterfly Counting

For the (2,2) case, there are 6 distinct motif types (not just grouped by edge count). The transformation matrix is 6×6, and the formula reduces to:

$$f_{BTF}(G') = \frac{\mu^4 n'_{\mathcal{M}_6} - \mu^3 n'_{\mathcal{M}_5} + \mu^2(n'_{\mathcal{M}_4} + n'_{\mathcal{M}_3}) - \mu n'_{\mathcal{M}_2} + n'_{\mathcal{M}_1}}{(\mu-1)^4}$$

This matches the formula in `non-interactive.tex`.

### Case 2: (2,q) - General Case

For (2,q), the formula is:

$$\hat{n}_{\mathcal{B}_{2q}} = \frac{1}{(1-\mu)^{2q}} \sum_{i=0}^{2q} (-\mu)^i \cdot m_i$$

where m_i counts motifs with i edges between two upper vertices and q lower vertices.

This is implemented in `biclique.cpp:3254-3258`.

## Key Challenges in Lean Formalization

1. **Combinatorial Complexity**: For large (p,q), there are 2^(pq) possible edge patterns. We group by edge count to reduce to (pq+1) types.

2. **Transformation Probability Calculation**: The formula involves binomial coefficients and sums over all valid flipping patterns. Need to prove correctness.

3. **Matrix Inversion**: Proving that T⁻¹ exists and has the claimed closed-form structure requires:
   - Showing T is invertible (determinant ≠ 0)
   - Computing the inverse symbolically
   - Verifying the last column matches the formula

4. **Expectation Calculations**: Working with probability distributions and expectations in Lean requires the ProbabilityTheory library.

## Implementation Roadmap

### Phase 1: Small Cases (2-3 weeks)
- [ ] Define motifs and edge counting for (2,2) case
- [ ] Compute transformation probabilities manually
- [ ] Build 5×5 transformation matrix (for 0-4 edges)
- [ ] Verify inverse matrix structure
- [ ] Prove unbiasedness for (2,2)

### Phase 2: General (2,q) (3-4 weeks)
- [ ] Generalize to arbitrary q
- [ ] Prove transformation probability formula
- [ ] Verify closed-form inverse
- [ ] Complete unbiasedness proof

### Phase 3: General (p,q) (4-6 weeks)
- [ ] Extend to arbitrary p
- [ ] Handle combinatorial complexity
- [ ] Prove general formula
- [ ] Optimize for computational efficiency

### Phase 4: Integration (2-3 weeks)
- [ ] Connect Lean proofs with C++ implementation
- [ ] Generate certified code from proofs
- [ ] Document the complete derivation

## Example: Computing Transformation Matrix for (2,2)

For (2,2) with μ = e^ε, the transformation matrix T (when normalized by 1/(1+μ)⁴) is:

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

The inverse T⁻¹ has the property that the last column (for estimating n_{M₆} = butterfly count) is:

T⁻¹[·][6] = [1, -μ, μ², μ², -μ³, μ⁴] / (μ-1)⁴

This gives the butterfly formula.

## Benefits of Lean Formalization

1. **Mathematical Rigor**: Prove correctness, not just empirical validation
2. **Automated Verification**: Catch errors in combinatorial calculations
3. **Reproducibility**: Complete derivation that can be checked by anyone
4. **Generalization**: Once proven for (2,2), extend to (2,q) and (p,q) systematically
5. **Research Value**: Formal proofs are stronger than empirical results

## References

- Paper: `overleaf-paper/3oneroundapproach.tex` - General (p,q) formula
- Paper: `overleaf-paper/non-interactive.tex` - (2,2) butterfly case
- Code: `biclique.cpp` - Implementation
- Lean File: `BicliqueTransformationMatrix.lean` - Formalization skeleton
