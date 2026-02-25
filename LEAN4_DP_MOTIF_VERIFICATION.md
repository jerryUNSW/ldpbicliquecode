# Automated Unbiased Estimator Derivation and Verification for DP-Protected Subgraph Counting

## Project Objective

To automate the derivation and formal verification of unbiased estimator formulas for subgraph (motif) counting in graphs protected by Differential Privacy (DP). This system bypasses human error in complex combinatorial derivations (especially for motifs with $\ge 5$ nodes) by using Lean 4 as a ground-truth verifier for LLM-generated conjectures.

## Core Architectural Components

### A. The Motif Specification Engine (Input)

**Representation:**
- Represent motifs using **Adjacency Matrices** or **Graph Isomorphism types** in Lean 4.
- Goal: Provide the LLM with the precise structure and its **Automorphism Group size** ($|Aut(G')|$), which is critical for correcting counting bias.

**Context Injection:**
- Feed the LLM existing papers/manuscripts as "Few-Shot" examples, showing how a 3-node triangle or 4-node clique was derived.

### B. The Inclusion-Exclusion Logic (Combinatorial Layer)

**Problem:**
- Counting motifs in DP involves "Edge-Induced" vs. "Node-Induced" subgraphs.
- When edges are added/removed (DP sensitivity), the change in one motif affects others.

**Lean's Role:**
- Use Lean's `Finset` (finite sets) and `Fintype` libraries to formally define the overlapping structures.
- Automation: The LLM generates the Inclusion-Exclusion terms; Lean verifies if the mapping between the "Private Observed Count" and the "True Count" is exhaustive and covers all isomorphic overlaps.

### C. The DP Noise & Unbiasedness Layer (Statistical Layer)

**Noise Modeling:**
- Define the DP mechanism (e.g., Laplace or Geometric) as a probability distribution in Lean.

**Unbiasedness Theorem:**
- Define the core theorem for the LLM to prove:
  $$\text{theorem unbiased\_motif\_est (m : Motif) : } \mathbb{E}[\hat{f}(\text{Data}_{private})] = \text{Count}_{true}(m)$$

**Expectation Arithmetic:**
- Use Lean's `ProbabilityTheory` library to prove that the added noise cancels out in the expectation of your estimator formula.

### D. The Iterative Refinement Loop (LLM + Lean)

1. **Conjecture:** The LLM (via Cursor) proposes a Lean 4 proof script for a specific motif (e.g., 5-clique).
2. **Verify:** You run the Lean compiler.
3. **Feedback:** If it fails (e.g., a "type mismatch" or "failed to unify"), copy the exact compiler error back into Cursor.
4. **Fix:** The LLM interprets the logical error (not just a syntax error) to adjust coefficients or combinatorial terms.

## Detailed Workflow for Cursor

When you start in Cursor, you should provide a System Prompt or a `.cursorrules` file with these instructions:

### Step 1: Define the Environment

```
"We are working in Lean 4. Use Mathlib4. We need to define a graph G and a motif M. 
We are applying a Differential Privacy mechanism where each edge presence is a Bernoulli 
variable or noise is added to the final count."
```

### Step 2: Formalize Inclusion-Exclusion

```
"Write a Lean 4 function that calculates the number of ways two k-node motifs can overlap. 
Use this to construct a system of linear equations that relates the observed noisy counts 
to the true counts."
```

### Step 3: Prove Unbiasedness

```
"For a 5-node motif M, given the sensitivity Δ, prove that the estimator 
M̂ = (M_obs - Bias) / Scale satisfies the property that its expected value equals 
the ground truth count."
```

## Why This is "High Value" (The Research Edge)

1. **Beyond Simulation:** Most DP papers use Python simulations to "show" an estimator is unbiased. You are proving it mathematically for any arbitrary motif.

2. **Machine-Certified Privacy:** By using the SampCert (or similar) framework in Lean, you ensure that your privacy budget $\epsilon$ is not just a claim, but a verified property of the code.

3. **Handling Complexity:** Manually deriving the 5-node motif constants is a nightmare. Your system treats it as a Search Problem—the LLM searches for the proof, Lean confirms the find.

## Next Steps to Feed into Cursor

### Prompt 1: Basic Setup
```
"Please set up a Lean 4 file that defines a simple undirected graph and a function 
to count triangles (3-cliques)."
```

### Prompt 2: Add Noise
```
"Now, add a Laplace noise term to this count. Define the expectation of this noisy 
count using Lean's ProbabilityTheory library."
```

### Prompt 3: Prove Unbiasedness
```
"Identify the bias introduced by the noise and propose an adjusted formula to make 
the estimator unbiased. Try to prove E[Estimator] = TrueCount."
```

## Technical Prerequisites

- **Lean 4** with **Mathlib4**
- **ProbabilityTheory** library in Lean
- **Finset** and **Fintype** libraries for combinatorial structures
- Understanding of:
  - Differential Privacy mechanisms (Laplace, Geometric)
  - Graph theory (motifs, cliques, automorphism groups)
  - Inclusion-Exclusion principle
  - Unbiased estimation theory

## Key Mathematical Components

### Sensitivity Calculation
For a motif $M$ in a graph $G$, the sensitivity $\Delta$ is the maximum change in motif count when a single edge is added or removed.

### Unbiased Estimator Form
For a noisy count $C_{obs} = C_{true} + \eta$ where $\eta \sim \text{Laplace}(0, \Delta/\epsilon)$, the unbiased estimator is:
$$\hat{C} = C_{obs} - \mathbb{E}[\eta] = C_{obs}$$

However, for more complex cases with overlapping motifs, the inclusion-exclusion principle must be applied to correct for double-counting.

### Inclusion-Exclusion for Overlapping Motifs
When motifs overlap, we need to account for:
- Direct counts of each motif type
- Overlaps between motif pairs
- Higher-order overlaps (triples, quadruples, etc.)

The general form:
$$C_{true}(M) = \sum_{S \subseteq \text{Overlaps}} (-1)^{|S|+1} \cdot C_{obs}(S)$$

## Implementation Roadmap

1. **Phase 1:** Set up Lean 4 environment with graph and motif definitions
2. **Phase 2:** Implement basic triangle counting with DP noise
3. **Phase 3:** Prove unbiasedness for 3-node motifs
4. **Phase 4:** Extend to 4-node motifs (cliques, paths, stars)
5. **Phase 5:** Automate 5+ node motif derivations using LLM-assisted proof search
6. **Phase 6:** Integrate with existing DP graph counting codebase

## References and Context

- Existing papers/manuscripts in `overleaf-paper/` directory
- Current implementation in `exactcounting/` directory
- Experimental results in various `jester-p-*` directories

## Notes

- This approach combines automated theorem proving with LLM-assisted conjecture generation
- The iterative refinement loop is key to handling complex combinatorial structures
- Formal verification provides mathematical guarantees beyond empirical validation
