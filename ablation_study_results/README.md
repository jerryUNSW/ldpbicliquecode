# Ablation Study for MRCN/MRCN++ (R2.O3)

## Overview
This ablation study demonstrates the incremental contribution of each component in MRCN++.

## Components Being Evaluated

Your implementation has two key optimizations that can be toggled:

1. **Multi-estimator optimization** (`multi_estimator_switch`)
   - Aggregates estimates from each p-tuple
   - Provides variance reduction through multiple perspectives

2. **Two noisy graph construction** (`two_noisy_graph_switch`)
   - Constructs independent noisy graphs from upper and lower vertices
   - Averages results for improved accuracy

## Algorithm Variants

| Algorithm | Code | Multi-estimator | Two Noisy Graphs | Description |
|-----------|------|----------------|------------------|-------------|
| ADV       | alg2 | ❌ | ❌ | Base MRCN (baseline) |
| ADV+      | alg3 | ✅ | ❌ | Base + multi-estimator |
| ADV++     | alg4 | ✅ | ✅ | Base + both optimizations |

## Running the Ablation Study

### Step 1: Run Experiments
```bash
./run_ablation_study.sh
```

This will:
- Run all three algorithm variants (ADV, ADV+, ADV++)
- Test on multiple datasets (rmwiki, lrcwiki, csbwiki, co, unicode, librec-filmtrust-ratings)
- Use consistent parameters (ε=1.0, p=2, q=2, 50 rounds)
- Save results to `ablation_study_results/`

### Step 2: Visualize Results
```bash
python3 plot_ablation_study.py
```

This generates:
- `ablation_comparison.pdf` - Bar charts comparing relative error, runtime, and improvement
- `ablation_table.tex` - LaTeX table for the paper
- Console summary showing percentage improvements

## Expected Outcomes

The ablation study should show:

1. **ADV → ADV+**: Improvement from multi-estimator optimization
2. **ADV+ → ADV++**: Additional improvement from two noisy graph construction
3. **Trade-offs**: Potential increase in runtime for accuracy gains

## Output Files

```
ablation_study_results/
├── ablation_study_TIMESTAMP.log           # Execution log
├── ablation_summary.csv                   # Raw results (CSV)
├── DATASET_alg2_ablation.txt             # Detailed output for each run
├── DATASET_alg3_ablation.txt
├── DATASET_alg4_ablation.txt
├── ablation_comparison.pdf                # Visualization
├── ablation_comparison.png
└── ablation_table.tex                     # LaTeX table for paper
```

## Customization

You can modify parameters in `run_ablation_study.sh`:
- `EPSILON`: Privacy budget
- `NUM_ROUNDS`: Number of iterations for averaging
- `P`, `Q`: Biclique size parameters
- `ALPHA0`, `ALPHA1`, `ALPHA2`: Privacy budget allocation ratios

## For the Paper

Include:
1. The generated LaTeX table showing comparative results
2. The ablation comparison figure showing component contributions
3. Discussion explaining:
   - Why multi-estimator helps (variance reduction)
   - Why two noisy graphs help (independent sampling, averaging)
   - Trade-offs between accuracy and computational cost

## Addressing Reviewer Comments

**R2.O3**: "I don't see an ablation study of MRCN/MRCN++."

This ablation study directly addresses this by:
- Systematically evaluating each component's contribution
- Showing incremental improvements from base → ADV+ → ADV++
- Quantifying both accuracy gains and computational costs
- Providing clear visualization of component impacts

