# Larger PQ Experiment

This experiment tests various (p,q) combinations with different algorithms and epsilon values to evaluate the performance of biclique counting algorithms.

## Configuration

- **Epsilon values**: 1, 2, 3
- **P values**: 2, 3, 4  
- **Q values**: 4, 5, 6, 7, 8
- **Algorithms**: 
  - 0: Naive
  - 1: ADV
  - 2: ADV+
  - 3: ADV++
- **Datasets**: co, unicode
- **Budget allocation**: alpha0=0.05, alpha1=0.6, alpha2=0.35
- **Probability filtering**: Disabled (0)

## Total Experiments

- **Per dataset**: 3 epsilons × 3 P values × 5 Q values × 4 algorithms = **180 experiments**
- **Total**: 2 datasets × 180 = **360 experiments**

## Files

### Scripts
- `run_larger_pq_experiment.sh` - Main experiment script
- `test_larger_pq_setup.sh` - Quick test to verify setup
- `analyze_larger_pq_experiment.py` - Results analysis and plotting

### Output Directories
- `larger_pq_experiment_logs/` - Individual experiment logs
- `larger_pq_experiment_results/` - Analysis results and plots

## Usage

### 1. Test Setup
```bash
./test_larger_pq_setup.sh
```

### 2. Run Full Experiment
```bash
./run_larger_pq_experiment.sh
```

### 3. Analyze Results
```bash
python3 analyze_larger_pq_experiment.py
```

## Expected Runtime

- **Per experiment**: ~1-30 seconds (depending on dataset size and parameters)
- **Total estimated time**: 1-3 hours for all 360 experiments
- **Timeout**: 30 minutes per experiment (will skip if exceeded)

## Output Analysis

The analysis script generates:

1. **Relative Error Plots**:
   - By algorithm and epsilon (averaged over p,q)
   - By (p,q) combinations
   - Heatmaps for each dataset

2. **Summary Tables**:
   - `summary_by_algorithm.csv` - Overall performance by algorithm
   - `detailed_results.csv` - All individual results

3. **Visualizations**:
   - PNG and PDF formats
   - Log-scale y-axis for relative errors
   - Color-coded heatmaps

## Monitoring Progress

During execution, you can monitor progress:

```bash
# Count completed experiments
ls larger_pq_experiment_logs/*.log | wc -l

# Check for errors
grep -l "Failed\|Error" larger_pq_experiment_logs/*.log

# View recent results
tail -f larger_pq_experiment_logs/co_eps1_p2_q4_alg2.log
```

## Expected Results

Based on previous experiments:
- **Naive**: Baseline performance, may struggle with larger p,q
- **ADV**: Basic differential privacy approach
- **ADV+**: Improved with multi-source estimator
- **ADV++**: Best performance with optimized budget allocation

The experiment will help identify:
1. Which algorithms perform best for different (p,q) combinations
2. How epsilon affects accuracy across algorithms
3. Optimal parameter ranges for each dataset
4. Computational efficiency trade-offs

## Troubleshooting

### Common Issues
1. **Dataset not found**: Ensure `../bidata/co.e` and `../bidata/unicode.e` exist
2. **Compilation errors**: Run `make clean && make` to rebuild
3. **Timeout issues**: Some large (p,q) combinations may be computationally expensive
4. **Memory issues**: Monitor system resources during execution

### Log Analysis
```bash
# Find failed experiments
grep -l "exit code" larger_pq_experiment_logs/*.log

# Check timeout issues  
grep -l "Timeout" larger_pq_experiment_logs/*.log

# View relative errors
grep "relative error" larger_pq_experiment_logs/*.log | head -10
```
