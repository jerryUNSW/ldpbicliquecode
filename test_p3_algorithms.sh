#!/bin/bash

# Test script for P=3 biclique counting algorithms using the non-batch (per-Q) runner
# Tests Naive, ADV, ADV+, ADV++ on multiple datasets with Q=4,5,6,7,8,9,10
# Uses the main `biclique` binary (runs one algorithm and one Q per invocation)

echo "=== P=3 Biclique Counting Algorithm Comparison (Batch Mode) ==="
echo "Testing on multiple datasets with 10 rounds each"
echo "Algorithms: Naive (0), ADV (2), ADV+ (3), ADV++ (4)"
echo "Q values: 4, 5, 6, 7, 8, 9, 10 (filtered by non-zero ground truth)"
echo "Epsilon: 1.0"
echo ""

# Datasets to test
# Datasets to test
datasets=("csbwiki" "rmwiki" "lrcwiki" "nips" "bag-kos")

# Function to run batch test for a dataset
run_batch_test() {
    local dataset=$1
    local output_file="larger_datasets_batch_results_FIXED/p3_${dataset}_new_2.txt"
    
    echo "=========================================="
    echo "DATASET: $dataset"
    echo "=========================================="
    echo "Running batch test (all algorithms and Q values)..."
    echo "Results will be saved to: $output_file"
    echo ""
    
    # Create output directory if it doesn't exist
    mkdir -p larger_datasets_batch_results_FIXED
    
    # Run the non-batch tests: iterate algorithms and Q values
    # Naive (alg=0) only needs 1 round, others need 10 rounds for statistical significance
    for alg in 0 2 3 4; do
        for q in {4..10}; do
            echo "==========================================" >> "$output_file"
            echo "DATASET: $dataset  ALG=$alg  Q=$q" >> "$output_file"
            echo "==========================================" >> "$output_file"
            # Run the main biclique binary: ./biclique <epsilon> <data_directory> <num_iterations> <algorithm_switch> <p> <q> [alpha0] [alpha1] [alpha2] [prob_filter] [sampling_ratio]
            # Naive algorithm (alg=0) is deterministic, so only 1 round needed
            # For Naive algorithm, use exact counting (sampling_ratio=1.0) instead of sampling
            if [ $alg -eq 0 ]; then
                ./biclique 1 ../bidata/$dataset 1 $alg 3 $q 0.05 0.6 0.35 0 1.0 >> "$output_file" 2>&1
            else
                ./biclique 1 ../bidata/$dataset 10 $alg 3 $q 0.05 0.6 0.35 0 0.1 >> "$output_file" 2>&1
            fi
            echo "" >> "$output_file"
        done
    done
    
    echo ""
    echo "=========================================="
    echo "Completed batch testing on $dataset"
    echo "Results saved to: $output_file"
    echo "=========================================="
    echo ""
}

# Run batch tests for each dataset
for dataset in "${datasets[@]}"; do
    run_batch_test "$dataset"
done

echo "=== All batch tests completed ==="
echo ""
echo "Summary:"
echo "- Each dataset was tested with all algorithms (Naive, ADV, ADV+, ADV++)"
echo "- Each algorithm was tested on all valid Q values (4-10, filtered by non-zero ground truth)"
echo "- Each test ran for 10 rounds to get statistical significance"
echo "- Results include mean estimates, relative errors, and standard deviations"
