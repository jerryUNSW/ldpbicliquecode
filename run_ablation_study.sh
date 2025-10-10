#!/bin/bash

# Ablation Study Script for MRCN/MRCN++ Components
# This script compares ADV (alg2), ADV+ (alg3), and ADV++ (alg4)
# to show the incremental contribution of each component

DATASETS=("rmwiki" "lrcwiki" "csbwiki" "co" "unicode" "librec-filmtrust-ratings")
EPSILON=1.0
NUM_ROUNDS=50
P=2
Q=2

# Budget allocation (using default: 0.05, 0.6, 0.35)
ALPHA0=0.05
ALPHA1=0.6
ALPHA2=0.35

OUTPUT_DIR="ablation_study_results"
mkdir -p "$OUTPUT_DIR"

LOG_FILE="$OUTPUT_DIR/ablation_study_$(date +%Y%m%d_%H%M%S).log"

echo "Starting Ablation Study: $(date)" | tee -a "$LOG_FILE"
echo "Comparing ADV (alg2), ADV+ (alg3), ADV++ (alg4)" | tee -a "$LOG_FILE"
echo "Parameters: epsilon=$EPSILON, rounds=$NUM_ROUNDS, p=$P, q=$Q" | tee -a "$LOG_FILE"
echo "Budget allocation: alpha0=$ALPHA0, alpha1=$ALPHA1, alpha2=$ALPHA2" | tee -a "$LOG_FILE"
echo "============================================" | tee -a "$LOG_FILE"

for dataset in "${DATASETS[@]}"; do
    echo "" | tee -a "$LOG_FILE"
    echo "Processing dataset: $dataset" | tee -a "$LOG_FILE"
    echo "--------------------------------------------" | tee -a "$LOG_FILE"
    
    for alg in 2 3 4; do
        ALG_NAME=""
        if [ $alg -eq 2 ]; then
            ALG_NAME="ADV (base MRCN)"
        elif [ $alg -eq 3 ]; then
            ALG_NAME="ADV+ (+ multi-estimator)"
        elif [ $alg -eq 4 ]; then
            ALG_NAME="ADV++ (+ two noisy graphs)"
        fi
        
        echo "  Running $ALG_NAME..." | tee -a "$LOG_FILE"
        
        OUTPUT_FILE="$OUTPUT_DIR/${dataset}_alg${alg}_ablation.txt"
        
        ./main $EPSILON $dataset $NUM_ROUNDS $alg $P $Q $ALPHA0 $ALPHA1 $ALPHA2 > "$OUTPUT_FILE" 2>&1
        
        # Extract key metrics
        MEAN=$(grep "Mean = " "$OUTPUT_FILE" | awk '{print $NF}')
        REL_ERR=$(grep "adv rel err = " "$OUTPUT_FILE" | awk '{print $NF}')
        TIME=$(grep "time:" "$OUTPUT_FILE" | awk -F: '{print $2}')
        
        echo "    Mean estimate: $MEAN" | tee -a "$LOG_FILE"
        echo "    Rel error: $REL_ERR" | tee -a "$LOG_FILE"
        echo "    Time: $TIME seconds" | tee -a "$LOG_FILE"
        
        # Save summary
        echo "$dataset,alg$alg,$ALG_NAME,$MEAN,$REL_ERR,$TIME" >> "$OUTPUT_DIR/ablation_summary.csv"
    done
done

echo "" | tee -a "$LOG_FILE"
echo "============================================" | tee -a "$LOG_FILE"
echo "Ablation study completed: $(date)" | tee -a "$LOG_FILE"
echo "Results saved in $OUTPUT_DIR/" | tee -a "$LOG_FILE"
echo "Summary CSV: $OUTPUT_DIR/ablation_summary.csv" | tee -a "$LOG_FILE"

