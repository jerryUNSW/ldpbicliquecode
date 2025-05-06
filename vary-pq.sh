#!/bin/bash

# Define experiment parameters
round=10
# DBpath="/data/yizhangh/biclq_counts.db"  # Define the path to your SQLite database (uncomment if needed)

# Define datasets (you can uncomment additional datasets as needed)
datasets=(
    # "to"
    # "co"
    # "big"
    # "crime"
    # "M_PL_030"
    "unicode"
    # "lrcwiki"
    # "librec-filmtrust-ratings"
    # "rmwiki"
    # "amazon-ratings"
    # "edit-iowiktionary"
    # "opsahl-collaboration"
    # "csbwiki"
    # "bag-kos"
    # "bpywiki"
    # "nips"
    # "lastfm_band"
    # "discogs_lstyle_lstyle"
    # "digg-votes"
    # "movielens-10m_rating"
)

# Define P, Q pairs as an array of space-separated strings
# p_q_pairs=("2 3" "2 4" "3 2" "3 3" "3 4" "4 2" "4 3" "4 4")
p_q_pairs=("2 2" )

# Loop over each dataset
for dataset in "${datasets[@]}"
do
    echo "Processing dataset: $dataset"

    # Loop over each P, Q pair
    for pq in "${p_q_pairs[@]}"
    do
        # Split the P, Q pair into separate variables
        read p q <<< "$pq"
        echo "  Running for P=$p, Q=$q"

        # Loop over epsilon values from 0.5 to 3.0 in increments of 0.5 (innermost loop)
        for epsilon in $(seq 0.5 0.5 3.0)
        do
            echo "    Epsilon: $epsilon"

            # Run the third algorithm (mode 3) with the current epsilon, dataset, and P, Q pair
            ./ldp-btf $epsilon ../bidata/$dataset $round 3 $p $q | tail | grep adv
        done
        echo " "
    done

    echo " "
done

# Optional: Clean up temporary file (uncomment if needed)
# rm tmp.txt