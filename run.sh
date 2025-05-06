#!/bin/bash

epsilon=2
round=100
# DBpath="/data/yizhangh/biclq_counts.db"  # Define the path to the SQLite database

datasets=(
	# "to"
	# "co"
	# "big"
	# "crime"
    # "M_PL_030"
	# "unicode"
	"lrcwiki"
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

## compare their performance difference when epsilon = 1 
## small datasets 
for dataset in "${datasets[@]}"
do
	echo "$dataset"

    # one-round algorithm:
    # ./ldp-btf $epsilon ../bidata/$dataset $round 1 btf 

	# wedget-based algorithm 
	# ./ldp-btf $epsilon ../bidata/$dataset $round 3 btf 

	# weighted one-side sampling 
	# ./ldp-btf $epsilon ../bidata/$dataset $round 1 btf | grep adv
	
	# adv:
	./ldp-btf $epsilon ../bidata/$dataset $round 3 $1 $2 
	# > tmp.txt

    # estimate_lines=($(grep -oP 'estimate = \K[0-9\.]+' tmp.txt))
    # error_lines=($(grep -oP 'relative error = \K[0-9\.]+' tmp.txt))

	# p=$1
	# q=$2
	# algorithm_name='adv'
    # for ((run=0; run<$rounds; run++)); do
    #     estimate=${estimate_lines[$run]}
    #     relative_error=${error_lines[$run]}

    #     # Insert the data into the SQLite database
    #     sqlite3 $DBpath "INSERT INTO approximation_results (datasetname, p, q, algorithm_name, epsilon, relative_error, estimate) 
    #                      VALUES ('$dataset', $p, $q, '$algorithm_name', $epsilon, $relative_error, $estimate);"
    # done


	echo " "
done
