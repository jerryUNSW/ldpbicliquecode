#!/bin/bash

epsilon=1
round=1
# DBpath="/data/yizhangh/biclq_counts.db"  # Define the path to the SQLite database


# need to input $1 and $2 as p and q. 

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

## compare their performance difference when epsilon = 1 
for dataset in "${datasets[@]}"
do
	echo "$dataset"

	# naive algorithm: 
	./biclique $epsilon ../bidata/$dataset $round 0 $1 $2 
	# | grep adv

    # # one-round algorithm
    # ./biclique $epsilon ../bidata/$dataset $round 1 $1 $2  | grep adv

	# # adv:
	# ./biclique $epsilon ../bidata/$dataset $round 3 $1 $2 | grep adv

	# # adv++:
	# ./biclique $epsilon ../bidata/$dataset $round 4 $1 $2 | grep adv

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
	echo " "
done
