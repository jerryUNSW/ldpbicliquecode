#!/bin/bash

epsilon=2.5
# issue:
# when epsilon is too large, naive becomes too good

round=10
# DBpath="/data/yizhangh/biclq_counts.db"  # Define the path to the SQLite database


# need to input $1 and $2 as p and q. 
datasets=(
	# "to"
	"co"
	# "big"
	# "crime"
    "M_PL_030"
	"unicode"
	"lrcwiki"
	"librec-filmtrust-ratings"
	"rmwiki"
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

	# naive algorithm: for naive running one round is enough
	./biclique $epsilon ../bidata/$dataset 1 0 $1 $2 | grep adv

    # # one-round algorithm, this is only feasible on To and Co.
    # ./biclique $epsilon ../bidata/$dataset $round 1 $1 $2  | grep adv
	
	# adv-base
	./biclique $epsilon ../bidata/$dataset $round 2 $1 $2 | grep adv

	# adv+
	./biclique $epsilon ../bidata/$dataset $round 3 $1 $2 | grep adv

	# adv++:
	./biclique $epsilon ../bidata/$dataset $round 4 $1 $2 | grep adv

	echo " "
	echo " "
done
