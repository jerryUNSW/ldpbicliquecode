#!/bin/bash

epsilon=1
round=10

datasets=(
	"unicode"
	"lrcwiki"
	"librec-filmtrust-ratings"
	"rmwiki"
	"opsahl-collaboration"
	"csbwiki"
	"bag-kos"
	"bpywiki"
	"nips"
	"lastfm_band"
    "discogs_lstyle_lstyle"
	"digg-votes"
	"movielens-10m_rating"
)

## compare their performance difference when epsilon = 1 
## small datasets 
for dataset in "${datasets[@]}"
do
	echo "$dataset"
	f1=extension/"$dataset-dbe-$epsilon.txt"
	f2=extension/"$dataset-multi-$epsilon.txt"
	f3=extension/"$dataset-wedge-$epsilon.txt"

    # state of the art is adv
	./ldp-btf $epsilon ../bidata/$dataset $round 1 btf  >> $f1
	./ldp-btf $epsilon ../bidata/$dataset $round 2 btf  >> $f2
	./ldp-btf $epsilon ../bidata/$dataset $round 3 btf  >> $f3

	echo " "
done