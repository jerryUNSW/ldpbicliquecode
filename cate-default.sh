#!/bin/bash

epsilon=2
round=10

datasets=(
	"unicode"
	"lrcwiki"
	"librec-filmtrust-ratings"
	"rmwiki"
	"opsahl-collaboration"
	"csbwiki"
	"bag-kos"
)
big=(
	"bpywiki"
	"nips"
	"lastfm_band"
	"discogs_lstyle_lstyle"
	"digg-votes"
	"movielens-10m_rating"
    # "netflix"
)
# for netflix we apply sampling p = 0.1
## small datasets 
for dataset in "${datasets[@]}"
do
	echo "$dataset"
	f1=JOC/"$dataset"-default-1-new.txt
	f2=JOC/"$dataset"-default-2-new.txt
	touch "$f1"
	touch "$f2"
	./ldp-btf $epsilon ../bidata/$dataset $round 1 cate  >> $f1
	./ldp-btf $epsilon ../bidata/$dataset $round 2 cate  >> $f2
	echo " "
done

## the big datasets 
for dataset in "${big[@]}"
do
	echo "$dataset"
	f1=JOC/"$dataset"-default-1.txt
	f2=JOC/"$dataset"-default-2.txt
	touch "$f1"
	touch "$f2"
	./ldp-btf $epsilon ../bidata/$dataset $round 1 cate  >> $f1
	./ldp-btf $epsilon ../bidata/$dataset $round 2 cate  >> $f2
	echo " "
done