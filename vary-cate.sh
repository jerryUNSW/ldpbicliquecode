#!/bin/bash

round=10

# Define the range of values 
start=1
end=3
step=0.5

datasets=(
	"unicode"
	"lrcwiki"
	"librec-filmtrust-ratings"
	"rmwiki"
  "csbwiki"
)
big=(
	"bag-kos"
	"bpywiki"
	"nips"
	"lastfm_band"
  "discogs_lstyle_lstyle"
	"digg-votes"
	"movielens-10m_rating"
)

## for small datasets 
# for dataset in "${datasets[@]}"
# do
#   echo "$dataset"
#   input="$start"
#   while (( $(bc <<< "$input <= $end") ))
#   do
#       echo "input = $input"

#       f1=JOC/"$dataset-vary-cate1.txt"
#       f2=JOC/"$dataset-vary-cate2.txt"

#       ./ldp-btf "$input" "../bidata/$dataset" "$round" 1 cate >> $f1
#       ./ldp-btf "$input" "../bidata/$dataset" "$round" 2 cate >> $f2

#       echo " "
#       input=$(bc <<< "$input + $step")
#   done
# done

## for big datasets 
for dataset in "${big[@]}"
do
  echo "$dataset"
  input="$start"
  while (( $(bc <<< "$input <= $end") ))
  do
      echo "input = $input"

      f1=JOC/"$dataset-vary-cate1.txt"
      f2=JOC/"$dataset-vary-cate2.txt"
      echo "" > "$f1"
      echo "" > "$f2"
      ./ldp-btf "$input" "../bidata/$dataset" "$round" 1 cate >> $f1
      ./ldp-btf "$input" "../bidata/$dataset" "$round" 2 cate >> $f2

      echo " "
      input=$(bc <<< "$input + $step")
  done
done