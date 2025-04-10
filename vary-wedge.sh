#!/bin/bash

round=10

# Define the range of values 
start=0.2
end=2.4
step=0.2

datasets=(
    # "M_PL_030"
	# "unicode"
	"lrcwiki"
	# "librec-filmtrust-ratings"
	# "rmwiki"
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

for dataset in "${datasets[@]}"
do
  echo "$dataset"
  input="$start"
  while (( $(bc <<< "$input <= $end") ))
  do
      echo "epsilon = $input"

      f1=extension/"$dataset-dbe-vary-$input.txt"
      f2=extension/"$dataset-multi-vary-$input.txt"
      f3=extension/"$dataset-wedge-vary-$input.txt"

      ./ldp-btf "$input" "../bidata/$dataset" "$round" 1 btf >> $f1
      ./ldp-btf "$input" "../bidata/$dataset" "$round" 2 btf >> $f2
      ./ldp-btf "$input" "../bidata/$dataset" "$round" 3 btf >> $f3

      echo " "
      input=$(bc <<< "$input + $step")
  done
done