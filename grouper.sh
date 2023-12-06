#!/bin/bash

META=$1
RDIR=$2

declare -a arr=("m1" "m2" "m3" "m4" "m5" "m6" "m7" "m8" "p1" "p2" "p3" "p4" "p5" "p6")

for i in "${arr[@]}"; do
  SAMPLES=$(awk -F'\t' '{print $1}' <(grep -w $i $META))
  mkdir -p $RDIR
  touch $RDIR/$i.txt
  for j in $SAMPLES; do
    echo $j >> $RDIR/$i.txt
  done
done
