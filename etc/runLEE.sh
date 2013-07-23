#!/bin/bash

ntoys=10
njobs=20
echo "Generating toys..."
for ijob in $(seq 1 $njobs); do
  iseed=1000
  (( iseed += $ijob ))
  command="nice combine -M GenerateOnly --saveToys -t $ntoys -s $iseed CombSplitAll_7P8TeV_135.0.txt"
  echo "Running: $command"
  date
  $command >& genonly.$ijob.log
  date
done
echo "==========================================\nStarting p-value jobs..."

chmod +x lxbatch_LEE.sh

#for i in CombSplitAll_7P8TeV_*.txt; do
for i in CombSplitAll_7P8TeV_135.0.txt CombSplitAll_7P8TeV_115.0.txt CombSplitAll_7P8TeV_145.0.txt; do
    [[ -e "$i" ]] || continue

  massVal=${#i}
  ((massVal -= 9))
  massVal=${i:$massVal:5}

  for ijob in $(seq 1 $njobs); do
    echo "Running on $i job number $ijob"
    command="bsub -q cmscaf1nd -o out_LEE_"$massVal"_$ijob lxbatch_LEE.sh $i $ijob"
    echo "Running: $command"
    $command
  done

done


