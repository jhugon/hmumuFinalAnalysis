#!/bin/bash

ntoys=5
ijobStart=6
njobs=50

#chmod +x lxbatch_LEE.sh
#chmod +x lxbatch_LEE_prep.sh



#echo "Generating toys..."
#for ijob in $(seq $ijobStart $njobs); do
#  command="bsub -q cmscaf1nh -o outGenOnly$ijob lxbatch_LEE_prep.sh CombSplitAll_7P8TeV_135.0.txt $ijob $ntoys"
#  echo "Running: $command"
#  $command
#done





echo "==========================================\nStarting p-value jobs..."
for i in CombSplitAll_7P8TeV_*.txt; do
    [[ -e "$i" ]] || continue

  massVal=${#i}
  ((massVal -= 9))
  massVal=${i:$massVal:5}

  for ijob in $(seq $ijobStart $njobs); do
    echo "Running on $i job number $ijob"
    command="bsub -q cmscaf1nd -o out_LEE_"$massVal"_$ijob -M 10000000 lxbatch_LEE.sh $i $ijob $ntoys"
    echo "Running: $command"
    $command
  done

done


