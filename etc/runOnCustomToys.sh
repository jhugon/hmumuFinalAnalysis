#!/bin/bash

ntoys=10

chmod +x lxbatch_customToys.sh

for i in *.txt; do
    [[ -e "$i" ]] || continue

  massVal=${#i}
  ((massVal -= 9))
  massVal=${i:$massVal:5}

  echo "Running on $i massVal: $massVal"
  command="bsub -q 1nd -o out_ToyJob_"$i" lxbatch_customToys.sh $i $ntoys"
  echo "Running: $command"
  $command
done


