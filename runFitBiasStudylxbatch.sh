#!/bin/bash

firstJobN=1
lastJobN=50
for iJob in $(seq $firstJobN $lastJobN); do
  command="bsub -q 1nd -o outFitBiasJob$iJob etc/lxbatch_fitBias.sh $iJob"
  echo "running: "$command
  #$command
done
