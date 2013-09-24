#!/bin/bash

fullDir=`pwd`

# must go inside etc/ or bsub won't be able to find executable
cd etc/

firstJobN=1
lastJobN=10
for iJob in $(seq $firstJobN $lastJobN); do
  #command="bsub -q 1nh -o $fullDir/outFitBiasJob$iJob lxbatch_fitBias.sh $iJob"
  command="bsub -q 1nd -o /dev/null lxbatch_fitBias.sh $iJob"
  echo "running: "$command
  $command
done

cd ..
