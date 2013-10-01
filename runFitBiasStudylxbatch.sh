#!/bin/bash

fullDir=`pwd`

# must go inside etc/ or bsub won't be able to find executable
cd etc/

firstJobN=1
lastJobN=100
for iJob in $(seq $firstJobN $lastJobN); do
  command="bsub -q cmscaf1nh -o $fullDir/jobOutFitBiasJob$iJob lxbatch_fitBias.sh $iJob"
  #command="bsub -q 1nd -o /dev/null lxbatch_fitBias.sh $iJob"
  echo "running: "$command
  $command
done

cd ..
