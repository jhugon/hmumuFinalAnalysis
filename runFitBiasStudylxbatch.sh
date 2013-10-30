#!/bin/bash

fullDir=`pwd`

# must go inside etc/ or bsub won't be able to find executable
cd etc/

###############################################
# To run every category in each job
#
#firstJobN=1
#lastJobN=100
#for iJob in $(seq $firstJobN $lastJobN); do
#  command="bsub -q 1nd -o $fullDir/jobOutFitBiasJob$iJob lxbatch_fitBias.sh $iJob"
#  echo "running: "$command
#  $command
#done
#
###############################################
# To run each cateogory as a seperate job

firstJobN=1
lastJobN=100
categories01Tight="Jets01PassPtG10BB Jets01PassPtG10BO Jets01PassPtG10BE Jets01PassPtG10OO Jets01PassPtG10OE Jets01PassPtG10EE"
categories01Loose="Jets01FailPtG10BB Jets01FailPtG10BO Jets01FailPtG10BE Jets01FailPtG10OO Jets01FailPtG10OE Jets01FailPtG10EE"
categories2="Jet2CutsVBFPass Jet2CutsGFPass Jet2CutsFailVBFGF"
categoriesAll=$categories01Tight" "$categories01Loose" "$categories2
categoriesImportant="Jets01PassPtG10BB Jets01PassPtG10BO Jets01PassPtG10BE Jet2CutsVBFPass Jet2CutsGFPass Jet2CutsFailVBFGF"
categoriesImportant01="Jets01PassPtG10BB Jets01PassPtG10BO Jets01PassPtG10BE"

for iJob in $(seq $firstJobN $lastJobN); do
  #for category in $categoriesAll; do
  #for category in $categories01Tight; do
  #for category in $categories2; do
  for category in $categoriesImportant; do
    command="bsub -q 1nd -o $fullDir/jobOutFitBiasJob$iJob$category lxbatch_fitBias.sh $iJob $category"
    echo "running: "$command
    $command
  done
done

###############################################

cd ..
