#!/bin/bash
echo "Sourcing cmsset_default.sh"
cd /afs/cern.ch/cms/sw
source cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc462
echo "SCRAM_ARCH is $SCRAM_ARCH"
cd $LS_SUBCWD
echo "In Directory: "
pwd
eval `scramv1 runtime -sh`
echo "cmsenv success!"
date

###########################################
###########################################
## Run with arguments: input_file_name.txt job_number number_of_toys
###########################################
###########################################


FILENAME=$1
echo "input filename: "$FILENAME
higgsmass=${FILENAME##*_}
higgsmass=${FILENAME%.txt}
echo "higgs mass: "$higgsmass
runname="_"${FILENAME%_*}
echo "runname: "$runname

ntoys=$3
echo "ntoys: "$ntoys
iseed=1000
((iseed += $2))
echo "iseed = $iseed"

command="nice combine -M GenerateOnly --saveToys -t $ntoys -s $iseed -m $higgsmass -n $runname $FILENAME"

echo $command
$command >& logGenOnly$2

echo "done"
date

