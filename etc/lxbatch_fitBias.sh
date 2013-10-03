#!/bin/bash
echo "Sourcing cmsset_default.sh"
cd /afs/cern.ch/cms/sw
source cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc462
echo "SCRAM_ARCH is $SCRAM_ARCH"
cd $LS_SUBCWD
# getting out of annoying etc/ directory!
cd ..
echo "In Directory: "
pwd
eval `scramv1 runtime -sh`
echo "cmsenv success!"
date

ijob=$1
echo "ijob = $ijob"
category=$2

command="./fitBiasStudy.py $ijob $category"
echo $command
$command >& logFitBiasJob$ijob$category

echo "done"
date

