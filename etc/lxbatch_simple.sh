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

TXTSUFFIX=".txt"
FILENAME=$1
DIRNAME="Dir"$1"Dir"
ROOTFILENAME=${1%$TXTSUFFIX}.root

mkdir $DIRNAME
cp $FILENAME $DIRNAME/
cp $ROOTFILENAME $DIRNAME/
cd $DIRNAME

echo "executing combine -M Asymptotic $FILENAME >& $FILENAME.out"

combine -M Asymptotic $FILENAME >& $FILENAME.out

cp $FILENAME.out ..

echo "done"
date

