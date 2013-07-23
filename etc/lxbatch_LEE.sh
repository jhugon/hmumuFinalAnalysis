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
DIRNAME="Dir"$1_"$2"_"Dir"
ROOTFILENAME=${1%$TXTSUFFIX}.root

mkdir $DIRNAME
cp $FILENAME $DIRNAME/
cp $ROOTFILENAME $DIRNAME/
cd $DIRNAME

ntoys=10
iseed=1000
((iseed += $2))
echo "iseed = $iseed"

toysfile="higgsCombineTest.GenerateOnly.mH120.$iseed.root"
echo "toysfile: $toysfile"
linkcommand="ln -s -T ../$toysfile $toysfile"
echo "running: $linkcommand"
$linkcommand

command="combine -M ProfileLikelihood --signif -t $ntoys --toysFile=$toysfile -v 2 $FILENAME"
echo $command
$command >& $FILENAME."$2".sigToys

echo "contents of dir:"
ls -lh

cp $FILENAME."$2".sigToys ..
cp higgs*ProfileLikelihood*.root ../$FILENAME."$2".sigToys.root

#cd ..
#rm -r $DIRNAME

echo "done"
date

