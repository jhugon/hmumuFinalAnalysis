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

## run like:
# bsub -q 1nd -o outMuToys0 lxbatch_MuBias.sh CombSplitAll_7P8TeV_125.0.txt 0 5
#
# or 
#
# for i in $(seq 1 5); do bsub -q 1nd -o outMuToys$i lxbatch_MuBias.sh CombSplitAll_7P8TeV_125.0.txt $i 5; done

FILENAME=$1
ijob=$2
ntoys=$3

TXTSUFFIX=".txt"
DIRNAME="Dir_"$FILENAME"_"$ijob"_Dir"
ROOTFILENAME=${FILENAME%$TXTSUFFIX}.root

mkdir $DIRNAME
cp $FILENAME $DIRNAME/
cp $ROOTFILENAME $DIRNAME/
cd $DIRNAME

iseed=1000
((iseed += $ijob))
echo "iseed = $iseed"

echo "input filename: "$FILENAME
echo "input root filename: "$ROOTFILENAME
echo "directory name: "$DIRNAME
higgsmass=$FILENAME
higgsmass=${higgsmass##*_}
higgsmass=${higgsmass%.txt}
higgsmass=${higgsmass%.*}
echo "higgs mass: "$higgsmass
runname="_"${FILENAME%_*}
echo "runname: "$runname

toysdir=/afs/cern.ch/user/j/jhugon/work/private/stats/CMSSW_6_1_1/finalAnalysisVoigtExp110160/statsCards/
toysfile=$toysdir/"higgsCombine$runname.GenerateOnly.mH$higgsmass.$iseed.root"
echo "toysfile: $toysfile"
linkcommand="ln -s -T $toysfile toysFile.root"
echo "running: $linkcommand"
$linkcommand

command="combine -M MaxLikelihoodFit -t $ntoys --rMin=-20 --rMax=20 --toysFile=toysFile.root -v 2 -s $iseed -m $higgsmass $FILENAME"
echo $command
$command >& $FILENAME."$ijob".muToys

echo "contents of dir:"
ls -lh

cp $FILENAME."$ijob".muToys ..
cp higgs*.root ../$FILENAME."$ijob".muToys.root

#cd ..
#rm -r $DIRNAME

echo "done"
date

