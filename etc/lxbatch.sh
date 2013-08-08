#!/bin/bash
echo "Sourcing cmsset_default.sh"
cd /afs/cern.ch/cms
source cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc472
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

echo "create directory $DIRNAME"
mkdir -p $DIRNAME
# check the directory exist
if [ ! -d $DIRNAME ]; then
    echo "folder $DIRNAME does not exist"
    echo "exiting..."
    exit
fi;

cp $FILENAME $DIRNAME/
cp $ROOTFILENAME $DIRNAME/
cd $DIRNAME

echo "executing combine -M Asymptotic $FILENAME >& $FILENAME.out"

combine -M Asymptotic $FILENAME >& $FILENAME.out

echo "executing combine -M ProfileLikelihood -d $FILENAME --signif >& $FILENAME.sig"

combine -M ProfileLikelihood -d $FILENAME --signif >& $FILENAME.sig
rm -f roostats*
rm -f higgsCombineTest*.root

echo "executing combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig"

combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 >& $FILENAME.expsig
#combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig
rm -f roostats*
rm -f higgsCombineTest*.root

echo "executing combine -M MaxLikelihoodFit --plots --saveNormalizations $FILENAME >& $FILENAME.mu"

#combine -M MaxLikelihoodFit --rMin -50 --rMax 50 --plots --saveNormalizations $FILENAME >& $FILENAME.mu
combine -M MaxLikelihoodFit --rMin -50 --rMax 50 --saveNormalizations  --expectSignal=1 -t 100 $FILENAME >& $FILENAME.mu
mv higgsCombineTest*.root $FILENAME.mu.root
rm -f roostats*
rm -f higgsCombineTest*.root

# only positive mu... not really good choice. Commenting it out
#combine -M ChannelCompatibilityCheck --saveFitResult --rMax 150 $FILENAME >> logCCC
#mv higgsCombineTest.ChannelCompatibilityCheck.*.root ../$FILENAME.CCC.root

#rm -f roostats*
#rm -f higgsCombineTest*.root

combine -M ChannelCompatibilityCheck --rMin -50 --rMax 50 --saveFitResult $FILENAME >> logCCC2
mv higgsCombineTest.ChannelCompatibilityCheck.*.root ../$FILENAME.CCC2.root

rm -f roostats*
rm -f higgsCombineTest*.root
cp $FILENAME.* ..
cp mlfit.root ../$FILENAME.root
#for subname in *_fit_s.png; do
#  cp $subname ../${FILENAME%$TXTSUFFIX}_$subname
#done
#cp $FILENAME.expsig ..

echo "done"
date

