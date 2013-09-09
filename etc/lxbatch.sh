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

RMAX=150
RMIN=-150

echo "executing combine -M Asymptotic --rMax $RMAX $FILENAME >& $FILENAME.out"

combine -M Asymptotic --rMax $RMAX $FILENAME >& $FILENAME.out

echo "executing combine -M ProfileLikelihood --rMax $RMAX -d $FILENAME --signif >& $FILENAME.sig"

combine -M ProfileLikelihood --rMax $RMAX -d $FILENAME --signif >& $FILENAME.sig
rm -f roostats*
rm -f higgsCombineTest*.root

#echo "executing combine -M ProfileLikelihood --rMax $RMAX -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig"
#
#combine -M ProfileLikelihood --rMax $RMAX -d $FILENAME --signif --expectSignal=1 -t -1 >& $FILENAME.expsig
##combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig
#rm -f roostats*
#rm -f higgsCombineTest*.root

#echo "executing combine -M MaxLikelihoodFit --rMin $RMIN --rMax $RMAX --plots --saveNormalizations $FILENAME >& $FILENAME.mu"
#
#combine -M MaxLikelihoodFit --rMin $RMIN --rMax $RMAX --plots --saveNormalizations $FILENAME >& $FILENAME.mu
#rm -f roostats*
#rm -f higgsCombineTest*.root

#combine -M ChannelCompatibilityCheck --saveFitResult --rMin $RMIN --rMax $RMAX $FILENAME >> logCCC2
#mv higgsCombineTest.ChannelCompatibilityCheck.*.root ../$FILENAME.CCC2.root
#
#rm -f roostats*
#rm -f higgsCombineTest*.root

####New Coupligs and LH Scan Stuff
#
#mv $FILENAME hmm.txt
#
#text2workspace.py hmm.txt -m 125 -D data_obs -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingXSHiggs --PO modes=ggH,qqH --PO ggHRange=-20:20 --PO qqHRange=-40:40 -o wsQQvGG.root
#
#combine -M MultiDimFit --algo=singles --cl=0.68 wsQQvGG.root >& $FILENAME.profileQQvGG
#mv higgsCombineTest*.root $FILENAME.profileQQvGG.root
#
#combine -M MultiDimFit --algo=grid --points=5000 --fastScan wsQQvGG.root >& $FILENAME.lhGridQQvGG
#mv higgsCombineTest*.root $FILENAME.lhGridQQvGG.root
#
#combine -M MultiDimFit --algo=grid --points=1000 --fastScan --rMin=-30 --rMax=30 hmm.txt >& $FILENAME.lhGridR
#mv higgsCombineTest*.root $FILENAME.lhGridR.root
#
##### What I think is supposed to happen to get ggH and VBF only limits profiling other modes
##### But I'm not sure if it works at all!!
#
#mv $FILENAME hmm.txt
#text2workspace.py hmm.txt -m 125 -D data_obs -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingXSHiggs --PO modes=ggH --PO ggHRange=-50:50 -o wsGGOnly.root
#text2workspace.py hmm.txt -m 125 -D data_obs -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingXSHiggs --PO modes=qqH --PO qqHRange=-50:50 -o wsQQOnly.root
#combine -M Asymptotic wsGGOnly.root >& $FILENAME.out.GGOnly
#combine -M Asymptotic wsQQOnly.root >& $FILENAME.out.QQOnly

#########################################################3

cp $FILENAME.* ..
cp mlfit.root ../$FILENAME.root

echo "done"
date

