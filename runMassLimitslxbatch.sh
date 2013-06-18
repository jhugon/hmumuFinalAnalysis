#!/bin/bash

REMOTEDIR=/afs/cern.ch/work/d/digiovan/code/CMSSW_6_1_1/batchBaselinePP/
if ! `ssh lxplus5 touch $REMOTEDIR/touchfile.txt`; then
  echo "Error: Remote directory $REMOTEDIR doesn't exist or isn't writable"
  exit
fi

#rm -f shapes/*
#rm -f statsCards/*
rm -f statsInput/*
#rm -f statsOutput/*

#for i in "115" "120" "123" "124" "124.5" "125" "125.5" "126" "126.5" "127" "130" "135"; do
#for i in "115" "120" "125" "130" "135"; do
#for i in "120" "125" "130"; do
#  nice ./makeCards.py -m $i
#done
echo "Removing files in lxplus5:$REMOTEDIR"
ssh lxplus5 "cd /tmp/digiovan/; rm -rf $REMOTEDIR/*;echo \"Contents of dir: \`ls $REMOTEDIR \`\""
echo "Copying input files to lxplus5..."
scp statsCards/* lxplus5:$REMOTEDIR/.
echo "Running combine on lxplus5..."
ssh lxplus5 "cd $REMOTEDIR; eval \`scramv1 runtime -sh\`;bash run.sh; bash getStatus2.sh"
echo "Copying output files from lxplus5..."
scp lxplus5:$REMOTEDIR/*.out statsInput/.
scp lxplus5:$REMOTEDIR/*.sig statsInput/.
scp lxplus5:$REMOTEDIR/*.mu statsInput/.
scp lxplus5:$REMOTEDIR/*.txt.root statsInput/.
scp lxplus5:$REMOTEDIR/*.png statsInput/.
scp lxplus5:$REMOTEDIR/*CCC*.root statsInput/.

#nice ./makeLimitPlots.py -m
#nice ./makeSigMuPlots.py -m
#nice ./makeCCCPlots.py
