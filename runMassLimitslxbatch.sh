#!/bin/bash

#REMOTEDIR=/afs/cern.ch/work/d/digiovan/code/CMSSW_6_1_1/batchBaselinePP/
REMOTEDIR=/afs/cern.ch/user/j/jhugon/work/private/stats/CMSSW_6_1_1/stats/
echo "Checking for remotedir to exist..."
if ! `ssh lxplus5 touch $REMOTEDIR/touchfile.txt`; then
  echo "Error: Remote directory $REMOTEDIR doesn't exist or isn't writable"
  exit
fi

echo "deleting local files..."
rm -f shapes/*
rm -f statsCards/*
rm -f statsInput/*
rm -f statsOutput/*

echo "making cards..."
#for i in $(seq 1245 1265); do i=${i:0:3}"."${i:3:4};    # For Coupling Combination
for i in $(seq 115 1 155); do
  nice ./makeCards.py -m $i
done
echo "Removing files in lxplus5:$REMOTEDIR"
ssh lxplus5 "cd /tmp/; rm -rf $REMOTEDIR/*;echo \"Contents of dir: \`ls $REMOTEDIR \`\""
echo "Copying input files to lxplus5..."
scp statsCards/* lxplus5:$REMOTEDIR/.
echo "Running combine on lxplus5..."
ssh lxplus5 "cd $REMOTEDIR; eval \`scramv1 runtime -sh\`;bash run.sh; bash getStatus2.sh"
echo "Copying output files from lxplus5..."
scp lxplus5:$REMOTEDIR/*.txt.* statsInput/.

nice ./makeLimitPlots.py -m
nice ./makeSigMuPlots.py -m
#nice ./makeCCCPlots.py
echo "done."
