#!/bin/bash

REMOTEDIR=/afs/cern.ch/user/j/jhugon/work/private/stats/CMSSW_5_2_5/stats/
if ! `ssh lxplus touch $REMOTEDIR/touchfile.txt`; then
  echo "Error: Remote directory $REMOTEDIR doesn't exist or isn't writable"
  exit
fi

rm -f shapes/*
rm -f statsCards/*
rm -f statsInput/*
rm -f statsOutput/*

#for i in "115" "120" "123" "124" "124.5" "125" "125.5" "126" "126.5" "127" "130" "135"; do
#for i in "115" "120" "125" "130" "135"; do
for i in "120" "125" "130"; do
  nice ./makeCards.py -m $i
done
echo "Removing files in lxplus:$REMOTEDIR"
ssh lxplus "cd /tmp/jhugon/; rm -rf $REMOTEDIR/*;echo \"Contents of dir: \`ls $REMOTEDIR \`\""
echo "Copying input files to lxplus..."
scp statsCards/* lxplus:/afs/cern.ch/user/j/jhugon/work/private/stats/CMSSW_5_2_5/stats/.
echo "Running combine on lxplus..."
ssh lxplus "cd $REMOTEDIR; eval \`scramv1 runtime -sh\`;bash run.sh; bash getStatus2.sh"
echo "Copying output files from lxplus..."
scp lxplus:$REMOTEDIR/*.out statsInput/.
scp lxplus:$REMOTEDIR/*.sig statsInput/.
scp lxplus:$REMOTEDIR/*.mu statsInput/.
scp lxplus:$REMOTEDIR/*.txt.root statsInput/.
scp lxplus:$REMOTEDIR/*.png statsInput/.
scp lxplus:$REMOTEDIR/*CCC*.root statsInput/.

nice ./makeLimitPlots.py -m
nice ./makeSigMuPlots.py -m
nice ./makeCCCPlots.py
