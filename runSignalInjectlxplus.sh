#!/bin/bash

SIGNALSTRENGTH=1.0

REMOTEDIR=/afs/cern.ch/user/j/jhugon/work/private/stats/CMSSW_6_1_1/stats/
if ! `ssh lxplus5 touch $REMOTEDIR/touchfile.txt`; then
  echo "Error: Remote directory $REMOTEDIR doesn't exist or isn't writable"
  exit
fi

rm -f shapes/*
rm -f statsCards/*
rm -f statsInput/*
rm -f statsOutput/*

nice ./makeCards.py --signalInject $SIGNALSTRENGTH #--toyData
echo "Removing files in lxplus5:$REMOTEDIR"
ssh lxplus5 "cd /tmp/jhugon/; rm -rf $REMOTEDIR/*;echo \"Contents of dir: \`ls $REMOTEDIR \`\""
echo "Copying input files to lxplus5..."
scp statsCards/* lxplus5:$REMOTEDIR/.
echo "Running combine on lxplus5..."
ssh lxplus5 "cd $REMOTEDIR; eval \`scramv1 runtime -sh\`;nice bash notlxbatch.sh;"
echo "Copying output files from lxplus5..."
scp lxplus5:$REMOTEDIR/*.out statsInput/.
scp lxplus5:$REMOTEDIR/*.sig statsInput/.
scp lxplus5:$REMOTEDIR/*.expsig statsInput/.
scp lxplus5:$REMOTEDIR/*.mu statsInput/.
scp lxplus5:$REMOTEDIR/*.txt.root statsInput/.
scp lxplus5:$REMOTEDIR/*.png statsInput/.

nice ./makeShapePlots.py --signalInject $SIGNALSTRENGTH --plotSignalStrength $SIGNALSTRENGTH
nice ./makeLimitPlots.py --signalInject $SIGNALSTRENGTH
nice ./makeSigMuPlots.py --signalInject $SIGNALSTRENGTH 
