#!/bin/bash

REMOTEDIR=/scratch/hpc/jhugon/stats/CMSSW_6_1_1/stats/
NODE=submit.hpc.ufl.edu
if ! `ssh $NODE touch $REMOTEDIR/touchfile.txt`; then
  echo "Error: Remote directory $REMOTEDIR doesn't exist or isn't writable"
  exit
fi

rm -f shapes/*
rm -f statsCards/*
rm -f statsInput/*
rm -f statsOutput/*

nice ./makeCards.py --bdtCut
echo "Removing files in $NODE:$REMOTEDIR"
ssh $NODE "rm -rf $REMOTEDIR/*;echo \"Contents of dir: \`ls $REMOTEDIR \`\""
echo "Copying input files to $NODE..."
scp statsCards/* $NODE:$REMOTEDIR/.
echo "Running combine on $NODE..."
ssh $NODE "cd $REMOTEDIR; bash runHPC.sh; bash getStatus2.sh "
echo "Copying output files from $NODE..."
rsync -az -e ssh  $NODE:$REMOTEDIR/*.out statsInput/.
rsync -az -e ssh  $NODE:$REMOTEDIR/*.sig statsInput/.
rsync -az -e ssh  $NODE:$REMOTEDIR/*.sigSM statsInput/.
rsync -az -e ssh $NODE:$REMOTEDIR/*.mu statsInput/.
rsync -az -e ssh $NODE:$REMOTEDIR/*.txt.root statsInput/.
rsync -az -e ssh $NODE:$REMOTEDIR/*.png statsInput/.
rsync -az -e ssh $NODE:$REMOTEDIR/*CCC*.root statsInput/.

#nice ./makeShapePlots.py
nice ./makeLimitPlots.py --bdtCut
#nice ./makeSigMuPlots.py --bdtCut
