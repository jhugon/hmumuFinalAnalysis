#!/bin/bash

REMOTEDIR=/scratch/hpc/jhugon/stats/CMSSW_5_2_5/stats/
NODE=submit.hpc.ufl.edu
if ! `ssh $NODE touch $REMOTEDIR/touchfile.txt`; then
  echo "Error: Remote directory $NODE:$REMOTEDIR doesn't exist or isn't writable"
  exit
fi

rm -f shapes/*
rm -f statsCards/*
rm -f statsInput/*
rm -f statsOutput/*

for i in "115" "120" "123" "124" "124.5" "125" "125.5" "126" "126.5" "127" "130" "135"; do
#for i in "115" "120" "125" "130" "135"; do
#for i in "120" "125" "130"; do
#for i in "125"; do
  nice ./makeCards.py -m $i
done
echo "Removing files in :$REMOTEDIR"
ssh $NODE "rm -rf $REMOTEDIR/*;echo \"Contents of dir: \`ls $REMOTEDIR \`\""
echo "Copying input files to $NODE..."
rsync -az -e ssh statsCards/* $NODE:$REMOTEDIR/.
echo "Running combine on $NODE..."
ssh $NODE "cd $REMOTEDIR; eval \`scramv1 runtime -sh\`;bash runHPC.sh; bash getStatus2.sh"
echo "Copying output files from $NODE..."
rsync -az -e ssh  $NODE:$REMOTEDIR/*.out statsInput/.
rsync -az -e ssh  $NODE:$REMOTEDIR/*.sig statsInput/.
rsync -az -e ssh  $NODE:$REMOTEDIR/*.sigSM statsInput/.
rsync -az -e ssh $NODE:$REMOTEDIR/*.mu statsInput/.
rsync -az -e ssh $NODE:$REMOTEDIR/*.txt.root statsInput/.
rsync -az -e ssh $NODE:$REMOTEDIR/*.png statsInput/.
rsync -az -e ssh $NODE:$REMOTEDIR/*CCC*.root statsInput/.

nice ./makeLimitPlots.py -m
nice ./makeSigMuPlots.py -m
nice ./makeCCCPlots.py
