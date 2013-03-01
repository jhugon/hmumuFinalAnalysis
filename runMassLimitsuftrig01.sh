#!/bin/bash

REMOTEDIR=/data/uftrig01b/jhugon/hmumu/stats/CMSSW_5_2_5/stats
if ! `ssh uftrig01 touch $REMOTEDIR/touchfile.txt`; then
  echo "Error: Remote directory $REMOTEDIR doesn't exist or isn't writable"
  exit
fi

rm -f shapes/*
rm -f statsCards/*
rm -f statsInput/*
rm -f statsOutput/*

for i in "115" "120" "123" "124" "124.5" "125" "125.5" "126" "126.5" "127" "130" "135"; do
#for i in "115" "120" "125" "130" "135"; do
  nice ./makeCards.py -m $i
done
echo "Removing files in uftrig01:$REMOTEDIR"
ssh uftrig01 "rm -rf $REMOTEDIR/*;echo \"Contents of dir: \`ls $REMOTEDIR \`\""
echo "Copying input files to uftrig01..."
scp statsCards/* uftrig01:$REMOTEDIR
echo "Running combine on uftrig01..."
ssh uftrig01 "cd $REMOTEDIR; echo \"In dir: \`pwd\`\"; eval \`scramv1 runtime -sh\`;nice bash notlxbatch.sh;"
echo "Copying output files from uftrig01..."
scp uftrig01:$REMOTEDIR/*.out statsInput/.
scp uftrig01:$REMOTEDIR/*.sig statsInput/.
scp uftrig01:$REMOTEDIR/*.mu statsInput/.
scp uftrig01:$REMOTEDIR/*.txt.root statsInput/.
scp uftrig01:$REMOTEDIR/*.png statsInput/.

nice ./makeLimitPlots.py -m
nice ./makeSigMuPlots.py -m
