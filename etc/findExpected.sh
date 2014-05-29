#!/bin/bash
echo "running findExpected.sh"
date
for i in *.txt; do
    [[ -e "$i" ]] || continue
FILENAME=$i

RMAX=150
RMIN=-150

echo "executing combine -M Asymptotic --rMax $RMAX  --run "expected" $FILENAME >& $FILENAME.explimit"

nice combine -M Asymptotic --rMax $RMAX --run "expected" $FILENAME >& $FILENAME.explimit
rm -f roostats*
rm -f higgsCombineTest*.root

echo "executing combine -M ProfileLikelihood --rMax $RMAX -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig"

nice combine -M ProfileLikelihood --rMax $RMAX -d $FILENAME --signif --expectSignal=1 -t -1 >& $FILENAME.expsig
rm -f roostats*
rm -f higgsCombineTest*.root

done

date
echo "done"

