#!/bin/bash
echo "running notlxbatch.sh"
date
for i in *.txt; do
    [[ -e "$i" ]] || continue
FILENAME=$i

RMAX=150
RMIN=-150

echo "executing combine -M Asymptotic --rMax $RMAX $FILENAME >& $FILENAME.out"

nice combine -M Asymptotic --rMax $RMAX $FILENAME >& $FILENAME.out
rm -f roostats*
rm -f higgsCombineTest*.root

echo "executing combine -M ProfileLikelihood --rMax $RMAX -d $FILENAME --signif >& $FILENAME.sig"

nice combine -M ProfileLikelihood --rMax $RMAX -d $FILENAME --signif >& $FILENAME.sig
rm -f roostats*
rm -f higgsCombineTest*.root

#echo "executing combine -M ProfileLikelihood --rMax $RMAX -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig"
#
#nice combine -M ProfileLikelihood --rMax $RMAX -d $FILENAME --signif --expectSignal=1 -t -1 >& $FILENAME.expsig
##combine -M ProfileLikelihood --rMax $RMAX -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig
#rm -f roostats*
#rm -f higgsCombineTest*.root

#echo "executing combine -M MaxLikelihoodFit --rMin $RMIN --rMax $RMAX --plots --saveNormalizations $FILENAME >& $FILENAME.mu"
#
#nice combine -M MaxLikelihoodFit --rMin $RMIN --rMax $RMAX --plots --saveNormalizations $FILENAME >& $FILENAME.mu
#rm -f roostats*
#rm -f higgsCombineTest*.root
#cp mlfit.root $FILENAME.root
#for subname in *_fit_s.png; do
#  cp $subname ${FILENAME%$TXTSUFFIX}_$subname
#done
#
#nice combine -M ChannelCompatibilityCheck --saveFitResult --rMin $RMIN --rMax $RMAX $FILENAME >> logCCC2
#mv higgsCombineTest.ChannelCompatibilityCheck.*.root ../$FILENAME.CCC2.root
#
#rm -f roostats*
#rm -f higgsCombineTest*.root

done

date
echo "done"

