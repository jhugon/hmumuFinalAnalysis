#!/bin/bash
echo "running uftrig01.sh"
date
for i in *.txt; do
    [[ -e "$i" ]] || continue
FILENAME=$i
echo "executing combine -M Asymptotic $FILENAME >& $FILENAME.out"

combine -M Asymptotic $FILENAME >& $FILENAME.out &

echo "executing combine -M ProfileLikelihood -d $FILENAME --signif --usePLC >& $FILENAME.sig"

combine -M ProfileLikelihood -d $FILENAME --signif --usePLC >& $FILENAME.sig &

echo "executing combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig"

combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 >& $FILENAME.expsig
#combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig

echo "executing combine -M MaxLikelihoodFit --rMin -50 --rMax 50 --plots --saveNormalizations $FILENAME >& $FILENAME.mu"

combine -M MaxLikelihoodFit --rMin -50 --rMax 50 --plots --saveNormalizations $FILENAME >& $FILENAME.mu
cp mlfit.root $FILENAME.root
for subname in *_fit_s.png; do
  cp $subname ${FILENAME%$TXTSUFFIX}_$subname
done

combine -M ChannelCompatibilityCheck --saveFitResult --rMax 50 $FILENAME >> logCCC
mv higgsCombineTest.ChannelCompatibilityCheck.*.root ../$FILENAME.CCC.root

rm -f roostats*
rm -f higgsCombineTest*.root


combine -M ChannelCompatibilityCheck --saveFitResult --rMin -150 --rMax 150 $FILENAME >> logCCC2
mv higgsCombineTest.ChannelCompatibilityCheck.*.root ../$FILENAME.CCC2.root

rm -f roostats*
rm -f higgsCombineTest*.root

wait

rm -f roostats*
rm -f higgsCombineTest*.root

done

date
echo "done"

