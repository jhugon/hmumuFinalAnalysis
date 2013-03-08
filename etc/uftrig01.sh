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

echo "executing combine -M MaxLikelihoodFit --rMax 50 --plots --saveNormalizations $FILENAME >& $FILENAME.mu"

combine -M MaxLikelihoodFit --rMax 50 --plots --saveNormalizations $FILENAME >& $FILENAME.mu
cp mlfit.root $FILENAME.root
for subname in *_fit_s.png; do
  cp $subname ${FILENAME%$TXTSUFFIX}_$subname
done

wait

rm -f roostats*
rm -f higgsCombineTest*.root

done

date
echo "done"

