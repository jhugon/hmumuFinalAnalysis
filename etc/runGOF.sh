#!/bin/bash

# Takes 30 min of UFTrig01 w/ 1 processor for 
#   BDTCutCatVBFBDTOnly_7P8TeV_125.0 1000 toys

STARTTIME=`date +%s`
echo "running runGOF.sh"
date
for i in *.txt; do
    [[ -e "$i" ]] || continue
FILENAME=$i

iSeed="0"
while true; do
  combine -M GoodnessOfFit --algorithm saturated --fixedSignalStrength=0.0 $FILENAME -t 100 -s $(( 124389 + $iSeed )) >> logGOF
  iSeed=$(( $iSeed + 1 ))
  if [ "$iSeed" -gt 10 ]; then
    break
  fi
done

hadd -f -k $FILENAME.GOF-Toys.root higgsCombineTest.GoodnessOfFit.*.root

rm -f roostats*
rm -f higgsCombineTest*.root

combine -M GoodnessOfFit --algorithm saturated --fixedSignalStrength=0.0 $FILENAME >> logGOF
mv higgsCombineTest.GoodnessOfFit.*.root $FILENAME.GOF.root

rm -f roostats*
rm -f higgsCombineTest*.root

done

date
echo "done"

ENDTIME=`date +%s`
echo "Took $(( $ENDTIME - $STARTTIME )) seconds"
