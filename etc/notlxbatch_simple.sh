#!/bin/bash
echo "running notlxbatch.sh Simple Style"
date
for i in *.txt; do
    [[ -e "$i" ]] || continue
FILENAME=$i
echo "executing combine -M Asymptotic $FILENAME >& $FILENAME.out"

nice combine -M Asymptotic --rMax 150 $FILENAME >& $FILENAME.out
rm -f roostats*
rm -f higgsCombineTest*.root

done

date
echo "done"

