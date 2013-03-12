#!/bin/bash
echo "running notlxbatch.sh Simple Style"
date
for i in *.txt; do
    [[ -e "$i" ]] || continue
FILENAME=$i
echo "executing combine -M Asymptotic $FILENAME >& $FILENAME.out"

combine -M Asymptotic -v 2 $FILENAME >& $FILENAME.out
rm -f roostats*
rm -f higgsCombineTest*.root

done

date
echo "done"
