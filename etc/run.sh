#!/bin/bash

chmod +x lxbatch.sh

for i in *.txt; do
    [[ -e "$i" ]] || continue
echo "Running on "$i
bsub -q cmscaf1nh -o /dev/null lxbatch.sh $i
#bsub -q 1nh -o /dev/null lxbatch.sh $i
done

