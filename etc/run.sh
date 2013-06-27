#!/bin/bash

chmod +x lxbatch.sh

for i in *.txt; do
    [[ -e "$i" ]] || continue
echo "Running on "$i

string=`echo $i`;
if [[ $string =~ .*Comb.* ]]
then
    bsub -q cmscaf1nd -o /dev/null lxbatch.sh $i
else
    bsub -q cmscaf1nh -o /dev/null lxbatch.sh $i
fi

done


