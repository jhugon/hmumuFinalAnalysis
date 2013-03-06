#!/bin/bash

NJOBS=`ls *.txt | wc -l`
sed "s/YAYYAYYAY/$NJOBS/" gofHPC_Template.sh > hpc.sh

chmod +x hpc.sh

qsub hpc.sh

echo "Submitted $NJOBS jobs"
