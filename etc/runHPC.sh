#!/bin/bash

NJOBS=`ls *.txt | wc -l`
sed "s/YAYYAYYAY/$NJOBS/" hpcTemplate.sh > hpc.sh

chmod +x hpc.sh

qsub hpc.sh

echo "Submitted $NJOBS jobs"
