#!/bin/bash

NJOBS=`ls *.txt | wc -l`
echo "$NJOBS Seperate Cards"
TOYJOBS=3
echo "$TOYJOBS Toy Jobs per Card"
NJOBS=$(( $NJOBS * $TOYJOBS ))
sed "s/YAYYAYYAY/$NJOBS/" compatHPC_Template.sh > hpcTmp.sh
sed "s/WOWOWOWO/$TOYJOBS/" hpcTmp.sh > hpc.sh

chmod +x hpc.sh

qsub hpc.sh

echo "Submitted $NJOBS jobs"
