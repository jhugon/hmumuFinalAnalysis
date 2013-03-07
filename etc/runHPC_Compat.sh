#!/bin/bash

NJOBS=`ls *.txt | wc -l`
NFILES=$NJOBS
echo "$NJOBS Seperate Cards" 1>&2
TOYJOBS=5
echo "$TOYJOBS Toy Jobs per Card" 1>&2
NJOBS=$(( $NJOBS * $TOYJOBS ))
NJOBSTORPL=$(( $NJOBS - 1 ))
sed "s/YAYYAYYAY/$NJOBSTORPL/" compatHPC_Template.sh > hpcTmp.sh
sed "s/WOWOWOWO/$TOYJOBS/" hpcTmp.sh > hpc.sh

chmod +x hpc.sh

qsub hpc.sh 1>&2

echo "Submitted $NJOBS jobs"  1>&2

MULTIPLIER=$(( $TOYJOBS-1))
echo $(( $NFILES * $MULTIPLIER ))
