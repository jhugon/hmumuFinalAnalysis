#!/bin/bash

#Vars to be replaced: TYPE, NJOBS, TIME, PREFIX

##########################################################

PREFIX="CombSplitAll"
TYPE="Comb"
TIME="0:30:00"

NJOBS=$(ls "$PREFIX"*.txt | wc -l)
sed "s/TYPE/$TYPE/;s/NJOBS/$NJOBS/;s/TIME/$TIME/;s/PREFIX/$PREFIX/" hpcTemplate.sh > hpc$TYPE.sh

chmod +x hpc$TYPE.sh

qsub hpc$TYPE.sh

echo "Submitted $NJOBS $TYPE jobs with prefix: '$PREFIX' and walltime $TIME"

##########################################################

PREFIX="Jet2"
TYPE="2Jet"
TIME="0:5:00"

NJOBS=$(ls "$PREFIX"*.txt | wc -l)
sed "s/TYPE/$TYPE/;s/NJOBS/$NJOBS/;s/TIME/$TIME/;s/PREFIX/$PREFIX/" hpcTemplate.sh > hpc$TYPE.sh

chmod +x hpc$TYPE.sh

qsub hpc$TYPE.sh

echo "Submitted $NJOBS $TYPE jobs with prefix: '$PREFIX' and walltime $TIME"

##########################################################

PREFIX="Jets01Fail"
TYPE="01Fail"
TIME="0:10:00"

NJOBS=$(ls "$PREFIX"*.txt | wc -l)
sed "s/TYPE/$TYPE/;s/NJOBS/$NJOBS/;s/TIME/$TIME/;s/PREFIX/$PREFIX/" hpcTemplate.sh > hpc$TYPE.sh

chmod +x hpc$TYPE.sh

qsub hpc$TYPE.sh

echo "Submitted $NJOBS $TYPE jobs with prefix: '$PREFIX' and walltime $TIME"

##########################################################

PREFIX="Jets01Pass"
TYPE="01Pass"
TIME="0:10:00"

NJOBS=$(ls "$PREFIX"*.txt | wc -l)
sed "s/TYPE/$TYPE/;s/NJOBS/$NJOBS/;s/TIME/$TIME/;s/PREFIX/$PREFIX/" hpcTemplate.sh > hpc$TYPE.sh

chmod +x hpc$TYPE.sh

qsub hpc$TYPE.sh

echo "Submitted $NJOBS $TYPE jobs with prefix: '$PREFIX' and walltime $TIME"

##########################################################

PREFIX="Jets01Split"
TYPE="01Split"
TIME="0:10:00"

NJOBS=$(ls "$PREFIX"*.txt | wc -l)
sed "s/TYPE/$TYPE/;s/NJOBS/$NJOBS/;s/TIME/$TIME/;s/PREFIX/$PREFIX/" hpcTemplate.sh > hpc$TYPE.sh

chmod +x hpc$TYPE.sh

qsub hpc$TYPE.sh

echo "Submitted $NJOBS $TYPE jobs with prefix: '$PREFIX' and walltime $TIME"

