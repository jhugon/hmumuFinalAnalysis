#!/bin/bash

#Vars to be replaced: TYPE, NJOBS, TIME, PREFIX

##########################################################

PREFIX="CombSplitAll"
TYPE="Comb"
TIME="7:00:00"

NJOBS=$(ls "$PREFIX"*.txt | wc -l)
sed "s/TYPE/$TYPE/;s/NJOBS/$NJOBS/;s/TIME/$TIME/;s/PREFIX/$PREFIX/" hpcTemplate.sh > hpc$TYPE.sh

chmod +x hpc$TYPE.sh

qsub hpc$TYPE.sh

echo "Submitted $NJOBS $TYPE jobs with prefix: '$PREFIX' and walltime $TIME"

##########################################################

PREFIX="Jet2"
TYPE="2Jet"
TIME="0:30:00"

NJOBS=$(ls "$PREFIX"*.txt | wc -l)
sed "s/TYPE/$TYPE/;s/NJOBS/$NJOBS/;s/TIME/$TIME/;s/PREFIX/$PREFIX/" hpcTemplate.sh > hpc$TYPE.sh

chmod +x hpc$TYPE.sh

qsub hpc$TYPE.sh

echo "Submitted $NJOBS $TYPE jobs with prefix: '$PREFIX' and walltime $TIME"

##########################################################

PREFIX="Jets01Fail"
TYPE="01Fail"
TIME="1:00:00"

NJOBS=$(ls "$PREFIX"*.txt | wc -l)
sed "s/TYPE/$TYPE/;s/NJOBS/$NJOBS/;s/TIME/$TIME/;s/PREFIX/$PREFIX/" hpcTemplate.sh > hpc$TYPE.sh

chmod +x hpc$TYPE.sh

qsub hpc$TYPE.sh

echo "Submitted $NJOBS $TYPE jobs with prefix: '$PREFIX' and walltime $TIME"

##########################################################

PREFIX="Jets01Pass"
TYPE="01Pass"
TIME="1:30:00"

NJOBS=$(ls "$PREFIX"*.txt | wc -l)
sed "s/TYPE/$TYPE/;s/NJOBS/$NJOBS/;s/TIME/$TIME/;s/PREFIX/$PREFIX/" hpcTemplate.sh > hpc$TYPE.sh

chmod +x hpc$TYPE.sh

qsub hpc$TYPE.sh

echo "Submitted $NJOBS $TYPE jobs with prefix: '$PREFIX' and walltime $TIME"

##########################################################

PREFIX="Jets01Split"
TYPE="01Split"
TIME="3:00:00"

NJOBS=$(ls "$PREFIX"*.txt | wc -l)
sed "s/TYPE/$TYPE/;s/NJOBS/$NJOBS/;s/TIME/$TIME/;s/PREFIX/$PREFIX/" hpcTemplate.sh > hpc$TYPE.sh

chmod +x hpc$TYPE.sh

qsub hpc$TYPE.sh

echo "Submitted $NJOBS $TYPE jobs with prefix: '$PREFIX' and walltime $TIME"

##########################################################

PREFIX="CombJets01PassBBBOJet2Split"
TYPE="BBBO2J"
TIME="1:00:00"

NJOBS=$(ls "$PREFIX"*.txt | wc -l)
sed "s/TYPE/$TYPE/;s/NJOBS/$NJOBS/;s/TIME/$TIME/;s/PREFIX/$PREFIX/" hpcTemplate.sh > hpc$TYPE.sh

chmod +x hpc$TYPE.sh

qsub hpc$TYPE.sh

echo "Submitted $NJOBS $TYPE jobs with prefix: '$PREFIX' and walltime $TIME"

##########################################################

PREFIX="CombJets01PassBBJet2Split"
TYPE="BB2J"
TIME="1:00:00"

NJOBS=$(ls "$PREFIX"*.txt | wc -l)
sed "s/TYPE/$TYPE/;s/NJOBS/$NJOBS/;s/TIME/$TIME/;s/PREFIX/$PREFIX/" hpcTemplate.sh > hpc$TYPE.sh

chmod +x hpc$TYPE.sh

qsub hpc$TYPE.sh

echo "Submitted $NJOBS $TYPE jobs with prefix: '$PREFIX' and walltime $TIME"

##########################################################

PREFIX="CombJets01PassBBJet2VBFPass"
TYPE="BBVBF"
TIME="1:00:00"

NJOBS=$(ls "$PREFIX"*.txt | wc -l)
sed "s/TYPE/$TYPE/;s/NJOBS/$NJOBS/;s/TIME/$TIME/;s/PREFIX/$PREFIX/" hpcTemplate.sh > hpc$TYPE.sh

chmod +x hpc$TYPE.sh

qsub hpc$TYPE.sh

echo "Submitted $NJOBS $TYPE jobs with prefix: '$PREFIX' and walltime $TIME"

##########################################################

PREFIX="CombJets01PassJet2Split"
TYPE="01JPass2J"
TIME="7:00:00"

NJOBS=$(ls "$PREFIX"*.txt | wc -l)
sed "s/TYPE/$TYPE/;s/NJOBS/$NJOBS/;s/TIME/$TIME/;s/PREFIX/$PREFIX/" hpcTemplate.sh > hpc$TYPE.sh

chmod +x hpc$TYPE.sh

qsub hpc$TYPE.sh

echo "Submitted $NJOBS $TYPE jobs with prefix: '$PREFIX' and walltime $TIME"



