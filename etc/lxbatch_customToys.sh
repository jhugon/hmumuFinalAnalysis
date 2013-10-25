#!/bin/bash
echo "Sourcing cmsset_default.sh"
cd /afs/cern.ch/cms/sw
source cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc462
echo "SCRAM_ARCH is $SCRAM_ARCH"
cd $LS_SUBCWD
echo "In Directory: "
pwd
eval `scramv1 runtime -sh`
echo "cmsenv success!"
date

TXTSUFFIX=".txt"
FILENAME=$1
DIRNAME="Dir"$1_"$2"_"Dir"
ROOTFILENAME=${1%$TXTSUFFIX}.root
TOYFILENAME="toyData_"${1%$TXTSUFFIX}.root

mkdir $DIRNAME
cp $FILENAME $DIRNAME/
cp $ROOTFILENAME $DIRNAME/
cd $DIRNAME

ntoys=$2

echo "TOYFILENAME: $TOYFILENAME"
linkcommand="ln -s -T ../$TOYFILENAME $TOYFILENAME"
echo "running: $linkcommand"
$linkcommand

command="combine -M Asymptotic -t $ntoys --toysFile=$TOYFILENAME $FILENAME"
echo $command
$command >& $FILENAME.outToys
mv higgs*Asymptotic*.root ../$FILENAME.outToys.root

command="combine -M Asymptotic"
echo $command
$command >& $FILENAME.out
mv higgs*Asymptotic*.root ../$FILENAME.out.root


echo "contents of dir:"
ls -lh

cp *.txt.* ../.

#cd ..
#rm -r $DIRNAME

echo "done"
date

