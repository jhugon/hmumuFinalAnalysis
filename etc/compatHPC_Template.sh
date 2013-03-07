#! /bin/bash
#
#PBS -r n

##Job settings
#PBS -N statsJobCompat
#PBS -o outstatsJobCompat
#PBS -e errstatsJobCompat
#PBS -m a
#PBS -M jhugon@phys.ufl.edu

#Multiple Job Submission:
#Jobs will have a variable called $PBS_ARRAYID
#that will be one of the following numbers
#PBS -t 1-YAYYAYYAY

##Job Resources
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=2500mb

# initialize environment for worker
STARTTIME=`date +%s`

export SCRAM_ARCH=slc5_amd64_gcc462
export OSG_APP=/osg/app
export VO_CMS_SW_DIR=${OSG_APP}/cmssoft/cms
export CMS_PATH=${VO_CMS_SW_DIR}
cd $CMS_PATH
source cmsset_default.sh;

echo "Job running on `hostname` at `date`"
echo "ArrayId: $PBS_ARRAYID"

# Vars
DIR=$PBS_O_WORKDIR

# enter working area
cd $DIR
echo "In Dir:"
pwd

eval `scram runtime -sh`

#####################################
#####Begin Real Job##################
#####################################

ifile="1"
for f in `ls *.txt`; do
  if [ "$ifile" -eq "$PBS_ARRAYID" ]; then
    FILENAME=$f
  fi
  ifile=$(( $ifile + 1 ))
done

echo "Input Filename: $FILENAME"

if [ -z "$FILENAME" ]; then
  echo "Error: Couldn't find filename for PBS_ARRAYID!! Exiting."
  exit "1"
fi

TXTSUFFIX=".txt"
DIRNAME="Dir"$FILENAME"Dir"
ROOTFILENAME=${FILENAME%$TXTSUFFIX}.root

mkdir $DIRNAME
cp $FILENAME $DIRNAME/
cp $ROOTFILENAME $DIRNAME/
cd $DIRNAME

STARTTIME=`date +%s`
echo "running runGOF.sh"
date

iSeed="0"
while true; do
  combine -M ChannelCompatibilityCheck --saveFitResult --rMax 50 $FILENAME -t 10 -s $(( 124389 + $iSeed )) >> logCCCToys
  iSeed=$(( $iSeed + 1 ))
  if [ "$iSeed" -gt 10 ]; then
    break
  fi
done

hadd -f -k $FILENAME.CCC-Toys.root higgsCombineTest.ChannelCompatibilityCheck.*.root

rm -f roostats*
rm -f higgsCombineTest*.root

combine -M ChannelCompatibilityCheck --saveFitResult --rMax 50 $FILENAME >> logCCC
mv higgsCombineTest.ChannelCompatibilityCheck.*.root $FILENAME.CCC.root

rm -f roostats*
rm -f higgsCombineTest*.root

cp *CCC*.root ..

date
echo "done"

ENDTIME=`date +%s`
echo "Took $(( $ENDTIME - $STARTTIME )) seconds"
