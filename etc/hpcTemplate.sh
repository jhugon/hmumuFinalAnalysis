#! /bin/bash
#
#PBS -r n

##Job settings
#PBS -N statsJob
#PBS -o outstatsJob
#PBS -e errstatsJob
#PBS -m a
#PBS -M jhugon@phys.ufl.edu

#Multiple Job Submission:
#Jobs will have a variable called $PBS_ARRAYID
#that will be one of the following numbers
#PBS -t 1-YAYYAYYAY

##Job Resources
#PBS -l walltime=00:10:00
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

echo "executing combine -M Asymptotic $FILENAME >& $FILENAME.out"

combine -M Asymptotic $FILENAME >& $FILENAME.out

echo "executing combine -M ProfileLikelihood -d $FILENAME --signif --usePLC >& $FILENAME.sig"

combine -M ProfileLikelihood -d $FILENAME --signif --usePLC >& $FILENAME.sig
rm -f roostats*
rm -f higgsCombineTest*.root

echo "executing combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig"

combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 >& $FILENAME.expsig
##combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig
rm -f roostats*
rm -f higgsCombineTest*.root

echo "executing combine -M MaxLikelihoodFit --rMax 50 --plots --saveNormalizations $FILENAME >& $FILENAME.mu"

combine -M MaxLikelihoodFit --rMax 50 --plots --saveNormalizations $FILENAME >& $FILENAME.mu
rm -f roostats*
rm -f higgsCombineTest*.root

cp $FILENAME.out ..
cp $FILENAME.mu ..
cp $FILENAME.sig ..
cp mlfit.root ../$FILENAME.root
for subname in *_fit_s.png; do
  cp $subname ../${FILENAME%$TXTSUFFIX}_$subname
done
#cp $FILENAME.expsig ..

echo "done"
ENDTIME=`date +%s`
echo "Took $(( $ENDTIME - $STARTTIME )) seconds"
date

