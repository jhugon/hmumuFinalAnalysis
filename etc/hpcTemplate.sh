#! /bin/bash
#
#PBS -r n

##Job settings
#PBS -N statsTYPE
#PBS -o joboutStatsJobTYPE
#PBS -e joberrStatsJobTYPE

#Multiple Job Submission:
#Jobs will have a variable called $PBS_ARRAYID
#that will be one of the following numbers
#PBS -t 1-NJOBS

##Job Resources
#PBS -l walltime=TIME
#PBS -l nodes=1:ppn=1
#PBS -l pmem=5000mb

# initialize environment for worker
starttime=$(date +%s)

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
for f in `ls PREFIX*.txt`; do
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

RMAX=150
RMIN=-150

echo "executing combine -M Asymptotic --rMax $RMAX $FILENAME >& $FILENAME.out"

combine -M Asymptotic --rMax $RMAX $FILENAME >& $FILENAME.out

echo "executing combine -M ProfileLikelihood --rMax $RMAX -d $FILENAME --signif >& $FILENAME.sig"

combine -M ProfileLikelihood --rMax $RMAX -d $FILENAME --signif >& $FILENAME.sig
rm -f roostats*
rm -f higgsCombineTest*.root

#echo "executing combine -M ProfileLikelihood --rMax $RMAX -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig"
#
#combine -M ProfileLikelihood --rMax $RMAX -d $FILENAME --signif --expectSignal=1 -t -1 >& $FILENAME.expsig
##combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig
#rm -f roostats*
#rm -f higgsCombineTest*.root

#echo "executing combine -M MaxLikelihoodFit --rMin $RMIN --rMax $RMAX --plots --saveNormalizations $FILENAME >& $FILENAME.mu"
#
#combine -M MaxLikelihoodFit --rMin $RMIN --rMax $RMAX --plots --saveNormalizations $FILENAME >& $FILENAME.mu
#rm -f roostats*
#rm -f higgsCombineTest*.root

#combine -M ChannelCompatibilityCheck --saveFitResult --rMin $RMIN --rMax $RMAX $FILENAME >> logCCC2
#mv higgsCombineTest.ChannelCompatibilityCheck.*.root ../$FILENAME.CCC2.root
#
#rm -f roostats*
#rm -f higgsCombineTest*.root

cp $FILENAME.* ..
cp mlfit.root ../$FILENAME.root

date
endtime=$(date +%s)
echo "Job took "$((endtime - starttime))" seconds"
echo "Job took "$(echo "scale = 2; ($endtime - $starttime)/60." | bc)" minutes"
echo "Job took "$(echo "scale = 2; ($endtime - $starttime)/3600." | bc)" hours"
echo "done"

