#!/usr/bin/env python

import optparse
parser = optparse.OptionParser(description="Merge two StatsInput Into One")
parser.add_option("-p","--pvalue", help="Select PValues Independently Using the Max Between the Two",action="store_true",default=False)
args, fakeargs = parser.parse_args()

import ROOT as root
import glob
import re
import os.path
import cPickle
from copy import deepcopy
from array import array
import sys

##############################################
# Only use Matplotlib if not in CMSSW_6*
cmsswVersion = ""
mplGood = True
if os.environ.has_key("CMSSW_VERSION"):
  cmsswVersion = os.environ["CMSSW_VERSION"]
#def ComparePlotTable(*arg):
#  #Dummy so that there are no errors
#  pass
if "CMSSW_6" in cmsswVersion:
  mplGood = False
if mplGood:
  from etc.evilMatPlotlibFunctions import *
##############################################


def selMaxObsLimit(fname,dirNameOne,dirNameTwo):

  obsLimitOne = -999.
  obsLimitTwo = -999.

  if os.path.exists(dirNameOne+fname):
    dataOne = getData(dirNameOne+fname)
    obsLimitOne = float(dataOne[1])

  if os.path.exists(dirNameTwo+fname):
    dataTwo = getData(dirNameTwo+fname)
    obsLimitTwo = float(dataTwo[1])
    
  if ( obsLimitOne >= obsLimitTwo ):
    return dirNameOne+fname
  else:
    return dirNameTwo+fname


def selMinSign(fname,dirNameOne,dirNameTwo):

  pValueOne = +999.
  pValueTwo = +999.

  if os.path.exists(dirNameOne+fname):
    dataOne = getDataSig(dirNameOne+fname)
    pValueOne = float(dataOne[1])
    
  if os.path.exists(dirNameTwo+fname):
    dataTwo = getDataSig(dirNameTwo+fname)
    pValueTwo = float(dataTwo[1])
    
  if ( pValueTwo >= pValueOne ):
    return dirNameOne+fname
  else:
    return dirNameTwo+fname



# Do match and don't match w/o extensions.  \.expsig and \.sig are added automatically
def getDataSig(fname,matchString=r"_([-\d.]+)\.txt",dontMatchStrings=[],doSort=True,getPValue=False,xMax=1e20):
  def sortfun(s):
    match = re.search(matchString,s)
    result = 1e12
    if match:
      result = float(match.group(1))
    return result

  obs = +100.0
  xNum = +100.0

  tmpF = open(fname)
  match = re.search(matchString+r"\.sig",fname)
  if match:
    xNum = match.group(1)
  
  for line in tmpF:
    obsMatch = None
    if getPValue:
      obsMatch = re.search(r"p-value = ([.\deE]+)",line)
    else:
      obsMatch = re.search(r"^Significance:[\s]+([.\deE]+)",line)
    if obsMatch:
      obs = obsMatch.group(1)

  thisPoint = [xNum,obs]
  return thisPoint

      
# edited getData to return the content from one file 
def getData(fname,matchString=r"_([-\d.]+)\.txt\.out",dontMatchStrings=[],doSort=True,xMax=1.0e20):
  def sortfun(s):
    match = re.search(matchString,s)
    result = 1e12
    if match:
      result = float(match.group(1))
    return result

  obs = -10.0
  median = -10.0
  low2sig = -10.0
  low1sig = -10.0
  high1sig = -10.0
  high2sig = -10.0
  xNum = -10.0

  tmpF = open(fname)
  match = re.search(matchString,fname)
  if match:
    xNum = match.group(1)

  for line in tmpF:
    obsMatch = re.search(r"Observed[\s]Limit:[^.\d]*< ([.\deE]+)",line)
    low2sigMatch = re.search(r"Expected.*2\.5.:[^.\d]*< ([.\deE]+)",line)
    low1sigMatch = re.search(r"Expected.*16.0[^.\d]*< ([.\deE]+)",line)
    medianMatch = re.search(r"Expected.*50\.0.*< ([.\deE]+)",line)
    high1sigMatch = re.search(r"Expected.*84.0.*< ([.\deE]+)",line)
    high2sigMatch = re.search(r"Expected.*97.5.*< ([.\deE]+)",line)
    if obsMatch:
      obs = obsMatch.group(1)
    if low2sigMatch:
      low2sig = low2sigMatch.group(1)
    if low1sigMatch:
      low1sig = low1sigMatch.group(1)
    if medianMatch:
      median = medianMatch.group(1)
    if high1sigMatch:
      high1sig = high1sigMatch.group(1)
    if high2sigMatch:
      high2sig = high2sigMatch.group(1)

  thisPoint = [xNum,obs,low2sig,low1sig,median,high1sig,high2sig]
  return thisPoint



if __name__ == "__main__":

  dirName          = "/data/uftrig01b/digiovan/baselinePP/m110to160_pixelLumi/hmumuFinalAnalysis/statsInput/"
  dirNameToCompare = "/data/uftrig01b/digiovan/baselinePP/m110to160_pixelLumi_voigtexp_zmass_zSig_Constr/hmumuFinalAnalysis_voigtexp_zmass_zSig_Constr/statsInput/"
  outDir = "statsInput/"
  
  #root.gErrorIgnoreLevel = root.kWarning
  #root.gROOT.SetBatch(True)

  # file to keep the record
  listOfLimits = open(outDir + 'OBSLIMITSREADME', 'w')
  listOfSignificances  = open(outDir + 'SIGNIFICANCEREADME', 'w')
 
  #for period in ["7TeV","8TeV","14TeV","7P8TeV"]:
  for period in ["7TeV","8TeV","7P8TeV"]:
    fnToGlob = dirName+"*_"+period+"_*.txt.out"
    allfiles = glob.glob(fnToGlob)
    
    for file in sorted(allfiles): 

      #if "CombSplitAll" not in file:
      #  continue
      #if "7P8TeV" not in period:
      #  continue

      # this will be useful for the comparison
      thefilename = file.replace(dirName,"")
   
      selectedLimit = selMaxObsLimit(thefilename,dirName,dirNameToCompare)
      os.system( "cp %s %s" % (selectedLimit,outDir) )
      listOfLimits.write("%s\n" % selectedLimit)

      if (args.pvalue):
        selectedSign = selMinSign(thefilename.replace(".out",".sig"),dirName,dirNameToCompare)
        os.system( "cp %s %s" % (selectedSign,outDir) )
        listOfSignificances.write("%s\n" % selectedSign)
      else:
        # just pick the one which corresponds to the max obs limit
        os.system( "cp %s %s" % (selectedLimit.replace(".out",".sig"),outDir) )
        listOfSignificances.write("%s\n" % selectedLimit)
        

  listOfLimits.close()
  if (args.pvalue):
    listOfSignificances.close()
