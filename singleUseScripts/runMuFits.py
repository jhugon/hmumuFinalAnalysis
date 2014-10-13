#!/usr/bin/env python

import sys
import os
import os.path
import subprocess

# catName : [lowerLimit,upperLimit] # best fit
limitMap = {}

limitMap["7TeV"] = {
#  "CombSplitAll"           : [  -1.0,  14.8], #    7.2
#  "Jets01SplitCatAll"      : [ -10.3,  10.6], #    0.2
#  "Jet2SplitCutsGFSplit"   : [   2.6,  32.4], #   18.8
#  "Jets01PassCatAll"       : [  -8.6,  13.0], #    2.3
#  "Jets01FailCatAll"       : [ -68.5,  13.5], #  -28.7
  "Jets01PassPtG10BB"      : [  -9.8,  25.2], #    8.0
  "Jets01PassPtG10BO"      : [ -26.2, 20.0], #   -5.9
  "Jets01PassPtG10BE"      : [ -23.6,  44.5], #   11.2
  "Jets01PassPtG10OO"      : [ -17.0,  35.4], #    9.5
  "Jets01PassPtG10OE"      : [ -64.6,  33.0], #  -15.8
  "Jets01PassPtG10EE"      : [-140.0, -59.2], # -104.6
  "Jets01FailPtG10BB"      : [-150.0, -42.5], #  -96.3
  "Jets01FailPtG10BO"      : [ -93.8,  35.4], #  -30.3
  #"Jets01FailPtG10BE"      : [-300.0, 300.0], #  150.0
  "Jets01FailPtG10BE"      : [-100.0, 350.0], #  150.0
  "Jets01FailPtG10OO"      : [-138.0, 119.0], #   -9.3
  "Jets01FailPtG10OE"      : [-160.0, 200.0], #   24.4
  "Jets01FailPtG10EE"      : [-380.0, 200.0], # -108.0
  "Jet2CutsVBFPass"        : [ -15.0,  20.0], #    8.4
  "Jet2CutsGFPass"         : [ -15.0,  60.0], #   24.0
  "Jet2CutsFailVBFGF"      : [  0.0,  80.0], #   41.2
}

limitMap["8TeV"] = {
#  "CombSplitAll"           : [  -4.6,   3.0], #   -0.8
#  "Jets01SplitCatAll"      : [   1.0,  11.6], #    6.4
#  "Jet2SplitCutsGFSplit"   : [ -12.0,  -2.8], #   -7.2
#  "Jets01PassCatAll"       : [   1.8,  12.9], #    7.5
#  "Jets01FailCatAll"       : [ -35.3,  11.3], #  -12.4
  "Jets01PassPtG10BB"      : [   4.9,  21.1], #   13.3
  "Jets01PassPtG10BO"      : [  -3.5,  17.6], #    7.2
  "Jets01PassPtG10BE"      : [ -19.1,  19.4], #    0.1
  "Jets01PassPtG10OO"      : [ -17.8,  10.1], #   -3.9
  "Jets01PassPtG10OE"      : [ -26.7,  58.4], #   16.2
  "Jets01PassPtG10EE"      : [-100.4,   9.8], #  -45.7
  "Jets01FailPtG10BB"      : [ -47.2,  22.2], #  -13.0
  "Jets01FailPtG10BO"      : [ -49.8,  35.6], #   -7.4
  "Jets01FailPtG10BE"      : [ -35.0, 155.0], #   62.5
  "Jets01FailPtG10OO"      : [-110.0,  10.0], #  -48.8
  "Jets01FailPtG10OE"      : [ -75.0, 150.0], #   36.5
  "Jets01FailPtG10EE"      : [-180.0,  20.0], #  -83.4
  "Jet2CutsVBFPass"        : [ -12.2,  -1.9], #   -6.8
  "Jet2CutsGFPass"         : [ -14.9,  16.8], #    1.0
  "Jet2CutsFailVBFGF"      : [ -33.1,  -2.8], #  -18.2
}

dirName = "statsCards/"
os.chdir(dirName)

#for period in ["7TeV","8TeV"]:
for period in ["8TeV"]:
  for cat in limitMap[period]:
    rLimits = limitMap[period][cat]
    logFileName =  "{0}_{1}_125.0.txt.mu2".format(cat,period)
    argString = "combine -M MultiDimFit --rMin={2} --rMax={3} --algo=singles --cl=0.68 {0}_{1}_125.0.txt".format(cat,period,rLimits[0],rLimits[1])
    if period == "8TeV":
      if "FailPtG10BE" in cat or "FailPtG10OO" in cat or "FailPtG10EE" in cat:
        argString = "combine -M MultiDimFit --rMin={2} --rMax={3} --robustFit=1 --stepSize=0.01 --algo=singles --cl=0.68 {0}_{1}_125.0.txt".format(cat,period,rLimits[0],rLimits[1])
    logFile = open(logFileName,"w")
    argList = argString.split(" ")
    print argString
    #print argList
    print subprocess.call(argList,stdout=logFile,stderr=logFile)
    logFile.close()
    
