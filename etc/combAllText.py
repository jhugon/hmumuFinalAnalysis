#!/usr/bin/env python

import sys
import os
import subprocess
import glob
import re

searchFor = {
  "Jets01PassCatAll" : "Jets01PassPtG10",
  "Jets01FailCatAll" : "Jets01FailPtG10",
  "Jets01SplitCatAll": "Jets01",
  "Jet2SplitCutsGFSplit" : "Jet2",
  "CombSplitAll" : "Jet",
}
nCards = {
  "CombSplitAll" : 15,
  "Jet2SplitCutsGFSplit" : 3,
  "Jets01SplitCatAll": 12,
  "Jets01PassCatAll" : 6,
  "Jets01FailCatAll" : 6,
}
combinations = [
  "CombSplitAll",
  "Jet2SplitCutsGFSplit",
  "Jets01SplitCatAll",
  "Jets01PassCatAll" ,
  "Jets01FailCatAll" ,
]
energies = ["7TeV","8TeV"]

allSingleCats = [
  "Jets01PassPtG10BB",
  "Jets01PassPtG10BO",
  "Jets01PassPtG10BE",
  "Jets01PassPtG10OO",
  "Jets01PassPtG10OE",
  "Jets01PassPtG10EE",
  "Jets01FailPtG10BB",
  "Jets01FailPtG10BO",
  "Jets01FailPtG10BE",
  "Jets01FailPtG10OO",
  "Jets01FailPtG10OE",
  "Jets01FailPtG10EE",
  "Jet2CutsVBFPass",
  "Jet2CutsGFPass",
  "Jet2CutsFailVBFGF",
]

nullFile = open(os.devnull,'r')

for mass in range(120,151):
  if len(glob.glob("*{0}.0.txt".format(mass)))==0:
      continue
  for combination in combinations:
    files7TeV = []
    files8TeV = []
    labels7TeV = []
    labels8TeV = []
    for cat in allSingleCats:
      if re.match("^"+searchFor[combination],cat):
        files7TeV.append("{0}_7TeV_{1}.0.txt".format(cat,mass))
        files8TeV.append("{0}_8TeV_{1}.0.txt".format(cat,mass))
        labels7TeV.append("{0}7TeV".format(cat))
        labels8TeV.append("{0}8TeV".format(cat))

    outName = combination+"_7P8TeV_{0:.1f}.txt".format(float(mass))
    outFile = open(outName,'w')
    callArgs = ["combineCards.py"]
    for infile,label in zip(files7TeV,labels7TeV):
      callArgs.append(label+"="+infile)
    for infile,label in zip(files8TeV,labels8TeV):
      callArgs.append(label+"="+infile)
    try:
      subprocess.check_call(callArgs,stdout=outFile,stderr=nullFile)
      #subprocess.check_call(callArgs,stdout=outFile)
      outFile.close()
    except:
      print "Error: Couldn't make ",outName
      outFile.close()
      os.remove(outName)
    
    outName = combination+"_7TeV_{0:.1f}.txt".format(float(mass))
    outFile = open(outName,'w')
    callArgs = ["combineCards.py"]
    for infile,label in zip(files7TeV,labels7TeV):
      callArgs.append(label+"="+infile)
    try:
      subprocess.check_call(callArgs,stdout=outFile,stderr=nullFile)
      #subprocess.check_call(callArgs,stdout=outFile)
      outFile.close()
    except:
      print "Error: Couldn't make ",outName
      outFile.close()
      os.remove(outName)
    
    outName = combination+"_8TeV_{0:.1f}.txt".format(float(mass))
    outFile = open(outName,'w')
    callArgs = ["combineCards.py"]
    for infile,label in zip(files8TeV,labels8TeV):
      callArgs.append(label+"="+infile)
    try:
      subprocess.check_call(callArgs,stdout=outFile,stderr=nullFile)
      #subprocess.check_call(callArgs,stdout=outFile)
      outFile.close()
    except:
      print "Error: Couldn't make ",outName
      outFile.close()
      os.remove(outName)

  # combined energy single channels
  for cat in allSingleCats:
    outName = cat+"_7P8TeV_{0:.1f}.txt".format(float(mass))

    callArgs = ["combineCards.py","{0}7TeV={0}_7TeV_{1}.0.txt".format(cat,mass),"{0}8TeV={0}_8TeV_{1}.0.txt".format(cat,mass)]
    outFile = open(outName,'w')
    try:
      subprocess.check_call(callArgs,stdout=outFile,stderr=nullFile)
      #subprocess.check_call(callArgs,stdout=outFile)
      outFile.close()
    except:
      print "Error: Couldn't make ",outName
      outFile.close()
      os.remove(outName)
      
nullFile.close()
