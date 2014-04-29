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

for combination in combinations:
  #for mass in range(120,151):
  for mass in range(125,126):
    outName = combination+"_7P8TeV_{0:.1f}.txt".format(float(mass))
    searchString = "{0}*{1}.0.txt".format(searchFor[combination],mass)
    files = sorted(glob.glob(searchString))
    nFiles = len(glob.glob(searchString))
    if mass == 125:
      print combination,searchString, nFiles
    if nFiles == 2*nCards[combination]:
      outFile = open(outName,'w')
      callArgs = ["combineCards.py"]
      for f in files:
        match = re.match(r"(.+)_([78P]+TeV)_[0-9.]+\.txt",f)
        assert(match)
        label=match.group(1)+match.group(2)
        callArgs.append(label+"="+f)
      #print callArgs
      #print outFile
      subprocess.call(callArgs,stdout=outFile)
      outFile.close()
    for energy in energies:
      outName = combination+"_{0}_{1:.1f}.txt".format(energy,float(mass))
      searchString = "{0}*_{1}_{2}.0.txt".format(searchFor[combination],energy,mass)
      files = sorted(glob.glob(searchString))
      nFiles = len(glob.glob(searchString))
      if mass == 125:
        print combination,searchString, nFiles
      if nFiles == nCards[combination]:
        outFile = open(outName,'w')
        callArgs = ["combineCards.py"]
        for f in files:
          match = re.match(r"(.+)_([78P]+TeV)_[0-9.]+\.txt",f)
          assert(match)
          label=match.group(1)+match.group(2)
          callArgs.append(label+"="+f)
        #print callArgs
        #print outFile
        subprocess.call(callArgs,stdout=outFile)
        outFile.close()
