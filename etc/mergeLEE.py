#!/usr/bin/env python

import re, glob, sys, os, os.path, string
from subprocess import call
import ROOT as root
root.gErrorIgnoreLevel = root.kWarning

massList = range(115,155+1)

minToys = 5

##########################################

filesDict = {}
for mass in massList:
  filesDict[str(mass)] = []
for iSeed in range(100):
  iSeedStr = str(iSeed)
  if len(glob.glob("CombSplitAll_7P8TeV_*.0.txt."+iSeedStr+".sigToys.root"))==0:
    continue
  badMassList = []
  for mass in massList:
    fn="CombSplitAll_7P8TeV_"+str(mass)+".0.txt."+iSeedStr+".sigToys.root"
    if not os.path.exists(fn):
      badMassList.append( mass)
      continue
    try:
      rf = root.TFile(fn)
      limit = rf.Get("limit")
      if limit.GetEntries() < minToys:
        badMassList.append( mass)
        continue
    except:
      #print("Root exception on file: %s" % fn)
      badMassList.append( mass)
      continue
  if len(badMassList) > 0:
    print("Bad masses for Seed %i " % iSeed)
    for mass in badMassList:
      print("  %i " % mass)
  else:
    print("Seed %i is done" % iSeed)
    for mass in massList:
      filesDict[str(mass)].append(iSeed)
for mass in massList:
  seedList = filesDict[str(mass)]
  if len(seedList) == 0:
    continue
  inputFns = []
  outputFn = "Added"+str(mass)+".root"
  for iSeed in seedList:
    iSeedStr = str(iSeed)
    fn="CombSplitAll_7P8TeV_"+str(mass)+".0.txt."+iSeedStr+".sigToys.root"
    inputFns.append(fn)
  command = ["hadd","-f",outputFn]+inputFns
  print string.join(command)
  call(command)
