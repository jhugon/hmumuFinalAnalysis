#!/usr/bin/env python

from xsec import *
from helpers import *
import ROOT as root
import os

dataDir = "input/lowPtCuts/"

CUTS="dimuonMass < 160. && dimuonMass > 110."

meNames = ["sigME","bakME","sigMEPdf","bakMEPdf"]
runPeriods = ["7TeV","8TeV"]
print("MENormDict = {}")
for RUNPERIOD in runPeriods:
   print("MENormDict['%s'] = {}" % (RUNPERIOD))
for RUNPERIOD in runPeriods:
  print ("#Running %s" % RUNPERIOD)
  meMaxes = {}
  for name in meNames:
    meMaxes[name] = 0.0
  for datasetName in ["ggHmumu125","DYJetsToLL"]:
    print ("#  Running %s" % datasetName)
    filename = dataDir+datasetName+"_"+RUNPERIOD+".root"
    rootFile = root.TFile(filename)
    tree = rootFile.Get("outtree")
    nEntries = tree.GetEntries()
    #toPrint = nEntries/10
    #toPrint /= 1000
    #toPrint *= 1000
    for i in range(nEntries):
      #if (i % toPrint) == 0:
      #  print ("#    Event: %i" % i)
      tree.GetEntry(i)
      if tree.dimuonMass < 110.:
        continue
      for name in meNames:
        meMaxes[name] = max(meMaxes[name],getattr(tree,name))
  for name in meNames:
    print("MENormDict['%s']['%s'] = %.4g" % (RUNPERIOD,name,1./meMaxes[name]))
