#! /usr/bin/env python

import datetime
import sys
import os
import os.path
import re
import math
import cPickle
import glob
import ROOT as root
root.gROOT.SetBatch(True)

from helpers import *

import scipy.stats
from numpy import mean, median, corrcoef, percentile
from numpy import std as stddev

#root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(1) #For MINUIT

TITLEMAP = {
  "Jets01PassPtG10BB": "0,1-Jet Tight BB",
  "Jets01PassPtG10BO": "0,1-Jet Tight BO",
  "Jets01PassPtG10BE": "0,1-Jet Tight BE",
  "Jets01PassPtG10OO": "0,1-Jet Tight OO",
  "Jets01PassPtG10OE": "0,1-Jet Tight OE",
  "Jets01PassPtG10EE": "0,1-Jet Tight EE",
  "Jets01PassCatAll" : "0,1-Jet Tight Combination",
                        
  "Jets01FailPtG10BB": "0,1-Jet Loose BB",
  "Jets01FailPtG10BO": "0,1-Jet Loose BO",
  "Jets01FailPtG10BE": "0,1-Jet Loose BE",
  "Jets01FailPtG10OO": "0,1-Jet Loose OO",
  "Jets01FailPtG10OE": "0,1-Jet Loose OE",
  "Jets01FailPtG10EE": "0,1-Jet Loose EE",
  "Jets01FailCatAll" : "0,1-Jet Loose Combination",

  "Jet2CutsVBFPass":"2-Jet VBF Tight",
  "Jet2CutsGFPass":"2-Jet GF Tight",
  "Jet2CutsFailVBFGF":"2-Jet Loose",
}

PDFTITLEMAP = {
    "ExpLog":"Exp(p_{1}m^{2}+p_{2}m+p_{3}ln(m))",
    "MOverSq":"#frac{m}{(m-p_{1})^{2}}",
    "Old":"Voigtian+Exp",
    "ExpMOverSq":"#frac{Exp(p_{1}m)}{(m-p_{2})^{2}}",
    "Bernstein":"Bernstein",
}


################################################################################################
################################################################################################
################################################################################################

def getData(refDir,altDir):
  refFns = glob.glob(refDir+"*.muToys.root")
  altFns = glob.glob(altDir+"*.muToys.root")
  refFns = sorted(refFns)
  altFns = sorted(altFns)
  refJobNs = []
  prefix="CombSplitAll_7P8TeV_125.0.txt"
  for refFn in refFns:
    refMatch = re.search(prefix+r"\.([0-9]+)\.muToys\.root",refFn)
    if refMatch:
      refJobNs.append(refMatch.group(1))
    else:
      refJobNs.append(-1)
  altJobNs = []
  for altFn in altFns:
    altMatch = re.search(prefix+r"\.([0-9]+)\.muToys\.root",altFn)
    if altMatch:
      altJobNs.append(altMatch.group(1))
    else:
      altJobNs.append(-1)

  altJobNsSet = set(altJobNs)
  refJobNsSet = set(refJobNs)

  refDifSet = refJobNsSet.difference(altJobNsSet)
  altDifSet = altJobNsSet.difference(refJobNsSet)

  for i in reversed(range(len(altFns))):
    n = altJobNs[i]
    if n == -1 or n in altDifSet:
      altFns.pop(i)
      altJobNs.pop(i)

  for i in reversed(range(len(refFns))):
    n = refJobNs[i]
    if n == -1 or n in refDifSet:
      refFns.pop(i)
      refJobNs.pop(i)

  assert(len(altFns)==len(refFns))
  assert(len(altFns)==len(altJobNs))
  assert(len(refFns)==len(refJobNs))

  result = []

  for i in range(len(altFns)):
    altFn = altFns[i]
    refFn = refFns[i]
    altN = altJobNs[i]
    refN = refJobNs[i]
    assert(altN==refN)
    altTF = root.TFile(altFn)
    refTF = root.TFile(refFn)
    altTree = altTF.Get("limit")
    refTree = refTF.Get("limit")

    # Figure out poper iToys
    iToySetAlt = set()
    for iLimit in range(3,altTree.GetEntries(),4):
      altTree.GetEntry(iLimit)
      iToySetAlt.add(altTree.iToy)
    iToySetRef = set()
    for iLimit in range(3,refTree.GetEntries(),4):
      refTree.GetEntry(iLimit)
      iToySetRef.add(refTree.iToy)
    iToySetBoth = iToySetAlt.intersection(iToySetRef)
    for iToy in iToySetBoth:
      for iLimit in range(3,refTree.GetEntries(),4):
        refTree.GetEntry(iLimit)
        if refTree.iToy == iToy:
            break
      for iLimit in range(3,altTree.GetEntries(),4):
        altTree.GetEntry(iLimit)
        if altTree.iToy == iToy:
            break
      assert(altTree.iToy==refTree.iToy)
      assert(altTree.iSeed==refTree.iSeed)
      point = {}
      point['muRef'] = refTree.limit
      point['muAlt'] = altTree.limit
      point['muErrRef'] = refTree.limitErr
      point['muErrAlt'] = altTree.limitErr
      result.append(point)
  return result
      

if __name__ == "__main__":
  inputRefDir= "/afs/cern.ch/user/j/jhugon/work/private/stats/CMSSW_6_1_1/finalAnalysisVoigtExp110160/statsCards/"
  inputAltDir= "/afs/cern.ch/user/j/jhugon/work/private/stats/CMSSW_6_1_1/finalAnalysis/statsCards/"
  data = getData(inputRefDir,inputAltDir)

  refMus = [toy['muRef'] for toy in data]
  altMus = [toy['muAlt'] for toy in data]
  refErrMus = [toy['muErrRef'] for toy in data]
  altErrMus = [toy['muErrAlt'] for toy in data]
  refZs = [toy['muRef']/toy['muErrRef'] for toy in data]
  altZs = [toy['muAlt']/toy['muErrAlt'] for toy in data]

  pulls = [(toy['muAlt']-toy['muRef'])/toy['muErrRef'] for toy in data]

  print "\n#########################################\n"
  print "NToys: ",len(data)
  print "\n#########################################\n"
  print "Reference Median Mu: ",median(refMus)
  print "Reference Mean Mu: ",mean(refMus)
  print "Reference Sigma Mu: ",stddev(refMus)
  print "Alternate Median Mu: ",median(altMus)
  print "Alternate Mean Mu: ",mean(altMus)
  print "Alternate Sigma Mu: ",stddev(altMus)
  print "\n#########################################\n"
  print "Reference Median ErrMu: ",median(refErrMus)
  print "Reference Mean ErrMu: ",mean(refErrMus)
  print "Reference Sigma ErrMu: ",stddev(refErrMus)
  print "Alternate Median ErrMu: ",median(altErrMus)
  print "Alternate Mean ErrMu: ",mean(altErrMus)
  print "Alternate Sigma ErrMu: ",stddev(altErrMus)
  print "\n#########################################\n"
  print "Reference Median Z: ",median(refZs)
  print "Reference Mean Z: ",mean(refZs)
  print "Reference Sigma Z: ",stddev(refZs)
  print "Alternate Median Z: ",median(altZs)
  print "Alternate Mean Z: ",mean(altZs)
  print "Alternate Sigma Z: ",stddev(altZs)
  print "\n#########################################\n"
  print "Median Bias: ",median(pulls)
  print "\n#########################################\n"
  
