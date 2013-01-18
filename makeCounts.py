#!/usr/bin/env python

import math
from math import sqrt
import ROOT as root
from helpers import *
import matplotlib.pyplot as mpl
import numpy
import glob
import os.path
import re

from xsec import *

root.gErrorIgnoreLevel = root.kWarning

channelNameMap = {
  "AllCat":"All Cat. Comb.",
  "IncCat":"Inc. Cat. Comb.",
  "VBFCat":"VBF Cat. Comb.",

  "Presel":"Presel. Comb.",
  "IncPresel":"Inc. Presel.",
  "VBFPresel":"VBF Presel.",

  "Pt0to30":"$p_{T}^{\mu\mu} \in [0,30]$",
  "Pt30to50":"$p_{T}^{\mu\mu} \in [30,50]$",
  "Pt50to125":"$p_{T}^{\mu\mu} \in [50,125]$",
  "Pt125to250":"$p_{T}^{\mu\mu} \in [125,250]$",
  "Pt250":"$p_{T}^{\mu\mu}>250$",

  "VBFLoose":"VBFL",
  "VBFMedium":"VBFM",
  "VBFTight":"VBFT",
  "VBFVeryTight":"VBFVT",

  "BDTCut":"BDT Comb.",
  "IncBDTCut":"Inc. BDT",
  "VBFBDTCut":"VBF BDT",

  "BDTCutCat":"BDT Res. Comb.",
  "IncBDTCutCat":"Inc. BDT Res.",
  "VBFBDTCutCat":"VBF BDT Res.",

  "PreselCat":"Presel. Res. Comb.",
  "IncPreselCat":"Inc. Res. Presel.",
  "VBFPreselCat":"VBF Res. Presel.",

  "IncBDTCutBB":"Inc. BDT BB",
  "IncBDTCutBO":"Inc. BDT BO",
  "IncBDTCutBE":"Inc. BDT BE",
  "IncBDTCutOO":"Inc. BDT OO",
  "IncBDTCutOE":"Inc. BDT OE",
  "IncBDTCutEE":"Inc. BDT EE",
  "IncBDTCutNotBB":"Inc. BDT !BB",
  "VBFBDTCutBB":"VBF BDT BB",
  "VBFBDTCutNotBB":"VBF BDT !BB",
  "IncPreselBB":"Inc. Presel. BB",
  "IncPreselBO":"Inc. Presel. BO",
  "IncPreselBE":"Inc. Presel. BE",
  "IncPreselOO":"Inc. Presel. OO",
  "IncPreselOE":"Inc. Presel. OE",
  "IncPreselEE":"Inc. Presel. EE",
  "IncPreselNotBB":"Inc. Presel. !BB",
  "VBFPreselBB":"VBF Presel. BB",
  "VBFPreselNotBB":"VBF Presel. !BB"
}

sampleNameMap = {
  "data_obs":"Data",
  "bak":"Background",
  "sig":"Signal"
}

class Counts:
  def __init__(self,filenames,categories,massBoudaries):
    data = {}
    self.filenames = filenames
    self.data = data
    self.categories = categories
    for f in filenames:
      rf = root.TFile(f)
      fName = os.path.basename(f)
      fNameNoExt = os.path.splitext(fName)[0]
      scaleBy = 1.0
      energyStr = None
      if "Run2012" in fNameNoExt:
        energyStr = "8TeV"
      elif "Run2011" in fNameNoExt:
        energyStr = "7TeV"
      else:
        energyStr = fNameNoExt.split('_')[1]
        if xsec.has_key(fNameNoExt) and nEventsMap.has_key(fNameNoExt):
          scaleBy *= xsec[fNameNoExt]/nEventsMap[fNameNoExt]*lumiDict[energyStr]
      data[fNameNoExt] = {}
      data[fNameNoExt]["misc"] = {
        "scaleBy":scaleBy,
        "energyStr": energyStr,
        "lumi": lumiDict[energyStr]
            }
      for i in categories:
        strToGet = i + '/mDiMu'
        strToGet = os.path.normpath(strToGet)
        if strToGet[0] == '/':
            strToGet = strToGet[1:]
        tmpHist = rf.Get(strToGet)
        nEvents = getIntegralAll(tmpHist,massBoudaries)
        nEvents *= scaleBy
        data[fNameNoExt][i] = nEvents
      
if __name__ == "__main__":
  
  filenames = glob.glob("input/vladEventCounts/*.root")
  categories = ["VBFPresel"]
  mBounds = [110.,160.]
  mBounds = [0.,1000.]
  c = Counts(filenames,categories,mBounds)

  print("=============================\nFor Vladimir: ({0})".format(categories[0]))
  fns = sorted(c.data.keys())
  maxFNLength = str(max([len(i) for i in fns])+2)
  for fn in fns:
    n = c.data[fn][categories[0]]
    toPrint = r"{0:<"+maxFNLength+r"} {1:<20.0f}"
    print(toPrint.format(fn,n))
  print("=============================")

 ######################################################


