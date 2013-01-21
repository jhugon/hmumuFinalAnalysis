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
  "IncPresel":"Inclusive Category",
  "VBFPresel":"VBF Category",

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
  def __init__(self,filenames,categories,massBoudaries,chopDir=True,chopExt=True):
    data = {}
    self.filenames = filenames
    self.data = data
    self.categories = categories
    for f in filenames:
      rf = root.TFile(f)
      fDir = os.path.dirname(f)
      fName = os.path.basename(f)
      fNameSplitExt = os.path.splitext(fName)
      fNameNoExt = fNameSplitExt[0]
      fNameExt = fNameSplitExt[1]
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
      fNameKey = fNameNoExt
      if not chopDir:
        fNameKey = fDir+'/'+fNameNoExt
      if not chopExt:
        fNameKey += fNameExt
      data[fNameKey] = {}
      data[fNameKey]["misc"] = {
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
        data[fNameKey][i] = nEvents

def compareDirs(dirNameDict,sigFileDict,bakFileNames,categories=["IncPresel","VBFPresel"],massBoundaries=[120,130]):
  data = {}
  dirNames = sorted(dirNameDict.keys())
  signals = sorted(sigFileDict.keys())
  #signals.append("bak")
  #sigFileDict["bak"] = "Background"
  for dirName in dirNames:
    data[dirName] = {}
    for sig in signals:
      data[dirName][sig] = {}
    data[dirName]["bak"] = {}
    data[dirName]["sig"] = {}
    data[dirName]["sob"] = {}

  energyStr = ""
  lumi = ""
  for dirName in dirNames:
    tmpCounts = Counts([dirName+fn for fn in signals+bakFileNames],categories,massBoundaries,False,False)
    for cat in categories:
      nbak = 0.0
      nsig = 0.0
      for sig in signals:
        energyStr = tmpCounts.data[dirName+sig]["misc"]["energyStr"]
        lumi = tmpCounts.data[dirName+sig]["misc"]["lumi"]
        tmp = tmpCounts.data[dirName+sig][cat]
        data[dirName][sig][cat] = tmp
        nsig += tmp
      for bak in bakFileNames:
        nbak += tmpCounts.data[dirName+bak][cat]
      data[dirName]["bak"][cat] = nbak
      data[dirName]["sig"][cat] = nsig
      if nbak > 0:
        data[dirName]["sob"][cat] = nsig/nbak
      else:
        data[dirName]["sob"][cat] = 0.0

  signals.append("bak")
  sigFileDict["bak"] = "Background"
  sigFileDict["sig"] = "Signal"
  sigFileDict["sob"] = "$S/B$"
  signals = ["sig","bak","sob"]

  ncols = len(categories)*len(signals)
  outString = ""
  # Header Line 1
  nsignals = len(signals)
  #outString += "Category &"
  outString += " &"
  for i in categories:
    outString += r" \multicolumn{"+str(nsignals)+r"}{|"+'c|'
    outString += r"}{"+channelNameMap[i]+"} &"
  outString = outString.rstrip(r"&")
  outString += r"\\ \hline" + '\n'

  #outString += "Sample &"
  outString += " &"
  for i in categories:
    for j in signals:
      outString += " "+ sigFileDict[j] + " &"
  outString = outString.rstrip(r"&")
  outString += r"\\ \hline" + '\n'

  for dirName in dirNames:
    outString += dirNameDict[dirName] + " &"
    for i in categories:
      for j in signals:
        #print("dir: {0:<18} sig: {1:<10} cat: {2:<10}".format(dirName,j,i))
        #print("dir: {0:<18} sig: {1:<10} cat: {2:<10} n: {3}".format(dirName,j,i,data[dirName][j][i]))
        if j == "bak":
          outString += " {0:<5.0f}".format(data[dirName][j][i]) + " &"
        else:
          outString += " {0:<5.2e}".format(data[dirName][j][i]) + " &"
    outString = outString.rstrip(r"&")
    outString += r"\\ \hline" + '\n'
  
  outString = r"\begin{tabular}{|l|"+'c|'*ncols+"} \hline"+"\n" + outString + r"\end{tabular}"+"\n"
  outString = outString.replace(r"%",r"\%")

  outString = re.sub(r"([-\d.+]+)e([-+])0([\d])",r"\1e\2\3",outString)
  outString = re.sub(r"([-\d.+]+)e([-\d+]+)",r"$\1 \\times 10^{\2}$",outString)

  lumi = "$\\mathcal{{L}}$ = {0:.0f} fb$^{{-1}}$".format(float(lumi))
  energyStr = "$\\sqrt{s}=$"+re.sub(r"TeV"," TeV",energyStr)
  outString += r"\\ "+lumi+", "+energyStr

  return '\n'+outString+'\n'
    
      
if __name__ == "__main__":
  
  """
  filenames = glob.glob("input/vladEventCounts/*.root")
  categories = ["VBFPresel"]
  mBounds = [110.,160.]
  c = Counts(filenames,categories,mBounds)

  print("=============================\nFor Vladimir: ({0})".format(categories[0]))
  fns = sorted(c.data.keys())
  maxFNLength = str(max([len(i) for i in fns])+2)
  for fn in fns:
    n = c.data[fn][categories[0]]
    toPrint = r"{0:<"+maxFNLength+r"} {1:<20.0f}"
    print(toPrint.format(fn,n))
  print("=============================")
  """

 ######################################################

  filenames = glob.glob("input/trk*/gg*.root")
  filenames += glob.glob("input/pf*/gg*.root")
  categories = ["VBFPresel"]
  mBounds = [110.,160.]

  dirs = {
    "input/trkLooseIso/":"Trk Loose Iso",
    "input/trkTightIso/":"Trk Tight Iso",
    "input/pfLooseIso/":"PF Loose Iso",
    "input/pfTightIso/":"PF Tight Iso"
  }
  sigFileNames = {
    "ggHmumu125_8TeV.root":"gg H",
    "vbfHmumu125_8TeV.root":"VBF H"
  }
  bakFileNames = [
    "SingleMuRun2012Av1.root",
    "SingleMuRun2012Bv1.root",
    "SingleMuRun2012Cv1.root",
    "SingleMuRun2012Cv2.root"
  ]

  print compareDirs(dirs,sigFileNames,bakFileNames)
