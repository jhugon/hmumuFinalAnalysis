#!/usr/bin/env python

import math
from math import sqrt
import ROOT as root
from helpers import *
import random
import glob
import os.path
import re

from xsec import *

import makeCards

root.gErrorIgnoreLevel = root.kWarning

GLOBALCOUNTER = 0

channelNameMap = {
  "AllCat":"All Cat. Comb.",
  "IncCat":"Non-VBF Cat. Comb.",
  "VBFCat":"VBF Cat. Comb.",

  "Presel":"Presel. Comb.",
  "IncPresel":"Non-VBF Presel",
  "VBFPresel":"VBF Presel",

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
  "IncBDTCut":"Non-VBF BDT",
  "VBFBDTCut":"VBF BDT",

  "BDTCutCat":"BDT Res. Comb.",
  "IncBDTCutCat":"Non-VBF BDT Res.",
  "VBFBDTCutCat":"VBF BDT Res.",

  "PreselCat":"Presel. Res. Comb.",
  "IncPreselCat":"Non-VBF Res. Presel.",
  "VBFPreselCat":"VBF Res. Presel.",

  "IncBDTCutBB":"Non-VBF BDT BB",
  "IncBDTCutBO":"Non-VBF BDT BO",
  "IncBDTCutBE":"Non-VBF BDT BE",
  "IncBDTCutOO":"Non-VBF BDT OO",
  "IncBDTCutOE":"Non-VBF BDT OE",
  "IncBDTCutEE":"Non-VBF BDT EE",
  "IncBDTCutNotBB":"Non-VBF BDT !BB",
  "VBFBDTCutBB":"VBF BDT BB",
  "VBFBDTCutNotBB":"VBF BDT Non-BB",
  "IncPreselBB":"Non-VBF Presel. BB",
  "IncPreselBO":"Non-VBF Presel. BO",
  "IncPreselBE":"Non-VBF Presel. BE",
  "IncPreselOO":"Non-VBF Presel. OO",
  "IncPreselOE":"Non-VBF Presel. OE",
  "IncPreselEE":"Non-VBF Presel. EE",
  "IncPreselNotBB":"Non-VBF Presel. !BB",
  "VBFPreselBB":"VBF Presel. BB",
  "VBFPreselNotBB":"VBF Presel. Non-BB",

  "IncPreselPtG10BB":"Non-VBF BB",
  "IncPreselPtG10BO":"Non-VBF BO",
  "IncPreselPtG10BE":"Non-VBF BE",
  "IncPreselPtG10OO":"Non-VBF OO",
  "IncPreselPtG10OE":"Non-VBF OE",
  "IncPreselPtG10EE":"Non-VBF EE",
  "IncPreselPtG10":"Non-VBF",
  "BDTCutCatVBFBDTOnly": "VBF & Non-VBF Combination",
  "all": "Inclusive",
}

sampleNameMap = {
  "data_obs":"Data",
  "data":"Data",
  "dataWide":"Data",
  "mcbak":"MC Background",
  "bak":"Background Prediction",
  "bakWide":"Data in Sidebands",
  "bakCheat":"Back Pred (Norm Cheat)",
  "bakErr":"Back Pred Diff",
  "width":"68% Signal Width (GeV)",
  "sig":"Signal",
  "sob":"$S/B$",
  "sosqrtb":"$S/\sqrt{B}$",
  "sosqrtsb":"$S/\sqrt{S+B}$",
  "ggHmumu125_8TeV":r"$gg \rightarrow H$",
  "vbfHmumu125_8TeV":r"VBF H",
  "wHmumu125_8TeV":r"WH",
  "zHmumu125_8TeV":r"ZH",
  "ggHmumu125_7TeV":r"$gg \rightarrow H$",
  "vbfHmumu125_7TeV":r"VBF H",

  "SingleMuRun2012Av1":"Run2012AReReco13Jul2012-v1",
  "SingleMuRun2012Bv1":"Run2012BReReco13Jul2012-v1",
  "SingleMuRun2012Cv1":"Run2012CReReco24Aug2012-v1",
  "SingleMuRun2012Cv2":"Run2012C-PromptReco-v2",
  "SingleMuRun2012D":"Run2012D-PromptReco-v1",
  "SingleMuRun2012Av1Recover":"Run2012A-recover-06Aug2012-v1",
}

def getFormatName(map,cat):
  if map.has_key(cat):
    return map[cat]
  elif "JetsL2" in cat:
    return cat.replace("JetsL2","$<2$ Jets ")
  elif "Jets2" in cat:
    return cat.replace("Jets2","$\geq2$ Jets ")
  elif "Jets1" in cat:
    return cat.replace("Jets1","1 Jet ")
  elif "Jets0" in cat:
    return cat.replace("Jets0","0 Jet ")
  else:
    return cat

def findMassBoundariesDG(fileNames,categories,cuts,level,useDG=True):
  global GLOBALCOUNTER
  massExtraCutString = " && dimuonMass > 110. && dimuonMass < 160."
  treeList = []
  rfList = []
  for fName in fileNames:
    rf = root.TFile(fName)
    rfList += [rf]
    tree = rf.Get("outtree")
    treeList += [tree]
  treesToUse = []
  for cat,cut in zip(categories,cuts):
    maxTree = None
    maxN = -1
    if cat != "":
      cut += massExtraCutString
    else:
      cut += massExtraCutString[3:]
    cutStr = treeCut(cat,cut)
    #print cutStr
    for tree in treeList:
      tmpHistName = 'hist'+str(GLOBALCOUNTER)
      GLOBALCOUNTER += 1
      varToDraw = "dimuonMass"
      drawStr = varToDraw+" >> "+tmpHistName
      #print drawStr
      tree.Draw(drawStr,cutStr)
      tmp = root.gDirectory.Get(tmpHistName)
      nEvents = getIntegralAll(tmp)
      if nEvents > maxN:
        maxTree = tree
        maxN = nEvents
    treesToUse += [maxTree]
  results = []
  for cat,cut,tree in zip(categories,cuts,treesToUse):
      tmpHistName = 'hist'+str(GLOBALCOUNTER)
      GLOBALCOUNTER += 1
      varToDraw = "dimuonMass"
      drawStr = varToDraw+" >> "+tmpHistName
      #print cat
      #print drawStr
      if cat != "":
        cut += massExtraCutString
      else:
        cut += massExtraCutString[3:]
      cutStr = treeCut(cat,cut)
      #print cutStr
      tree.Draw(drawStr,cutStr)
      tmp = root.gDirectory.Get(tmpHistName)
      quantiles = None
      if useDG:
        quantiles = fitDGFindQuantiles(tmp,level)
      else:
        quantiles = getMedianAndQuantileInterval(tmp,level)
      results += [[quantiles[2],quantiles[0]]]
  return results

class Counts:
  def __init__(self,filenames,categories,cuts,massBoundaries):
    global GLOBALCOUNTER
    assert(type(massBoundaries)==list and len(massBoundaries)>0)
    if type(massBoundaries[0]) != list:
      massBoundariesBak = massBoundaries
      massBoundaries = []
      for i in categories:
        massBoundaries += [massBoundariesBak]
    assert(len(massBoundaries)==len(categories))
    assert(len(categories)==len(cuts))
    data = {}
    self.filenames = filenames
    self.data = data
    self.categories = categories
    self.massBoundaries = massBoundaries
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
        if nEventsMap.has_key(fNameNoExt):
          scaleBy = xsec[fNameNoExt]/nEventsMap[fNameNoExt]*lumiDict[energyStr]*1000.
      fNameKey = fNameNoExt
      data[fNameKey] = {}
      data[fNameKey]["misc"] = {
        "scaleBy":scaleBy,
        "energyStr": energyStr,
        "lumi": lumiDict[energyStr]
            }
      tree = rf.Get("outtree")
      for i,j,mb in zip(categories,cuts,massBoundaries):
        tmpHistName = 'hist'+str(GLOBALCOUNTER)

        GLOBALCOUNTER += 1
        massExtraCutString = " && {0} < dimuonMass && {1} > dimuonMass".format(*mb)
        varToDraw = "dimuonMass"
        drawStr = varToDraw+" >> "+tmpHistName
        #print i
        #print drawStr
        if j != "":
          j += massExtraCutString
        else:
          j += massExtraCutString[3:]
        cutStr = treeCut(i,j)
        #print cutStr
        tree.Draw(drawStr,cutStr)
        tmp = root.gDirectory.Get(tmpHistName)
        nEvents = getIntegralAll(tmp)
        nEvents *= scaleBy
        if i == '':
          i = "all"
        data[fNameKey][i] = nEvents

  def getBoundariesDict(self):
    result = {}
    for cat,mb in zip(self.categories,self.massBoundaries):
      result[cat] = mb
    return result
  def getWidthDict(self):
    result = self.getBoundariesDict()
    for cat in result:
      w = result[cat][1]-result[cat][0]
      result[cat] = w
    return result

def getBakPrediction(filenames,categories,massBoundaries,cheat=False):
  assert(len(massBoundaries)==4)
  maxMass = massBoundaries[3]
  minMass = massBoundaries[0]
  mMuMu = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)
  #mMuMu.setRange("z",88,94)
  #mMuMu.setRange("verylow",controlRegionVeryLow[0],controlRegionVeryLow[1])
  mMuMu.setRange("low",minMass,massBoundaries[1])
  mMuMu.setRange("high",massBoundaries[2],maxMass)
  mMuMu.setRange("signal",massBoundaries[1],massBoundaries[2])
  data = {}
  files = []
  for f in filenames:
    rf = root.TFile(f)
    files.append(rf)
  for i in categories:
    sumHist = None
    for rf in files:
      strToGet = i + '/mDiMu'
      strToGet = os.path.normpath(strToGet)
      if strToGet[0] == '/':
          strToGet = strToGet[1:]
      tmpHist = rf.Get(strToGet)
      if sumHist == None:
        sumHist = tmpHist.Clone()
      else:
        sumHist.Add(tmpHist)
    if i == '':
      i = 'all'
    tmpName = i+str(random.randint(0,10000))
    class Stupid:
      pass
    histIntegral = getIntegralAll(sumHist,[minMass,maxMass])
    tmpParams, bakNormTup = makeCards.makePDFBak(tmpName,sumHist,mMuMu,
                                       minMass,maxMass,None)
    print "integral, estimate for all: {}, {}".format(histIntegral,bakNormTup[0]*bakNormTup[1])
    if cheat:
      nBak = (1.0-1.0/bakNormTup[1])*histIntegral
    else:
      nBak = (1.0-1.0/bakNormTup[1])*bakNormTup[0]*bakNormTup[1]
    data[i] = nBak
  return data

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
        elif j == "sob":
          outString += " {0:<5.2e}".format(data[dirName][j][i]) + " &"
        else:
          outString += " {0:<5.3g}".format(data[dirName][j][i]) + " &"
    outString = outString.rstrip(r"&")
    outString += r"\\ \hline" + '\n'
  
  outString = r"\begin{tabular}{|l|"+'c|'*ncols+"} \hline"+"\n" + outString + r"\end{tabular}"+"\n"
  outString = outString.replace(r"%",r"\%")

  outString = re.sub(r"([-\d.+]+)e([-+])0([\d])",r"\1e\2\3",outString)
  outString = re.sub(r"([-\d.+]+)e[+]*([-\d]+)",r"$\1 \\times 10^{\2}$",outString)

  lumi = "$\\mathcal{{L}}$ = {0:.0f} fb$^{{-1}}$".format(float(lumi))
  energyStr = "$\\sqrt{s}$ = "+re.sub(r"TeV"," TeV",energyStr)
  massStr = r"{0} GeV $< m_{{\mu\mu}} <$ {1} GeV".format(*massBoundaries)
  outString += r"\\ "+massStr+", "+energyStr + ", "+lumi+"\n"

  return '\n'+outString+'\n'

def cutFlow(sigFileNames,bakFileNames,channel,massBoundaries=[120.,130.]):
  data = {}
  data["bak"] = {}
  data["sig"] = {}
  data["sob"] = {}
  data["sosqrtb"] = {}
  data["sosqrtsb"] = {}

  categories = ['',"IncPresel","VBFPresel","VBFBDTCut","IncBDTCut"]

  energyStr = ""
  lumi = ""
  sigCounts = Counts(sigFileNames,categories,massBoundaries)
  bakCounts = Counts(bakFileNames,categories,massBoundaries)
  for sig in sigCounts.data:
    data[sig] = {}
  categories.remove('')
  categories.append("all")
  for cat in categories:
    nBak = 0.0
    for bak in bakCounts.data:
      nBak += bakCounts.data[bak][cat]
    data["bak"][cat] = nBak
    nSig = 0.0
    for sig in sigCounts.data:
      energyStr = sigCounts.data[sig]["misc"]["energyStr"]
      lumi = sigCounts.data[sig]["misc"]["lumi"]
      nSig += sigCounts.data[sig][cat]
      data[sig][cat] = sigCounts.data[sig][cat]
    data["sig"][cat] = nSig
    data["sob"][cat] = float(nSig)/nBak
    data["sosqrtb"][cat] = float(nSig)/sqrt(float(nBak))
    data["sosqrtsb"][cat] = float(nSig)/sqrt(float(nBak+nSig))

  samples = ["ggHmumu125_8TeV","vbfHmumu125_8TeV","bak","sob","sosqrtsb"]
  ncols = len(samples)

  outString = " &"
  for i in samples:
    outString += " "+sampleNameMap[i]+" &"
  outString = outString.rstrip(r"&")
  outString += r"\\ \hline" + '\n'

  outString += "No Cuts &"
  for s in samples:
    if s=="bak" or s=="sob" or s=="sosqrtb" or s=="sosqrtsb" :
        outString += "  &"
        continue
    outString += " {0:.2f} &".format(xsec[s]*lumiDict[energyStr]*1000.)
  outString = outString.rstrip(r"&")
  outString += r"\\ \hline" + '\n'

  outString += "Muon Cuts &"
  for s in samples:
    n = data[s]['all']
    if s=="bak" or s=="sob" or s=="sosqrtb" or s=="sosqrtsb" :
       outString += " {0:.2e} &".format(n)
       continue
    outString += " {0:.2f} &".format(n)
  outString = outString.rstrip(r"&")
  outString += r"\\ \hline" + '\n'

  if channel == "Inc":
    outString += "!VBF Presel. &"
  else:
    outString += "VBF Presel. &"
  for s in samples:
    n = data[s][channel+"Presel"]
    if s=="bak" or s=="sob" or s=="sosqrtb" or s=="sosqrtsb" :
       outString += " {0:.2e} &".format(n)
       continue
    outString += " {0:.2f} &".format(n)
  outString = outString.rstrip(r"&")
  outString += r"\\ \hline" + '\n'

  if channel == "Inc":
    outString += "!VBF BDT Cut &"
  else:
    outString += "VBF BDT Cut &"
  for s in samples:
    n = data[s][channel+"BDTCut"]
    if s=="bak" or s=="sob" or s=="sosqrtb" or s=="sosqrtsb" :
       outString += " {0:.2e} &".format(n)
       continue
    outString += " {0:.2f} &".format(n)
  outString = outString.rstrip(r"&")
  outString += r"\\ \hline" + '\n'

  outString = r"\begin{tabular}{|l|"+'c|'*ncols+"} \hline"+"\n" + outString + r"\end{tabular}"+"\n"
  outString = outString.replace(r"%",r"\%")

  outString = re.sub(r"([-\d.+]+)e([-+])0([\d])",r"\1e\2\3",outString)
  outString = re.sub(r"([-\d.+]+)e[+]*([-\d]+)",r"$\1 \\times 10^{\2}$",outString)

  lumi = "$\\mathcal{{L}}$ = {0:.0f} fb$^{{-1}}$".format(float(lumi))
  energyStr = "$\\sqrt{s}$ = "+re.sub(r"TeV"," TeV",energyStr)
  massStr = r"{0} GeV $< m_{{\mu\mu}} <$ {1} GeV".format(*massBoundaries)
  outString += r"\\ "+massStr+", "+energyStr + ", "+lumi+"\n"

  return outString

def inMassWindowTableMkr(sigFileNames,bakFileNames,datFileNames,categories,cuts,massBoundaries=None,samples=None):
  data = {}
  data["bak"] = {}
  data["bakCheat"] = {}
  data["bakErr"] = {}
  data["mcbak"] = {}
  data["sig"] = {}
  data["sob"] = {}
  data["sosqrtb"] = {}
  data["sosqrtsb"] = {}
  data["bakWide"] = {}
  data["data"] = {}
  data["dataWide"] = {}

  energyStr = ""
  lumi = ""
  if massBoundaries == None:
    massBoundaries = findMassBoundariesDG(sigFileNames,categories,cuts,0.68)
  sigCounts = Counts(sigFileNames,categories,cuts,massBoundaries)
  bakCounts = Counts(bakFileNames,categories,cuts,massBoundaries)
  datCounts = Counts(datFileNames,categories,cuts,massBoundaries)
  datCountsWide = Counts(datFileNames,categories,cuts,[110.,160.])
  for key in sigCounts.data.keys()+datCounts.data.keys():
    data[key] = {}
  for cat in categories:
    if cat == '':
        cat = "all"
    nBak = 0.0
    for bak in bakCounts.data:
      nBak += bakCounts.data[bak][cat]
    data["mcbak"][cat] = nBak
    nSig = 0.0
    for sig in sigCounts.data:
      energyStr = sigCounts.data[sig]["misc"]["energyStr"]
      lumi = sigCounts.data[sig]["misc"]["lumi"]
      nSig += sigCounts.data[sig][cat]
      data[sig][cat] = sigCounts.data[sig][cat]
    nDat = 0.0
    for dat in datCounts.data:
      nDat += datCounts.data[dat][cat]
      data[dat][cat] = datCounts.data[dat][cat]
    data["data"][cat] = nDat
    nDatWide = 0.0
    for dat in datCountsWide.data:
      nDatWide += datCountsWide.data[dat][cat]
    data["dataWide"][cat] = nDatWide

    data["sig"][cat] = nSig
    if nBak>0:
      data["sob"][cat] = float(nSig)/nBak
      data["sosqrtb"][cat] = float(nSig)/sqrt(float(nBak))
    else:
      data["sob"][cat] = float('nan')
      data["sosqrtb"][cat] = float('nan')
    if nSig+nBak>0:
      data["sosqrtsb"][cat] = float(nSig)/sqrt(float(nBak+nSig))
    else:
      data["sosqrtsb"][cat] = float('nan')

  #data["bak"] = getBakPrediction(datFileNames,categories,bakPredBounds)
  #data["bakCheat"] = getBakPrediction(datFileNames,categories,bakPredBounds,cheat=True)

  if categories.count('')==1:
    allIndex = categories.index('')
    categories[allIndex] = "all"
#  countsWide = Counts(bakFileNames,categories,[110.,160.])
#  for cat in categories:
#    nBakWide = 0.0
#    for bak in bakCounts.data:
#      nBakWide += countsWide.data[bak][cat]
#    data["bakWide"][cat] = nBakWide-data["bak"][cat]
#  for cat in categories:
#    data["bakErr"][cat] = (data["bak"][cat]-data["data"][cat])/data["data"][cat]

  if samples==None:
    samples = ["sig","mcbak","data","sob","width"]
    samples = ["width","sob"]
    #samples = sorted(datCounts.data.keys()) + ['data']
    #samples = ["ggHmumu125_8TeV","vbfHmumu125_8TeV","wHmumu125_8TeV","zHmumu125_8TeV","sig"]
    #samples += ["ggHmumu125_8TeV","vbfHmumu125_8TeV","wHmumu125_8TeV","zHmumu125_8TeV"]
    #samples += ["bakCheat"]
    #samples += ["bakErr"]
  ncols = len(samples)

  widths = sigCounts.getWidthDict()

  outString = ""
  if True:
    outString += " &"
    for i in samples:
      outString += " "+sampleNameMap[i]+" &"
    outString = outString.rstrip(r"&")
    outString += r"\\ \hline" + '\n'
  
    for cat in categories:
      outString += getFormatName(channelNameMap,cat) +" &"
      for s in samples:
        if s=="width":
          n = widths[cat]
          outString += " {0:.1f} &".format(n)
          continue
        n = data[s][cat]
        if s=="sob" or s=="sosqrtb" or s=="sosqrtsb":
           outString += " {0:.2e} &".format(n)
        elif s=="bakErr":
           outString += " {0:.2%} &".format(n)
        elif s=="data_obs" or s=="data" or s=="dataWide" or ("SingleMuRun" in s):
          outString += " {0:.0f} &".format(n)
        else:
          outString += " {0:.2f} &".format(n)
      outString = outString.rstrip(r"&")
      outString += r"\\ \hline" + '\n'

    outString = r"\begin{tabular}{|l|"+'c|'*ncols+"} \hline"+"\n" + outString + r"\end{tabular}"+"\n"
    outString = outString.replace(r"%",r"\%")
  
    outString = re.sub(r"([-\d.+]+)e([-+])0([\d])",r"\1e\2\3",outString)
    outString = re.sub(r"([-\d.+]+)e[+]*([-\d]+)",r"$\1 \\times 10^{\2}$",outString)
  
    lumi = "$\\mathcal{{L}}$ = {0:.0f} fb$^{{-1}}$".format(float(lumi))
    energyStr = "$\\sqrt{s}$ = "+re.sub(r"TeV"," TeV",energyStr)
    massStr = r"{0} GeV $< m_{{\mu\mu}} <$ {1} GeV".format(*massBoundaries)
    #outString += r"\\ "+massStr+", "+energyStr + ", "+lumi+"\n"
    outString += r"\\ "+energyStr + ", "+lumi+"\n"
  if True:
    firstLen = 4
    for i in samples:
      tmp = len(sampleNameMap[i])
      if firstLen < tmp:
        firstLen = tmp
    firstLen = str(firstLen+1)
    printString = ("{0:<"+firstLen+"} ").format("")
    colWidths = []
    for i in samples:
      colLen = len(sampleNameMap[i])+3
      if colLen < 6:
        colLen = 6
      colLen = str(colLen)
      printString += ("{0:^"+colLen+"}").format(sampleNameMap[i])
      colWidths.append(colLen)
    printString += "\n"
    for cat in categories:
      printString += ("{0:<"+firstLen+"} ").format(getFormatName(channelNameMap,cat))
      for s,colLen in zip(samples,colWidths):
        if s=="width":
          n = widths[cat]
          printString += ("{0:^"+colLen+".1f}").format(n)
          continue
        n = data[s][cat]
        if s=="sob" or s=="sosqrtb" or s=="sosqrtsb":
           printString += ("{0:^"+colLen+".2e}").format(n)
        elif s=="bakErr":
           printString += ("{0:"+colLen+".2%}").format(n)
        elif s=="data_obs" or s=="data" or s=="dataWide" or ("SingleMu" in s):
          printString += ("{0:^"+colLen+".0f}").format(n)
        else:
          printString += ("{0:^"+colLen+".2f}").format(n)
      printString += '\n'

    printString += '\n'

    printString += ("{0:<"+firstLen+"}").format("")
    printString += ("{0:^"+colLen+"}").format("Width (GeV)")
    printString += '\n'
    for cat in categories:
      printString += ("{0:<"+firstLen+"} ").format(getFormatName(channelNameMap,cat))
      printString += ("{0:^"+colLen+".2f}").format(widths[cat])
      printString += '\n'
    print(printString)

  return outString
      
if __name__ == "__main__":
  root.gROOT.SetBatch(True)
#  
#  filenames = glob.glob("input/vladEventCounts/*.root")
#  categories = ["VBFPresel"]
#  mBounds = [110.,160.]
#  c = Counts(filenames,categories,mBounds)
#
#  print("=============================\nFor Vladimir: ({0})".format(categories[0]))
#  fns = sorted(c.data.keys())
#  maxFNLength = str(max([len(i) for i in fns])+2)
#  for fn in fns:
#    n = c.data[fn][categories[0]]
#    toPrint = r"{0:<"+maxFNLength+r"} {1:<20.0f}"
#    print(toPrint.format(fn,n))
#  print("=============================")
#
#  ######################################################
#
#  filenames = glob.glob("input/trk*/gg*.root")
#  filenames += glob.glob("input/pf*/gg*.root")
#  categories = ["VBFPresel"]
#  #mBounds = [110.,160.]
#  mBounds = [120.,130.]
#
#  dirs = {
#    "input/trkLooseIso/":"Trk Loose Iso",
#    "input/trkTightIso/":"Trk Tight Iso",
#    "input/pfLooseIso/":"PF Loose Iso",
#    "input/pfTightIso/":"PF Tight Iso"
#  }
#  sigFileNames = {
#    "ggHmumu125_8TeV.root":"gg H",
#    "vbfHmumu125_8TeV.root":"VBF H"
#  }
#  bakFileNames = [
#    "SingleMuRun2012Av1.root",
#    "SingleMuRun2012Bv1.root",
#    "SingleMuRun2012Cv1.root",
#    "SingleMuRun2012Cv2.root"
#  ]
#
#  isoCountsTable = compareDirs(dirs,sigFileNames,bakFileNames,massBoundaries=mBounds)
#  isoTableFile = open("isoTable.tex",'w')
#  isoTableFile.write( r"""
#\documentclass[12pt,a4paper]{article}
#\usepackage{lscape}
#\begin{document}
#%\tiny
#\small
#\begin{center}
#
#"""
#  )
#  isoTableFile.write(isoCountsTable)
#  isoTableFile.write(     r"""
#
#\end{center}
#\end{document}
#"""
#  )
#  isoTableFile.close()
#
#
#  ######################################################
  signames = ["ggHmumu125_8TeV.root","vbfHmumu125_8TeV.root"]
  baknames = [
    "DYJetsToLL_8TeV.root",
    "ttbar_8TeV.root"
  ]
  datnames = [
    "SingleMuRun2012Av1.root",
    "SingleMuRun2012Av1Recover.root",
    "SingleMuRun2012Bv1.root",
    "SingleMuRun2012Cv1.root",
    "SingleMuRun2012Cv2.root",
    "SingleMuRun2012D.root"
  ]
  indir = "input/V00-01-10/"
  signames = [indir+i for i in signames]
  baknames = [indir+i for i in baknames]
  datnames = [indir+i for i in datnames]

#  cutTableFile = open("cutTable.tex",'w')
#  cutTableFile.write( r"""
#\documentclass[12pt,a4paper]{article}
#\usepackage{lscape}
#\begin{document}
#%\tiny
#\small
#\begin{center}
#
#"""
#  )
#  cutTable =  cutFlow(signames,baknames,"Inc")
#  cutTableFile.write(r" {\large Not-VBF Category} \\"+'\n')
#  cutTableFile.write(cutTable)
#  cutTable =  cutFlow(signames,baknames,"VBF")
#  cutTableFile.write(r"\vspace {4em}\\ {\large VBF Category} \\"+'\n')
#  cutTableFile.write(cutTable)
#  cutTableFile.write(     r"""
#
#\end{center}
#\end{document}
#"""
#  )
#  cutTableFile.close()
#
#  #####################################################

  inMassWindowTableFile = open("inMassWindowTable.tex",'w')
  inMassWindowTableFile.write( r"""
\documentclass[12pt,a4paper]{article}
\usepackage{lscape}
\begin{document}
%\tiny
\small
\begin{center}

"""
  )
  cats = [""]
  cuts = [""]

#  cats += ["Jets2"]
#  cuts += ["nJets>=2"]
#  cats += ["Jets2PtMissL100"]
#  cuts += ["nJets>=2 && ptMiss<100."]
#  cats += ["Jets2CutPass"]
#  cuts += ["nJets>=2 && (deltaEtaJets>3. && dijetMass>750.) && ptMiss<100."]
#  cats += ["Jets2CutFail"]
#  cuts += ["nJets>=2 && !(deltaEtaJets>3. && dijetMass>750.) && ptMiss<100."]

#  cats += ["Jets0"]
#  cuts += ["nJets==0"]
#  cats += ["Jets0Pass"]
#  cuts += ["nJets==0 && dimuonPt>10"]
#  cats += ["Jets0Fail"]
#  cuts += ["nJets==0 && !(dimuonPt>10)"]
#
#  cats += ["Jets1"]
#  cuts += ["nJets==1"]
#  cats += ["Jets1Pass"]
#  cuts += ["nJets==1 && dimuonPt>10"]
#  cats += ["Jets1Fail"]
#  cuts += ["nJets==1 && !(dimuonPt>10)"]

  #cats += ["JetsL2"]
  #cuts += ["nJets<2"]
  #cats += ["JetsL2BB"]
  #cuts += ["nJets<2"]
  #cats += ["JetsL2BO"]
  #cuts += ["nJets<2"]
  #cats += ["JetsL2BE"]
  #cuts += ["nJets<2"]
  #cats += ["JetsL2OO"]
  #cuts += ["nJets<2"]
  #cats += ["JetsL2OE"]
  #cuts += ["nJets<2"]
  #cats += ["JetsL2EE"]
  #cuts += ["nJets<2"]

  cats += ["Jets0"]
  cuts += ["nJets==0"]
  cats += ["Jets0BB"]
  cuts += ["nJets==0"]
  cats += ["Jets0BO"]
  cuts += ["nJets==0"]
  cats += ["Jets0BE"]
  cuts += ["nJets==0"]
  cats += ["Jets0OO"]
  cuts += ["nJets==0"]
  cats += ["Jets0OE"]
  cuts += ["nJets==0"]
  cats += ["Jets0EE"]
  cuts += ["nJets==0"]

  #cats += ["Jets1"]
  #cuts += ["nJets==1"]
  #cats += ["Jets1BB"]
  #cuts += ["nJets==1"]
  #cats += ["Jets1BO"]
  #cuts += ["nJets==1"]
  #cats += ["Jets1BE"]
  #cuts += ["nJets==1"]
  #cats += ["Jets1OO"]
  #cuts += ["nJets==1"]
  #cats += ["Jets1OE"]
  #cuts += ["nJets==1"]
  #cats += ["Jets1EE"]
  #cuts += ["nJets==1"]


  samples=['width','sob','dataWide']

  inMassWindowTable =  inMassWindowTableMkr(signames,baknames,datnames,cats,cuts,samples=samples)
  #inMassWindowTableFile.write(r" {\large Not-VBF BDT Cut} \\"+'\n')
  inMassWindowTableFile.write(inMassWindowTable)

  inMassWindowTableFile.write(     r"""

\end{center}
\end{document}
"""
  )
  inMassWindowTableFile.close()

