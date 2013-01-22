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
  "IncCat":"!VBF Cat. Comb.",
  "VBFCat":"VBF Cat. Comb.",

  "Presel":"Presel. Comb.",
  "IncPresel":"!VBF Category",
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
  "IncBDTCut":"!VBF BDT",
  "VBFBDTCut":"VBF BDT",

  "BDTCutCat":"BDT Res. Comb.",
  "IncBDTCutCat":"!VBF BDT Res.",
  "VBFBDTCutCat":"VBF BDT Res.",

  "PreselCat":"Presel. Res. Comb.",
  "IncPreselCat":"!VBF Res. Presel.",
  "VBFPreselCat":"VBF Res. Presel.",

  "IncBDTCutBB":"!VBF BDT BB",
  "IncBDTCutBO":"!VBF BDT BO",
  "IncBDTCutBE":"!VBF BDT BE",
  "IncBDTCutOO":"!VBF BDT OO",
  "IncBDTCutOE":"!VBF BDT OE",
  "IncBDTCutEE":"!VBF BDT EE",
  "IncBDTCutNotBB":"!VBF BDT !BB",
  "VBFBDTCutBB":"VBF BDT BB",
  "VBFBDTCutNotBB":"VBF BDT !BB",
  "IncPreselBB":"!VBF Presel. BB",
  "IncPreselBO":"!VBF Presel. BO",
  "IncPreselBE":"!VBF Presel. BE",
  "IncPreselOO":"!VBF Presel. OO",
  "IncPreselOE":"!VBF Presel. OE",
  "IncPreselEE":"!VBF Presel. EE",
  "IncPreselNotBB":"!VBF Presel. !BB",
  "VBFPreselBB":"VBF Presel. BB",
  "VBFPreselNotBB":"VBF Presel. !BB"
}

sampleNameMap = {
  "data_obs":"Data",
  "bak":"Background",
  "bakWide":"Background in Sidebands",
  "sig":"Signal",
  "sob":"$S/B$",
  "sosqrtb":"$S/\sqrt{B}$",
  "sosqrtsb":"$S/\sqrt{S+B}$",
  "ggHmumu125_8TeV":r"$gg \rightarrow H$",
  "vbfHmumu125_8TeV":r"VBF H",
  "ggHmumu125_7TeV":r"$gg \rightarrow H$",
  "vbfHmumu125_7TeV":r"VBF H"
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
          scaleBy = xsec[fNameNoExt]/nEventsMap[fNameNoExt]*lumiDict[energyStr]*1000.
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
        if i == '':
          i = "all"
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

def inMassWindowTableMkr(sigFileNames,bakFileNames,channel,massBoundaries=[120.,130.]):
  data = {}
  data["bak"] = {}
  data["sig"] = {}
  data["sob"] = {}
  data["sosqrtb"] = {}
  data["sosqrtsb"] = {}

  energyStr = ""
  lumi = ""
  categories = [channel]
  sigCounts = Counts(sigFileNames,categories,massBoundaries)
  bakCounts = Counts(bakFileNames,categories,massBoundaries)
  for sig in sigCounts.data:
    data[sig] = {}
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

  countsWide = Counts(bakFileNames,categories,[110.,160.])
  nBakWide = 0.0
  for bak in bakCounts.data:
    nBakWide += countsWide.data[bak][channel]
  data["bakWide"] = {channel:nBakWide-data["bak"][channel]}
  print nBakWide

  samples = ["ggHmumu125_8TeV","vbfHmumu125_8TeV","bak","bakWide"]
  ncols = len(samples)

  outString = " &"
  for i in samples:
    outString += " "+sampleNameMap[i]+" &"
  outString = outString.rstrip(r"&")
  outString += r"\\ \hline" + '\n'

  outString += "Events &"
  for s in samples:
    n = data[s][channel]
    if s=="bak" or s=="sob" or s=="sosqrtb" or s=="sosqrtsb" or s=="bakWide" :
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
  #mBounds = [110.,160.]
  mBounds = [120.,130.]

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

  isoCountsTable = compareDirs(dirs,sigFileNames,bakFileNames,massBoundaries=mBounds)
  isoTableFile = open("isoTable.tex",'w')
  isoTableFile.write( r"""
\documentclass[12pt,a4paper]{article}
\usepackage{lscape}
\begin{document}
%\tiny
\small
\begin{center}

"""
  )
  isoTableFile.write(isoCountsTable)
  isoTableFile.write(     r"""

\end{center}
\end{document}
"""
  )
  isoTableFile.close()


  ######################################################
  signames = ["ggHmumu125_8TeV.root","vbfHmumu125_8TeV.root"]
  baknames = [
    "SingleMuRun2012Av1.root",
    "SingleMuRun2012Bv1.root",
    "SingleMuRun2012Cv1.root",
    "SingleMuRun2012Cv2.root"
  ]
  signames = ["input/notblind/"+i for i in signames]
  baknames = ["input/notblind/"+i for i in baknames]

  cutTableFile = open("cutTable.tex",'w')
  cutTableFile.write( r"""
\documentclass[12pt,a4paper]{article}
\usepackage{lscape}
\begin{document}
%\tiny
\small
\begin{center}

"""
  )
  cutTable =  cutFlow(signames,baknames,"Inc")
  cutTableFile.write(r" {\large Not-VBF Category} \\"+'\n')
  cutTableFile.write(cutTable)
  cutTable =  cutFlow(signames,baknames,"VBF")
  cutTableFile.write(r"\vspace {4em}\\ {\large VBF Category} \\"+'\n')
  cutTableFile.write(cutTable)
  cutTableFile.write(     r"""

\end{center}
\end{document}
"""
  )
  cutTableFile.close()

  #####################################################

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
  inMassWindowTable =  inMassWindowTableMkr(signames,baknames,"IncBDTCut")
  inMassWindowTableFile.write(r" {\large Not-VBF BDT Cut} \\"+'\n')
  inMassWindowTableFile.write(inMassWindowTable)
  inMassWindowTable =  inMassWindowTableMkr(signames,baknames,"VBFBDTCut")
  inMassWindowTableFile.write(r"\vspace {4em}\\ {\large VBF BDT Cut} \\"+'\n')
  inMassWindowTableFile.write(inMassWindowTable)
  inMassWindowTableFile.write(     r"""
\end{center}
\end{document}
"""
  )
  inMassWindowTableFile.close()

