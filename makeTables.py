#!/usr/bin/env python

import math
from math import sqrt
import ROOT as root
from helpers import *
import matplotlib.pyplot as mpl
import numpy
import glob
import re

from xsec import *

from ROOT import gSystem
gSystem.Load('libRooFit')

root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

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
errNameMap = {
  "expParam":"Exponential Param.",
  "voitSig":"Voigtian $\sigma$",
  "voitmZ":"Voigtian Mass",
  "mixParam":"V/E Ratio",
  "TotalSyst":"Total Systematic",
  "Stat":"Statistical",

  "lumi":"Luminosity",
  "br_Hmm":r"BR($H\rightarrow\mu\mu$)",
  "xs_ggH":r"$\sigma(gg\rightarrow H)$",
  "xs_vbfH":r"$\sigma(gg \rightarrow H)$",
  "xs_wH":r"$\sigma(WH)$",
  "xs_zH":r"$\sigma(ZH)$",

  "ResASig":r"Muon Smearing $\sigma$",
  "ResSig":r"Muon Smearing Mixing",
}

class TableMaker:
  def __init__(self,filename,outNameBase,titleMap,sampleNameMap,errNameMap,xRange=[120.,130.]):
    assert(len(xRange) ==2)
    self.textFileName = os.path.splitext(filename)[0]+".txt"
    self.processNameMap, self.params = getattr(self,"readCard")(self.textFileName)
    self.titleMap = titleMap
    self.sampleNameMap = sampleNameMap
    self.errNameMap = errNameMap
    self.lumi = -1
    self.lumiStr = ""
    self.energyStr = ""
    tmpMatch = re.search(r"([\w]*)_(.+)_([.0-9]+)\.root",filename)
    if tmpMatch:
      self.lumi = int(float(tmpMatch.group(3)))
      self.lumiStr = "$\\mathcal{{L}}$ = {0} fb$^{{-1}}$".format(self.lumi)
      self.energyStr = "$\\sqrt{s}$="+re.sub(r"TeV"," TeV",tmpMatch.group(2))

    rooXRange = root.RooFit.Range("tableRange")
    rooAllNormRange = root.RooFit.NormRange("all")
    #rooNormSet = root.RooFit.NormSet(rooAllNormRange)

    #print("==================================================")
    #print(filename)
    self.f = root.TFile(filename)
    dataDict2 = {}
    dataDict = {}
    for channelKey in self.f.GetListOfKeys():
      if channelKey.GetClassName() != "RooWorkspace":
        continue
      channelName = channelKey.GetName()
      dataDict[channelName] = {}
      dataDict2[channelName] = {}
      channelWS = channelKey.ReadObj()
      mMuMu = channelWS.var("mMuMu")
      mMuMu.setRange("tableRange",xRange[0],xRange[1])
      mMuMuArgSet = root.RooArgSet(mMuMu)
      rooNormSet = root.RooFit.NormSet(mMuMuArgSet)
      #mMuMu.Print("verbose")
      bakPDF = channelWS.pdf("bak")
      data_obs = channelWS.data("data_obs")

      nData = data_obs.sumEntries("mMuMu > {0} && mMuMu < {1}".format(*xRange))
      nDataTotal = data_obs.sumEntries("mMuMu > 0.0")
      bakInt = bakPDF.createIntegral(mMuMuArgSet,rooXRange,rooNormSet)
      nBak = bakInt.getVal()*nDataTotal

      pdfParams = bakPDF.getParameters(data_obs)
      itr = pdfParams.createIterator()

      nBakErrUpList = []
      nBakErrDownList = []
      for iParam in range(pdfParams.getSize()):
        param = itr.Next()
        paramName = param.GetName()
        if self.params.has_key(paramName):
          nominal,err = self.params[paramName]
          nominal = float(nominal)
          err = float(err)
          param.setVal(nominal+err)
          nBakErrUp = bakInt.getVal()*nDataTotal
          param.setVal(nominal-err)
          nBakErrDown = bakInt.getVal()*nDataTotal
          nBakErrUp -= nBak
          nBakErrDown -= nBak
          if nBakErrUp > 0.0:
            nBakErrUpList.append(nBakErrUp**2)
          else:
            nBakErrDownList.append(nBakErrUp**2)
          if nBakErrDown > 0.0:
            nBakErrUpList.append(nBakErrDown**2)
          else:
            nBakErrDownList.append(nBakErrDown**2)
          dataDict2[channelName][paramName] = {"nom":nominal,"err":err}
      nBakErrUp = 0.
      nBakErrDown = 0.
      #print [sqrt(x) for x in nBakErrUpList]
      #print [sqrt(x) for x in nBakErrDownList]
#      if len(nBakErrUpList)>0:
#        nBakErrUp = sqrt(reduce(lambda x, y: x+y,nBakErrUpList))
#      if len(nBakErrDownList)>0:
#        nBakErrDown = sqrt(reduce(lambda x, y: x+y,nBakErrDownList))
      if len(nBakErrUpList)>0:
        nBakErrUp = sqrt(max(nBakErrUpList))/nBak
      if len(nBakErrDownList)>0:
        nBakErrDown = sqrt(max(nBakErrDownList))/nBak
      #print("{0:<20} Data: {1:<10} Predict: {2:<10.1f} +{3:<10.2%} -{4:<10.2%}".format(channelName,nData,nBak,nBakErrUp,nBakErrDown))
      dataDict[channelName]["nData"] = nData
      dataDict[channelName]["nBak"] = nBak
      dataDict[channelName]["bakUpErr"] = nBakErrUp
      dataDict[channelName]["bakDownErr"] = nBakErrDown
#      #Templates Time
#      for processName in self.processNameMap[channelName]:
#        if processName == "bak":
#          continue
#        template = channelWS.data(processName+"_Template")
#        nTemplateTotal = data_obs.sumEntries("mMuMu > 0.0")
#        eff = float(nTemplateTotal)/nEventsMap[processName]
    self.dataDict = dataDict
    self.dataDict2 = dataDict2


  def readCard(self,fn):
    f = open(fn)
    foundBin = False
    binList = []
    processList = []
    rateList = []
    paramMap = {}
    for line in f:
      if re.search("^bin",line):
        if foundBin:
          m =  re.findall("[\s]+[\w]+",line)
          binList.extend([i for i in m])
        else:
          foundBin = True
      if re.search("^process[\s]+[a-zA-Z]+",line):
          m =  re.findall("[\s]+[\w]+",line)
          processList.extend([i for i in m])
      if re.search("^rate[\s]+[-+eE.0-9]+",line):
          m =  re.findall("[\s]+[-+eE.0-9]+",line)
          rateList.extend([float(i) for i in m])
      paramMatch = re.search(r"([a-zA-Z0-9_]+)[\s]+param[\s]+([-.+eE0-9]+)[\s]+([-.+eE0-9]+)",line)
      if paramMatch:
        gs = paramMatch.groups()
        paramMap[gs[0]] = [gs[1],gs[2]]
    binList = [x.replace(r" ","") for x in binList]
    processList = [x.replace(r" ","") for x in processList]
    result = {}
    for i in binList:
      if not result.has_key(i):
        result[i] = {}
    for b,p,r in zip(binList,processList,rateList):
      result[b][p] = r
    return result, paramMap

  def printBakPredict(self):
    dataDict = self.dataDict
    titleMap = self.titleMap
    outString = r"Category & $N_{Data}$ & $N_{Predicted}$ & Background Error \\ \hline \hline"+'\n'
    maxChannelNameLength = max([len(i) for i in dataDict])
    numberLength = 8
    channels = dataDict.keys()
    channels.sort()
    channels.reverse()
    for channelName in channels:
      errList = []
      errNamesString = ""
      iErrName = 0
      errList.append(titleMap[channelName])
      errNamesString += "{"+str(iErrName)+":<"+str(maxChannelNameLength)+"} & "
      iErrName +=1
      errList.append(dataDict[channelName]["nData"])
      errNamesString += "{"+str(iErrName)+":<"+str(numberLength)+".0f} & "
      iErrName +=1
      errList.append(dataDict[channelName]["nBak"])
      errNamesString += "{"+str(iErrName)+":<"+str(numberLength)+".1f} & "
      iErrName +=1
      errList.append(dataDict[channelName]["bakUpErr"])
      errNamesString += "+{"+str(iErrName)+":<"+str(numberLength)+".2%} "
      iErrName +=1
      errList.append(dataDict[channelName]["bakDownErr"])
      errNamesString += "-{"+str(iErrName)+":<"+str(numberLength)+".2%} "
      iErrName +=1
      errNamesString += r"\\ \hline"+"\n"
      outString += errNamesString.format(*errList)

    outString = r"\begin{tabular}{|l|c|c|c|c|} \hline"+"\n" + outString + r"\end{tabular}"+"\n"
    outString += r"\\ "+self.lumiStr+", "+self.energyStr
    outString = outString.replace(r"%",r"\%")
    print
    print outString

  def printBakErrors(self):
    dataDict = self.dataDict2
    titleMap = self.titleMap
    errNameMap = self.errNameMap
    outString = r"Category & "
    maxErrNameLengthList = []
    channelNames = dataDict.keys()
    channelNames.sort()
    for channelName in channelNames:
      channel = dataDict[channelName]
      errNames = channel.keys()
      errNames.sort()
      errNames.reverse()
      for errName in errNames:
        maxErrNameLengthList.append(errName)
      break
    maxErrNameLength = max([len(i) for i in maxErrNameLengthList])
    errNameList = []
    errNamesString = ""
    iErrName = 0
    for channelName in channelNames:
      channel = dataDict[channelName]
      errNames = channel.keys()
      errNames.sort()
      errNames.reverse()
      for errName in errNames:
        errNamesString += " {"+str(iErrName)+":<"+str(maxErrNameLength)+"} &"
        errNameList.append(errNameMap[re.sub(r".*_","",errName)])
        iErrName +=1
      break
    errNamesString = errNamesString.rstrip("&")
    outString += errNamesString.format(*errNameList)

    outString += r" \\ \hline \hline"+'\n'
    maxChannelNameLength = max([len(i) for i in dataDict])
    numberLength = 8
    for channelName in channelNames:
      errList = []
      errNamesString = ""
      iErrName = 0
      errList.append(titleMap[channelName])
      errNamesString += "{"+str(iErrName)+":<"+str(maxChannelNameLength)+"} & "
      iErrName += 1
      channel = dataDict[channelName]
      errNames = channel.keys()
      errNames.sort()
      errNames.reverse()
      for errName in errNames:
        errDict = channel[errName]
        errList.append(errDict["nom"])
        errList.append(errDict["err"])
        errNamesString += " {"+str(iErrName)+":<"+str(numberLength)+".3f} $\pm$ {"+str(iErrName+1)+":<"+str(numberLength)+".3f} &"
        iErrName += 2
      errNamesString = errNamesString.rstrip("&")
      errNamesString = errNamesString.format(*errList)
      errNamesString += r"\\ \hline"+"\n"
      outString += errNamesString

    
    outString = r"\begin{tabular}{|l|"+'c|'*len(maxErrNameLengthList)+"} \hline"+"\n" + outString + r"\end{tabular}"+"\n"
    outString += r"\\ "+self.lumiStr+", "+self.energyStr
    outString = outString.replace(r"%",r"\%")
    print
    print outString

    
if __name__ == "__main__":
 import subprocess
 import os
  
 ggFileName = "input/ggHmumu125_8TeV.root"
 vbfFileName = "input/vbfHmumu125_8TeV.root"
 dataDir = "statsCards/"
 #filenames = ["statsCards/BDTCut_8TeV_20.root"]
 #filenames = ["statsCards/BDTCutCat_8TeV_20.root"]
 filenames = glob.glob(dataDir+"*.root")

 for fn in filenames:
    tm = TableMaker(fn,"sillyTables",channelNameMap,sampleNameMap,errNameMap)
    print("===========================")
    print(fn)
    print
    #tm.printBakPredict()
    tm.printBakErrors()


