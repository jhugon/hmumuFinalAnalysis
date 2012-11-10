#!/usr/bin/env python

import math
from math import sqrt
import ROOT as root
from helpers import *
import matplotlib.pyplot as mpl
import numpy
import glob

from xsec import *

import makeShapePlots

channelNameMap = {
  "AllCat":"All Cat. Comb.",
  "IncCat":"Inc. Cat. Comb.",
  "VBFCat":"VBF Cat. Comb.",

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

  "BDTSig80":"BDT Cut Comb.",
  "IncBDTSig80":"Inc. BDT Cut",
  "VBFBDTSig80":"VBF BDT Cut"
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
  "Stat":"Background Stat.",
}

def convertHistToCounts(dataDict,xlow,xhigh):
  outDict = {}
  for channelName in dataDict:
    tmpDict = {}
    channel = dataDict[channelName]
    lowBin = 0
    highBin = 0
    for histName in channel:
      hist = channel[histName]
      val = getIntegralAll(hist,[xlow,xhigh])
      tmpDict[histName] = val
    outDict[channelName] = tmpDict
  return outDict
  
def getShapeErrorsFromCounts(data):
  outDict = {}
  for channelName in data:
    tmpDict = {}
    channel = data[channelName]
    for name in channel:
      if name == "data_obs":
        tmpDict["data_obs"] = {"nom":channel[name]}
        continue
      matchUp = re.match(r"(.+)_(.+)Up",name)
      matchDown = re.match(r"(.+)_(.+)Down",name)
      if matchUp:
        histName = matchUp.group(1)
        parName = matchUp.group(2)
        if tmpDict.has_key(histName):
          if tmpDict[histName].has_key(parName):
            tmpDict[histName][parName]["Up"] = channel[name]
          else:
            tmpDict[histName][parName] = {}
            tmpDict[histName][parName]["Up"] = channel[name]
        else:
          tmpDict[histName] = {}
          tmpDict[histName][parName] = {}
          tmpDict[histName][parName]["Up"] = channel[name]
      elif matchDown:
        histName = matchDown.group(1)
        parName = matchDown.group(2)
        if tmpDict.has_key(histName):
          if tmpDict[histName].has_key(parName):
            tmpDict[histName][parName]["Down"] = channel[name]
          else:
            tmpDict[histName][parName] = {}
            tmpDict[histName][parName]["Down"] = channel[name]
        else:
          tmpDict[histName] = {}
          tmpDict[histName][parName] = {}
          tmpDict[histName][parName]["Down"] = channel[name]
      else:
        if tmpDict.has_key(name):
          tmpDict[name]["nom"] = channel[name]
        else:
          tmpDict[name] = {"nom":channel[name]}
    ##################################
    # Now Process into Error Fractions
    for histName in tmpDict:
      nomVal = tmpDict[histName]["nom"]
      upSum2 = 0.0
      downSum2 = 0.0
      for errName in tmpDict[histName]:
        if errName == "nom":
            continue
        upVal = tmpDict[histName][errName]["Up"]
        downVal = tmpDict[histName][errName]["Down"]
        if upVal < downVal:
          tmpVal = upVal
          upVal = downVal
          downVal = tmpVal
        upErr = (upVal-nomVal)/nomVal
        downErr = (nomVal-downVal)/nomVal
        tmpDict[histName][errName]["Up"] = upErr
        tmpDict[histName][errName]["Down"] = downErr
        upSum2 += upErr**2
        downSum2 += downErr**2
      if histName == "bak":
        tmpDict[histName]["TotalSyst"] = {"Up":sqrt(upSum2),"Down":sqrt(downSum2)}
        tmpDict[histName]["Stat"] = {"Up":sqrt(nomVal)/nomVal,"Down":sqrt(nomVal)/nomVal}

    for histName in tmpDict:
        print(histName+"  "+str(tmpDict[histName].keys()))
        
    outDict[channelName] = tmpDict
      
      
  return outDict

def writeErrorTable(data,latex,niceTitles):
  outString = ""
  # Get Widths
  maxChannelWidth = 0
  maxSampleWidth = 0
  maxErrWidth = 0
  for channelName in data:
     channel = data[channelName]
     if niceTitles:
        channelName = channelNameMap[channelName]
     if len(channelName) > maxChannelWidth:
        maxChannelWidth = len(channelName)
     for sampleName in channel:
       sample = channel[sampleName]
       if niceTitles:
         sampleName = sampleNameMap[sampleName]
       if len(sampleName) > maxSampleWidth:
          maxSampleWidth = len(sampleName)
       for errName in sample:
         if errName == "nom":
            continue
         if niceTitles:
           errName = errNameMap[errName]
         if len(errName) > maxErrWidth:
            maxErrWidth = len(errName)
  maxChannelWidth = str(maxChannelWidth+2)
  maxSampleWidth = str(maxSampleWidth+2)
  maxErrWidth = str(maxErrWidth+2)
  # Get Err Names
  errNames = []
  errNamesString = " "*int(maxChannelWidth)
  if latex:
    errNamesString += "&"
  iErrName = 0
  for channelName in data:
     for sampleName in data[channelName]:
       for errName in data[channelName][sampleName]:
         if errName == "nom":
            continue
         if niceTitles:
           errName = errNameMap[errName]
         errNames.append(errName)
         errNamesString += "{"+str(iErrName)+":^"+str(len(errName)+2)+"}"
         if latex:
           errNamesString += "&"
         iErrName += 1
     break
  if latex:
    errNamesString = errNamesString.rstrip("&")
    errNamesString += r"\\ \hline \hline"
  errNamesString = errNamesString.format(*errNames)
  outString += "\n"+errNamesString+"\n"
  for channelName in data:
    errVals = [channelName]
    if niceTitles:
      errVals = [channelNameMap[channelName]]
    errVals2 = [""]
    errValsString = "{0:<"+maxChannelWidth+"}"
    errValsString2 = "{0:<"+maxChannelWidth+"}"
    if latex:
      errValsString += r"&"
      errValsString2 += r"&"
    iErrVal = 1
    for sampleName in data[channelName]:
      for errName in data[channelName][sampleName]:
        if errName == "nom":
            continue
        errVals.append("+{0:.3%}".format(data[channelName][sampleName][errName]["Up"]))
        errVals2.append("-{0:.3%}".format(data[channelName][sampleName][errName]["Down"]))
        if niceTitles:
            errName = errNameMap[errName]
        errValsString += "{"+str(iErrVal)+":^"+str(len(errName)+2)+"}"
        errValsString2 += "{"+str(iErrVal)+":^"+str(len(errName)+2)+"}"
        if latex:
           errValsString += "&"
           errValsString2 += "&"
        iErrVal += 1
    errValsString = errValsString.format(*errVals)
    errValsString2 = errValsString2.format(*errVals2)
    if latex:
      errValsString = errValsString.rstrip("&")
      errValsString2 = errValsString2.rstrip("&")
      errValsString += r"\\"
      errValsString2 += r"\\ \hline"
    errValsString += "\n"+errValsString2+"\n"
    outString += errValsString

  if latex:
    columnFormat = "|l|" + "c|"*len(errNames)
    outString = r"\begin{tabular}{"+columnFormat+"} \hline" + outString + r"\end{tabular}"+"\n"
    outString = outString.replace(r"%",r"\%")

  return outString
    
if __name__ == "__main__":
 import subprocess
 import os
  
 dataDir = "statsCards/"
 filenames = ["statsCards/BDTSig80_8TeV_20.root"]
 #filenames = glob.glob(dataDir+"20*.root")

 for fn in filenames:
  sPlotter = makeShapePlots.ShapePlotter(fn,makeShapePlots.titleMap)
    
  data = getShapeErrorsFromCounts(convertHistToCounts(sPlotter.data,123.,127.))

  f = open("tableShapeErrors.tex","w")
  f.write(writeErrorTable(data,True,True))
  f.close()


  f = open("tableShapeErrorsTest.tex","w")
  f.write(r"""
\documentclass[12pt,a4paper]{article}
\usepackage{lscape}
\begin{document}
\begin{landscape}
%\tiny
\small
\input{shapeErrors}
\end{landscape}
\end{document}
         """)
  f.close()
  subprocess.call(["latex","tableShapeErrorsTest.tex"])
  subprocess.call(["dvipdf","tableShapeErrorsTest.dvi"])
  os.remove("tableShapeErrorsTest.aux")
  os.remove("tableShapeErrorsTest.log")
  os.remove("tableShapeErrorsTest.dvi")

