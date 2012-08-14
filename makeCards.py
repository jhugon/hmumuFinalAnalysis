#!/usr/bin/python

import math
import ROOT as root
from helpers import *
import datetime

from xsec import *

class Analysis:
  def __init__(self,directory,signalName,backgroundNames,analysis,massRange=[123,127]):
    self.sigFile = root.TFile(directory+signalName+".root")
    self.name = analysis
    self.bakNames = backgroundNames
    self.sigHist = self.sigFile.Get("mDiMu"+analysis)
    self.sigName = signalName

    self.bakFiles = []
    self.bakHists = []
    for name in backgroundNames:
      tmpF = root.TFile(directory+name+".root")
      tmpH = tmpF.Get("mDiMu"+analysis)
      self.bakFiles.append(tmpF)
      self.bakHists.append(tmpH)

    effMap = {}
    xsecMap = {}
    axis = self.sigHist.GetXaxis()
    lowBin = axis.FindBin(massRange[0])
    assert(axis.GetBinLowEdge(lowBin)==massRange[0])
    highBin = axis.FindBin(massRange[1])
    highBin -= 1
    assert(axis.GetBinUpEdge(highBin)==massRange[1])

    countsSig = self.sigHist.Integral(lowBin,highBin)
    effSig = countsSig/nEventsMap[signalName]
    xsecSig = effSig*xsec[signalName]
    self.effSig = effSig
    self.xsecSig = xsecSig

    self.xsecBakTotal = 0.0
    self.xsecBakList = []
    self.effBakList = []
    for h,name in zip(self.bakHists,backgroundNames):
      countsBak = h.Integral(lowBin,highBin)
      effBak = countsBak/nEventsMap[name]
      xsecBak = effBak*xsec[name]
      self.xsecBakTotal += xsecBak
      self.xsecBakList.append(xsecBak)
      self.effBakList.append(effBak)

  def getSigEff(self):
    return self.effSig
  def getSigXSec(self):
    return self.xsecSig
  def getBakXSecTotal(self):
    return self.xsecBakTotal
  def getBakXSec(self,bakName):
    result = -1.0
    if self.bakNames.count(bakName)>0:
        i = self.bakNames.index(bakName)
        result = self.xsecBakList[i]
    return result

class DataCardMaker:
  def __init__(self,directory,signalAnalysisPairs,backgroundNames,massRange=[123,127],nuisanceMap=None):
  #def __init__(self,directory,signalName,backgroundNames,analysisList,massRange=[123,127]):
    channels = []
    self.channelNames = []
    for pair in signalAnalysisPairs:
      tmp = Analysis(directory,pair[0],backgroundNames,pair[1],massRange=massRange)
      channels.append(tmp)
      self.channelNames.append(pair[1])
    self.channels = channels

    if nuisanceMap == None:
      self.nuisance = {}
      self.nuisance["lumi"] = (0.044,["vbfHmumu125","ggHmumu125"])
      self.nuisance["xs_ggH"] = (0.147,["ggHmumu125"])
      self.nuisance["xs_vbfH"] = (0.03,["vbfHmumu125"])
      self.nuisance["br_Hmm"] = (0.06,["ggHmumu125","vbfHmumu125"])
      self.nuisance["bg_dy"] = (0.05,["DYJetsToLL"])
      self.nuisance["bg_tt"] = (0.05,["ttbar"])
    else:
      self.nuisance = nuisanceMap

    self.largestChannelName = 0
    for name in self.channelNames:
        if len(name)>self.largestChannelName:
          self.largestChannelName = len(name)
    for channel in channels:
      if len(channel.sigName)>self.largestChannelName:
          self.largestChannelName = len(channel.sigName)
      for name in channel.bakNames:
        if len(name)>self.largestChannelName:
          self.largestChannelName = len(name)
    if self.largestChannelName < 8:
        self.largestChannelName = 8
    self.largestChannelName += 2
    self.bakNames = backgroundNames

    if self.channelNames.count("")>0:
      i = self.channelNames.index("")
      self.channelNames[i] = "Inc"

  def write(self,outfilename,lumi):
    print("Writing Card: {0}".format(outfilename))
    lumi *= 1000.0
    nuisance = self.nuisance
    outfile = open(outfilename,"w")
    outfile.write("# Hmumu combine datacard produced by makeTables.py\n")
    now = datetime.datetime.now().replace(microsecond=0).isoformat(' ')
    outfile.write("# {0}\n".format(now))
    outfile.write("############################### \n")
    outfile.write("############################### \n")
    outfile.write("imax {0}\n".format(len(self.channels)))
    #outfile.write("jmax {0}\n".format(len(backgroundNames)))
    outfile.write("jmax {0}\n".format("*"))
    outfile.write("kmax {0}\n".format(len(nuisance)))
    outfile.write("------------\n")
    outfile.write("# Channels, observed N events:\n")
    # Make Channels String
    binFormatString = "bin           "
    observationFormatString = "observation  "
    binFormatList = self.channelNames
    observationFormatList = []
    iParam = 0
    for channel,channelName in zip(self.channels,self.channelNames):
      binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
      binFormatList.append(channelName)
      observationFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
      observationFormatList.append(int(channel.getBakXSecTotal()*lumi))
      iParam += 1
    binFormatString+= "\n"
    observationFormatString+= "\n"
    outfile.write(binFormatString.format(*binFormatList))
    outfile.write(observationFormatString.format(*observationFormatList))
    outfile.write("------------\n")
    outfile.write("# Expected N events:\n")

    binFormatString = "bin           "
    proc1FormatString = "process       "
    proc2FormatString = "process       "
    rateFormatString = "rate          "
    binFormatList = []
    proc1FormatList = []
    proc2FormatList = []
    rateFormatList = []
    iParam = 0
    for channel,channelName in zip(self.channels,self.channelNames):
        iProc = 0
        binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
        binFormatList.append(channelName)
  
        proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
        proc1FormatList.append(channel.sigName)
  
        proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
        proc2FormatList.append(iProc)
  
        expNum = channel.getSigXSec()*lumi
        decimals = ".4f"
        if expNum>1000.0:
            decimals = ".4e"
        rateFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+decimals+"} "
        rateFormatList.append(expNum)
        iParam += 1
        iProc += 1
  
        for bakName in self.bakNames:
          binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          binFormatList.append(channelName)
  
          proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc1FormatList.append(bakName)
  
          proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc2FormatList.append(iProc)
  
          expNum = channel.getBakXSec(bakName)*lumi
          decimals = ".4f"
          if expNum>1000.0:
            decimals = ".4e"
          rateFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+decimals+"} "
          rateFormatList.append(expNum)
  
          iParam += 1
          iProc += 1
    binFormatString+= "\n"
    proc1FormatString+= "\n"
    proc2FormatString+= "\n"
    rateFormatString+= "\n"
    outfile.write(binFormatString.format(*binFormatList))
    outfile.write(proc1FormatString.format(*proc1FormatList))
    outfile.write(proc2FormatString.format(*proc2FormatList))
    outfile.write(rateFormatString.format(*rateFormatList))
    outfile.write("------------\n")
    outfile.write("# Uncertainties:\n")

    for nu in nuisance:
      thisNu = nuisance[nu]
      formatString = "{0:<8} {1:^4} "
      formatList = [nu,"lnN"]
      iParam = 2
      for channel,channelName in zip(self.channels,self.channelNames):
          iProc = 2
          formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          value = "-"
          if thisNu[1].count(channel.sigName)>0:
            value = thisNu[0]+1.0
          formatList.append(value)
          iParam += 1
    
          for bakName in self.bakNames:
            formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            value = "-"
            if thisNu[1].count(bakName)>0:
              value = thisNu[0]+1.0
            formatList.append(value)
            iParam += 1
      formatString += "\n"
      #print formatString
      #print formatList
      outfile.write(formatString.format(*formatList))
    outfile.close()

if __name__ == "__main__":
  print "Started makeCards.py"

  directory = "input/open/"
  outDir = "statsCards/"
  analysisList = [
        ["ggHmumu125","PtL30"],
        ["ggHmumu125","Pt30to50"],
        ["ggHmumu125","Pt50to75"],
        ["ggHmumu125","Pt75to125"],
        ["ggHmumu125","Pt125"],
        ["vbfHmumu125","VBFL"],
        ["vbfHmumu125","VBFM"],
        ["vbfHmumu125","VBFT"]
  ]
  backgroundNames= ["DYJetsToLL","ttbar"]
  lumiList = [5,10,15,20,25,30,40,50,75,100,200,500,1000]

  dataCard = DataCardMaker(directory,analysisList,backgroundNames)
  for i in lumiList:
    dataCard.write(outDir+"combined_"+str(i)+".txt",i)

  for i in analysisList:
    title = i[1]
    if title=="":
        title="Inc"
    dataCard = DataCardMaker(directory,[i],backgroundNames)
    for j in lumiList:
      dataCard.write(outDir+title+"_"+str(j)+".txt",j)

  for i in range(0,20):
    amount = 0.25+i*0.25
    dataCard = DataCardMaker(directory,analysisList,backgroundNames,massRange=[125.0-amount,125.0+amount])
    title = "combine125PM"+str(amount)
    dataCard.write(outDir+title+"_"+str(20)+".txt",20)

  runFile = open(outDir+"run.sh","w")
  runFile.write("#!/bin/bash\n")
  runFile.write("\n")
  runFile.write("for i in *.txt; do\n")
  runFile.write('    [[ -e "$i" ]] || continue\n')
  runFile.write('echo "Running on "$i\n')
  runFile.write('combine -M Asymptotic $i > $i.out\n')
  runFile.write('done\n')
  runFile.close()
