#!/usr/bin/env python

import math
import ROOT as root
from helpers import *
import datetime
import sys
import os.path
import copy
import multiprocessing
import time
myThread = multiprocessing.Process

from ROOT import gSystem
gSystem.Load('libRooFit')

#root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
#PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

class Param:
  def __init__(self,name,nominal,lowErr,highErr):
    self.name = name
    self.nominal = nominal
    self.lowErr = lowErr
    self.highErr = highErr
  def __str__(self):
    return "{0}: {1:.3g} +{2:.3g} -{3:.3g}".format(self.name,self.nominal,self.lowErr,self.highErr)
  def __repr__(self):
    return str(self)

class Analysis:
  def __init__(self,name,lumi,sigInject=0.0):
    massVeryLowRange = [80,95]
    massLowRange = [95,120]
    massHighRange = [130,160]
 
    rooParamList = []
    
    maxMass = massHighRange[1]
    minMass = massVeryLowRange[0]
    mMuMu = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)
    mMuMu.setRange("z",88,94)
    mMuMu.setRange("verylow",massVeryLowRange[0],massVeryLowRange[1])
    mMuMu.setRange("low",massLowRange[0],massLowRange[1])
    mMuMu.setRange("high",massHighRange[0],massHighRange[1])
    mMuMu.setRange("signal",massLowRange[1],massHighRange[0])

    #Background PDF
    bwLambda = root.RooRealVar("bwLambda","bwLambda",-1e-03,-1e-01,-1e-05)
    expMmumu = root.RooExponential("expMmumu","expMmumu",mMuMu,bwLambda)

    mean = root.RooRealVar("mean","mean",125.0,100.0,150.0)
    width = root.RooRealVar("width","width",8.0,2.0,20.0)
    gausMmumu = root.RooGaussian("gausMmumu","gausMmumu",mMuMu,mean,width)

    rooParamList.append(bwLambda)
    rooParamList.append(mean)
    rooParamList.append(width)

    N = 10
    dataset = expMmumu.generate(root.RooArgSet(mMuMu),N)
    dataHist = dataset.binnedClone()

    expMmumu.fitTo(dataHist)

    NSig = 2
    datasetSig = gausMmumu.generate(root.RooArgSet(mMuMu),NSig)
    dataHistSig = datasetSig.binnedClone()

    gausMmumu.fitTo(dataHistSig)

    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    N = int(N*lumi)
    print "N = "
    print N
    dataset = expMmumu.generate(root.RooArgSet(mMuMu),N)
    dataHist = dataset.binnedClone()

    dataHist.SetName("data_obs")
    gausMmumu.SetName("sig")
    expMmumu.SetName("bak")

    w = root.RooWorkspace(name)
    wimport = getattr(w,"import")
    #wimport(mMuMu)
    wimport(gausMmumu)
    wimport(expMmumu)
    wimport(dataHist)
    #for i in rooParamList:
    #    wimport(i)


    self.w = w
    self.paramList = paramList
    self.obsNorm = N
    self.bakNorm = N
    self.sigNorm = NSig

class DataCardMaker:
  def __init__(self):
    pass

  def write(self,outfilename,lumi,sigInject=0.0):
    self.channels = [Analysis("yay",lumi,sigInject)]
    self.channelNames = ["yay"]
    self.largestChannelName = 10

    print("Writing Card: {0}".format(outfilename))
    rootfileName = re.sub(r"\.txt",".root",outfilename)
    outfile = open(outfilename,"w")
    outfile.write("# Hmumu combine datacard produced by makeTables.py\n")
    now = datetime.datetime.now().replace(microsecond=0).isoformat(' ')
    outfile.write("# {0}\n".format(now))
    outfile.write("############################### \n")
    outfile.write("############################### \n")
    outfile.write("imax {0}\n".format(len(self.channels)))
    #outfile.write("jmax {0}\n".format(len(backgroundNames)))
    outfile.write("jmax {0}\n".format("*"))
    outfile.write("kmax {0}\n".format("*"))
    outfile.write("------------\n")
    outfile.write("shapes * * {0} $CHANNEL:$PROCESS\n".format(rootfileName))
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
      print("Writing Pretend Data Counts")
      counts = channel.obsNorm
      observationFormatList.append(int(counts))
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
        #for sigName in self.sigNames:
        if True:
          sigName = "sig"
          binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          binFormatList.append(channelName)
  
          proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc1FormatList.append(sigName)
  
          proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc2FormatList.append(iProc)
  
          expNum = channel.sigNorm*lumi
          decimals = ".4f"
          if expNum>1000.0:
            decimals = ".4e"
          rateFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+decimals+"} "
          rateFormatList.append(expNum)
  
          iParam += 1
          iProc += 1
        #for bakName in self.bakNames:
        if True:
          bakName = "bak"
          binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          binFormatList.append(channelName)
  
          proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc1FormatList.append(bakName)
  
          proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc2FormatList.append(iProc)
  
          expNum = channel.bakNorm*lumi
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

    """
    for nu in nuisance:
      thisNu = nuisance[nu]
      formatString = "{0:<8} {1:^4} "
      formatList = [nu,"lnN"]
      iParam = 2
      for channel,channelName in zip(self.channels,self.channelNames):
          for sigName in self.sigNames:
            formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            value = "-"
            if thisNu.has_key(sigName):
              value = thisNu[sigName]+1.0
            formatList.append(value)
            iParam += 1
          for bakName in self.bakNames:
            formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            value = "-"
            if thisNu.has_key(bakName):
              value = thisNu[bakName]+1.0
            formatList.append(value)
            iParam += 1
      formatString += "\n"
      #print formatString
      #print formatList
      outfile.write(formatString.format(*formatList))
    """
    # Bak Shape Uncertainties (All Correlated)
    for channel,channelName in zip(self.channels,self.channelNames):
      for nu in channel.paramList:
        nuisanceName = nu.name
        formatString = "{0:<8} {1:^5} {2:^5.2g} {3:^5.3g}"
        formatList = [nuisanceName,"param",nu.nominal,nu.highErr]
        formatString += "\n"
        #print formatString
        #print formatList
        outfile.write(formatString.format(*formatList))
      break


    outfile.write("#################################\n")
    outfile.close()

    rootfile = root.TFile(rootfileName,"RECREATE")
    for channel in self.channels:
      channel.w.Write()
    rootfile.Close()

if __name__ == "__main__":
  dc = DataCardMaker()
  dc.write("testCard.txt",20.0)
