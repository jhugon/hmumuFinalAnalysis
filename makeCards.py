#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description="Makes cards for use in the CMS Combine tool.")
parser.add_argument("--signalInject", help="Inject Signal with Strength into data_obs",type=float,default=0.0)
parser.add_argument("--toyData", help="Make Toy Data from PDFs for data_obs",action="store_true",default=False)
args = parser.parse_args()

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

root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

NPROCS = 1

BAKUNC = 0.1

BAKUNCON = True
SIGUNCON = False

SIGGAUS = True

from xsec import *

if scaleHiggsBy != 1.0:
  print("Error: higgs xsec is scaled!!! Return to 1. Exiting.")
  sys.exit(1)
def vetoOutOfBoundsEvents(hist,boundaries=[]):
  xbinLow = None
  xbinHigh = None
  if len(boundaries)==2:
    xbinLow, xbinHigh = getXbinsHighLow(hist,boundaries[0],boundaries[1])
  else:
    print("Error: vetoOutOfBoundsEvents: boundaries must be length 2, exiting.")
    sys.exit(1)
  for i in range(0,xbinLow):
    hist.SetBinContent(i,0.0)
    hist.SetBinError(i,0.0)
  for i in range(xbinHigh+1,hist.GetNbinsX()+2):
    hist.SetBinContent(i,0.0)
    hist.SetBinError(i,0.0)

def getRooVars(directory,signalNames,histNameBase,analysis):
    hist = None
    is2D = False
    for name in signalNames:
      filename = directory+name+".root"
      histName = histNameBase+analysis
      #print("file name: {0}".format(filename))
      #print("hist name: {0}".format(histName))
      tmpF = root.TFile(filename)
      hist = tmpF.Get(histName)
      break
    if hist.InheritsFrom("TH2"):
      is2D = True

    x = root.RooRealVar('mMuMu','mMuMu',
                    hist.GetXaxis().GetXmin(),
                    hist.GetXaxis().GetXmax()
                    )
    if is2D:
      y = root.RooRealVar('mva','mva',
                    hist.GetYaxis().GetXmin(),
                    hist.GetYaxis().GetXmax()
                    )
    if is2D:
      return [x,y]
    else:
      return [x]

###################################################################################

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
  def getErrString(self):
    if self.lowErr == self.highErr:
        return "{0:.5g}".format(self.lowErr)
    else:
        return "-{0:.5g}/+{1:.5g}".format(self.lowErr,self.highErr)

def makePDFBak(name,hist,mMuMu,minMass,maxMass,workspaceImportFn):
    debug = ""
    debug += "### makePDFBak: "+name+"\n"

    channelName = name

    voitWidth = root.RooRealVar(channelName+"_voitWidth","voitWidth",2.4952)
    voitmZ = root.RooRealVar(channelName+"_voitmZ","voitmZ",85,95)
    voitSig = root.RooRealVar(channelName+"_voitSig","voitSig",0.0,30.0)
    voitMmumu = root.RooVoigtian("bak_voitMmumu","voitMmumu",mMuMu,voitmZ,voitWidth,voitSig)

    expParam = root.RooRealVar(channelName+"_expParam","expParam",-1,0)
    expMmumu = root.RooExponential("bak_expMmumu","expMmumu",mMuMu,expParam)

    mixParam = root.RooRealVar(channelName+"_mixParam","mixParam",0,1)

    pdfMmumu = root.RooAddPdf("bak","bak",root.RooArgList(voitMmumu,expMmumu),root.RooArgList(mixParam))
    
    mMuMuRooDataHist = root.RooDataHist("bak_Template","bak_Template",root.RooArgList(mMuMu),hist)
    hist.Print()
    mMuMuRooDataHist.Print()

    voitMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("z"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    voitmZ.setConstant(True)
    voitSig.setConstant(True)

    expMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("high"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    expParam.setConstant(True)
    
    pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    chi2 = pdfMmumu.createChi2(mMuMuRooDataHist)

    ## Error time

    rooParamList = [voitmZ,voitSig,expParam,mixParam]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    workspaceImportFn(pdfMmumu)
    workspaceImportFn(mMuMuRooDataHist)

    #mMuMuRooDataHist2 = mMuMuRooDataHist.reduce(root.RooFit.CutRange("low,signal,high"))
    #mMuMuRooDataHist2.SetName("bak_TemplateNoVeryLow")
    #workspaceImportFn(mMuMuRooDataHist2)

    return paramList

def makePDFSig(name,hist,mMuMu,minMass,maxMass,workspaceImportFn,channelName):

    debug = ""
    debug += "### makePDFSig: "+channelName+": "+name+"\n"

    mean = root.RooRealVar(channelName+"_"+name+"_Mean",channelName+"_"+name+"_Mean",125.,100.,150.)
    width = root.RooRealVar(channelName+"_"+name+"_Width",channelName+"_"+name+"_Width",5.0,0.5,20.0)
    pdfMmumu = root.RooGaussian(name,name,mMuMu,mean,width)
    
    mMuMuRooDataHist = root.RooDataHist(name+"_Template",channelName+"_"+name+"_Template",root.RooArgList(mMuMu),hist)

    pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.SumW2Error(False),PRINTLEVEL)

    ## Error time

    rooParamList = [mean,width]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]
    for i in rooParamList:
       i.setConstant(True)

    workspaceImportFn(pdfMmumu)
    workspaceImportFn(mMuMuRooDataHist)

###################################################################################

class Analysis:
  def __init__(self,directory,signalNames,backgroundNames,dataNames,analysis,lumi,controlRegionVeryLow,controlRegionLow,controlRegionHigh,histNameBase="mDiMu",rebin=[],histNameSuffix="",toyData=False,sigInject=0.0):
    self.sigNames = signalNames
    self.bakNames = backgroundNames
    self.datNames = dataNames
    self.controlRegionVeryLow = controlRegionVeryLow
    self.controlRegionLow = controlRegionLow
    self.controlRegionHigh = controlRegionHigh
    self.analysis = analysis
    self.params = []

    self.workspace = root.RooWorkspace(analysis)
    wImport = getattr(self.workspace,"import")

    maxMass = controlRegionHigh[1]
    minMass = controlRegionVeryLow[0]
    mMuMu = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)
    minMass = controlRegionLow[0]
    mMuMu.setRange("z",88,94)
    mMuMu.setRange("verylow",controlRegionVeryLow[0],controlRegionVeryLow[1])
    mMuMu.setRange("low",controlRegionLow[0],controlRegionLow[1])
    mMuMu.setRange("high",controlRegionHigh[0],controlRegionHigh[1])
    mMuMu.setRange("signal",controlRegionLow[1],controlRegionHigh[0])

    self.sigFiles = []
    self.sigHistsRaw = []
    for name in signalNames:
      tmpF = root.TFile(directory+name+".root")
      tmpH = tmpF.Get(histNameBase+analysis+histNameSuffix)
      self.sigFiles.append(tmpF)
      self.sigHistsRaw.append(tmpH)

    self.bakFiles = []
    self.bakHistsRaw = []
    for name in backgroundNames:
      tmpF = root.TFile(directory+name+".root")
      tmpH = tmpF.Get(histNameBase+analysis+histNameSuffix)
      self.bakFiles.append(tmpF)
      self.bakHistsRaw.append(tmpH)

    self.datFiles = []
    self.datHists = []
    for name in dataNames:
      tmpF = root.TFile(directory+name+".root")
      tmpH = tmpF.Get(histNameBase+analysis+histNameSuffix)
      self.datFiles.append(tmpF)
      self.datHists.append(tmpH)

    #Rebin
    rb = rebin
    if type(rb) != list:
      print("Error: Analysis.rebin: argument must be a list!!  Exiting.")
      sys.exit(1)
    if len(rb) == 2 and False:
        for hist in self.sigHistsRaw:
          hist.Rebin2D(*rb)
        for hist in self.bakHistsRaw:
          hist.Rebin2D(*rb)
        for hist in self.datHists:
          hist.Rebin2D(*rb)
    elif len(rb) == 1 and True:
        for hist in self.sigHistsRaw:
          hist.Rebin(*rb)
        for hist in self.bakHistsRaw:
          hist.Rebin(*rb)
        for hist in self.datHists:
          hist.Rebin(*rb)
    elif len(rb) == 0:
      pass
    else:
      print("Error: Analysis.rebin: argument must be len 0, 1, or 2 list!!  Exiting.")
      print("  Must also be same length as dimension of hist, if not 0.")
      sys.exit(1)

    effMap = {}
    xsecMap = {}
    lowBin = 0
    highBin = self.sigHistsRaw[0].GetNbinsX()+1
    #massBounds = [controlRegionLow[0],controlRegionHigh[1]]
    #massBounds = [controlRegionVeryLow[0],controlRegionHigh[1]]
    massBounds = [controlRegionVeryLow[0],controlRegionHigh[1]]
    self.massBounds = massBounds

    self.xsecBakTotal = 0.0
    self.xsecBakList = []
    self.effBakList = []
    self.bakHists = []
    self.bakHistTotal = None
    for h,name in zip(self.bakHistsRaw,backgroundNames):
      counts = getIntegralAll(h,boundaries=massBounds)
      eff = counts/nEventsMap[name]*efficiencyMap[getPeriod(name)]
      xs = eff*xsec[name]
      self.xsecBakTotal += xs
      self.xsecBakList.append(xs)
      self.effBakList.append(eff)
      h.Scale(xsec[name]/nEventsMap[name]*efficiencyMap[getPeriod(name)])
      self.bakHists.append(h)
      if self.bakHistTotal == None:
        self.bakHistTotal = h.Clone("bak")
      else:
        self.bakHistTotal.Add(h)

    self.bakHistTotal.Scale(lumi)
    self.bakHistTotalReal = self.bakHistTotal.Clone("data_obs")

    self.dataCountsTotal = None
    self.datHistTotal = None
    for h,name in zip(self.datHists,dataNames):
      counts = getIntegralAll(h,boundaries=massBounds)
      if self.dataCountsTotal == None:
        self.dataCountsTotal = counts
      else:
        self.dataCountsTotal += counts
      if self.datHistTotal == None:
        self.datHistTotal = h.Clone("data_obs")
      else:
        self.datHistTotal.Add(h)

    histToUseForBak = self.datHistTotal
    if histToUseForBak is None:
        histToUseForBak = self.bakHistTotal

    self.params.extend(makePDFBak(analysis,histToUseForBak,
                                mMuMu, minMass,maxMass,wImport
                                 )
                       )

    for name, hist in zip(signalNames,self.sigHistsRaw):
        makePDFSig(name,hist,mMuMu,minMass,maxMass,wImport,analysis)

    self.xsecSigTotal = 0.0
    self.xsecSigList = []
    self.effSigList = []
    self.sigHists = []
    for h,name in zip(self.sigHistsRaw,signalNames):
      counts = getIntegralAll(h,boundaries=massBounds)
      eff = counts/nEventsMap[name]*efficiencyMap[getPeriod(name)]
      xs = eff*xsec[name]
      self.xsecSigTotal += xs
      self.xsecSigList.append(xs)
      self.effSigList.append(eff)

    self.countsSigTotal = self.xsecSigTotal*lumi
    self.countsBakTotal = self.xsecBakTotal*lumi
    self.countsSigList = [x*lumi for x in self.xsecSigList]
    self.countsBakList = [x*lumi for x in self.xsecBakList]

    if toyData:
      bakPDF = self.workspace.pdf("bak")
      sigPDFList = [self.workspace.pdf(i) for i in signalNames]
      self.dataCountsTotal = int(self.countsBakTotal)
      toyDataset = bakPDF.generate(root.RooArgSet(mMuMu),self.dataCountsTotal)
      if sigInject>0.0:
        for name,counts in zip(signalNames,self.countsSigList):
          counts = int(sigInect*counts)
          if counts < 1:
            continue
          tmpSigPDF = self.workspace.pdf(name)
          tmpSigDataset = tmpSigPDF.generate(root.RooArgSet(mMuMu),counts)
          toyDataset.append(tmpSigdataset)
          self.dataCountsTotal += counts

      toyDataHist = toyDataSet.binnedClone("data_obs","Toy Data")
      wImport(toyDataHist)
    elif self.dataCountsTotal is None:
      self.dataCountsTotal = int(self.countsBakTotal)
      bakData = self.workspace.data("bak_Template")
      obsData = bakData.Clone("data_obs")
      obsData.SetTitle("MC Full-Sim Data")
      wImport(obsData)
      self.bakHistTotal.SetName("data_obs_"+analysis)
      degubf = root.TFile("debug.root","recreate")
      self.bakHistTotal.Write()
      degubf.Close()
    else:
      realDataHist = root.RooDataHist("data_obs","Real Observed Data",root.RooArgList(mMuMu),self.datHistTotal)
      wImport(realDataHist)
      #realDataHistNotVeryLow = realDataHist.reduce(root.RooFit.CutRange("low,signal,high"))
      #wImport(realDataHistNotVeryLow)

  def getSigEff(self,name):
    result = -1.0
    if self.sigNames.count(name)>0:
        i = self.sigNames.index(name)
        result = self.effSigList[i]
    return result
  def getSigXSec(self,name):
    result = -1.0
    if self.sigNames.count(name)>0:
        i = self.sigNames.index(name)
        result = self.xsecSigList[i]
    return result
  def getSigXSecTotal(self):
    return self.xsecSigTotal
  def getBakXSecTotal(self):
    return self.xsecBakTotal
  def getBakXSec(self,bakName):
    result = -1.0
    if self.bakNames.count(bakName)>0:
        i = self.bakNames.index(bakName)
        result = self.xsecBakList[i]
    return result
  def getDataTotal(self):
    return self.dataCountsTotal
  def getParamList(self):
    return self.params
  def getSigCounts(self,name):
    result = -1.0
    if self.sigNames.count(name)>0:
        i = self.sigNames.index(name)
        result = self.countsSigList[i]
    return result
  def getBakCounts(self,bakName):
    result = -1.0
    if self.bakNames.count(bakName)>0:
        i = self.bakNames.index(bakName)
        result = self.countsBakList[i]
    return result
  def getSigCountsTotal(self):
    return self.countsSigTotal
  def getBakCountsTotal(self):
    return self.countsBakTotal

###################################################################################

class DataCardMaker:
  def __init__(self,directory,analysisNames,signalNames,backgroundNames,dataNames,outfilename,lumi,nuisanceMap=None,histNameBase="",controlRegionLow=[110.,115],controlRegionHigh=[135,150],controlRegionVeryLow=[80.,110.],rebin=[],histNameSuffix="",sigInject=0.0,toyData=False):

    ########################
    ## Setup

    channels = []
    self.channelNames = copy.deepcopy(analysisNames)
    self.is2D = False

    lumi *= 1000.0

    for analysis in analysisNames:
      tmp = Analysis(directory,signalNames,backgroundNames,dataNames,analysis,lumi,controlRegionVeryLow,controlRegionLow,controlRegionHigh,histNameBase=histNameBase,rebin=rebin,histNameSuffix=histNameSuffix,toyData=toyData,sigInject=sigInject)
      channels.append(tmp)
    self.channels = channels

    if nuisanceMap == None:
      self.nuisance = {}
    else:
      self.nuisance = nuisanceMap
    nuisance = self.nuisance


    self.largestChannelName = 0
    for name in self.channelNames:
        if len(name)>self.largestChannelName:
          self.largestChannelName = len(name)
    for channel in channels:
      for name in channel.sigNames:
        if len(name)>self.largestChannelName:
          self.largestChannelName = len(name)
      for name in channel.bakNames:
        if len(name)>self.largestChannelName:
          self.largestChannelName = len(name)
    if self.largestChannelName < 8:
        self.largestChannelName = 8
    self.largestChannelName += 2
    self.sigNames = signalNames
    self.bakNames = backgroundNames

    if self.channelNames.count("")>0:
      i = self.channelNames.index("")
      self.channelNames[i] = "Inc"

    self.controlRegionHigh = controlRegionHigh
    self.controlRegionLow = controlRegionLow
    self.controlRegionVeryLow = controlRegionVeryLow
    self.toyData = toyData

    ### ROOT Part
    ##########################################################
    outRootFilename = re.sub(r"\.txt",r".root",outfilename)
    outRootFile = root.TFile(outRootFilename, "RECREATE")
    outRootFile.cd()

    rootDebugString = ""

    for channel in self.channels:
        channel.workspace.Write()

    outRootFile.Close()

    ### Text Part
    ##########################################################

    print("Writing Card: {0} & {1}".format(outfilename,outRootFilename))
    outfile = open(outfilename,"w")
    outfile.write("# Hmumu shape combine datacard produced by makeCards.py\n")
    now = datetime.datetime.now().replace(microsecond=0).isoformat(' ')
    outfile.write("# {0}\n".format(now))
    outfile.write("############################### \n")
    outfile.write("############################### \n")
    outfile.write("imax {0}\n".format(len(self.channels)))
    #outfile.write("jmax {0}\n".format(len(backgroundNames)))
    outfile.write("jmax {0}\n".format("*"))
    outfile.write("kmax {0}\n".format("*"))
    outfile.write("------------\n")
    outfile.write("shapes * * {0} $CHANNEL:$PROCESS\n".format( os.path.basename(outRootFilename)))
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
      observedNumber = channel.getDataTotal()
      observationFormatList.append(observedNumber)
      #print("text Observed {0}: {1}".format(channelName,observedNumber))
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
        ## Signal Time

        iProc = -len(channel.sigNames)+1
        for sigName in self.sigNames:
          binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          binFormatList.append(channelName)
  
          proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc1FormatList.append(sigName)
  
          proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc2FormatList.append(iProc)
  
          expNum = channel.getSigCounts(sigName)
          decimals = ".4f"
          if expNum>1000.0:
            decimals = ".4e"
          rateFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+decimals+"} "
          rateFormatList.append(expNum)
  
          iParam += 1
          iProc += 1

        ## Background Time

        binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
        binFormatList.append(channelName)
  
        proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
        proc1FormatList.append("bak")
  
        proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
        proc2FormatList.append(iProc)

        expNum = channel.getBakCountsTotal()
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

    # lnN Uncertainties
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
          if True:
              bakName="bak"
              formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
              value = "-"
              if thisNu.has_key(bakName):
                value = thisNu[bakName]+1.0
              formatList.append(value)
              iParam += 1
          else:
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

    # Parameter Uncertainties
    for channel,channelName in zip(self.channels,self.channelNames):
      for nu in channel.params:
        nuisanceName = nu.name
        formatString = "{0:<25} {1:<6} {2:<10.5g} {3:<10}"
        formatList = [nuisanceName,"param",nu.nominal,nu.getErrString()]
        formatString += "\n"
        #print formatString
        #print formatList
        outfile.write(formatString.format(*formatList))

    #Debugging
    outfile.write("#################################\n")
    for channel,channelName in zip(self.channels,self.channelNames):
        outfile.write("#\n")
        outfile.write("#info: channel {0}: \n".format(channelName))
    outfile.write(rootDebugString)
    outfile.close()

class ThreadedCardMaker(myThread):
  def __init__(self,*args,**dictArgs):
    myThread.__init__(self)
    self.args = args
    self.dictArgs = dictArgs
    self.started = False
  def run(self):
    #try:
      self.started = True
      dataCardMassShape = None
      dataCardMassShape = DataCardMaker(*(self.args),**(self.dictArgs))
    #except Exception as e:
    #  print("Error: Exception: {0}".format(e))
    #  print("  Thread Arguments: {0}, {1}".format(self.args,self.dictArgs))

###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################

if __name__ == "__main__":
  print "Started makeCards.py"
  root.gROOT.SetBatch(True)

  directory = "input/"
  outDir = "statsCards/"
  periods = ["7TeV","8TeV"]
  periods = ["8TeV"]
  analysesInc = ["IncPresel","IncBDTCut"]
  analysesVBF = ["VBFPresel","VBFBDTCut"]
  analysesInc = ["IncPresel"]
  analysesVBF = ["VBFPresel"]
  analyses = analysesInc + analysesVBF
  analyses = ["IncPresel"]
  categoriesInc = ["BB","BO","BE","OO","OE","EE"]
  categoriesVBF = ["BB","NotBB"]
  tmpList = []
  for a in analysesInc:
    for c in categoriesInc:
        tmpList.append(a+c)
#  analyses += tmpList
  tmpList = []
  for a in analysesVBF:
    for c in categoriesVBF:
        tmpList.append(a+c)
#  analyses += tmpList
  combinations = []
  combinationsLong = []
#  combinations.append((
#        ["IncBDTCut"+x for x in categoriesInc],"IncBDTCutCat"
#  ))
#  combinations.append((
#        ["VBFBDTCut"+x for x in categoriesVBF],"VBFBDTCutCat"
#  ))
#  combinations.append((
#        ["IncPresel"+x for x in categoriesInc],"IncPreselCat"
#  ))
#  combinations.append((
#        ["VBFPresel"+x for x in categoriesVBF],"VBFPreselCat"
#  ))
#  combinations.append((
#        ["IncBDTCut","VBFBDTCut"],"BDTCut"
#  ))
#  combinations.append((
#        ["IncPresel","VBFPresel"],"Presel"
#  ))
#  combinations.append((
#        ["VBFPresel"+x for x in categoriesVBF]+["IncPresel"+x for x in categoriesInc],"PreselCat"
#  ))
#  combinations.append((
#        ["VBFBDTCut"+x for x in categoriesVBF]+["IncBDTCut"+x for x in categoriesInc],"BDTCutCat"
#  ))
#  combinationsLong.append((
#        ["IncBDTCut","VBFBDTCut"],"BDTCut"
#  ))
#  combinationsLong.append((
#        ["VBFBDTCut"+x for x in categoriesVBF]+["IncBDTCut"+x for x in categoriesInc],"BDTCutCat"
#  ))
#  combinationsLong.append((
#        ["VBFPresel"+x for x in categoriesVBF]+["IncPresel"+x for x in categoriesInc],"PreselCat"
#  ))
  histPostFix="/mDiMu"
  #analyses = ["mDiMu"]
  #histPostFix=""
  signalNames=["ggHmumu125","vbfHmumu125","wHmumu125","zHmumu125"]
  backgroundNames= ["DYJetsToLL","ttbar"]
  dataDict = {}
  dataDict["8TeV"] = [
    #"SingleMuRun2012Av1",
    #"SingleMuRun2012Bv1",
    #"SingleMuRun2012Cv1",
    #"SingleMuRun2012Cv2"
  ]
  dataDict["7TeV"] = [
    #"SingleMuRun2011Av1",
    #"SingleMuRun2011Bv1"
  ]
  dataDict["14TeV"] = []
  lumiListLong = [5,10,15,20,25,30,40,50,75,100,200,500,1000,2000,5000]
  lumiListLong = [20,30,50,100,500,1000,5000]
  lumiList = [lumiDict["8TeV"],20,25,30]
  lumiList = [20]
  #lumiListLong = lumiList

  MassRebin = 1 # 4 Bins per GeV originally
  controlRegionVeryLow=[60,110]
  controlRegionLow=[110,120]
  controlRegionHigh=[130,160]

  shape=True
  toyData=args.toyData

  print("Simple Analyses to run:")
  for a in analyses:
    print("  {0}".format(a))
  print("Combination Analyses to run:")
  for c in combinations:
    print("  {0}".format(c[1]))
    for a in c[0]:
      print("    {0}".format(a))

  print("Creating Threads...")
  threads = []
  for p in periods:
    for i in lumiList:
      if p == "7TeV":
        i = lumiDict[p]
      for ana in analyses:
        tmp = ThreadedCardMaker(
          #__init__ args:
          directory,[ana],
          appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),dataDict[p],
          rebin=[MassRebin],
          controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix=histPostFix,
          controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,
          #write args:
          outfilename=outDir+ana+"_"+p+"_"+str(i)+".txt",lumi=i
          )
        threads.append(tmp)
      for comb in combinations:
       threads.append(
        ThreadedCardMaker(
          #__init__ args:
          directory,
          comb[0],
          appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),dataDict[p],
          rebin=[MassRebin],
          controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix=histPostFix,
          controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,
          #write args:
          outfilename=outDir+comb[1]+"_"+p+"_"+str(i)+".txt",lumi=i
        )
       )
      if p == "7TeV":
        break

  for i in lumiListLong:
    p = "14TeV"
    for comb in combinationsLong:
       threads.append(
        ThreadedCardMaker(
          #__init__ args:
          directory,
          comb[0],
          appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),dataDict[p],
          rebin=[MassRebin],
          controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix=histPostFix,
          controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,
          #write args:
          outfilename=outDir+comb[1]+"_"+p+"_"+str(i)+".txt",lumi=i
        )
       )

  nThreads = len(threads)
  print("nProcs: {0}".format(NPROCS))
  print("nCards: {0}".format(nThreads))

  threadsNotStarted = copy.copy(threads)
  threadsRunning = []
  threadsDone = []
  while True:
    iThread = 0
    while iThread < len(threadsRunning):
        alive = threadsRunning[iThread].is_alive()
        if not alive:
          tmp = threadsRunning.pop(iThread)
          threadsDone.append(tmp)
        else:
          iThread += 1

    nRunning = len(threadsRunning)
    if nRunning < NPROCS and len(threadsNotStarted) > 0:
        tmp = threadsNotStarted.pop()
        tmp.start()
        threadsRunning.append(tmp)

    nRunning = len(threadsRunning)
    nNotStarted = len(threadsNotStarted)
    if nRunning == 0 and nNotStarted == 0:
        break

    time.sleep(0.1)
      

  runFile = open(outDir+"run.sh","w")
  batchString = \
"""#!/bin/bash

chmod +x lxbatch.sh

for i in *.txt; do
    [[ -e "$i" ]] || continue
echo "Running on "$i
bsub lxbatch.sh $i
#bsub -q 1nh lxbatch.sh $i
done
"""
  runFile.write(batchString)
  runFile.close()

  runFile = open(outDir+"lxbatch.sh","w")
  batchString = \
"""#!/bin/bash
echo "Sourcing cmsset_default.sh"
cd /afs/cern.ch/cms/sw
source cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc462
echo "SCRAM_ARCH is $SCRAM_ARCH"
cd $LS_SUBCWD
echo "In Directory: "
pwd
eval `scramv1 runtime -sh`
echo "cmsenv success!"
date

TXTSUFFIX=".txt"
FILENAME=$1
DIRNAME="Dir"$1"Dir"
ROOTFILENAME=${1%$TXTSUFFIX}.root

mkdir $DIRNAME
cp $FILENAME $DIRNAME/
cp $ROOTFILENAME $DIRNAME/
cd $DIRNAME

echo "executing combine -M Asymptotic $FILENAME >& $FILENAME.out"

combine -M Asymptotic $FILENAME >& $FILENAME.out

cp $FILENAME.out ..


echo "done"
date
"""
  runFile.write(batchString)
  runFile.close()

  runFile = open(outDir+"notlxbatch.sh","w")
  batchString = \
"""#!/bin/bash
echo "running notlxbatch.sh"
date
for i in *.txt; do
    [[ -e "$i" ]] || continue
FILENAME=$i
echo "executing combine -M Asymptotic $FILENAME >& $FILENAME.out"

combine -M Asymptotic $FILENAME >& $FILENAME.out
rm -f roostats*
rm -f higgsCombineTest*.root

echo "executing combine -M ProfileLikelihood -d $FILENAME --signif >& $FILENAME.sig"

combine -M ProfileLikelihood -d $FILENAME --signif >& $FILENAME.sig
rm -f roostats*
rm -f higgsCombineTest*.root

echo "executing combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig"

combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 >& $FILENAME.expsig
#combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig
rm -f roostats*
rm -f higgsCombineTest*.root

echo "executing combine -M MaxLikelihoodFit $FILENAME >& $FILENAME.mu"

combine -M MaxLikelihoodFit $FILENAME >& $FILENAME.mu
rm -f roostats*
rm -f higgsCombineTest*.root

done

date
echo "done"
"""
  runFile.write(batchString)
  runFile.close()

  runFile = open(outDir+"getStatus.sh","w")
  batchString = \
"""#!/bin/bash

echo "==========================="
echo "Files Found: `ls *.out | wc -l` of `ls *.txt | wc -l`"
echo "==========================="
for i in *.out; do wc $i; done
echo "==========================="
"""
  runFile.write(batchString)
  runFile.close()
