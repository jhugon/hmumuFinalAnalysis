#!/usr/bin/python

import math
import ROOT as root
from helpers import *
import datetime
import sys
import os.path
import copy

from ROOT import gSystem
gSystem.Load('libRooFit')

from xsec import *

if scaleHiggsBy != 1.0:
  print("Error: higgs xsec is scaled!!! Return to 1. Exiting.")
  sys.exit(1)

def getIntegralAll(hist):
  if hist.InheritsFrom("TH2"):
    nBinsX = hist.GetXaxis().GetNbins()
    nBinsY = hist.GetYaxis().GetNbins()
    return hist.Integral(0,nBinsX+1,0,nBinsY+1)
  elif hist.InheritsFrom("TH1"):
    nBinsX = hist.GetXaxis().GetNbins()
    return hist.Integral(0,nBinsX+1)
  else:
    return -1

class MVAvMassPDF:
  def __init__(self,name,canvas,hist2D,massLowRange,massHighRange):
    maxMass = massHighRange[1]
    minMass = massLowRange[0]
    mMuMu = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)
    mMuMu.setRange("low",massLowRange[0],massLowRange[1])
    mMuMu.setRange("high",massHighRange[0],massHighRange[1])
    mMuMu.setRange("signal",massLowRange[1],massLowRange[0])
    mva = root.RooRealVar("mva","mva",-1,1)
    
    a0 = root.RooRealVar("a0","a0",-10.0,10.0)
    a1 = root.RooRealVar("a1","a1",-10.0,10.0)
    a2 = root.RooRealVar("a2","a2",-10.0,10.0)
    
    polyArgsMmumu = root.RooArgList(a0,a1,a2)
    polyMmumu = root.RooPolynomial("polyFunc","polyFunc",mMuMu,polyArgsMmumu)
    
    tmpAxis = hist2D.GetXaxis()
    lowBin = tmpAxis.FindBin(minMass)
    highBin = tmpAxis.FindBin(maxMass)
    mMuMuHist = hist2D.ProjectionX("_mMuMuHist",lowBin,highBin)
    mMuMuRooDataHist = root.RooDataHist("mMuMuRooDataHist","mMuMuRooDataHist",root.RooArgList(mMuMu),mMuMuHist)
    
    pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("low,high"))

    ###################
    
    mvaHistLow = hist2D.ProjectionY("_mvaLow",lowBin,tmpAxis.FindBin(massLowRange[1]))
    mvaHistHigh = hist2D.ProjectionY("_mvaHigh",tmpAxis.FindBin(massHighRange[0]),highBin)
    
    mvaHist = mvaHistLow.Clone()
    mvaHist.Add(mvaHistHigh)
    
    mvaRooDataHist = root.RooDataHist("mvaRooDataHist","mvaRooDataHist",root.RooArgList(mva),mvaHist)
    
    pdfMva = root.RooHistPdf("pdfMva","pdfMva",root.RooArgSet(mva),mvaRooDataHist)
    
    pdfMva.fitTo(mvaRooDataHist)

    ###################

    pdf2d = root.RooProdPdf("pdf2d","pdf2d",RooArgList(pdfMmumu,pdfMva))

    ###################

    plotMmumu = mMuMu.frame()
    plotMmva = mva.frame()

    mMuMuRooDataHist.plotOn(plotMmumu)
    pdfMmumu.plotOn(plotMmumu)

    mvaRooDataHist.plotOn(plotMva)
    pdfMva.plotOn(plotMva)

    plotMmumu.Draw()
    saveAs(canvas,name+"_mMuMu")
    plotMva.Draw()
    saveAs(canvas,name+"_mva")
    
    mMuMuBinning = root.RooFit.Binning(mMuMuHist.GetNbinsX(),minMass,maxMass)
    mvaBinning = root.RooFit.Binning(mvaHist.GetNbinsX(),-1,1)
    pdf2dHist = pdf2d.createHistogram("pdf2dHist",mMuMu,mMuMuBinning,root.RooFit.YVar(mva,mvaBinning))

    canvas.Clear()
    canvas.Divide(2,2)

    canvas.cd(1)
    pdf2dHist.SetLineColor(root.kRed)
    pdf2dHist.Draw("surf")
    hist2d.GetXaxis().SetRangeUser(minMass,maxMass)
    hist2d.GetYaxis().SetRangeUser(-1,1)
    hist2d.Draw("surf same")

    canvas.cd(2)
    hist2d.Draw("colz")
    canvas.cd(3)
    pdf2dHist.Draw("colz")

    saveAs(canvas,name+"_2d")

    #########################

    self.maxMass = maxMass

    self.maxMass = maxMass
    self.minMass = minMass
    self.mMuMu = mMuMu
    self.mva = mva
    self.a0 = a0
    self.a1 = a1
    self.a2 = a2
    self.polyMmumu = polyMmumu
    self.mMuMuHist = mMuMuHist
    self.mMuMuRooDataHist = mMuMuRooDataHist
    self.mvaHistLow = mvaHistLow
    self.mvaHistHigh = mvaHistHigh
    self.mvaHist = mvaHist
    self.mvaRooDataHist = mvaRooDataHist
    self.pdfMva = pdfMva
    self.pdf2d = pdf2d
    self.plotMmumu = plotMmumu
    self.plotMmva = plotMmva
    self.pdf2dHist = pdf2dHist

class Analysis:
  def __init__(self,directory,signalNames,backgroundNames,analysis,histNameBase="mDiMu"):
    self.sigNames = signalNames
    self.bakNames = backgroundNames

    self.sigFiles = []
    self.sigHistsRaw = []
    for name in signalNames:
      tmpF = root.TFile(directory+name+".root")
      tmpH = tmpF.Get(histNameBase+analysis)
      self.sigFiles.append(tmpF)
      self.sigHistsRaw.append(tmpH)

    self.bakFiles = []
    self.bakHistsRaw = []
    for name in backgroundNames:
      tmpF = root.TFile(directory+name+".root")
      tmpH = tmpF.Get(histNameBase+analysis)
      self.bakFiles.append(tmpF)
      self.bakHistsRaw.append(tmpH)

    effMap = {}
    xsecMap = {}
    lowBin = 0
    highBin = self.sigHistsRaw[0].GetNbinsX()+1

    self.xsecSigTotal = 0.0
    self.xsecSigList = []
    self.effSigList = []
    self.sigHists = []
    for h,name in zip(self.sigHistsRaw,signalNames):
      counts = getIntegralAll(h)
      eff = counts/nEventsMap[name]
      xs = eff*xsec[name]
      self.xsecSigTotal += xs
      self.xsecSigList.append(xs)
      self.effSigList.append(eff)
      h.Scale(xsec[name]/nEventsMap[name])
      self.sigHists.append(h)

    self.xsecBakTotal = 0.0
    self.xsecBakList = []
    self.effBakList = []
    self.bakHists = []
    for h,name in zip(self.bakHistsRaw,backgroundNames):
      #counts = h.Integral(lowBin,highBin)
      counts = h.Integral()
      eff = counts/nEventsMap[name]
      xs = eff*xsec[name]
      self.xsecBakTotal += xs
      self.xsecBakList.append(xs)
      self.effBakList.append(eff)
      h.Scale(xsec[name]/nEventsMap[name])
      self.bakHists.append(h)

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
  def getSigHist(self,sigName):
    result = -1.0
    if self.sigNames.count(sigName)>0:
        i = self.sigNames.index(sigName)
        result = self.sigHists[i]
    return result
  def getBakHist(self,bakName):
    result = -1.0
    if self.bakNames.count(bakName)>0:
        i = self.bakNames.index(bakName)
        result = self.bakHists[i]
    return result

class DataCardMaker:
  def __init__(self,directory,analysisNames,signalNames,backgroundNames,nuisanceMap=None,histNameBase="mDiMu"):
    channels = []
    self.channelNames = copy.deepcopy(analysisNames)
    for analysis in analysisNames:
      tmp = Analysis(directory,signalNames,backgroundNames,analysis,histNameBase=histNameBase)
      channels.append(tmp)
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
        iProc = -len(channel.sigNames)+1
        for sigName in self.sigNames:
          binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          binFormatList.append(channelName)
  
          proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc1FormatList.append(sigName)
  
          proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc2FormatList.append(iProc)
  
          expNum = channel.getSigXSec(sigName)*lumi
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
          for sigName in self.sigNames:
            formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            value = "-"
            if thisNu[1].count(sigName)>0:
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

class ShapeDataCardMaker(DataCardMaker):
  def __init__(self,directory,analysisNames,signalNames,backgroundNames,nuisanceMap=None,histNameBase="",rebin=[],useTH1=False):
    DataCardMaker.__init__(self,directory,analysisNames,signalNames,backgroundNames,nuisanceMap,histNameBase)
    if len(rebin) == 2:
      for channel in self.channels:
        for hist in channel.sigHists:
          hist.Rebin2D(*rebin)
        for hist in channel.bakHists:
          hist.Rebin2D(*rebin)
    elif len(rebin) == 1:
      for channel in self.channels:
        for hist in channel.sigHists:
          hist.Rebin(*rebin)
        for hist in channel.bakHists:
          hist.Rebin(*rebin)

    self.is2D = False
    self.useTH1 = useTH1

    for channel in self.channels:
      for hist in channel.sigHists:
        self.x = RooRealVar('x','x',hist.GetXaxis().GetXmin(),hist.GetXaxis().GetXmax())
        if hist.InheritsFrom("TH2"):
          self.y = RooRealVar('y','y',hist.GetYaxis().GetXmin(),hist.GetYaxis().GetXmax())
          self.is2D = True

  def MakeRFHistWrite(self,hist):
    if self.useTH1:
      hist.Write()
    else:
      rfHist = None
      if self.is2D:
        rfHist = RooDataHist(hist.GetName(),hist.GetName(),RooArgList(RooArgSet(self.x,self.y)),hist)
      else:
        rfHist = RooDataHist(hist.GetName(),hist.GetName(),RooArgList(RooArgSet(self.x)),hist)
      rfHist.Write()

  def histFloor(self,hist,integral=False):
    nBinsX = hist.GetNbinsX()
    if integral:
      integral = getIntegralAll(hist)
      desiredIntegral = math.floor(integral)
      hist.Scale(desiredIntegral/integral)
    else:
      if hist.InheritsFrom("TH2"):
        nBinsY = hist.GetNbinsY()
        for i in range(0,nBinsX+2):
         for j in range(0,nBinsY+2):
          tmp = hist.GetBinContent(i,j)
          tmp = math.floor(tmp)
          hist.SetBinContent(i,j,tmp)
      else:
       for i in range(0,nBinsX+2):
        tmp = hist.GetBinContent(i)
        tmp = math.floor(tmp)
        hist.SetBinContent(i,tmp)

  def write(self,outfilename,lumi):
    outRootFilename = re.sub(r"\.txt",r".root",outfilename)
    print("Writing Card: {0} & {1}".format(outfilename,outRootFilename))
    lumi *= 1000.0
    nuisance = self.nuisance
    outfile = open(outfilename,"w")
    outRootFile = root.TFile(outRootFilename, "RECREATE")
    outRootFile.cd()
    outfile.write("# Hmumu shape combine datacard produced by makeTables.py\n")
    now = datetime.datetime.now().replace(microsecond=0).isoformat(' ')
    outfile.write("# {0}\n".format(now))
    outfile.write("############################### \n")
    outfile.write("############################### \n")
    outfile.write("imax {0}\n".format(len(self.channels)))
    #outfile.write("jmax {0}\n".format(len(backgroundNames)))
    outfile.write("jmax {0}\n".format("*"))
    outfile.write("kmax {0}\n".format(len(nuisance)))
    outfile.write("------------\n")
    outfile.write("shapes * * {0} $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC\n".format( os.path.basename(outRootFilename)))
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
      print("text Observed {}: {}".format(channelName,int((channel.getSigXSecTotal()+channel.getBakXSecTotal())*lumi)))
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
        iProc = -len(channel.sigNames)+1
        sumAllMCHist = None
        tmpDir = outRootFile.mkdir(channelName)
        tmpDir.cd()
        for sigName in self.sigNames:
          binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          binFormatList.append(channelName)
  
          proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc1FormatList.append(sigName)
  
          proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc2FormatList.append(iProc)
  
          expNum = channel.getSigXSec(sigName)*lumi
          decimals = ".4f"
          if expNum>1000.0:
            decimals = ".4e"
          rateFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+decimals+"} "
          rateFormatList.append(expNum)

          tmpHist = channel.getSigHist(sigName).Clone(sigName)
          tmpHist.Scale(lumi)
          self.MakeRFHistWrite(tmpHist)
          
          if sumAllMCHist == None:
            sumAllMCHist = tmpHist.Clone("data_obs")
          else:
            sumAllMCHist.Add(tmpHist)
  
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

          tmpHist = channel.getBakHist(bakName).Clone(bakName)
          tmpHist.Scale(lumi)
          self.MakeRFHistWrite(tmpHist)
          if sumAllMCHist == None:
            sumAllMCHist = tmpHist.Clone("data_obs")
          else:
            sumAllMCHist.Add(tmpHist)
  
          iParam += 1
          iProc += 1
        self.histFloor(sumAllMCHist,integral=True)
        print("hist Observed generic integral: {}".format(sumAllMCHist.Integral()))
        print("hist Observed include all integral: {}".format(getIntegralAll(sumAllMCHist)))
        self.MakeRFHistWrite(sumAllMCHist)
        outRootFile.cd()
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
          for sigName in self.sigNames:
            formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            value = "-"
            if thisNu[1].count(sigName)>0:
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

    outRootFile.Close()

if __name__ == "__main__":
  print "Started makeCards.py"

  directory = "input/"
  outDir = "statsCards/"
  analyses = ["mDiMu"]
  analyses2D = ["likelihoodHistMuonOnlyVMass","likelihoodHistVBFVMass","BDTHistMuonOnlyVMass","BDTHistVBFVMass"]
  signalNames=["ggHmumu125","vbfHmumu125"]
  backgroundNames= ["DYJetsToLL","ttbar"]
  lumiList = [5,10,15,20,25,30,40,50,75,100,200,500,1000]
  lumiList = [10,20,30,100]

  for ana in analyses:
    dataCardMassShape = ShapeDataCardMaker(directory,[ana],signalNames,backgroundNames,rebin=[4])
    for i in lumiList:
      dataCardMassShape.write(outDir+ana+"_"+str(i)+".txt",i)

  for ana in analyses2D:
    dataCardMassShape = ShapeDataCardMaker(directory,[ana],signalNames,backgroundNames,rebin=[4,40])
    for i in lumiList:
      dataCardMassShape.write(outDir+ana+"_"+str(i)+".txt",i)

  dataCardBDTComb = ShapeDataCardMaker(directory,["BDTHistMuonOnlyVMass","BDTHistVBFVMass"],signalNames,backgroundNames,rebin=[4,40])
  for i in lumiList:
      dataCardBDTComb.write(outDir+"BDTComb"+"_"+str(i)+".txt",i)

  dataCardLHComb = ShapeDataCardMaker(directory,["likelihoodHistMuonOnlyVMass","likelihoodHistVBFVMass"],signalNames,backgroundNames,rebin=[4,40])
  for i in lumiList:
      dataCardLHComb.write(outDir+"LHComb"+"_"+str(i)+".txt",i)

  runFile = open(outDir+"run.sh","w")
  runFile.write("#!/bin/bash\n")
  runFile.write("\n")
  runFile.write("for i in *.txt; do\n")
  runFile.write('    [[ -e "$i" ]] || continue\n')
  runFile.write('echo "Running on "$i\n')
  runFile.write('nice combine -M Asymptotic $i > $i.out\n')
  runFile.write('done\n')
  runFile.close()
