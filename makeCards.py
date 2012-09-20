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

###################################################################################

class MVAvMassPDFBak:
  def __init__(self,name,hist2D,massLowRange,massHighRange,smooth=False):
    hist2DSmooth = hist2D.Clone(hist2D.GetName()+"_smoothed")
    if smooth:
      hist2DSmooth.Smooth()

    maxMass = massHighRange[1]
    minMass = massLowRange[0]
    mMuMu = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)
    mMuMu.setRange("low",massLowRange[0],massLowRange[1])
    mMuMu.setRange("high",massHighRange[0],massHighRange[1])
    mMuMu.setRange("signal",massLowRange[1],massLowRange[0])
    mva = root.RooRealVar("mva","mva",-1,1)
    
    bwWidth = root.RooRealVar("bwWidth","bwWidth",0.0,30.0)
    bwmZ = root.RooRealVar("bwmZ","bwmZ",85,95)
    pdfMmumu = root.RooBreitWigner("pdfMmumu","pdfMmumu",mMuMu,bwmZ,bwWidth)
    
    tmpAxis = hist2D.GetXaxis()
    lowBin = tmpAxis.FindBin(minMass)
    highBin = tmpAxis.FindBin(maxMass)

    mMuMuHist = hist2D.ProjectionX("_mMuMuHist")
    mMuMuRooDataHist = root.RooDataHist("mMuMuRooDataHist","mMuMuRooDataHist",root.RooArgList(mMuMu),mMuMuHist)
    mMuMuHistSmooth = hist2DSmooth.ProjectionX("_mMuMuHist")
    mMuMuRooDataHistSmooth = root.RooDataHist("mMuMuRooDataHistSmooth","mMuMuRooDataHistSmooth",root.RooArgList(mMuMu),mMuMuHistSmooth)
    
    pdfMmumu.fitTo(mMuMuRooDataHistSmooth,root.RooFit.Range("low,high"))

    ###################
    
    mvaHistLow = hist2D.ProjectionY("_mvaLow",lowBin,tmpAxis.FindBin(massLowRange[1]))
    mvaHistHigh = hist2D.ProjectionY("_mvaHigh",tmpAxis.FindBin(massHighRange[0]),highBin)
    mvaHist = mvaHistLow.Clone()
    mvaHist.Add(mvaHistHigh)
    mvaRooDataHist = root.RooDataHist("mvaRooDataHist","mvaRooDataHist",root.RooArgList(mva),mvaHist)


    mvaHistLowSmooth = hist2DSmooth.ProjectionY("_mvaLow",lowBin,tmpAxis.FindBin(massLowRange[1]))
    mvaHistHighSmooth = hist2DSmooth.ProjectionY("_mvaHigh",tmpAxis.FindBin(massHighRange[0]),highBin)
    mvaHistSmooth = mvaHistLowSmooth.Clone()
    mvaHistSmooth.Add(mvaHistHighSmooth)
    mvaRooDataHistSmooth = root.RooDataHist("mvaRooDataHistSmooth","mvaRooDataHistSmooth",root.RooArgList(mva),mvaHistSmooth)
    
    templatePdfMva = root.RooHistPdf("templatePdfMva","templatePdfMva",root.RooArgSet(mva),mvaRooDataHistSmooth)

    gMvaWidth = root.RooRealVar("gMvaWidth","gMvaWidth",0.0,10)
    gMvaMean = root.RooRealVar("gMvaMean","gMvaMean",-1,1)
    gausPdfMva = root.RooGaussian("gausPdfMva","gausPdfMva",mva,gMvaMean,gMvaWidth)

    pdfMva = root.RooFFTConvPdf("pdfMva","pdfMva",mva,templatePdfMva,gausPdfMva)

    pdfMva.fitTo(mvaRooDataHistSmooth)

    ###################

    pdf2d = root.RooProdPdf("pdf2d","pdf2d",root.RooArgList(pdfMmumu,pdfMva))

    ###################

    plotMmumu = mMuMu.frame()
    plotMva = mva.frame()

    mMuMuRooDataHist.plotOn(plotMmumu)
    pdfMmumu.plotOn(plotMmumu)

    mvaRooDataHist.plotOn(plotMva)
    pdfMva.plotOn(plotMva)

    mMuMuBinning = root.RooFit.Binning(mMuMuHist.GetNbinsX(),minMass,maxMass)
    mvaBinning = root.RooFit.Binning(mvaHist.GetNbinsX(),-1,1)
    pdf2dHist = pdf2d.createHistogram("pdf2dHist",mMuMu,mMuMuBinning,root.RooFit.YVar(mva,mvaBinning))

    #########################

    self.hist2D = hist2D
    self.hist2DSmooth = hist2DSmooth
    self.mMuMuHistSmooth = mMuMuHistSmooth
    self.mMuMuRooDataHistSmooth = mMuMuRooDataHistSmooth
    self.mvaHistLowSmooth = mvaHistLowSmooth
    self.mvaHistHighSmooth = mvaHistHighSmooth
    self.mvaHistSmooth = mvaHistSmooth
    self.mvaRooDataHistSmooth = mvaRooDataHistSmooth

    self.maxMass = maxMass
    self.minMass = minMass
    self.mMuMu = mMuMu
    self.mva = mva
    self.bwWidth = bwWidth
    self.bwmZ = bwmZ
    self.pdfMmumu = pdfMmumu
    self.mMuMuHist = mMuMuHist
    self.mMuMuRooDataHist = mMuMuRooDataHist
    self.mvaHistLow = mvaHistLow
    self.mvaHistHigh = mvaHistHigh
    self.mvaHist = mvaHist
    self.mvaRooDataHist = mvaRooDataHist
    self.pdfMva = pdfMva
    self.pdf2d = pdf2d
    self.plotMmumu = plotMmumu
    self.plotMva = plotMva
    self.pdf2dHist = pdf2dHist
    self.name = name

    self.templatePdfMva = templatePdfMva
    self.gMvaWidth = gMvaWidth
    self.gMvaMean = gMvaMean
    self.gausPdfMva = gausPdfMva

  def writeDebugHists(self,canvas,compareHist=None):
    canvas.cd()
    canvas.SetLogy(1)
    self.plotMmumu.Draw()
    saveAs(canvas,self.name+"_mMuMu")
    self.plotMva.Draw()
    saveAs(canvas,self.name+"_mva")

    canvas.Clear()
    canvas.Divide(2,2)

    canvas.cd(1)
    self.pdf2dHist.SetLineColor(root.kRed)
    self.pdf2dHist.Draw("surf")
    self.hist2D.GetXaxis().SetRangeUser(self.minMass,self.maxMass)
    self.hist2D.GetYaxis().SetRangeUser(-1,1)
    self.hist2D.Draw("surf same")

    canvas.cd(2)
    self.hist2D.Draw("colz")
    canvas.cd(4)
    self.pdf2dHist.Draw("colz")

    if compareHist != None:
      compareHist = compareHist.Clone("mySig")
      compareHist.Scale(self.pdf2dHist.Integral()/compareHist.Integral())
      canvas.cd(3)
      self.pdf2dHist.Draw("colz")
      compareHist.SetFillStyle(0)
      compareHist.SetFillColor(0)
      compareHist.SetLineStyle(1)
      compareHist.SetLineColor(1)
      compareHist.Draw("box same")

    saveAs(canvas,self.name+"_2d")

    debugFileName = "pdfHists_debug_"+self.name+".root"
    tmpFile = root.TFile(debugFileName,"RECREATE")
    tmpFile.cd()
    self.pdf2dHist.Write()
    self.hist2D.Write()
    if compareHist != None:
      compareHist.Write()
    self.plotMmumu.Write()
    self.plotMva.Write()
    self.pdfMmumu.Write()
    self.pdfMva.Write()
    self.pdf2d.Write()
    tmpFile.Close()
    
    canvas.SetLogy(0)


###################################################################################

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

###################################################################################

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
      #self.nuisance["bg_dy"] = (0.05,["DYJetsToLL"])
      #self.nuisance["bg_tt"] = (0.05,["ttbar"])
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

###################################################################################

class ShapeDataCardMaker(DataCardMaker):
  def __init__(self,directory,analysisNames,signalNames,backgroundNames,nuisanceMap=None,histNameBase="",rebin=[],useTH1=False,controlRegionLow=[110,115],controlRegionHigh=[135,150]):
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
    self.controlRegionHigh = controlRegionHigh
    self.controlRegionLow = controlRegionLow

    for channel in self.channels:
      for hist in channel.sigHists:
        self.x = root.RooRealVar('x','x',hist.GetXaxis().GetXmin(),hist.GetXaxis().GetXmax())
        if hist.InheritsFrom("TH2"):
          self.y = root.RooRealVar('y','y',hist.GetYaxis().GetXmin(),hist.GetYaxis().GetXmax())
          self.is2D = True

  def MakeRFHistWrite(self,hist,thisDir=None):
    if thisDir != None:
      thisDir.cd()
    if self.useTH1:
      hist.Write()
    else:
      rfHist = None
      if self.is2D:
        rfHist = root.RooDataHist(hist.GetName(),hist.GetName(),root.RooArgList(root.RooArgSet(self.x,self.y)),hist)
      else:
        rfHist = root.RooDataHist(hist.GetName(),hist.GetName(),root.RooArgList(root.RooArgSet(self.x)),hist)
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

  def write(self,outfilename,lumi,sumAllBak=True,writeBakPDF=True,smooth=False,includeSigInAllMC=False):
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
      observedNumber = int(channel.getBakXSecTotal()*lumi)
      if includeSigInAllMC:
        observedNumber = int((channel.getSigXSecTotal()+channel.getBakXSecTotal())*lumi)
      observationFormatList.append(observedNumber)
      print("text Observed {}: {}".format(channelName,observedNumber))
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
        sumAllSigMCHist = None
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
          self.MakeRFHistWrite(tmpHist,tmpDir)
          
          if includeSigInAllMC:
            if sumAllMCHist == None:
              sumAllMCHist = tmpHist.Clone("data_obs")
            else:
              sumAllMCHist.Add(tmpHist)
          if sumAllMCHist == None:
            sumAllSigMCHist = tmpHist.Clone("sig")
          else:
            sumAllSigMCHist.Add(tmpHist)
  
          iParam += 1
          iProc += 1

        if sumAllBak:
          sumAllBakMCHist = None
          for bakName in self.bakNames:
            tmpHist = channel.getBakHist(bakName).Clone(bakName)
            tmpHist.Scale(lumi)
            if sumAllBakMCHist == None:
              sumAllBakMCHist = tmpHist.Clone("bak")
            else:
              sumAllBakMCHist.Add(tmpHist)

          binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          binFormatList.append(channelName)
    
          proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc1FormatList.append("bak")
    
          proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc2FormatList.append(iProc)

          expNum = channel.getBakXSecTotal()*lumi
          decimals = ".4f"
          if expNum>1000.0:
            decimals = ".4e"
          rateFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+decimals+"} "
          rateFormatList.append(expNum)

          if writeBakPDF:
            c1 = root.TCanvas("c1")
            bakPDFMaker = MVAvMassPDFBak("pdfHists_"+channelName,sumAllBakMCHist,
                                self.controlRegionLow, self.controlRegionHigh,smooth=smooth)
            bakPDFMaker.writeDebugHists(c1,compareHist=sumAllSigMCHist)
            bakPDFMaker.pdf2d.SetName("bak")
            tmpDir.cd()
            bakPDFMaker.pdf2d.Write()

            outfile.write("# Background Debug: Breit-Wigner mZ    = {0:.4g}\n".format(bakPDFMaker.bwmZ.getVal()))
            outfile.write("# Background Debug: Breit-Wigner width = {0:.4g}\n".format(bakPDFMaker.bwWidth.getVal()))
          else:
            self.MakeRFHistWrite(sumAllBakMCHist,tmpDir)

          if sumAllMCHist == None:
            sumAllMCHist = sumAllBakMCHist.Clone("data_obs")
          else:
            sumAllMCHist.Add(sumAllBakMCHist)
      
          iParam += 1
          iProc += 1
        else:
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
            self.MakeRFHistWrite(tmpHist,tmpDir)
            if sumAllMCHist == None:
              sumAllMCHist = tmpHist.Clone("data_obs")
            else:
              sumAllMCHist.Add(tmpHist)
    
            iParam += 1
            iProc += 1
        self.histFloor(sumAllMCHist,integral=True)
        print("hist Observed generic integral: {}".format(sumAllMCHist.Integral()))
        print("hist Observed include all integral: {}".format(getIntegralAll(sumAllMCHist)))
        self.MakeRFHistWrite(sumAllMCHist,tmpDir)
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
          if sumAllBak:
              bakName="bak"
              formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
              value = "-"
              if thisNu[1].count(bakName)>0:
                value = thisNu[0]+1.0
              formatList.append(value)
              iParam += 1
          else:
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
  analyses2D = ["likelihoodHistMuonOnlyVMass","likelihoodHistVBFVMass","BDTHistMuonOnlyVMass","BDTHistVBFVMass"]
  signalNames=["ggHmumu125","vbfHmumu125"]
  backgroundNames= ["DYJetsToLL","ttbar"]
  lumiList = [5,10,15,20,25,30,40,50,75,100,200,500,1000]
  lumiList = [10,20,30,100]

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
