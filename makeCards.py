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

root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

NPROCS = 1

BAKUNC = 0.1

BAKUNCON = True

from xsec import *

if scaleHiggsBy != 1.0:
  print("Error: higgs xsec is scaled!!! Return to 1. Exiting.")
  sys.exit(1)

def getXbinsHighLow(hist,low,high):
  axis = hist.GetXaxis()
  xbinLow = axis.FindBin(low)
  xbinHigh = axis.FindBin(high)
  #print("xbinhigh: {}, {}, {}".format(xbinHigh,axis.GetBinLowEdge(xbinHigh),float(high)))
  if axis.GetBinLowEdge(xbinHigh)==float(high):
    xbinHigh -= 1
  return xbinLow, xbinHigh

def getIntegralAll(hist,boundaries=[]):
  xbinLow = None
  xbinHigh = None
  if len(boundaries)==0:
    xbinLow = 0
    xbinHigh = hist.GetXaxis().GetNbins()
  elif len(boundaries)==2:
    xbinLow, xbinHigh = getXbinsHighLow(hist,boundaries[0],boundaries[1])
  else:
    return -1
  if hist.InheritsFrom("TH2"):
    nBinsY = hist.GetYaxis().GetNbins()
    return hist.Integral(xbinLow,xbinHigh,0,nBinsY+1)
  elif hist.InheritsFrom("TH1"):
    return hist.Integral(xbinLow,xbinHigh)
  else:
    return -1

def getIntegralLowHigh(hist,lowBoundaries,highBoundaries):
  lowInt = getIntegralAll(hist,lowBoundaries)
  highInt = getIntegralAll(hist,highBoundaries)
  return lowInt+highInt

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
      tmpF = root.TFile(directory+name+".root")
      hist = tmpF.Get(histNameBase+analysis)
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

class MassPDFBakMSSM:
  def __init__(self,name,hist,massLowRange,massHighRange,massVeryLowRange,rooVars=None,smooth=False,hack=True):
    if rooVars == None:
        print("Error: MVAvMassPDFBakMSSM requires rooVars list of variables, exiting.")
        sys.exit(1)
    if len(massVeryLowRange)==0:
        print("Error: MVAvMassPDFBakMSSM requires verylow mass range, exiting.")
        sys.exit(1)

    self.debug = ""
    self.debug += "### MassPDFBakMSSM: "+name+"\n"

    maxMass = massHighRange[1]
    minMass = massVeryLowRange[0]
    mMuMu = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)
    mMuMu.setRange("z",88,94)
    mMuMu.setRange("verylow",massVeryLowRange[0],massVeryLowRange[1])
    mMuMu.setRange("low",massLowRange[0],massLowRange[1])
    mMuMu.setRange("high",massHighRange[0],massHighRange[1])
    mMuMu.setRange("signal",massLowRange[1],massHighRange[0])

    bwmZ = root.RooRealVar("bwmZ","bwmZ",85,95)
    bwSig = root.RooRealVar("bwSig","bwSig",0.0,30.0)
    bwLambda = root.RooRealVar("bwLambda","bwLambda",-1e-03,-1e-01,-1e-04)
    bwMmumu = root.RooGenericPdf("bwMmumu","exp(bwLambda*mMuMu)*bwSig/(((mMuMu-bwmZ)*(mMuMu-bwmZ) + bwSig*bwSig*0.25))",root.RooArgList( mMuMu, bwLambda, bwSig, bwmZ))

    phoMmumu = root.RooGenericPdf("phoMmumu","exp(bwLambda*mMuMu)/pow(mMuMu,2)",root.RooArgList(mMuMu,bwLambda))

    mixParam = root.RooRealVar("mixParam","mixParam",0.5,0,1)

    pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",root.RooArgList(bwMmumu,phoMmumu),root.RooArgList(mixParam))
    
    tmpAxis = hist.GetXaxis()
    lowBin = tmpAxis.FindBin(minMass)
    highBin = tmpAxis.FindBin(maxMass)
    lowBin, highBin = getXbinsHighLow(hist,minMass,maxMass)
    nBinsX = highBin - lowBin + 1
    normalization = getIntegralAll(hist,[massLowRange[0],massHighRange[1]])
    self.debug += "# bak Hist no bounds: {:.3g}\n".format(hist.Integral())
    self.debug += "# bak Hist bounds:    {:.3g}\n".format(normalization)

    mMuMuRooDataHist = root.RooDataHist(name+"DataHist",name+"DataHist",root.RooArgList(mMuMu),hist)

    bwMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("z"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    bwmZ.setConstant(True)
    bwSig.setConstant(True)

    #phoMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("high"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    #bwLambda.setConstant(True)
    
    pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    chi2 = pdfMmumu.createChi2(mMuMuRooDataHist)

    plotMmumu = mMuMu.frame()

    mMuMuRooDataHist.plotOn(plotMmumu)
    pdfMmumu.plotOn(plotMmumu,root.RooFit.Range("low,high"))
    pdfMmumu.plotOn(plotMmumu,root.RooFit.LineStyle(2),root.RooFit.Range(minMass,maxMass))
    pdfMmumu.plotOn(plotMmumu,root.RooFit.Components("phoMmumu"),root.RooFit.LineStyle(2),root.RooFit.Range(minMass,maxMass),root.RooFit.LineColor(root.kGreen+1))

    mMuMuBinning = root.RooFit.Binning(nBinsX,minMass,maxMass)
    nominalHist = pdfMmumu.createHistogram("pdf2dHist",mMuMu,mMuMuBinning)
    nominalHist.Scale(normalization/getIntegralAll(nominalHist,[massLowRange[0],massHighRange[1]]))

    self.name = name
    self.hist = hist
    self.mMuMuRooDataHist = mMuMuRooDataHist
    self.lowBin = lowBin
    self.highBin = highBin
    self.nBinsX = nBinsX
    self.bwmZ = bwmZ
    self.bwSig = bwSig
    self.bwLambda = bwLambda
    self.bwMmumu = bwMmumu
    self.phoMmumu = phoMmumu
    self.mixParam = mixParam
    self.pdfMmumu = pdfMmumu
    self.mMuMuBinning = mMuMuBinning
    self.nominalHist = nominalHist
    self.maxMass = maxMass
    self.minMass = minMass
    self.mMuMu = mMuMu
    self.plotMmumu = plotMmumu
    self.chi2 = chi2

    self.debug += "# nominal Integral: {0:.3g}\n".format(getIntegralAll(nominalHist))
    self.debug += "# BW Sigma: {0:.3g} +/- {1:.3g}\n".format(bwSig.getVal(),bwSig.getError())
    self.debug += "# BW mZ:    {0:.3g} +/- {1:.3g}\n".format(bwmZ.getVal(),bwmZ.getError())
    self.debug += "# Lam:      {0:.3g} +/- {1:.3g}\n".format(bwLambda.getVal(),bwLambda.getError())
    self.debug += "# BW Coef:  {0:.3g} +/- {1:.3g}\n".format(mixParam.getVal(),mixParam.getError())
    self.debug += "# chi2/ndf: {0:.3g}\n".format(chi2.getVal()/(nBinsX-1))

    ## Error time

    self.errNames = []
    self.errHists = {}
    if BAKUNCON:
      for errVar in [bwmZ,bwSig,bwLambda,mixParam]:
        val = errVar.getVal()
        err = errVar.getError()
        varName = errVar.GetName()

        errVar.setVal(val+err)
        pdfMmumu.plotOn(plotMmumu,root.RooFit.LineColor(root.kGreen-1),root.RooFit.Range(110,150),root.RooFit.LineStyle(3))
        upHist = pdfMmumu.createHistogram(varName+"Up",mMuMu,mMuMuBinning)

        errVar.setVal(val-err)
        pdfMmumu.plotOn(plotMmumu,root.RooFit.LineColor(root.kGreen-1),root.RooFit.Range(110,150),root.RooFit.LineStyle(3))
        downHist = pdfMmumu.createHistogram(varName+"Down",mMuMu,mMuMuBinning)

        errVar.setVal(val)

        upHist.Scale(normalization/getIntegralAll(upHist,[massLowRange[0],massHighRange[1]]))
        downHist.Scale(normalization/getIntegralAll(downHist,[massLowRange[0],massHighRange[1]]))

        setattr(self,varName+"UpHist",upHist)
        setattr(self,varName+"DownHist",downHist)

        self.errNames.append(varName)
        self.errHists[varName+"Up"] = upHist
        self.errHists[varName+"Down"] = downHist

  def writeDebugHistsToCurrTDir(self,compareHist=None):
    canvas = root.TCanvas("canvas")
    canvas.cd()
    #canvas.SetLogy(1)
    self.plotMmumu.Draw()
    canvas.SetName(self.name+"Canvas")
    canvas.Write()

class MassPDFBak:
  def __init__(self,name,hist,massLowRange,massHighRange,massVeryLowRange,rooVars=None,smooth=False,hack=True):
    if rooVars == None:
        print("Error: MVAvMassPDFBak requires rooVars list of variables, exiting.")
        sys.exit(1)
    if len(massVeryLowRange)==0:
        print("Error: MVAvMassPDFBak requires verylow mass range, exiting.")
        sys.exit(1)

    self.debug = ""
    self.debug += "### MassPDFBak: "+name+"\n"

    maxMass = massHighRange[1]
    minMass = massVeryLowRange[0]
    mMuMu = root.RooRealVar("mMuMu2","mMuMu",minMass,maxMass)
    mMuMu.setRange("z",88,94)
    mMuMu.setRange("verylow",massVeryLowRange[0],massVeryLowRange[1])
    mMuMu.setRange("low",massLowRange[0],massLowRange[1])
    mMuMu.setRange("high",massHighRange[0],massHighRange[1])
    mMuMu.setRange("signal",massLowRange[1],massHighRange[0])

    voitWidth = root.RooRealVar("voitWidth","voitWidth",2.4952)
    voitmZ = root.RooRealVar("voitmZ","voitmZ",85,95)
    voitSig = root.RooRealVar("voitSig","voitSig",0.0,30.0)
    voitMmumu = root.RooVoigtian("voitMmumu","voitMmumu",mMuMu,voitmZ,voitWidth,voitSig)

    expParam = root.RooRealVar("expParam","expParam",-1,0)
    expMmumu = root.RooExponential("expMmumu","expMmumu",mMuMu,expParam)

    mixParam = root.RooRealVar("mixParam","mixParam",0,1)

    pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",root.RooArgList(voitMmumu,expMmumu),root.RooArgList(mixParam))
    
    tmpAxis = hist.GetXaxis()
    lowBin = tmpAxis.FindBin(minMass)
    highBin = tmpAxis.FindBin(maxMass)
    lowBin, highBin = getXbinsHighLow(hist,minMass,maxMass)
    nBinsX = highBin - lowBin + 1
    normalization = getIntegralAll(hist,[massLowRange[0],massHighRange[1]])
    self.debug += "# bak Hist no bounds: {:.3g}\n".format(hist.Integral())
    self.debug += "# bak Hist bounds:    {:.3g}\n".format(normalization)

    mMuMuRooDataHist = root.RooDataHist(name+"DataHist",name+"DataHist",root.RooArgList(mMuMu),hist)

    voitMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("z"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    voitmZ.setConstant(True)
    voitSig.setConstant(True)

    expMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("high"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    expParam.setConstant(True)
    
    pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    chi2 = pdfMmumu.createChi2(mMuMuRooDataHist)

    plotMmumu = mMuMu.frame()

    mMuMuRooDataHist.plotOn(plotMmumu)
    pdfMmumu.plotOn(plotMmumu,root.RooFit.Range("low,high"))
    pdfMmumu.plotOn(plotMmumu,root.RooFit.LineStyle(2),root.RooFit.Range(minMass,maxMass))
    pdfMmumu.plotOn(plotMmumu,root.RooFit.Components("expMmumu"),root.RooFit.LineStyle(2),root.RooFit.Range(minMass,maxMass),root.RooFit.LineColor(root.kGreen+1))

    mMuMuBinning = root.RooFit.Binning(nBinsX,minMass,maxMass)
    nominalHist = pdfMmumu.createHistogram("pdf2dHist",mMuMu,mMuMuBinning)
    nominalHist.Scale(normalization/getIntegralAll(nominalHist,[massLowRange[0],massHighRange[1]]))

    self.name = name
    self.hist = hist
    self.mMuMuRooDataHist = mMuMuRooDataHist
    self.lowBin = lowBin
    self.highBin = highBin
    self.nBinsX = nBinsX
    self.voitWidth = voitWidth
    self.voitmZ = voitmZ
    self.voitSig = voitSig
    self.voitMmumu = voitMmumu
    self.expParam = expParam
    self.expMmumu = expMmumu
    self.mixParam = mixParam
    self.pdfMmumu = pdfMmumu
    self.mMuMuBinning = mMuMuBinning
    self.nominalHist = nominalHist
    self.maxMass = maxMass
    self.minMass = minMass
    self.mMuMu = mMuMu
    self.plotMmumu = plotMmumu
    self.chi2 = chi2

    self.debug += "# nominal Integral: {0:.3g}\n".format(getIntegralAll(nominalHist))
    self.debug += "# V Width: {0:.3g} +/- {1:.3g}\n".format(voitWidth.getVal(),voitWidth.getError())
    self.debug += "# V Sigma: {0:.3g} +/- {1:.3g}\n".format(voitSig.getVal(),voitSig.getError())
    self.debug += "# V mZ:    {0:.3g} +/- {1:.3g}\n".format(voitmZ.getVal(),voitmZ.getError())
    self.debug += "# Exp Par: {0:.3g} +/- {1:.3g}\n".format(expParam.getVal(),expParam.getError())
    self.debug += "# V Coef:  {0:.3g} +/- {1:.3g}\n".format(mixParam.getVal(),mixParam.getError())
    self.debug += "# chi2/ndf: {0:.3g}\n".format(chi2.getVal()/(nBinsX-1))

    ## Error time

    self.errNames = []
    self.errHists = {}
    if BAKUNCON:
      for errVar in [voitmZ,voitSig,expParam,mixParam]:
        val = errVar.getVal()
        err = errVar.getError()
        varName = errVar.GetName()

        errVar.setVal(val+err)
        pdfMmumu.plotOn(plotMmumu,root.RooFit.LineColor(root.kGreen-1),root.RooFit.Range(110,150),root.RooFit.LineStyle(3))
        upHist = pdfMmumu.createHistogram(varName+"Up",mMuMu,mMuMuBinning)

        errVar.setVal(val-err)
        pdfMmumu.plotOn(plotMmumu,root.RooFit.LineColor(root.kGreen-1),root.RooFit.Range(110,150),root.RooFit.LineStyle(3))
        downHist = pdfMmumu.createHistogram(varName+"Down",mMuMu,mMuMuBinning)

        errVar.setVal(val)

        upHist.Scale(normalization/getIntegralAll(upHist,[massLowRange[0],massHighRange[1]]))
        downHist.Scale(normalization/getIntegralAll(downHist,[massLowRange[0],massHighRange[1]]))

        setattr(self,varName+"UpHist",upHist)
        setattr(self,varName+"DownHist",downHist)

        self.errNames.append(varName)
        self.errHists[varName+"Up"] = upHist
        self.errHists[varName+"Down"] = downHist

  def writeDebugHistsToCurrTDir(self,compareHist=None):
    canvas = root.TCanvas("canvas")
    canvas.cd()
    #canvas.SetLogy(1)
    self.plotMmumu.Draw()
    canvas.SetName(self.name+"Canvas")
    canvas.Write()

class MVAvMassPDFBak:
  def __init__(self,name,hist2D,massLowRange,massHighRange,rooVars=None,smooth=False,hack=True):
    if rooVars == None:
        print("Error: MVAvMassPDFBak requires rooVars list of variables, exiting.")
        sys.exit(1)
    print("original hist bins X: {} Y: {}".format(hist2D.GetNbinsX(),hist2D.GetNbinsY()))

    hist2DSmooth = hist2D.Clone(name+hist2D.GetName()+"_smoothed")
    if smooth:
      hist2DSmooth.Smooth()

    maxMass = massHighRange[1]
    minMass = massLowRange[0]
    mMuMu = rooVars[0]
    mMuMu.setRange("low",massLowRange[0],massLowRange[1])
    mMuMu.setRange("high",massHighRange[0],massHighRange[1])
    mMuMu.setRange("signal",massLowRange[1],massHighRange[0])
    mva = rooVars[1]
    
    bwWidth = root.RooRealVar("bwWidth","bwWidth",2.4952)
    bwmZ = root.RooRealVar("bwmZ","bwmZ",85,95)
    #pdfMmumu = root.RooBreitWigner("pdfMmumu","pdfMmumu",mMuMu,bwmZ,bwWidth)
    voitSigma = root.RooRealVar("voitSigma","voitSigma",0.0,80.0)
    pdfMmumu = root.RooVoigtian("pdfMmumu","pdfMmumu",mMuMu,bwmZ,bwWidth,voitSigma)
    
    tmpAxis = hist2D.GetXaxis()
    lowBin, highBin = getXbinsHighLow(hist2D,minMass,maxMass)

    mMuMuHist = hist2D.ProjectionX("_mMuMuHist")
    mMuMuRooDataHist = root.RooDataHist(name+"mMuMuRooDataHist","mMuMuRooDataHist",root.RooArgList(mMuMu),mMuMuHist)
    mMuMuHistSmooth = hist2DSmooth.ProjectionX("_mMuMuHistSmooth")
    mMuMuRooDataHistSmooth = root.RooDataHist(name+"mMuMuRooDataHistSmooth","mMuMuRooDataHistSmooth",root.RooArgList(mMuMu),mMuMuHistSmooth)
    
    pdfMmumu.fitTo(mMuMuRooDataHistSmooth,root.RooFit.Range("low,high"),PRINTLEVEL)

    ###################
    
    silly, lowRangeHighBin = getXbinsHighLow(hist2D,minMass,massLowRange[1])
    mvaHistLow = hist2D.ProjectionY("_mvaLow",lowBin,lowRangeHighBin)
    mvaHistHigh = hist2D.ProjectionY("_mvaHigh",tmpAxis.FindBin(massHighRange[0]),highBin)
    mvaHist = mvaHistLow.Clone()
    mvaHist.Add(mvaHistHigh)
    mvaRooDataHist = root.RooDataHist(name+"mvaRooDataHist","mvaRooDataHist",root.RooArgList(mva),mvaHist)


    mvaHistLowSmooth = hist2DSmooth.ProjectionY("_mvaLowSmooth",lowBin,lowRangeHighBin)
    mvaHistHighSmooth = hist2DSmooth.ProjectionY("_mvaHighSmooth",tmpAxis.FindBin(massHighRange[0]),highBin)
    mvaHistSmooth = mvaHistLowSmooth.Clone()
    mvaHistSmooth.Add(mvaHistHighSmooth)
    mvaRooDataHistSmooth = root.RooDataHist(name+"mvaRooDataHistSmooth","mvaRooDataHistSmooth",root.RooArgList(mva),mvaHistSmooth)
    
    #templatePdfMva = root.RooHistPdf("templatePdfMva","templatePdfMva",root.RooArgSet(mva),mvaRooDataHistSmooth)
    pdfMva = root.RooHistPdf("pdfMva","pdfMva",root.RooArgSet(mva),mvaRooDataHistSmooth)

    #gMvaWidth = root.RooRealVar("gMvaWidth","gMvaWidth",0.0,10)
    #gMvaMean = root.RooRealVar("gMvaMean","gMvaMean",-1,1)
    #gausPdfMva = root.RooGaussian("gausPdfMva","gausPdfMva",mva,gMvaMean,gMvaWidth)

    #pdfMva = root.RooFFTConvPdf("pdfMva","pdfMva",mva,templatePdfMva,gausPdfMva)

    pdfMva.fitTo(mvaRooDataHistSmooth,PRINTLEVEL)

    ###################

    pdf2d = root.RooProdPdf("pdf2d","pdf2d",root.RooArgList(pdfMmumu,pdfMva))

    ###################

    plotMmumu = mMuMu.frame()
    plotMva = mva.frame()

    mMuMuRooDataHist.plotOn(plotMmumu)
    pdfMmumu.plotOn(plotMmumu)
    pdfMmumu.plotOn(plotMmumu,root.RooFit.LineStyle(2),root.RooFit.Range(minMass,maxMass))

    mvaRooDataHist.plotOn(plotMva)
    pdfMva.plotOn(plotMva)

    tmpLowBin,tmpHighBin = getXbinsHighLow(hist2D,minMass,maxMass)
    nBinsX = tmpHighBin - tmpLowBin

    mMuMuBinning = root.RooFit.Binning(nBinsX,minMass,maxMass-0.00000001)
    mvaBinning = root.RooFit.Binning(mvaHist.GetNbinsX(),-1,1)
    pdf2dHist = pdf2d.createHistogram("pdf2dHist",mMuMu,mMuMuBinning,root.RooFit.YVar(mva,mvaBinning))

    #####
    ## mMuMu Errors
    for var,col in zip([bwmZ,voitSigma],[root.kRed,root.kGreen]):
      original = var.getVal()
      err = var.getError()
      var.setVal(original+err)
      pdfMmumu.plotOn(plotMmumu,root.RooFit.LineStyle(2),root.RooFit.LineColor(col),root.RooFit.NormRange("low,high"),root.RooFit.Range(minMass,maxMass))

      var.setVal(original-err)
      pdfMmumu.plotOn(plotMmumu,root.RooFit.LineStyle(2),root.RooFit.LineColor(col),root.RooFit.NormRange("low,high"),root.RooFit.Range(minMass,maxMass))
      var.setVal(original)

    #########################

    self.hist2D = hist2D
    self.hist2DSmooth = hist2DSmooth
    self.mMuMuHistSmooth = mMuMuHistSmooth
    self.mMuMuRooDataHistSmooth = mMuMuRooDataHistSmooth
    self.mvaHistLowSmooth = mvaHistLowSmooth
    self.mvaHistHighSmooth = mvaHistHighSmooth
    self.mvaHistSmooth = mvaHistSmooth
    self.mvaRooDataHistSmooth = mvaRooDataHistSmooth

    self.smooth = smooth
    self.massLowRange = massLowRange
    self.massHighRange = massHighRange
    self.maxMass = maxMass
    self.minMass = minMass
    self.mMuMu = mMuMu
    self.mva = mva
    self.bwWidth = bwWidth
    self.bwmZ = bwmZ
    self.voitSigma = voitSigma
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
    self.debug = ""

    #self.templatePdfMva = templatePdfMva
    #self.gMvaWidth = gMvaWidth
    #self.gMvaMean = gMvaMean
    #self.gausPdfMva = gausPdfMva

    histForErrs = pdf2dHist
    self.hackHist = None
    if hack:
      hackHist = pdf2dHist.Clone("hackHist")
      hackHist.Scale(hist2D.Integral()/hackHist.Integral())
      self.hackHist = hackHist
      self.pdf2d = root.RooDataHist(hist2D.GetName(),hist2D.GetName(),root.RooArgList(mMuMu,mva),hackHist)
      histForErrs = hackHist

    lowPertBin = pdf2dHist.GetXaxis().FindBin(120.0)
    highPertBin = pdf2dHist.GetXaxis().FindBin(131.0)
    lowPertBin, highPertBin = getIntegralAll(pdf2dHist,120.0,131.0)
    print("pdf2dHist bins X: {} Y: {}".format(pdf2dHist.GetNbinsX(),pdf2dHist.GetNbinsY()))
    self.nuisanceNames = []
    self.pdf2dHistErrs = {}
    bakPlus = pdf2dHist.Clone("bakUnc"+"Up")
    bakMinus = pdf2dHist.Clone("bakUnc"+"Down")
    errIntNominal = 0.
    for xBin in range(lowPertBin,highPertBin):
      for yBin in range(0,pdf2dHist.GetNbinsY()+2):
        orig = pdf2dHist.GetBinContent(xBin,yBin)
        err = pdf2dHist.GetBinError(xBin,yBin)*BAKUNC
        bakPlus.SetBinContent(xBin,yBin,orig+err)
        bakMinus.SetBinContent(xBin,yBin,orig-err)
        errIntNominal += pdf2dHist.GetBinError(xBin,yBin)**2
    errIntNominal = sqrt(errIntNominal)
    bakMinus.Scale((pdf2dHist.GetSumOfWeights()-errIntNominal)/bakMinus.GetSumOfWeights())
    bakPlus.Scale((pdf2dHist.GetSumOfWeights()+errIntNominal)/bakMinus.GetSumOfWeights())
    self.pdf2dHistErrs[bakPlus.GetName()] = bakPlus
    self.pdf2dHistErrs[bakMinus.GetName()] = bakMinus
    self.nuisanceNames.append("bakUnc")

  def dump(self):
    print("#####################################")
    print("MVAvMassPDFBak Dump: ")
    print("name: {0}".format(self.name))
    print("hist2D: {0}".format(self.hist2D))
    print("hist2DSmooth: {0}".format(self.hist2DSmooth))
    print("mMuMuHistSmooth: {0}".format(self.mMuMuHistSmooth))
    print("mMuMuRooDataHistSmooth: {0}".format(self.mMuMuRooDataHistSmooth))
    print("mvaHistLowSmooth: {0}".format(self.mvaHistLowSmooth))
    print("mvaHistHighSmooth: {0}".format(self.mvaHistHighSmooth))
    print("mvaHistSmooth: {0}".format(self.mvaHistSmooth))
    print("mvaRooDataHistSmooth: {0}".format(self.mvaRooDataHistSmooth))
    print("smooth: {0}".format(self.smooth))
    print("massLowRange: {0}".format(self.massLowRange))
    print("massHighRange: {0}".format(self.massHighRange))
    print("maxMass: {0}".format(self.maxMass))
    print("minMass: {0}".format(self.minMass))
    print("mMuMu: {0}".format(self.mMuMu))
    print("mva: {0}".format(self.mva))
    print("bwWidth: {0}".format(self.bwWidth))
    print("bwmZ: {0}".format(self.bwmZ))
    print("voitSigma: {0}".format(self.voitSigma))
    print("pdfMmumu: {0}".format(self.pdfMmumu))
    print("mMuMuHist: {0}".format(self.mMuMuHist))
    print("mMuMuRooDataHist: {0}".format(self.mMuMuRooDataHist))
    print("mvaHistLow: {0}".format(self.mvaHistLow))
    print("mvaHistHigh: {0}".format(self.mvaHistHigh))
    print("mvaHist: {0}".format(self.mvaHist))
    print("mvaRooDataHist: {0}".format(self.mvaRooDataHist))
    print("pdfMva: {0}".format(self.pdfMva))
    print("pdf2d: {0}".format(self.pdf2d))
    print("plotMmumu: {0}".format(self.plotMmumu))
    print("plotMva: {0}".format(self.plotMva))
    print("pdf2dHist: {0}".format(self.pdf2dHist))
    print("hackHist: {0}".format(self.hackHist))
    print("#####################################")

  def writeDebugHistsToCurrTDir(self,compareHist=None):

    plottingMassRange = [110,150]
    print("justin writeDebugHistsToCurrTDir")
    canvas = root.TCanvas("canvas")
    canvas.cd()
    #canvas.SetLogy(1)
    self.plotMmumu.Draw()
    canvas.SetName("plotMmumuCanvas")
    canvas.Write()
    self.plotMva.Draw()
    canvas.SetName("plotMvaCanvas")
    canvas.Write()

    canvas.Clear()
    canvas.Divide(2,2)

    canvas.cd(1)
    self.hist2D.SetTitle("Original 2D Hist")
    self.hist2D.SetStats(False)
    self.hist2D.GetXaxis().SetTitle("mMuMu")
    self.hist2D.GetYaxis().SetTitle("MVA")
    self.hist2D.Draw("colz")
    canvas.cd(3)
    self.pdf2dHist.SetTitle("")
    self.pdf2dHist.SetStats(False)

    tmp2 = self.pdf2dHist.Clone("tmp2Hist")
    tmp2.SetTitle("PDF Assuming M & MVA Uncorrelated")
    tmp2.GetXaxis().SetRangeUser(110,self.maxMass)
    tmp2.GetYaxis().SetRangeUser(-1,1)
    tmp2.Draw("colz")

    tmp3 = None
    if compareHist != None:
      compareHist = compareHist.Clone("mySig")
      compareHist.Scale(self.pdf2dHist.Integral()/compareHist.Integral())
      canvas.cd(4)
      
      tmp3 = self.pdf2dHist.Clone("tmp3Hist")
      tmp3.SetTitle("PDF, Compare Signal in Black Boxes")
      tmp3.GetXaxis().SetRangeUser(*plottingMassRange)
      tmp3.GetYaxis().SetRangeUser(-1,1)
      tmp3.Draw("colz")
      compareHist.SetFillStyle(0)
      compareHist.SetFillColor(0)
      compareHist.SetLineStyle(1)
      compareHist.SetLineColor(1)
      compareHist.Draw("box same")

    canvas.cd(2)
    tmp1 = self.pdf2dHist.Clone("tmp1Hist")
    tmp1.SetLineColor(root.kRed)
    tmp1.GetXaxis().SetRangeUser(*plottingMassRange)
    tmp1.GetYaxis().SetRangeUser(-1,1)
    tmp1.SetTitle("Blue Original TH2, Red Final 2D PDF")
    tmp1.Draw("surf")
    self.hist2D.GetXaxis().SetRangeUser(*plottingMassRange)
    self.hist2D.GetYaxis().SetRangeUser(-1,1)
    self.hist2D.Draw("surf same")

    canvas.SetName("2dCanvas")
    canvas.Write()

    tmpPave = root.TPaveText(0,0,1,1)
    tmpPave.SetFillColor(0)
    tmpPave.SetLineColor(1)
    tmpPave.AddText("HistName: "+self.hist2D.GetName())
    tmpPave.AddText("Smoothed: {}".format(self.smooth))
    tmpPave.AddText("Control Region Low: {}".format(self.massLowRange))
    tmpPave.AddText("Control Region High: {}".format(self.massHighRange))
    tmpPave.AddText("pdfMmumu:")
    tmpPave.AddText("  BW mZ: {0:.2g} +/- {1:.2g}".format(self.bwmZ.getVal(),self.bwmZ.getError()))
    tmpPave.AddText("  BW Width: {0:.2g} +/- {1:.2g}".format(self.bwWidth.getVal(),self.bwWidth.getError()))
    tmpPave.AddText("  Voit Sigma: {0:.2g} +/- {1:.2g}".format(self.voitSigma.getVal(),self.voitSigma.getError()))

    canvas.Clear()
    tmpPave.Draw()
    canvas.SetName("infoCanvas")
    canvas.Write()

    self.pdf2dHist.Write()
    self.hist2D.Write()
    if compareHist != None:
      compareHist.Write()
    self.plotMmumu.Write()
    self.plotMva.Write()
    self.pdfMmumu.Write()
    self.pdfMva.Write()
    self.pdf2d.Write()
    if self.hackHist != None:
        self.hackHist.SetName("hackHist")
        self.hackHist.Write()

    """
    first = True
    for key in self.pdf2dHistErrs:
      hist = self.pdf2dHistErrs[key]
      #hist.Add(self.pdf2dHist,-1.0)
      hist.Write()
    """
    #canvas.SetLogy(0)


###################################################################################

class Analysis:
  def __init__(self,directory,signalNames,backgroundNames,dataNames,analysis,x,y,controlRegionLow,controlRegionHigh,histNameBase="mDiMu",bakShape=False,rebin=[],histNameSuffix="",controlRegionVeryLow=[]):
    self.bakShape = bakShape
    self.sigNames = signalNames
    self.bakNames = backgroundNames
    self.datNames = dataNames
    self.controlRegionVeryLow = controlRegionVeryLow
    self.controlRegionLow = controlRegionLow
    self.controlRegionHigh = controlRegionHigh
    self.analysis = analysis

    self.is2D = False
    if y != None:
      self.is2D = True
    self.x = x
    self.y = y
    self.x1d = None

    self.sigFiles = []
    self.sigHistsRaw = []
    for name in signalNames:
      tmpF = root.TFile(directory+name+".root")
      tmpH = tmpF.Get(histNameBase+analysis+histNameSuffix)
      self.sigFiles.append(tmpF)
      self.sigHistsRaw.append(tmpH)
      if tmpH.InheritsFrom("TH2"):
        self.is2D = True

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
    if len(rb) == 2 and self.is2D:
        for hist in self.sigHistsRaw:
          hist.Rebin2D(*rb)
        for hist in self.bakHistsRaw:
          hist.Rebin2D(*rb)
        for hist in self.datHists:
          hist.Rebin2D(*rb)
    elif len(rb) == 1 and not self.is2D:
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
    massBounds = [controlRegionLow[0],controlRegionHigh[1]]
    self.massBounds = massBounds

    for hist in self.sigHistsRaw:
      vetoOutOfBoundsEvents(hist,boundaries=massBounds)

    self.xsecSigTotal = 0.0
    self.xsecSigList = []
    self.effSigList = []
    self.sigHists = []
    for h,name in zip(self.sigHistsRaw,signalNames):
      counts = getIntegralAll(h,boundaries=massBounds)
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
    self.bakHistTotal = None
    for h,name in zip(self.bakHistsRaw,backgroundNames):
      counts = getIntegralAll(h,boundaries=massBounds)
      eff = counts/nEventsMap[name]
      xs = eff*xsec[name]
      self.xsecBakTotal += xs
      self.xsecBakList.append(xs)
      self.effBakList.append(eff)
      h.Scale(xsec[name]/nEventsMap[name])
      self.bakHists.append(h)
      if self.bakHistTotal == None:
        self.bakHistTotal = h.Clone("bak")
      else:
        self.bakHistTotal.Add(h)

    self.bakHistTotalReal = self.bakHistTotal.Clone("data_obs")

    self.bakShape = bakShape
    self.bakShapeMkr = None

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
  def getBakHistTotal(self,lumi):
    analysis = self.analysis
    bakShape = self.bakShape
    if self.datHistTotal == None:
      self.bakHistTotal = self.bakHistTotalReal.Clone("stuff")
      self.bakHistTotal.Scale(lumi)
    else:
      self.bakHistTotal = self.datHistTotal.Clone("stuff")
    if bakShape and self.is2D:
      bakShapeMkr = MVAvMassPDFBak("pdfHists_"+analysis,
                                self.bakHistTotal,
                                self.controlRegionLow,self.controlRegionHigh,
                                rooVars = [self.x,self.y]
                                )
      self.bakShapeMkr = bakShapeMkr
      self.bakHistTotal = bakShapeMkr.hackHist
      #self.xsecBakTotal = getIntegralAll(self.bakHistTotal,boundaries=massBounds)
    elif bakShape:
      bakShapeMkr = MassPDFBak("pdfHists_"+analysis,
                                self.bakHistTotal,
                                self.controlRegionLow,self.controlRegionHigh,
                                self.controlRegionVeryLow,
                                rooVars = [self.x]
                                )
      self.bakShapeMkr = bakShapeMkr
      self.bakHistTotal = bakShapeMkr.nominalHist
      #self.xsecBakTotal = getIntegralAll(self.bakHistTotal,boundaries=self.massBounds)
    return self.bakHistTotal
  def dump(self):
    print("##########################################")
    print("SigXSecTotal: {0}".format(self.xsecSigTotal))
    print("BakXSecTotal: {0}".format(self.xsecBakTotal))
    print("signames: {0}".format(self.sigNames))
    print("baknames: {0}".format(self.bakNames))
    print("sighists: {0}".format(self.sigHists))
    for i in self.sigHists:
        i.Print()
    print("bakhists: {0}".format(self.bakHists))
    for i in self.bakHists:
        i.Print()
    print("bakHistTotal: {0}".format(self.bakHistTotal))
    self.bakHistTotal.Print()
    print("##########################################")

###################################################################################

class DataCardMaker:
  def __init__(self,directory,analysisNames,signalNames,backgroundNames,dataNames,nuisanceMap=None,histNameBase="",controlRegionLow=[80,115],controlRegionHigh=[135,150],controlRegionVeryLow=[],bakShape=False,rebin=[],histNameSuffix=""):
    channels = []
    self.channelNames = copy.deepcopy(analysisNames)
    self.is2D = False

    x = None
    y = None
    for analysis in analysisNames:
      tmpList = getRooVars(directory,signalNames,histNameBase,analysis+histNameSuffix)
      x = tmpList[0]
      if len(tmpList)==2:
        self.is2D = True
        y = tmpList[1]
    self.x = x
    self.y = y
    self.shape = bakShape

    for analysis in analysisNames:
      tmp = Analysis(directory,signalNames,backgroundNames,dataNames,analysis,x,y,controlRegionLow,controlRegionHigh,controlRegionVeryLow=controlRegionVeryLow,histNameBase=histNameBase,bakShape=bakShape,rebin=rebin,histNameSuffix=histNameSuffix)
      channels.append(tmp)
    self.channels = channels

    if nuisanceMap == None:
      self.nuisance = {}
      self.nuisance["lumi"] = (0.044,["vbfHmumu125","ggHmumu125","wHmumu125","zHmumu125"])
      self.nuisance["xs_ggH"] = (0.147,["ggHmumu125"])
      self.nuisance["xs_vbfH"] = (0.03,["vbfHmumu125"])
      self.nuisance["xs_wH"] = (0.041,["wHmumu125"])
      self.nuisance["xs_zH"] = (0.051,["zHmumu125"])
      self.nuisance["br_Hmm"] = (0.06,["ggHmumu125","vbfHmumu125","wHmumu125","zHmumu125"])
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
      if channel.dataCountsTotal == None:
        print("Writing Pretend Data Counts")
        observationFormatList.append(int(channel.getBakXSecTotal()*lumi))
      else:
        print("Writing Real Data Counts")
        observationFormatList.append(channel.dataCountsTotal)
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

    #Debugging
    outfile.write("#################################\n")
    for channel,channelName in zip(self.channels,self.channelNames):
        outfile.write("#\n")
        outfile.write("#info: channel {0}: \n".format(channelName))
        for sigName in self.sigNames:
          outfile.write("#  {0} XS, Eff: {1:.3g}, {2:.3%} \n".format(sigName,
                        channel.getSigXSec(sigName), channel.getSigEff(sigName)))
        outfile.write("#   Mass: {0} \n".format(channel.massBounds))
    outfile.close()

###################################################################################

class ShapeDataCardMaker(DataCardMaker):
  def __init__(self,directory,analysisNames,signalNames,backgroundNames,dataNames,nuisanceMap=None,histNameBase="",rebin=[],useTH1=False,controlRegionLow=[80,115],controlRegionHigh=[135,200],controlRegionVeryLow=[],bakShape=False,histNameSuffix="",toyData=False):
    DataCardMaker.__init__(self,directory,analysisNames,signalNames,backgroundNames,dataNames,nuisanceMap,histNameBase,controlRegionLow,controlRegionHigh,controlRegionVeryLow=controlRegionVeryLow,bakShape=bakShape,rebin=rebin,histNameSuffix=histNameSuffix)

    self.useTH1 = useTH1
    self.controlRegionHigh = controlRegionHigh
    self.controlRegionLow = controlRegionLow
    self.controlRegionVeryLow = controlRegionVeryLow
    self.toyData = toyData

  def makeRFHistWrite(self,channel,hist,thisDir,isData=True,compareHist=None,writeBakShape=False):
    thisDir.cd()
    is2D = hist.InheritsFrom("TH2")
    if self.useTH1 and not is2D:
        hist.Write()
        return
    if not is2D:
      hist = shrinkTH1(hist,self.controlRegionLow[0],self.controlRegionHigh[1])
    x = channel.x
    origHist = hist
    if is2D:
      hist = hist2to1(hist)
      hist.SetName(re.sub("_1d","",hist.GetName()))
      if channel.x1d == None:
        channel.x1d = root.RooRealVar("x1d","x1d",0,hist.GetNbinsX()+2)
      x = channel.x1d
    rfHist = root.RooDataHist(hist.GetName(),hist.GetName(),root.RooArgList(root.RooArgSet(x)),hist)
    rfHistPdf = root.RooHistPdf(hist.GetName(),hist.GetName(),root.RooArgSet(x),rfHist)
    if isData:
      rfHist.Write()
    else:
      rfHistPdf.Write()
    debugDir = thisDir.FindObject('debug')
    if  debugDir==None:
      debugDir = thisDir.mkdir("debug")
    debugDir.cd()
    hist.Write()
    if is2D:
      xBinning = root.RooFit.Binning(origHist.GetNbinsX())
      yBinning = root.RooFit.Binning(origHist.GetNbinsY())
      x = channel.x
      y = channel.y
      rfHistTH2 = rfHist.createHistogram(origHist.GetName()+"rfHist2d",x,xBinning,root.RooFit.YVar(y,yBinning))
      rfHistPdfTH2 = rfHistPdf.createHistogram(origHist.GetName()+"rfHistPdf2d",x,xBinning,root.RooFit.YVar(y,yBinning))
      rfHistTH2.Write()
      rfHistPdfTH2.Write()
      hist.Write()
      if channel.bakShape and compareHist != None:
        channel.bakShapeMkr.writeDebugHistsToCurrTDir(compareHist)
    else:
      plot = x.frame()
      rfHist.plotOn(plot)
      rfHistPdf.plotOn(plot)
      plot.SetName(hist.GetName())
      plot.Write()
      if channel.bakShape and writeBakShape:
        channel.bakShapeMkr.writeDebugHistsToCurrTDir()

  def write(self,outfilename,lumi,sumAllBak=True,includeSigInAllMC=False):
    lumi *= 1000.0
    nuisance = self.nuisance

    ### ROOT Part
    ##########################################################
    outRootFilename = re.sub(r"\.txt",r".root",outfilename)
    outRootFile = root.TFile(outRootFilename, "RECREATE")
    outRootFile.cd()

    rootDebugString = ""

    observedN = {}

    for channel,channelName in zip(self.channels,self.channelNames):
        tmpDir = outRootFile.mkdir(channelName)
        tmpDir.cd()
        sumAllMCHist = None
        sumAllSigMCHist = None
        sumAllBakMCHist = None
        rootDebugString += "# channelName: {0}\n".format(channelName)
        for sigName in self.sigNames:
          tmpHist = channel.getSigHist(sigName).Clone(sigName)
          tmpHist.Scale(lumi)
          self.makeRFHistWrite(channel,tmpHist,tmpDir)
          #rootDebugString += "#     {0}: {1}\n".format(sigName,getIntegralAll(tmpHist))
          
          if includeSigInAllMC:
            if sumAllMCHist == None:
              sumAllMCHist = tmpHist.Clone("data_obs")
            else:
              sumAllMCHist.Add(tmpHist)
          if sumAllSigMCHist == None:
            sumAllSigMCHist = tmpHist.Clone("sig")
          else:
            sumAllSigMCHist.Add(tmpHist)
  
        if sumAllBak:
          #channel.dump()
          sumAllBakMCHist = channel.getBakHistTotal(lumi).Clone("bak")
          #channel.bakShapeMkr.dump()
          if self.shape:
            rootDebugString += channel.bakShapeMkr.debug
            for nuName in channel.bakShapeMkr.errNames:
                origUp = channel.bakShapeMkr.errHists[nuName+"Up"]
                origDown = channel.bakShapeMkr.errHists[nuName+"Down"]
                tmpUp = origUp.Clone("bak_"+nuName+"Up")
                tmpDown = origDown.Clone("bak_"+nuName+"Down")
                self.makeRFHistWrite(channel,tmpUp,tmpDir)
                self.makeRFHistWrite(channel,tmpDown,tmpDir)

          # for simulated data_obs:
          sumAllBakMCHistReal = channel.bakHistTotalReal.Clone("data_obs")
          sumAllBakMCHistReal.Scale(lumi)

          if sumAllMCHist == None:
            sumAllMCHist = sumAllBakMCHistReal
          else:
            sumAllMCHist.Add(sumAllBakMCHistReal)
        else:
          for bakName in self.bakNames:
            tmpHist = channel.getBakHist(bakName).Clone(bakName)
            tmpHist.Scale(lumi)
            self.makeRFHistWrite(channel,tmpHist,tmpDir)
            #rootDebugString += "#     {0}: {1}\n".format(bakName,getIntegralAll(tmpHist,boundaries=massBounds))
            
            if sumAllMCHist == None:
                sumAllMCHist = tmpHist.Clone("data_obs")
            else:
                sumAllMCHist.Add(tmpHist)
        massLimits = [self.controlRegionLow[0],self.controlRegionHigh[1]]
        sumAllMCHist.Scale(int(getIntegralAll(sumAllMCHist,boundaries=massLimits))/getIntegralAll(sumAllMCHist,boundaries=massLimits)) # Make Integer
        if channel.datHistTotal == None:
          if self.toyData:
            print("Writing Toy Data Histogram")
            toy = sumAllBakMCHist.Clone("data_obs")
            toyHistogram(toy)
            observedN[channelName] = getIntegralAll(toy,boundaries=massLimits)
            self.makeRFHistWrite(channel,toy,tmpDir) #Pretend Toy Data
          else:
            print("Writing Pretend Data Histogram")
            observedN[channelName] = getIntegralAll(sumAllMCHist,boundaries=massLimits)
            self.makeRFHistWrite(channel,sumAllMCHist,tmpDir) #Pretend Data
        else:
          print("Writing Real Data Histogram")
          observedN[channelName] = getIntegralAll(channel.datHistTotal,boundaries=massLimits)
          self.makeRFHistWrite(channel,channel.datHistTotal,tmpDir) #Real Data
        self.makeRFHistWrite(channel,sumAllSigMCHist,tmpDir) #Pretend Signal
        self.makeRFHistWrite(channel,sumAllBakMCHist,tmpDir,compareHist=sumAllSigMCHist,writeBakShape=True) #Background Sum
        rootDebugString += "#######################\n"
        rootDebugString += "#     {0}\n".format(channelName)
        rootDebugString += "#     Pretend Obs: {0}\n".format(getIntegralAll(sumAllMCHist,boundaries=massLimits))
        rootDebugString += "#     All Signal:  {0}\n".format(getIntegralAll(sumAllSigMCHist,boundaries=massLimits))
        rootDebugString += "#     All Bak:     {0}\n".format(getIntegralAll(sumAllBakMCHist,boundaries=massLimits))
        #rootDebugString += "#     Pretend Obs: {0}\n".format(getIntegralAll(sumAllMCHist))
        #rootDebugString += "#     All Signal:  {0}\n".format(getIntegralAll(sumAllSigMCHist))
        #rootDebugString += "#     All Bak:     {0}\n".format(getIntegralAll(sumAllBakMCHist))

    outRootFile.Close()

    ### Text Part
    ##########################################################

    print("Writing Card: {0} & {1}".format(outfilename,outRootFilename))
    outfile = open(outfilename,"w")
    outfile.write("# Hmumu shape combine datacard produced by makeTables.py\n")
    now = datetime.datetime.now().replace(microsecond=0).isoformat(' ')
    outfile.write("# {0}\n".format(now))
    outfile.write("############################### \n")
    outfile.write("############################### \n")
    outfile.write("imax {0}\n".format(len(self.channels)))
    #outfile.write("jmax {0}\n".format(len(backgroundNames)))
    outfile.write("jmax {0}\n".format("*"))
    outfile.write("kmax {0}\n".format("*"))
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
      observedNumber = observedN[channelName]
      observationFormatList.append(observedNumber)
      #print("text Observed {}: {}".format(channelName,observedNumber))
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

        if sumAllBak:

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

    # Bak Shape Uncertainties (All Correlated)
    for channel,channelName in zip(self.channels,self.channelNames):
      for nuisanceName in channel.bakShapeMkr.errNames:
        formatString = "{0:<8} {1:^4} "
        formatList = [nuisanceName,"shape"]
        iParam = 2
        for channel2,channelName2 in zip(self.channels,self.channelNames):
          for sigName in self.sigNames:
            formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            value = "-"
            formatList.append(value)
            iParam += 1
          if sumAllBak:
              bakName="bak"
              formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
              value = "1"
              formatList.append(value)
              iParam += 1
          else:
            for bakName in self.bakNames:
              formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
              value = "-"
              formatList.append(value)
              iParam += 1
        formatString += "\n"
        #print formatString
        #print formatList
        outfile.write(formatString.format(*formatList))
      break


    #Debugging
    outfile.write("#################################\n")
    for channel,channelName in zip(self.channels,self.channelNames):
        outfile.write("#\n")
        outfile.write("#info: channel {0}: \n".format(channelName))
        outfile.write("#  x var name: {0} \n".format(channel.x.GetName()))
        outfile.write("#  x var range: [{0:.3g},{1:.3g}] \n".format(channel.x.getMin(),channel.x.getMax()))
        if channel.is2D:
          outfile.write("#  y var name: {0} \n".format(channel.y.GetName()))
          outfile.write("#  y var range: [{0:.3g},{1:.3g}] \n".format(channel.y.getMin(),channel.y.getMax()))
    outfile.write(rootDebugString)
    outfile.close()

class ThreadedCardMaker(myThread):
  def __init__(self,*args,**dictArgs):
    myThread.__init__(self)
    self.writeArgs = (dictArgs["outfilename"],dictArgs["lumi"])
    self.writeArgsDict = {}
    if dictArgs.has_key("sumAllBak"):
        self.writeArgsDict["sumAllBak"] = dictArgs["sumAllBak"]
    if dictArgs.has_key("includeSigInAllMC"):
        self.writeArgsDict["includeSigInAllMC"] = dictArgs["includeSigInAllMC"]
    self.shapeDataCardMaker = True
    if dictArgs.has_key("shapeDataCardMaker"):
        self.shapeDataCardMaker = dictArgs["shapeDataCardMaker"]
        dictArgs.pop("shapeDataCardMaker",None)
    self.args = args
    dictArgs.pop("sumAllBak",None)
    dictArgs.pop("includeSigInallMC",None)
    dictArgs.pop("outfilename",None)
    dictArgs.pop("lumi",None)
    self.dictArgs = dictArgs
    self.started = False
  def run(self):
    self.started = True
    dataCardMassShape = None
    if self.shapeDataCardMaker:
      dataCardMassShape = ShapeDataCardMaker(*(self.args),**(self.dictArgs))
    else:
      dataCardMassShape = DataCardMaker(*(self.args),**(self.dictArgs))
    dataCardMassShape.write(*(self.writeArgs),**(self.writeArgsDict))

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
  analyses = ["VBFPresel","IncPresel","VBFLoose","VBFMedium","VBFTight","VBFVeryTight","Pt0to30","Pt30to50","Pt50to125","Pt125to250","Pt250","IncBDTSig80","VBFBDTSig80"]
  histPostFix="/mDiMu"
  #analyses = ["mDiMu"]
  #histPostFix=""
  signalNames=["ggHmumu125","vbfHmumu125","wHmumu125","zHmumu125"]
  backgroundNames= ["DYToMuMu","ttbar","WW","WZ","ZZ"]
  dataNames=[]
  #dataNames=["SingleMuRun2012Av1.root","SingleMuRun2012Bv1.root","SingleMuRun2012Cv1.root"]
  #lumiList = [5,10,15,20,25,30,40,50,75,100,200,500,1000]
  lumiList = [10,20,30,100]
  lumiList = [20]

  MassRebin = 4 # 4 Bins per GeV originally
  controlRegionVeryLow=[80,110]
  controlRegionLow=[110,120]
  controlRegionHigh=[130,180]

  shape=True
  toyData=False

  print("Creating Threads...")
  threads = []
  for i in lumiList:
    threads.append(
      ThreadedCardMaker(
        #__init__ args:
        directory,["VBFLoose","VBFMedium","VBFTight","VBFVeryTight","Pt0to30","Pt30to50","Pt50to125","Pt125to250","Pt250"],signalNames,backgroundNames,dataNames,
        rebin=[MassRebin], bakShape=shape,
        controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix=histPostFix,
        controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,
        #write args:
        outfilename=outDir+"AllCat"+"_"+str(i)+".txt",lumi=i
      )
    )
    threads.append(
      ThreadedCardMaker(
        #__init__ args:
        directory,["VBFLoose","VBFMedium","VBFTight","VBFVeryTight"],signalNames,backgroundNames,dataNames,
        rebin=[MassRebin], bakShape=shape,
        controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix=histPostFix,
        controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,
        #write args:
        outfilename=outDir+"VBFCat"+"_"+str(i)+".txt",lumi=i
      )
    )
    threads.append(
      ThreadedCardMaker(
        #__init__ args:
        directory,["Pt0to30","Pt30to50","Pt50to125","Pt125to250","Pt250"],signalNames,backgroundNames,dataNames,
        rebin=[MassRebin], bakShape=shape,
        controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix=histPostFix,
        controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,
        #write args:
        outfilename=outDir+"IncCat"+"_"+str(i)+".txt",lumi=i
      )
    )
    threads.append(
      ThreadedCardMaker(
        #__init__ args:
        directory,["IncBDTSig80","VBFBDTSig80"],signalNames,backgroundNames,dataNames,
        rebin=[MassRebin], bakShape=shape,
        controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix=histPostFix,
        controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,
        #write args:
        outfilename=outDir+"BDTSig80"+"_"+str(i)+".txt",lumi=i
      )
    )
    for ana in analyses:
      tmp = ThreadedCardMaker(
        #__init__ args:
        directory,[ana],signalNames,backgroundNames,dataNames,rebin=[MassRebin],bakShape=shape,
        controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix=histPostFix,
        controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,
        #write args:
        outfilename=outDir+ana+"_"+str(i)+".txt",lumi=i
        )
      threads.append(tmp)
#    ## Cut and Count!!!
#    for ana in analyses:
#      tmp = ThreadedCardMaker(
#        #__init__ args:
#        directory,[ana],signalNames,backgroundNames,dataNames,rebin=[MassRebin],bakShape=shape,
#        controlRegionLow=[123.0,125.0],controlRegionHigh=[125.0,127.0],histNameSuffix=histPostFix,
#        controlRegionVeryLow=controlRegionVeryLow,
#        #write args:
#        outfilename=outDir+"CNC_"+ana+"_"+str(i)+".txt",lumi=i,
#
#        shapeDataCardMaker=False
#        )
#      threads.append(tmp)

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
