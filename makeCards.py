#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description="Makes cards for use in the CMS Combine tool.")
parser.add_argument("--signalInject", help="Inject Signal with Strength into data_obs",type=float,default=0.0)
parser.add_argument("--signalInjectMass", help="Mass For Injected Signal",type=float,default=125.0)
parser.add_argument("--toyData", help="Make Toy Data from PDFs for data_obs",action="store_true",default=False)
parser.add_argument("--bdtCut", help="Creates Cards with different BDT Cuts",action="store_true",default=False)
parser.add_argument("--gaussian", help="Use A Gaussian Signal Template with floating width",type=float,default=-1.0)
parser.add_argument("-m","--higgsMass", help="Use This Higgs Mass",type=float,default=-1.0)
args = parser.parse_args()

import math
import ROOT as root
from helpers import *
import datetime
import sys
import os.path
import copy
import random
import shutil
import multiprocessing
import time
import string
myThread = multiprocessing.Process

from ROOT import gSystem
gSystem.Load('libRooFit')

root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

NPROCS = 1

#Scaling Parameter for Bak norm uncertainty
BAKUNC = 1.0

BAKUNCON = True
SIGUNCON = False

FREEBAKPARAMS = False
LIMITTOSIGNALREGION = False

SIGNALFIT = [110.,140.]

if args.bdtCut:
  SIGUNCON = False

from xsec import *

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

def makePDFBakExpLog(name,hist,mMuMu,minMass,maxMass,workspaceImportFn):
    debug = ""
    debug += "### makePDFBakExpLog: "+name+"\n"

    channelName = name

    p1 = root.RooRealVar(channelName+"_p1","p1", 0.0, -1., 1.)
    p2 = root.RooRealVar(channelName+"_p2","p2", 0.0, -1., 1.)
    p3 = root.RooRealVar(channelName+"_p3","p3", 0.0, -1., 1.)
    pdfMmumu = root.RooGenericPdf("bak","TMath::Exp(@0*@0*@1 + @0*@2 + @3*TMath::Log(@0) )",root.RooArgList(mMuMu,p1,p2,p3))

    mMuMuRooDataHist = root.RooDataHist("bak_Template","bak_Template",root.RooArgList(mMuMu),hist)

    fr = pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    chi2 = pdfMmumu.createChi2(mMuMuRooDataHist)

    rooParamList = [p1,p2,p3]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    for param in rooParamList:
        param.setConstant(not FREEBAKPARAMS)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(mMuMuRooDataHist)
      workspaceImportFn(fr)

    #Norm Time
    bakNormTup = None
    if LIMITTOSIGNALREGION:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(mMuMu,"signal")
      getSidebandString = "mMuMu < {0} || mMuMu > {1}".format(*signalRangeList)
      getSignalString = "mMuMu > {0} && mMuMu < {1}".format(*signalRangeList)
      nSideband =  mMuMuRooDataHist.sumEntries(getSidebandString)
      nData =  mMuMuRooDataHist.sumEntries(getSignalString)
      signalFraction = signalIntegral.getVal()/wholeIntegral.getVal()
      bakNormTup = (nSideband,(signalFraction)/(1.0-signalFraction))
      print("N_side: {0:.2f}, alpha: {1:.2f}".format(bakNormTup[0],bakNormTup[1]))
      print("Signal Region: {0:.2f} Prediction: {1:.2f}".format(nData,bakNormTup[0]*bakNormTup[1]))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:%} true error: {2:%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
      mMuMu.Print()
    else:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(mMuMu,"signal")
      getSidebandString = "mMuMu < {0} || mMuMu > {1}".format(*signalRangeList)
      nSideband =  mMuMuRooDataHist.sumEntries(getSidebandString)
      nData =  mMuMuRooDataHist.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
      mMuMu.Print()
    #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #mMuMuRooDataHist2 = mMuMuRooDataHist.reduce(root.RooFit.CutRange("low,signal,high"))
    #mMuMuRooDataHist2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(mMuMuRooDataHist2)

    ## Debug Time
    frame = mMuMu.frame()
    frame.SetName("bak_Plot")
    mMuMuRooDataHist.plotOn(frame)
    pdfMmumu.plotOn(frame)
    canvas = root.TCanvas()
    frame.Draw()
    canvas.SaveAs("debug_"+name+channelName+".png")

    return paramList, bakNormTup


def makePDFBakExpMOverSq(name,hist,mMuMu,minMass,maxMass,workspaceImportFn):
    debug = ""
    debug += "### makePDFBakExpMOverSq: "+name+"\n"

    channelName = name

    InvPolMass = root.RooRealVar(channelName+"_InvPolMass","InvPolMass", 91., 70., 95.)
    ExpMass = root.RooRealVar(channelName+"_ExpMass","ExpMass", 0.0, -1., 1.)
    pdfMmumu = root.RooGenericPdf("bak","TMath::Exp(@0*@2)/(@0-@1)/(@0-@1)",root.RooArgList(mMuMu,InvPolMass,ExpMass))

    mMuMuRooDataHist = root.RooDataHist("bak_Template","bak_Template",root.RooArgList(mMuMu),hist)

    fr = pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    chi2 = pdfMmumu.createChi2(mMuMuRooDataHist)

    rooParamList = [InvPolMass,ExpMass]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    for param in rooParamList:
        param.setConstant(not FREEBAKPARAMS)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(mMuMuRooDataHist)
      workspaceImportFn(fr)

    #Norm Time
    bakNormTup = None
    if LIMITTOSIGNALREGION:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(mMuMu,"signal")
      getSidebandString = "mMuMu < {0} || mMuMu > {1}".format(*signalRangeList)
      getSignalString = "mMuMu > {0} && mMuMu < {1}".format(*signalRangeList)
      nSideband =  mMuMuRooDataHist.sumEntries(getSidebandString)
      nData =  mMuMuRooDataHist.sumEntries(getSignalString)
      signalFraction = signalIntegral.getVal()/wholeIntegral.getVal()
      bakNormTup = (nSideband,(signalFraction)/(1.0-signalFraction))
      print("N_side: {0:.2f}, alpha: {1:.2f}".format(bakNormTup[0],bakNormTup[1]))
      print("Signal Region: {0:.2f} Prediction: {1:.2f}".format(nData,bakNormTup[0]*bakNormTup[1]))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:%} true error: {2:%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
      mMuMu.Print()
    else:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(mMuMu,"signal")
      getSidebandString = "mMuMu < {0} || mMuMu > {1}".format(*signalRangeList)
      nSideband =  mMuMuRooDataHist.sumEntries(getSidebandString)
      nData =  mMuMuRooDataHist.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
      mMuMu.Print()
   #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #mMuMuRooDataHist2 = mMuMuRooDataHist.reduce(root.RooFit.CutRange("low,signal,high"))
    #mMuMuRooDataHist2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(mMuMuRooDataHist2)

    ## Debug Time
    frame = mMuMu.frame()
    frame.SetName("bak_Plot")
    mMuMuRooDataHist.plotOn(frame)
    pdfMmumu.plotOn(frame)
    canvas = root.TCanvas()
    frame.Draw()
    canvas.SaveAs("debug_"+name+channelName+".png")

    return paramList, bakNormTup

def makePDFBakMOverSq(name,hist,mMuMu,minMass,maxMass,workspaceImportFn):
    debug = ""
    debug += "### makePDFBakMOverSq: "+name+"\n"

    channelName = name

    InvPolMass = root.RooRealVar(channelName+"_InvPolMass","InvPolMass", 91., 70., 95.)
    pdfMmumu = root.RooGenericPdf("bak","@0/(@0-@1)/(@0-@1)",root.RooArgList(mMuMu,InvPolMass))

    mMuMuRooDataHist = root.RooDataHist("bak_Template","bak_Template",root.RooArgList(mMuMu),hist)

    fr = pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    chi2 = pdfMmumu.createChi2(mMuMuRooDataHist)

    rooParamList = [InvPolMass]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    for param in rooParamList:
        param.setConstant(False)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(mMuMuRooDataHist)
      workspaceImportFn(fr)

    #Norm Time
    bakNormTup = None
    if LIMITTOSIGNALREGION:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(mMuMu,"signal")
      getSidebandString = "mMuMu < {0} || mMuMu > {1}".format(*signalRangeList)
      getSignalString = "mMuMu > {0} && mMuMu < {1}".format(*signalRangeList)
      nSideband =  mMuMuRooDataHist.sumEntries(getSidebandString)
      nData =  mMuMuRooDataHist.sumEntries(getSignalString)
      signalFraction = signalIntegral.getVal()/wholeIntegral.getVal()
      bakNormTup = (nSideband,(signalFraction)/(1.0-signalFraction))
      print("N_side: {0:.2f}, alpha: {1:.2f}".format(bakNormTup[0],bakNormTup[1]))
      print("Signal Region: {0:.2f} Prediction: {1:.2f}".format(nData,bakNormTup[0]*bakNormTup[1]))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:%} true error: {2:%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
      mMuMu.Print()
    else:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(mMuMu,"signal")
      getSidebandString = "mMuMu < {0} || mMuMu > {1}".format(*signalRangeList)
      nSideband =  mMuMuRooDataHist.sumEntries(getSidebandString)
      nData =  mMuMuRooDataHist.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
      mMuMu.Print()
    #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #mMuMuRooDataHist2 = mMuMuRooDataHist.reduce(root.RooFit.CutRange("low,signal,high"))
    #mMuMuRooDataHist2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(mMuMuRooDataHist2)

    ## Debug Time
    frame = mMuMu.frame()
    frame.SetName("bak_Plot")
    mMuMuRooDataHist.plotOn(frame)
    pdfMmumu.plotOn(frame)
    canvas = root.TCanvas()
    frame.Draw()
    canvas.SaveAs("debug_"+name+channelName+".png")

    return paramList, bakNormTup



def makePDFBakOld(name,hist,mMuMu,minMass,maxMass,workspaceImportFn):
    debug = ""
    debug += "### makePDFBakOld: "+name+"\n"

    channelName = name

    voitWidth = root.RooRealVar(channelName+"_voitWidth","voitWidth",2.4952)
    voitmZ = root.RooRealVar(channelName+"_voitmZ","voitmZ",85,95)
    voitSig = root.RooRealVar(channelName+"_voitSig","voitSig",0.0,30.0)
    voitMmumu = root.RooVoigtian("bak_voitMmumu","voitMmumu",mMuMu,voitmZ,voitWidth,voitSig)

    expParam = root.RooRealVar(channelName+"_expParam","expParam",-1,0)
    expMmumu = root.RooExponential("bak_expMmumu","expMmumu",mMuMu,expParam)

    mixParam = root.RooRealVar(channelName+"_mixParam","mixParam",0,1)

    pdfMmumu = root.RooAddPdf("bak","bak",root.RooArgList(voitMmumu,expMmumu),root.RooArgList(mixParam))

    # Just For Z-Peak Part

    mMuMuZ = root.RooRealVar("mMuMu","mMuMu",88.0,94.0)
    voitMmumuZ = root.RooVoigtian("bak_voitMmumuZ","voitMmumuZ",mMuMuZ,voitmZ,voitWidth,voitSig)
    mMuMuZRooDataHist = root.RooDataHist("bak_TemplateZ","bak_TemplateZ",root.RooArgList(mMuMuZ),hist)

    voitMmumuZ.fitTo(mMuMuZRooDataHist,root.RooFit.SumW2Error(False),PRINTLEVEL)
    voitmZ.setConstant(True)
    voitSig.setConstant(True)

#    ## Debug Time
#    frameZ = mMuMuZ.frame()
#    frameZ.SetName("bak_PlotZ")
#    mMuMuZRooDataHist.plotOn(frameZ)
#    voitMmumuZ.plotOn(frameZ)
#    canvas = root.TCanvas()
#    frameZ.Draw()
#    saveAs(canvas,"debug_bakZ")

    # Back to everywhere else

    mMuMuRooDataHist = root.RooDataHist("bak_Template","bak_Template",root.RooArgList(mMuMu),hist)

    expMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("high"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    expParam.setConstant(True)
    
    fr = pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    chi2 = pdfMmumu.createChi2(mMuMuRooDataHist)

    rooParamList = [voitmZ,voitSig,expParam,mixParam]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    if FREEBAKPARAMS:
      for param in rooParamList:
        param.setConstant(False)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(mMuMuRooDataHist)
      workspaceImportFn(fr)

    #Norm Time
    bakNormTup = None
    if LIMITTOSIGNALREGION:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(mMuMu,"signal")
      getSidebandString = "mMuMu < {0} || mMuMu > {1}".format(*signalRangeList)
      getSignalString = "mMuMu > {0} && mMuMu < {1}".format(*signalRangeList)
      nSideband =  mMuMuRooDataHist.sumEntries(getSidebandString)
      nData =  mMuMuRooDataHist.sumEntries(getSignalString)
      signalFraction = signalIntegral.getVal()/wholeIntegral.getVal()
      bakNormTup = (nSideband,(signalFraction)/(1.0-signalFraction))
      print("N_side: {0:.2f}, alpha: {1:.2f}".format(bakNormTup[0],bakNormTup[1]))
      print("Signal Region: {0:.2f} Prediction: {1:.2f}".format(nData,bakNormTup[0]*bakNormTup[1]))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:%} true error: {2:%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
      mMuMu.Print()
    else:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(mMuMu,"signal")
      getSidebandString = "mMuMu < {0} || mMuMu > {1}".format(*signalRangeList)
      nSideband =  mMuMuRooDataHist.sumEntries(getSidebandString)
      nData =  mMuMuRooDataHist.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
      mMuMu.Print()
    #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #mMuMuRooDataHist2 = mMuMuRooDataHist.reduce(root.RooFit.CutRange("low,signal,high"))
    #mMuMuRooDataHist2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(mMuMuRooDataHist2)

    ## Debug Time
    frame = mMuMu.frame()
    frame.SetName("bak_Plot")
    mMuMuRooDataHist.plotOn(frame)
    pdfMmumu.plotOn(frame)
    canvas = root.TCanvas()
    frame.Draw()
    canvas.SaveAs("debug_"+name+channelName+".png")

    return paramList, bakNormTup

def makePDFSigCBPlusGaus(name,hist,mMuMu,minMass,maxMass,workspaceImportFn,channelName,forceMean=-1.,sigInject=0):

    debug = ""
    debug += "### makePDFSigCBPlusGaus: "+channelName+": "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,mMuMu.GetName(),maxMass)

    mean = root.RooRealVar(channelName+"_"+name+"_Mean",channelName+"_"+name+"_Mean",125.,100.,150.)
    width = root.RooRealVar(channelName+"_"+name+"_Width",channelName+"_"+name+"_Width",5.0,0.5,20.0)
    width2 = root.RooRealVar(channelName+"_"+name+"_Width2",channelName+"_"+name+"_Width2",5.0,0.1,20.0)
    alpha = root.RooRealVar(channelName+"_"+name+"_Alpha",channelName+"_"+name+"_Alpha",1.0,0.1,10.0)
    n = root.RooRealVar(channelName+"_"+name+"_n",channelName+"_"+name+"_n",1.0,0.1,10.0)
    mix = root.RooRealVar(channelName+"_"+name+"_mix",channelName+"_"+name+"_mix",0.5,0.0,1.0)
    cb = root.RooCBShape(name+"_CB",name+"_CB",mMuMu,mean,width,alpha,n)
    gaus = root.RooGaussian(name+"_Gaus",name+"_Gaus",mMuMu,mean,width2)
    pdfMmumu = root.RooAddPdf(name,name,cb,gaus,mix)
    
    mMuMuRooDataHist = root.RooDataHist(name+"_Template",channelName+"_"+name+"_Template",root.RooArgList(mMuMu),hist)

    fr = pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True),root.RooFit.Range("signalfit"))
    fr.SetName(name+"_fitResult")

    #mean.setVal(125.)
    #width.setVal(2.)
    #alpha.setVal(1.4)
    #n.setVal(2.2)

    if forceMean > 0.0:
        mean.setVal(forceMean)

    ## Error time

    rooParamList = [mean,width,width2,alpha,n,mix]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]
    for i in rooParamList:
       i.setConstant(True)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(mMuMuRooDataHist)
      workspaceImportFn(fr)

    ## Debug Time
#    frame = mMuMu.frame()
#    frame.SetName(name+"_Plot")
#    mMuMuRooDataHist.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    canvas.SaveAs("debug_"+name+".png")

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())

    sigInjectDataset = None
    if sigInject > 0:
      sigInjectDataset = pdfMmumu.generate(root.RooArgSet(mMuMu),sigInject,False)

    return paramList, debug, sigInjectDataset

def makePDFSigCB(name,hist,mMuMu,minMass,maxMass,workspaceImportFn,channelName,forceMean=-1.,sigInject=0):

    debug = ""
    debug += "### makePDFSigCB: "+channelName+": "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,mMuMu.GetName(),maxMass)

    mean = root.RooRealVar(channelName+"_"+name+"_Mean",channelName+"_"+name+"_Mean",125.,100.,150.)
    width = root.RooRealVar(channelName+"_"+name+"_Width",channelName+"_"+name+"_Width",5.0,0.5,20.0)
    alpha = root.RooRealVar(channelName+"_"+name+"_Alpha",channelName+"_"+name+"_Alpha",1.0,0.1,10.0)
    n = root.RooRealVar(channelName+"_"+name+"_n",channelName+"_"+name+"_n",1.0,0.1,10.0)
    #pdfMmumu = root.RooGaussian(name,name,mMuMu,mean,width)
    pdfMmumu = root.RooCBShape(name,name,mMuMu,mean,width,alpha,n)
    
    mMuMuRooDataHist = root.RooDataHist(name+"_Template",channelName+"_"+name+"_Template",root.RooArgList(mMuMu),hist)

    fr = pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True),root.RooFit.Range("signalfit"))
    fr.SetName(name+"_fitResult")

    #mean.setVal(125.)
    #width.setVal(2.)
    #alpha.setVal(1.4)
    #n.setVal(2.2)

    if forceMean > 0.0:
        mean.setVal(forceMean)

    ## Error time

    rooParamList = [mean,width,alpha,n]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]
    for i in rooParamList:
       i.setConstant(True)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(mMuMuRooDataHist)
      workspaceImportFn(fr)


    ## Debug Time
#    frame = mMuMu.frame()
#    frame.SetName(name+"_Plot")
#    mMuMuRooDataHist.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    canvas.SaveAs("debug_"+name+".png")

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())

    sigInjectDataset = None
    if sigInject > 0:
      sigInjectDataset = pdfMmumu.generate(root.RooArgSet(mMuMu),sigInject,False)

    return paramList, debug, sigInjectDataset

def makePDFSigGaus(name,hist,mMuMu,minMass,maxMass,workspaceImportFn,channelName,forceMean=-1.,forceWidth=-1.,sigInject=0):

    debug = ""
    debug += "### makePDFSigGaus: "+channelName+": "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,mMuMu.GetName(),maxMass)

    mean = root.RooRealVar(channelName+"_"+name+"_Mean",channelName+"_"+name+"_Mean",125.,100.,150.)
    width = root.RooRealVar(channelName+"_"+name+"_Width",channelName+"_"+name+"_Width",5.0,0.1,20.0)
    pdfMmumu = root.RooGaussian(name,name,mMuMu,mean,width)
    
    mMuMuRooDataHist = root.RooDataHist(name+"_Template",channelName+"_"+name+"_Template",root.RooArgList(mMuMu),hist)
    #debug += "#    DatasetInfo: \n"
    #debug += "#      sumEntries: {0:.2f} nEntries: {1}\n".format(mMuMuRooDataHist.sumEntries(),mMuMuRooDataHist.numEntries())
    #debug += "#      Mean: {0:.2f} Sigma: {1:.2f}\n".format(mMuMuRooDataHist.mean(mMuMu),mMuMuRooDataHist.sigma(mMuMu))

    fr = pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True),root.RooFit.Range("signalfit"))
    fr.SetName(name+"_fitResult")

    #mean.setVal(125.)
    #width.setVal(2.)

    #mMuMu.Print()
    #pdfMmumu.Print()
    #hist.Print()
    #print hist.GetNbinsX()

    if forceMean > 0.0:
        mean.setVal(forceMean)
    if forceWidth > 0.0:
        width.setVal(forceWidth)

    ## Error time

    rooParamList = [mean,width]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]
    for i in rooParamList:
       i.setConstant(True)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(mMuMuRooDataHist)
      workspaceImportFn(fr)

    ### Debug Time
#    frame = mMuMu.frame()
#    frame.SetName(name+"_Plot")
#    mMuMuRooDataHist.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    canvas.SaveAs("debug_"+name+".png")

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())

    sigInjectDataset = None
    if sigInject > 0:
      sigInjectDataset = pdfMmumu.generate(root.RooArgSet(mMuMu),sigInject,False)

    return paramList, debug, sigInjectDataset

makePDFSig = makePDFSigCBPlusGaus
makePDFBak = makePDFBakOld
#makePDFBak = makePDFBakMOverSq
#makePDFBak = makePDFBakExpMOverSq
#makePDFBak = makePDFBakExpLog

###################################################################################

class Analysis:
  def __init__(self,directory,signalNames,backgroundNames,dataNames,analysis,lumi,controlRegionVeryLow,controlRegionLow,controlRegionHigh,histNameBase="mDiMu",rebin=[],histNameSuffix="",toyData=False,sigInject=0.0,sigInjectMass=125.0,bdtCut=None,energyStr="8TeV"):
    getCutHist = getattr(self,"getCutHist")
    doSigInject = getattr(self,"doSigInject")

    higgsPeakMean = args.higgsMass - 0.3

    self.lumi = lumi
    self.sigNames = signalNames
    self.bakNames = backgroundNames
    self.datNames = dataNames
    self.controlRegionVeryLow = controlRegionVeryLow
    self.controlRegionLow = controlRegionLow
    self.controlRegionHigh = controlRegionHigh
    self.analysis = analysis
    self.params = []
    self.debug = ""
    self.debug += "#Nominal Higgs Mass: "+str(args.higgsMass) +"\n"
    self.debug += "#Peak Centered at: "+str(higgsPeakMean) +"\n"

    self.workspace = root.RooWorkspace(analysis+energyStr)
    self.workspaceName = analysis+energyStr
    wImport = getattr(self.workspace,"import")
    self.sigInjectWorkspaces = []

    maxMass = controlRegionHigh[1]
    minMass = controlRegionLow[0]
    mMuMu = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)
    #mMuMu.setRange("z",88,94)
    #mMuMu.setRange("verylow",controlRegionVeryLow[0],controlRegionVeryLow[1])
    mMuMu.setRange("low",controlRegionLow[0],controlRegionLow[1])
    mMuMu.setRange("high",controlRegionHigh[0],controlRegionHigh[1])
    mMuMu.setRange("signal",controlRegionLow[1],controlRegionHigh[0])
    mMuMu.setRange("signalfit",SIGNALFIT[0],SIGNALFIT[1])
    self.mMuMu = mMuMu

    print("{} {} {}".format("low",controlRegionLow[0],controlRegionLow[1]))
    print("{} {} {}".format("high",controlRegionHigh[0],controlRegionHigh[1]))
    print("{} {} {}".format("signal",controlRegionLow[1],controlRegionHigh[0]))
    print("{} {} {}".format("signalfit",SIGNALFIT[0],SIGNALFIT[1]))

    self.sigFiles = []
    self.sigHistsRaw = []
    for name in signalNames:
      tmpF = root.TFile(directory+name+".root")
      tmpHLoc = histNameBase+analysis+histNameSuffix
      if tmpHLoc[0] == '/':
        tmpHLoc = tmpHLoc[1:]
      tmpH = tmpF.Get(tmpHLoc)
      if bdtCut != None:
        tmpH1Loc = histNameBase+analysis+"/mDiMu"
        if tmpH1Loc[0] == '/':
          tmpH1Loc = tmpH1Loc[1:]
        tmpH1 = tmpF.Get(tmpH1Loc)
        tmpH1.Reset()
        tmpH = self.getCutHist(tmpH,tmpH1,bdtCut)
        tmpH = tmpH1
      self.sigFiles.append(tmpF)
      self.sigHistsRaw.append(tmpH)

    self.bakFiles = []
    self.bakHistsRaw = []
    for name in backgroundNames:
      tmpF = root.TFile(directory+name+".root")
      tmpHLoc = histNameBase+analysis+histNameSuffix
      if tmpHLoc[0] == '/':
        tmpHLoc = tmpHLoc[1:]
      tmpH = tmpF.Get(tmpHLoc)
      if bdtCut != None:
        tmpH1Loc = histNameBase+analysis+"/mDiMu"
        if tmpH1Loc[0] == '/':
          tmpH1Loc = tmpH1Loc[1:]
        tmpH1 = tmpF.Get(tmpH1Loc)
        tmpH1.Reset()
        tmpH = self.getCutHist(tmpH,tmpH1,bdtCut)
        tmpH = tmpH1
      self.bakFiles.append(tmpF)
      self.bakHistsRaw.append(tmpH)

    self.datFiles = []
    self.datHists = []
    for name in dataNames:
      tmpF = root.TFile(directory+name+".root")
      tmpHLoc = histNameBase+analysis+histNameSuffix
      if tmpHLoc[0] == '/':
        tmpHLoc = tmpHLoc[1:]
      tmpH = tmpF.Get(tmpHLoc)
      if bdtCut != None:
        tmpH1Loc = histNameBase+analysis+"/mDiMu"
        if tmpH1Loc[0] == '/':
          tmpH1Loc = tmpH1Loc[1:]
        tmpH1 = tmpF.Get(tmpH1Loc)
        tmpH1.Reset()
        tmpH = self.getCutHist(tmpH,tmpH1,bdtCut)
        tmpH = tmpH1
      self.datFiles.append(tmpF)
      self.datHists.append(tmpH)

    # Signal Shape systematics
    self.sigErrHistsMap = {}
    if SIGUNCON:
      for f in self.sigFiles:
        name = histNameBase+analysis+"/mDiMu"
        tmpDir = None
        if name[0] == '/':
            tmpDir = f
        else:
          name = name.split("/")
          tmpDirKey = f.GetKey(name[0]) #Will break in main dir
          tmpDir = tmpDirKey.ReadObj()
        #tmpDir.Print()
        for key in tmpDir.GetListOfKeys():
          matchUp = re.match('mDiMu'+".+Up",key.GetName())
          matchDown = re.match('mDiMu'+".+Down",key.GetName())
          if matchUp or matchDown:
            self.sigErrHistsMap[re.sub('mDiMu',"",key.GetName())] = []
        break
      for f in self.sigFiles:
        for sysName in self.sigErrHistsMap:
          tmpHist = f.Get(histNameBase+analysis+histNameSuffix+sysName)
          self.sigErrHistsMap[sysName].append(tmpHist)

    #Rebin
    rb = rebin
    if type(rb) != list:
      print("Error: Analysis.rebin: argument must be a list!!  Exiting.")
      sys.exit(1)
    try:
      if len(rb) == 2 and False:
          for hist in self.sigHistsRaw:
            hist.Rebin2D(*rb)
          for hist in self.bakHistsRaw:
            hist.Rebin2D(*rb)
          for hist in self.datHists:
            hist.Rebin2D(*rb)
          for key in self.sigErrHistsMap:
            for hist in self.sigErrHistsMap[key]:
              hist.Rebin2D(*rb)
      elif len(rb) == 1 and True:
          for hist in self.sigHistsRaw:
            hist.Rebin(*rb)
          for hist in self.bakHistsRaw:
            hist.Rebin(*rb)
          for hist in self.datHists:
            hist.Rebin(*rb)
          for key in self.sigErrHistsMap:
            for hist in self.sigErrHistsMap[key]:
              hist.Rebin(*rb)
      elif len(rb) == 0:
        pass
      else:
        print("Error: Analysis.rebin: argument must be len 0, 1, or 2 list!!  Exiting.")
        print("  Must also be same length as dimension of hist, if not 0.")
        sys.exit(1)
    except Exception as e:
      print("Error: Analysis threw exception while rebinning data: {0}".format(e))
      print("  Analysis Name: {0}, directory: {1}, Rebin: {2}".format(analysis,directory,rb))
      print("  hist name: {0}".format(histNameBase+analysis+histNameSuffix))
      if bdtCut != None:
        print("  bdtCut hist name: {0}".format(histNameBase+analysis+"/mDiMu"))
      print("  signals: {0}".format(signalNames))
      print("  backgrounds: {0}".format(backgroundNames))
      print("  data: {0}".format(dataNames))
      sys.exit(1)
    effMap = {}
    xsecMap = {}
    lowBin = 0
    highBin = self.sigHistsRaw[0].GetNbinsX()+1
    #massBounds = [controlRegionLow[0],controlRegionHigh[1]]
    #massBounds = [controlRegionVeryLow[0],controlRegionHigh[1]]
    massBounds = [controlRegionLow[0],controlRegionHigh[1]]
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
        doSigInject(self.bakHistTotal,sigInject,sigInjectMass)
        histToUseForBak = self.bakHistTotal
    else:
        self.dataCountsTotal += doSigInject(self.datHistTotal,sigInject,sigInjectMass)

    tmpBakParams, self.bakNormTup = makePDFBak(analysis+energyStr,histToUseForBak,
                                        mMuMu, minMass,maxMass,wImport
                                       )
    self.params.extend(tmpBakParams)
    self.countsBakTotal = self.bakNormTup[0]*self.bakNormTup[1]

    if LIMITTOSIGNALREGION:
      minMass = controlRegionLow[1]
      maxMass = controlRegionHigh[0]
      if self.dataCountsTotal == None:
        self.bakHistTotal = shrinkTH1(self.bakHistTotal,minMass,maxMass)
      else:
        self.datHistTotal = shrinkTH1(self.datHistTotal,minMass,maxMass)
        self.dataCountsTotal = getIntegralAll(self.datHistTotal)
      mMuMu.setRange(minMass,maxMass)
      for hist in self.sigHistsRaw:
        vetoOutOfBoundsEvents(hist,[minMass,maxMass])
      for iSignal in range(len(signalNames)):
        for errName in self.sigErrHistsMap:
          hist = self.sigErrHistsMap[errName][iSignal]
          vetoOutOfBoundsEvents(hist,[minMass,maxMass])

    self.sigParamListList = []
    if args.gaussian > 0.:
      for name, hist in zip(signalNames,self.sigHistsRaw):
        sigParams, sigDebug, tmpDS = makePDFSigGaus(name,hist,mMuMu,minMass,maxMass,wImport,analysis+energyStr,forceMean=higgsPeakMean,forceWidth=args.gaussian)
        self.sigParamListList.append(sigParams)
        self.debug += sigDebug
    else:
      for name, hist in zip(signalNames,self.sigHistsRaw):
        sigParams, sigDebug, tmpDS = makePDFSig(name,hist,mMuMu,minMass,maxMass,wImport,analysis+energyStr,forceMean=higgsPeakMean)
        self.sigParamListList.append(sigParams)
        self.debug += sigDebug
        
    if args.gaussian > 0.:
      for i in self.sigParamListList:
        for j in i:
          if "Width" in j.name:
            j.lowErr = max(j.nominal/2.0,0.75)
            j.highErr = max(j.nominal/2.0,0.75)
            print("Adding sig param to list:")
            print j
            #self.params.append(j)
    elif SIGUNCON:
      for name, iSignal, paramsNoErr in zip(signalNames,
            range(len(signalNames)),self.sigParamListList):
        firstHist = True
        paramErrList = []
        for errName in self.sigErrHistsMap:
          nameNew = name+"_"+errName
          hist = self.sigErrHistsMap[errName][iSignal]
          sigParams, sigDebug, tmpDS = makePDFSig(nameNew,hist,mMuMu,minMass,maxMass,None,
                                                            analysis+energyStr,
                                                            forceMean=higgsPeakMean
                                                              )
          self.debug += sigDebug
          for curErr, nominal, i in zip(sigParams,paramsNoErr,range(len(sigParams))):
                val = curErr.nominal
                nomVal = nominal.nominal
                err = abs(val-nomVal)
                if firstHist:
                  paramErrList.append(err)
                else:
                  if err > paramErrList[i]:
                    paramErrList[i] = err
          firstHist = False
        for i in range(len(paramErrList)):
            paramsNoErr[i].highErr = paramErrList[i]
            paramsNoErr[i].lowErr = paramErrList[i]
      for i in self.sigParamListList:
        for j in i:
          if "Width2" in j.name:
            self.params.append(j)
      # To Make sure nothing got messed up!
      for name, hist in zip(signalNames,self.sigHistsRaw):
        sigParams, sigDebug, tmpDS = makePDFSig(name,hist,mMuMu,minMass,maxMass,None,analysis+energyStr)

    self.xsecSigTotal = 0.0
    self.xsecSigList = []
    self.effSigList = []
    self.sigHists = []
    for h,name in zip(self.sigHistsRaw,signalNames):
      counts = getIntegralAll(h,boundaries=massBounds)
      eff = counts/nEventsMap[name]*efficiencyMap[getPeriod(name)]
      xs = eff*xsec[name]
      if args.higgsMass > 0.0:
        prec = "0"
        if args.higgsMass % 1 > 0.0:
          prec = "1"
        higgsMassString =("{0:."+prec+"f}").format(args.higgsMass)
        tmpName = name.replace("125",higgsMassString)
        xs = eff*xsec[tmpName]
      self.xsecSigTotal += xs
      self.xsecSigList.append(xs)
      self.effSigList.append(eff)

    self.countsSigTotal = self.xsecSigTotal*lumi
    #self.countsBakTotal = self.xsecBakTotal*lumi
    self.countsSigList = [x*lumi for x in self.xsecSigList]
    self.countsBakList = [x*lumi for x in self.xsecBakList]

    if self.dataCountsTotal is None:
      bakDataTH1 = self.bakHistTotal.Clone("bak_Template"+str(random.randint(0,10000)))
      bakDataTH1 = shrinkTH1(bakDataTH1,minMass,maxMass)
      for i in range(bakDataTH1.GetNbinsX()+2):
        binCenter = bakDataTH1.GetXaxis().GetBinCenter(i)
        contentInt = int("{0:.0f}".format(bakDataTH1.GetBinContent(i)))
        bakDataTH1.SetBinContent(i,0.0)
        bakDataTH1.SetBinError(i,0.0)
        for j in range(contentInt):
          bakDataTH1.Fill(binCenter)
      self.dataCountsTotal = int(bakDataTH1.Integral())
      obsData = root.RooDataHist("data_obs","MC Full-Sim Data",root.RooArgList(mMuMu),bakDataTH1)
      print "counts: {} obsData: {}".format(self.dataCountsTotal,obsData.sumEntries())
      wImport(obsData)
      self.bakHistTotal.SetName("data_obs_"+analysis)
      #degubf = root.TFile("debug.root","recreate")
      #self.bakHistTotal.Write()
      #degubf.Close()
    elif toyData:
      bakPDF = self.workspace.pdf("bak")
      sigPDFList = [self.workspace.pdf(i) for i in signalNames]
      toyDataset = bakPDF.generate(root.RooArgSet(mMuMu),int(self.dataCountsTotal))
      doSigInject(toyDataset,sigInject,sigInjectMass)
      toyDataHist = toyDataset.binnedClone("data_obs","Toy Data")
      self.dataCountsTotal = int(toyDataHist.sumEntries())
      wImport(toyDataHist)
    else:
      realDataHist = root.RooDataHist("data_obs","Real Observed Data",root.RooArgList(mMuMu),self.datHistTotal)
      wImport(realDataHist)
      #realDataHistNotVeryLow = realDataHist.reduce(root.RooFit.CutRange("low,signal,high"))
      #wImport(realDataHistNotVeryLow)

  def getCutHist(self,inHist,outHist,cut):
    if not inHist.InheritsFrom("TH1"):
      print("Error: Analysis.getCutHist(): inHist is not a hist exiting.")
      sys.exit(1)
    if not inHist.InheritsFrom("TH2"):
      print("Error: Analysis.getCutHist(): inHist is not 2D!! exiting.")
      sys.exit(1)
    if not outHist.InheritsFrom("TH1"):
      print("Error: Analysis.getCutHist(): outHist is not a hist exiting.")
      sys.exit(1)
    if outHist.InheritsFrom("TH2"):
      print("Error: Analysis.getCutHist(): outHist is 2D!! exiting.")
      sys.exit(1)
    cutBin = inHist.GetYaxis().FindBin(cut)
    nBinsY = inHist.GetYaxis().GetNbins()
    nBinsX = inHist.GetXaxis().GetNbins()
    #print("max bin: {0}, cutBin: {1}".format(nBinsY,cutBin))
    for iX in range(0,nBinsX+2):
      mySum = 0.
      for iY in range(cutBin,nBinsY+2):
        mySum += inHist.GetBinContent(iX,iY)
      outHist.SetBinContent(iX,mySum)

  def doSigInject(self,dataHist,sigStrength,sigMass):
    if sigStrength <= 0.0:
        return 0
    mMuMu = self.mMuMu
    w = root.RooWorkspace("tmpWorkspace"+str(random.randint(0,10000)))
    self.sigInjectWorkspaces.append(w)
    wImport = getattr(w,"import")
    countsList = []
    for h,name in zip(self.sigHistsRaw,self.sigNames):
      counts = getIntegralAll(h,boundaries=self.massBounds)
      eff = counts/nEventsMap[name]*efficiencyMap[getPeriod(name)]
      xs = eff*xsec[name]
      if sigMass > 0.0:
        prec = "0"
        if sigMass % 1 > 0.0:
          prec = "1"
        higgsMassString =("{0:."+prec+"f}").format(sigMass)
        tmpName = name.replace("125",higgsMassString)
        xs = eff*xsec[tmpName]
      countsList.append(int(xs*self.lumi*sigStrength))
    for name, hist, counts in zip(self.sigNames,self.sigHistsRaw,countsList):
      if counts == 0:
        continue
      tmpName = name+"TmpSigInject"+str(random.randint(0,10000))
      sigParams, sigDebug, rooData = makePDFSig(tmpName,hist,mMuMu,sigMass-20,sigMass+20,wImport,tmpName,forceMean=sigMass-0.3,sigInject=counts)
      if dataHist.InheritsFrom("RooAbsData"):
        dataHist.append(rooData)
      else:
        tmpHist = dataHist.Clone("sigInjectHistTmp"+str(random.randint(0,10000)))
        tmpHist.Reset()
        rooData.fillHistogram(tmpHist,root.RooArgList(mMuMu))
        dataHist.Add(tmpHist)
    return sum(countsList)

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
  def __init__(self,directory,analysisNames,signalNames,backgroundNames,dataNames,outfilename,lumi,nuisanceMap=None,histNameBase="",controlRegionLow=[110.,115],controlRegionHigh=[135,150],controlRegionVeryLow=[80.,110.],rebin=[],histNameSuffix="",sigInject=0.0,sigInjectMass=125.0,toyData=False,bdtCut=None,energyStr="8TeV"):

    ########################
    ## Setup

    channels = []
    self.is2D = False

    if type(energyStr) == list:
      for analysis in analysisNames:
        for es,sn,bn,dn,lu in zip(energyStr,signalNames,backgroundNames,dataNames,lumi):
          lu *= 1000.0
          tmp = Analysis(directory,sn,bn,dn,analysis,lu,controlRegionVeryLow,controlRegionLow,controlRegionHigh,histNameBase=histNameBase,rebin=rebin,histNameSuffix=histNameSuffix,toyData=toyData,sigInject=sigInject,sigInjectMass=sigInjectMass,bdtCut=bdtCut,energyStr=es)
          channels.append(tmp)
      self.channels = channels
    else:
      lumi *= 1000.0
      for analysis in analysisNames:
        tmp = Analysis(directory,signalNames,backgroundNames,dataNames,analysis,lumi,controlRegionVeryLow,controlRegionLow,controlRegionHigh,histNameBase=histNameBase,rebin=rebin,histNameSuffix=histNameSuffix,toyData=toyData,sigInject=sigInject,sigInjectMass=sigInjectMass,bdtCut=bdtCut,energyStr=energyStr)
        channels.append(tmp)
      self.channels = channels

    self.channelNames = [i.workspaceName for i in channels]

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
        print(channel.workspace.data("data_obs").GetTitle())
        rootDebugString += "#"+channel.workspace.data("data_obs").GetTitle()+"\n"
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
        for sigName in channel.sigNames:
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
    for nu in nuisance.keys():
      formatString = "{0:<8} {1:^4} "
      formatList = [nu,"lnN"]
      iParam = 2
      for channel,channelName in zip(self.channels,self.channelNames):
          for sigName in channel.sigNames:
            formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            value = nuisance(nu,sigName)
            if value == None:
              value = "-"
            formatList.append(value)
            iParam += 1
          if True:
              bakName="bak"
              formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
              value = nuisance(nu,bakName)
              if value == None:
                value = "-"
              formatList.append(value)
              iParam += 1
          else:
            for bakName in self.bakNames:
              formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
              value = nuisance(nu,bakName)
              if value == None:
                value = "-"
              formatList.append(value)
              iParam += 1
      formatString += "\n"
      #print formatString
      #print formatList
      outfile.write(formatString.format(*formatList))

    # Bak norm uncertainty
    for channel1,channel1Name in zip(self.channels,self.channelNames):
      tmpTup = channel1.bakNormTup
      formatString = "{0:<8} {1:^3} {2:.0f}"
      formatList = ["bkN"+channel1Name,"gmN",tmpTup[0]]
      iParam = 3
      if FREEBAKPARAMS:
        formatString = "{0:<8} {1:^3}"
        formatList = ["bkN"+channel1Name,"lnU"]
        iParam = 2
      for channel in self.channels:
          for sigName in channel.sigNames:
            formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            value = "-"
            formatList.append(value)
            iParam += 1
          value = '-'
          tmpString = "{"+str(iParam)+":^"+str(self.largestChannelName)+"}"
          if channel == channel1:
            value = tmpTup[1]
            if FREEBAKPARAMS:
              value = 1.0/sqrt(tmpTup[1]*tmpTup[0])
              if value <= 0.01/20.:
                  value = 0.01
              else:
                  value *= 20
              value += 1.0
            tmpString = "{"+str(iParam)+":^"+str(self.largestChannelName)+".2f}"
          formatString += tmpString
          formatList.append(value)
          iParam += 1
      formatString += "\n"
      outfile.write(formatString.format(*formatList))

    if not FREEBAKPARAMS:
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
        outfile.write(channel.debug)
    outfile.write(rootDebugString)
    outfile.close()

class ThreadedCardMaker(myThread):
  def __init__(self,*args,**dictArgs):
    myThread.__init__(self)
    self.args = args
    self.dictArgs = dictArgs
    self.started = False
    #print("\nPos Args:\n {0} \nKeyword Args:\n {1}\n\n".format(self.args,self.dictArgs))
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

  directory = "input/preApproveSample/"
  outDir = "statsCards/"
  periods = ["7TeV","8TeV"]
  periods = ["8TeV"]
  analysesInc = ["IncPresel","IncBDTCut"]
  analysesVBF = ["VBFPresel","VBFBDTCut"]
  analyses = analysesInc + analysesVBF
  categoriesInc = ["BB","BO","BE","OO","OE","EE"]
  categoriesVBF = []
  tmpList = []
  for a in analysesInc:
    for c in categoriesInc:
        tmpList.append(a+c)
  analyses += tmpList
  tmpList = []
  for a in analysesVBF:
    for c in categoriesVBF:
        tmpList.append(a+c)
  analyses += tmpList
  analyses += ["IncPreselPtG10"+ x for x in categoriesInc]
  analyses = ["VBFBDTCut","IncPreselPtG10"]
  combinations = []
  combinationsLong = []
  combinations.append((
        ["IncPreselPtG10"+x for x in categoriesInc],"IncPreselCat"
  ))
  """
  combinations.append((
        ["VBFPresel"]+["IncPreselPtG10"],"BDTCutVBFBDTOnly"
  ))
  """
  combinations.append((
        ["VBFBDTCut"]+["IncPreselPtG10"+x for x in categoriesInc],"BDTCutCatVBFBDTOnly"
  ))

#  combinationsLong.append((
#        ["IncBDTCut","VBFBDTCut"],"BDTCut"
#  ))
#  combinationsLong.append((
#        ["VBFBDTCut"]+["IncBDTCut"+x for x in categoriesInc],"BDTCutCat"
#  ))
#  combinationsLong.append((
#        ["VBFPresel"]+["IncPresel"+x for x in categoriesInc],"PreselCat"
#  ))

  combinationsBDTCut = []
  #combinationsBDTCut.append((
  #  ["IncPresel"],"IncBDTCut",0.025,-0.7,-0.25,"BDTHistMuonOnlyVMass"
  #))
  #combinationsBDTCut.append((
  #  ["IncPresel"],"IncPtCut",1.0,0.0,20.0,"ptVmDiMu"
  #))
  combinationsBDTCut.append((
    ["VBFPresel"],"VBFBDTCut",0.01,-0.2,0.2,"BDTHistVBFVMass"
  ))
  #combinationsBDTCut.append((
  #  ["IncPresel"+x for x in categoriesInc],"IncBDTCutCat",0.025,-0.7,-0.35,"BDTHistMuonOnlyVMass"
  #))
  #combinationsBDTCut.append((
  #  ["VBFPresel"+x for x in categoriesVBF],"VBFBDTCutCat",0.025,-0.2,0.0,"BDTHistVBFVMass"
  #))

  histPostFix="/mDiMu"
  signalNames=["ggHmumu125","vbfHmumu125","wHmumu125","zHmumu125"]
  backgroundNames= ["DYJetsToLL","ttbar"]
  dataDict = {}
  dataDict["8TeV"] = [
    "SingleMuRun2012Av1",
    "SingleMuRun2012Av1Recover",
    "SingleMuRun2012Bv1",
    "SingleMuRun2012Cv1",
    "SingleMuRun2012Cv2",
    "SingleMuRun2012D",
  ]
  dataDict["7TeV"] = [
    "SingleMuRun2011Av1",
    "SingleMuRun2011Bv1"
  ]
  dataDict["14TeV"] = []
  lumiListLong = [5,10,15,20,25,30,40,50,75,100,200,500,1000,2000,5000]
  lumiListLong = [20,30,50,100,500,1000,5000]
  lumiList = [lumiDict["8TeV"],20,25,30]
  lumiList = [lumiDict["8TeV"]]
  #lumiListLong = lumiList

  MassRebin = 1 # 4 Bins per GeV originally
  controlRegionVeryLow=[60,110]
  controlRegionLow=[110,120]
  controlRegionHigh=[130,160]
  if args.higgsMass > 0.0:
    controlRegionLow=[110,args.higgsMass-5.0]
    controlRegionHigh=[args.higgsMass+5.0,160]

  shape=True
  toyData=args.toyData

  threads = []
  if not args.bdtCut:
    print("Simple Analyses to run:")
    for a in analyses:
      print("  {0}".format(a))
    print("Combination Analyses to run:")
    for c in combinations:
      print("  {0}".format(c[1]))
      for a in c[0]:
        print("    {0}".format(a))
  
    print("Creating Threads...")
    for p in periods:
      for i in lumiList:
        if p == "7TeV":
          i = lumiDict[p]
        name = i
        if args.higgsMass > 0.0:
          name = args.higgsMass
        for ana in analyses:
          tmp = ThreadedCardMaker(
            #__init__ args:
            directory,[ana],
            appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),dataDict[p],
            rebin=[MassRebin],
            controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix=histPostFix,
            controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,sigInjectMass=args.signalInjectMass,
            energyStr=p,
            #write args:
            outfilename=outDir+ana+"_"+p+"_"+str(name)+".txt",lumi=i
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
            controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,sigInjectMass=args.signalInjectMass,
            energyStr=p,
            #write args:
            outfilename=outDir+comb[1]+"_"+p+"_"+str(name)+".txt",lumi=i
          )
         )
        if p == "7TeV":
          break

    # 2011 + 2012 combination
    if len(periods)>1:
        filenameLumi = str(sum([float(lumiDict[p]) for p in periods]))
        if args.higgsMass > 0.0:
          filenameLumi = str(args.higgsMass)
        filenamePeriod = ""
        for p in periods:
          filenamePeriod += re.sub("TeV","P",p)
        filenamePeriod = filenamePeriod.rstrip("P")
        filenamePeriod += "TeV"
        for ana in analyses:
          tmp = ThreadedCardMaker(
            #__init__ args:
            directory,[ana],
            [appendPeriod(signalNames,p) for p in periods],
            [appendPeriod(backgroundNames,p) for p in periods],
            [dataDict[p] for p in periods],
            rebin=[MassRebin],
            controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix=histPostFix,
            controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,sigInjectMass=args.signalInjectMass,
            energyStr=periods,
            #write args:
            outfilename=outDir+ana+"_"+filenamePeriod+"_"+filenameLumi+".txt",lumi=[float(lumiDict[p]) for p in periods]
            )
          threads.append(tmp)
        for comb in combinations:
         threads.append(
          ThreadedCardMaker(
            #__init__ args:
            directory,
            comb[0],
            [appendPeriod(signalNames,p) for p in periods],
            [appendPeriod(backgroundNames,p) for p in periods],
            [dataDict[p] for p in periods],
            rebin=[MassRebin],
            controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix=histPostFix,
            controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,sigInjectMass=args.signalInjectMass,
            energyStr=periods,
            #write args:
            outfilename=outDir+comb[1]+"_"+filenamePeriod+"_"+filenameLumi+".txt",lumi=[float(lumiDict[p]) for p in periods]
          )
         )
  
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
            controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,sigInjectMass=args.signalInjectMass,
            energyStr=p,
            #write args:
            outfilename=outDir+comb[1]+"_"+p+"_"+str(i)+".txt",lumi=i
          )
         )
  else: #Do bdtCut Stuff
    for p in periods:
      i = lumiDict[p]
      for comb in combinationsBDTCut:
          bdtCutList = []
          tmpCut = float(comb[3])
          while True:
            if tmpCut > float(comb[4]):
                break
            bdtCutList.append(tmpCut)
            tmpCut += float(comb[2])
            if abs(tmpCut) < 1e-10:
              tmpCut = 0.0
          for bdtCutVal in bdtCutList:
            print("BDT Cut Card: {0} {1} {2} fb^-1 hist={3} cut={4}".format(p,comb[1],i,comb[5],bdtCutVal))
            threads.append(
             ThreadedCardMaker(
               #__init__ args:
               directory,
               comb[0],
               appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),[],
               rebin=[MassRebin],
               controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix="/"+comb[5],
               controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,sigInjectMass=args.signalInjectMass,
               energyStr=p,
               #write args:
               outfilename=outDir+comb[1]+str(i)+"_"+p+"_"+str(bdtCutVal)+".txt",lumi=i,
               bdtCut=bdtCutVal
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
bsub -o /dev/null lxbatch.sh $i
#bsub -q 1nh -o /dev/null lxbatch.sh $i
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

echo "executing combine -M ProfileLikelihood -d $FILENAME --signif >& $FILENAME.sig"

combine -M ProfileLikelihood -d $FILENAME --signif >& $FILENAME.sig
rm -f roostats*
rm -f higgsCombineTest*.root

echo "executing combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig"

combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 >& $FILENAME.expsig
##combine -M ProfileLikelihood -d $FILENAME --signif --expectSignal=1 -t -1 --toysFreq >& $FILENAME.expsig
rm -f roostats*
rm -f higgsCombineTest*.root

echo "executing combine -M MaxLikelihoodFit --plots --saveNormalizations $FILENAME >& $FILENAME.mu"

combine -M MaxLikelihoodFit --plots --saveNormalizations $FILENAME >& $FILENAME.mu
rm -f roostats*
rm -f higgsCombineTest*.root

cp $FILENAME.out ..
cp $FILENAME.mu ..
cp $FILENAME.sig ..
cp mlfit.root ../$FILENAME.root
cp data_fit_s.png ../$FILENAME.png
#cp $FILENAME.expsig ..

echo "done"
date
"""
  simpleBatchString = \
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
  #runFile.write(simpleBatchString)
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

echo "executing combine -M MaxLikelihoodFit --plots --saveNormalizations $FILENAME >& $FILENAME.mu"

combine -M MaxLikelihoodFit --plots --saveNormalizations $FILENAME >& $FILENAME.mu
rm -f roostats*
rm -f higgsCombineTest*.root
cp mlfit.root $FILENAME.root
cp data_fit_s.png $FILENAME.png

done

date
echo "done"
"""

  simpleBatchString = \
"""#!/bin/bash
echo "running notlxbatch.sh Simple Style"
date
for i in *.txt; do
    [[ -e "$i" ]] || continue
FILENAME=$i
echo "executing combine -M Asymptotic $FILENAME >& $FILENAME.out"

combine -M Asymptotic -v 2 $FILENAME >& $FILENAME.out
rm -f roostats*
rm -f higgsCombineTest*.root
done

date
echo "done"
"""

  runFile.write(batchString)
#  runFile.write(simpleBatchString)
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

  runFile = open(outDir+"getStatus2.sh","w")
  batchString = \
"""#!/bin/bash

echo "==========================="

NJOBS=`ls *.txt | wc -l`
STARTTIME=`date +%s`

while true; do
  NCOMPLETE="0"
  for i in *.out; do
    NTMP=`cat $i | wc -l`
    if [ "$NTMP" -gt "0" ]; then
        NCOMPLETE=$(( $NCOMPLETE + 1 ))
    fi  
  done
  NSTARTED=`ls *.out | wc -l`
  echo "`date --rfc-3339=seconds` Jobs: $NJOBS Started: $NSTARTED Complete: $NCOMPLETE"
  if [ "$NCOMPLETE" -eq "$NJOBS" ]; then
    ENDTIME=`date +%s`
    echo "Took $(( $ENDTIME - $STARTTIME )) seconds"
    echo "All Jobs Complete"
    echo "==========================="
    exit "0"
  fi
  sleep 10
done

"""
  runFile.write(batchString)
  runFile.close()
  shutil.copy("etc/nuisanceDiff.py",outDir+"diffNuisances.py")
  shutil.copy("etc/fitNormsToText_mlfit.py",outDir+"fitNormsToText_mlfit.py")
