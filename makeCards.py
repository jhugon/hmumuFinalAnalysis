#!/usr/bin/env python

import optparse
parser = optparse.OptionParser(description="Makes cards for use in the CMS Combine tool.")
parser.add_option("--signalInject", help="Inject Signal with Strength into data_obs",type=float,default=0.0)
parser.add_option("--signalInjectMass", help="Mass For Injected Signal",type=float,default=125.0)
parser.add_option("--toyData", help="Make Toy Data from PDFs for data_obs",action="store_true",default=False)
parser.add_option("--cutOpt", help="Creates Cards with Different Cut Values",action="store_true",default=False)
parser.add_option("--gaussian", help="Use A Gaussian Signal Template with floating width",type=float,default=-1.0)
parser.add_option("-m","--higgsMass", help="Use This Higgs Mass",type=float,default=-1.0)
parser.add_option("--combinationsOnly", help="Run Only Combinations of Channels",action="store_true",default=False)
args, fakeargs = parser.parse_args()

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

from signalfits import getRooFitSignalPars as sigFits
#from signalfitsNoMuScle import getRooFitSignalPars as sigFits
effReader = EfficiencyReader()

NPROCS = 8

#Scaling Parameter for Bak norm uncertainty
BAKUNC = 1.0

BAKUNCON = True
SIGUNCON = False

FREEBAKPARAMS = True

USEGPANNA = False
if args.cutOpt:
  USEGPANNA = False

SIGNALFIT = [110.,140.]

if args.cutOpt:
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

def convertSigName(name):
  if "ggH" in name:
    return "ggH"
  if "vbfH" in name:
    return "qqH"
  if "wH" in name:
    return "WH"
  if "zH" in name:
    return "ZH"

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

def makePDFBakExpLog(name,rooDataset,mMuMu,minMass,maxMass,workspaceImportFn,mMuMuZ=None,rooDatasetZ=None):
    debug = ""
    debug += "### makePDFBakExpLog: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,mMuMu.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    p1 = root.RooRealVar(channelName+"_p1","p1", 0.0, -1., 1.)
    p2 = root.RooRealVar(channelName+"_p2","p2", 0.0, -1., 1.)
    p3 = root.RooRealVar(channelName+"_p3","p3", 0.0, -1., 1.)
    pdfMmumu = root.RooGenericPdf("bak","TMath::Exp(@0*@0*@1 + @0*@2 + @3*TMath::Log(@0) )",root.RooArgList(mMuMu,p1,p2,p3))

    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    #chi2 = pdfMmumu.createChi2(rooDataset)

    rooParamList = [p1,p2,p3]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    for param in rooParamList:
        param.setConstant(not FREEBAKPARAMS)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(fr)

    #Norm Time
    bakNormTup = None
    if True:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(mMuMu,"signal")
      getSidebandString = "mMuMu < {0} || mMuMu > {1}".format(*signalRangeList)
      nSideband =  rooDataset.sumEntries(getSidebandString)
      nData =  rooDataset.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
    #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #rooDataset2 = rooDataset.reduce(root.RooFit.CutRange("low,signal,high"))
    #rooDataset2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(rooDataset2)

#    ## Debug Time
#    frame = mMuMu.frame()
#    frame.SetName("bak_Plot")
#    rooDataset.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    canvas.SaveAs("debug_"+name+channelName+".png")

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
    debug += "#    Bak Norm Tuple: {0:.2f} {1:.2f}\n".format(*bakNormTup)

    return paramList, bakNormTup, debug


def makePDFBakExpMOverSq(name,rooDataset,mMuMu,minMass,maxMass,workspaceImportFn,mMuMuZ=None,rooDatasetZ=None):
    debug = ""
    debug += "### makePDFBakExpMOverSq: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,mMuMu.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    InvPolMass = root.RooRealVar(channelName+"_InvPolMass","InvPolMass", 91., 70., 95.)
    ExpMass = root.RooRealVar(channelName+"_ExpMass","ExpMass", 0.0, -1., 1.)
    pdfMmumu = root.RooGenericPdf("bak","TMath::Exp(@0*@2)/(@0-@1)/(@0-@1)",root.RooArgList(mMuMu,InvPolMass,ExpMass))

    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    #chi2 = pdfMmumu.createChi2(rooDataset)

    rooParamList = [InvPolMass,ExpMass]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    for param in rooParamList:
        param.setConstant(not FREEBAKPARAMS)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(fr)

    #Norm Time
    bakNormTup = None
    if True:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(mMuMu,"signal")
      getSidebandString = "mMuMu < {0} || mMuMu > {1}".format(*signalRangeList)
      nSideband =  rooDataset.sumEntries(getSidebandString)
      nData =  rooDataset.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
   #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #rooDataset2 = rooDataset.reduce(root.RooFit.CutRange("low,signal,high"))
    #rooDataset2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(rooDataset2)

#    ## Debug Time
#    frame = mMuMu.frame()
#    frame.SetName("bak_Plot")
#    rooDataset.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    canvas.SaveAs("debug_"+name+channelName+".png")

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
    debug += "#    Bak Norm Tuple: {0:.2f} {1:.2f}\n".format(*bakNormTup)

    return paramList, bakNormTup, debug

def makePDFBakMOverSq(name,rooDataset,mMuMu,minMass,maxMass,workspaceImportFn,mMuMuZ=None,rooDatasetZ=None):
    debug = ""
    debug += "### makePDFBakMOverSq: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,mMuMu.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    InvPolMass = root.RooRealVar(channelName+"_InvPolMass","InvPolMass", 91., 70., 95.)
    pdfMmumu = root.RooGenericPdf("bak","@0/(@0-@1)/(@0-@1)",root.RooArgList(mMuMu,InvPolMass))

    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    #chi2 = pdfMmumu.createChi2(rooDataset)

    rooParamList = [InvPolMass]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    for param in rooParamList:
        param.setConstant(False)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(fr)

    #Norm Time
    bakNormTup = None
    if True:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(mMuMu,"signal")
      getSidebandString = "mMuMu < {0} || mMuMu > {1}".format(*signalRangeList)
      nSideband =  rooDataset.sumEntries(getSidebandString)
      nData =  rooDataset.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
    #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #rooDataset2 = rooDataset.reduce(root.RooFit.CutRange("low,signal,high"))
    #rooDataset2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(rooDataset2)

#    ## Debug Time
#    frame = mMuMu.frame()
#    frame.SetName("bak_Plot")
#    rooDataset.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    canvas.SaveAs("debug_"+name+channelName+".png")

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
    debug += "#    Bak Norm Tuple: {0:.2f} {1:.2f}\n".format(*bakNormTup)

    return paramList, bakNormTup, debug



def makePDFBakOld(name,rooDataset,mMuMu,minMass,maxMass,workspaceImportFn,mMuMuZ=None,rooDatasetZ=None):
    debug = ""
    debug += "### makePDFBakOld: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,mMuMu.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    voitWidth = root.RooRealVar(channelName+"_voitWidth","voitWidth",2.4952)
    voitWidth.setConstant(True)
    voitmZ = root.RooRealVar(channelName+"_voitmZ","voitmZ",85,95)
    voitSig = root.RooRealVar(channelName+"_voitSig","voitSig",0.0,30.0)
    voitMmumu = root.RooVoigtian("bak_voitMmumu","voitMmumu",mMuMu,voitmZ,voitWidth,voitSig)

    expParam = root.RooRealVar(channelName+"_expParam","expParam",-1,0)
    expMmumu = root.RooExponential("bak_expMmumu","expMmumu",mMuMu,expParam)

    mixParam = root.RooRealVar(channelName+"_mixParam","mixParam",0,1)

    pdfMmumu = root.RooAddPdf("bak","bak",root.RooArgList(voitMmumu,expMmumu),root.RooArgList(mixParam))

    # Just For Z-Peak Part
    assert(mMuMuZ != None)
    assert(rooDatasetZ != None)
    voitMmumuZ = root.RooVoigtian("bak_voitMmumuZ","voitMmumuZ",mMuMuZ,voitmZ,voitWidth,voitSig)

    voitMmumuZ.fitTo(rooDatasetZ,root.RooFit.SumW2Error(False),PRINTLEVEL)
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

    expMmumu.fitTo(rooDataset,root.RooFit.Range("high"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    expParam.setConstant(True)
    
    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    #chi2 = pdfMmumu.createChi2(rooDataset)

    rooParamList = [voitmZ,voitSig,expParam,mixParam]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    if FREEBAKPARAMS:
      for param in rooParamList:
        param.setConstant(False)
      voitWidth.setConstant(True)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(fr)

    #Norm Time
    bakNormTup = None
    if True:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(mMuMu),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(mMuMu,"signal")
      getSidebandString = "mMuMu < {0} || mMuMu > {1}".format(*signalRangeList)
      nSideband =  rooDataset.sumEntries(getSidebandString)
      nData =  rooDataset.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
    #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #rooDataset2 = rooDataset.reduce(root.RooFit.CutRange("low,signal,high"))
    #rooDataset2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(rooDataset2)

#    ## Debug Time
#    frame = mMuMu.frame()
#    frame.SetName("bak_Plot")
#    rooDataset.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    canvas.SaveAs("debug_"+name+channelName+".png")

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
    debug += "#    Bak Norm Tuple: {0:.2f} {1:.2f}\n".format(*bakNormTup)

    return paramList, bakNormTup, debug

def makePDFSigCBPlusGaus(name,rooDataset,mMuMu,minMass,maxMass,workspaceImportFn,channelName,forceMean=-1.,sigInject=0):

    debug = ""
    debug += "### makePDFSigCBPlusGaus: "+channelName+": "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,mMuMu.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    mean = root.RooRealVar(channelName+"_"+name+"_Mean",channelName+"_"+name+"_Mean",125.,100.,150.)
    width = root.RooRealVar(channelName+"_"+name+"_Width",channelName+"_"+name+"_Width",5.0,0.5,20.0)
    width2 = root.RooRealVar(channelName+"_"+name+"_Width2",channelName+"_"+name+"_Width2",5.0,0.1,20.0)
    alpha = root.RooRealVar(channelName+"_"+name+"_Alpha",channelName+"_"+name+"_Alpha",1.0,0.1,10.0)
    n = root.RooRealVar(channelName+"_"+name+"_n",channelName+"_"+name+"_n",1.0,0.1,10.0)
    mix = root.RooRealVar(channelName+"_"+name+"_mix",channelName+"_"+name+"_mix",0.5,0.0,1.0)
    cb = root.RooCBShape(name+"_CB",name+"_CB",mMuMu,mean,width,alpha,n)
    gaus = root.RooGaussian(name+"_Gaus",name+"_Gaus",mMuMu,mean,width2)
    pdfMmumu = root.RooAddPdf(convertSigName(name),convertSigName(name),cb,gaus,mix)
    
    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True),root.RooFit.Range("signalfit"))
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
      #workspaceImportFn(rooDataset)
      workspaceImportFn(fr)

    ## Debug Time
#    frame = mMuMu.frame()
#    frame.SetName(name+"_Plot")
#    rooDataset.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    now = datetime.datetime.now().isoformat()
#    canvas.SaveAs("debug_"+name+now+".png")

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())

    sigInjectDataset = None
    if sigInject > 0:
      sigInjectDataset = pdfMmumu.generate(root.RooArgSet(mMuMu),sigInject,False)

    return paramList, debug, sigInjectDataset

def makePDFSigDG(name,rooDataset,mMuMu,minMass,maxMass,workspaceImportFn,channelName,forceMean=-1.,sigInject=0):

    debug = ""
    debug += "### makePDFSigDG: "+channelName+": "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,mMuMu.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    rooParamList = []

    meanG1 = root.RooRealVar(channelName+"_"+name+"_MeanG1",
                             channelName+"_"+name+"_MeanG1", 
                             123.,100.,150.)
    meanG2 = root.RooRealVar(channelName+"_"+name+"_MeanG2",
                             channelName+"_"+name+"_MeanG2", 
                             125.,100.,150.)
    
    widthG1 = root.RooRealVar(channelName+"_"+name+"_WidthG1",
                             channelName+"_"+name+"_WidthG1", 
                             5.,2.,10.)
    widthG2 = root.RooRealVar(channelName+"_"+name+"_WidthG2",
                              channelName+"_"+name+"_WidthG2", 
                              1.,0.5,5.)
    
    mixGG = root.RooRealVar(channelName+"_"+name+"_mixGG",
                            channelName+"_"+name+"_mixGG", 
                            0.5,0.,1.)
    gaus1 = root.RooGaussian(channelName+"_"+name+"_gaus1",
                             channelName+"_"+name+"_gaus1",
                             mMuMu,meanG1,widthG1)
    gaus2 = root.RooGaussian(channelName+"_"+name+"_gaus2",
                             channelName+"_"+name+"_gaus2",
                             mMuMu,meanG2,widthG2)
    pdfMmumu = root.RooAddPdf(convertSigName(name),
                              name,
                              gaus1,gaus2,mixGG)
    rooParamList += [meanG1,meanG2,widthG1,widthG2,mixGG]
    
    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True),root.RooFit.Range("signalfit"))
    fr.SetName(name+"_fitResult")

    #mean.setVal(125.)
    #width.setVal(2.)
    #alpha.setVal(1.4)
    #n.setVal(2.2)

    if forceMean > 0.0:
        meanG2.setVal(forceMean)

    ## Error time

    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]
    for i in rooParamList:
       i.setConstant(True)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      #workspaceImportFn(rooDataset)
      workspaceImportFn(fr)

    ## Debug Time
#    frame = mMuMu.frame()
#    frame.SetName(name+"_Plot")
#    rooDataset.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    now = datetime.datetime.now().isoformat()
#    canvas.SaveAs("debug_"+name+now+".png")

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())

    sigInjectDataset = None
    if sigInject > 0:
      sigInjectDataset = pdfMmumu.generate(root.RooArgSet(mMuMu),sigInject,False)

    return paramList, debug, sigInjectDataset

def makePDFSigNew(channelName,name,mMuMu,mass,workspaceImportFn,useDG=True):

    debug = ""
    debug += "### makePDFSigNew: "+channelName+": "+name+"\n"
    debug += "#    useDG: {0}\n".format(useDG)
    matchChannel = re.match(r"(.*)([\d]TeV)",channelName)
    assert(matchChannel)
    category = matchChannel.group(1)
    energy = matchChannel.group(2)
    massStr = "{0:.0f}".format(mass)
    if mass*10 % 1 > 0.0:
      massStr = "{0:.1f}".format(mass)
    fitTypeString = "CBG"
    if useDG:
      fitTypeString = "DG"
    params = sigFits.parameters[energy][massStr][category][fitTypeString]

    rooParamList = []
    if useDG:
      if len(params.keys())>4:
        # define the Double Gaussian
        meanG1 = root.RooRealVar(channelName+"_"+name+"_MeanG1",
                                 channelName+"_"+name+"_MeanG1", 
                                 params['meanG1'])
        meanG2 = root.RooRealVar(channelName+"_"+name+"_MeanG2",
                                 channelName+"_"+name+"_MeanG2", 
                                 params['meanG2'])
        
        widthG1 = root.RooRealVar(channelName+"_"+name+"_WidthG1",
                                 channelName+"_"+name+"_WidthG1", 
                                 params['widthG1'])
        widthG2 = root.RooRealVar(channelName+"_"+name+"_WidthG2",
                                  channelName+"_"+name+"_WidthG2", 
                                  params['widthG2'])
        
        mixGG = root.RooRealVar(channelName+"_"+name+"_mixGG",
                                channelName+"_"+name+"_mixGG", 
                                params['mixGG'])
        gaus1 = root.RooGaussian(channelName+"_"+name+"_gaus1",
                                 channelName+"_"+name+"_gaus1",
                                 mMuMu,meanG1,widthG1)
        gaus2 = root.RooGaussian(channelName+"_"+name+"_gaus2",
                                 channelName+"_"+name+"_gaus2",
                                 mMuMu,meanG2,widthG2)
        pdfMmumu = root.RooAddPdf(convertSigName(name),
                                  name,
                                  gaus1,gaus2,mixGG)
        workspaceImportFn(pdfMmumu)
        rooParamList += [meanG1,meanG2,widthG1,widthG2,mixGG]
      else:
        # define the Single Gaussian for the EE category
        meanSG  = root.RooRealVar(channelName+"_"+name+"_MeanSG", 
                                  "MeanSG", 
                                  params['meanSG'])
        widthSG = root.RooRealVar(channelName+"_"+name+"_WidthSG",
                                  channelName+"_"+name+"_WidthSG", 
                                  params['widthSG'])
        pdfMmumu = root.RooGaussian(convertSigName(name),name,mMuMu,meanSG,widthSG)
        workspaceImportFn(pdfMmumu)
        rooParamList += [meanSG,widthSG]
        debug += "#    using Single Gaussian (Probably for EE)"
    else:
      mean = root.RooRealVar(channelName+"_"+name+"_Mean",
                             channelName+"_"+name+"_Mean",
                             params['mean'])
      width = root.RooRealVar(channelName+"_"+name+"_Width",
                              channelName+"_"+name+"_Width",
                              params['width1'])
      width2 = root.RooRealVar(channelName+"_"+name+"_Width2",
                               channelName+"_"+name+"_Width2",
                               params['width2'])
      alpha = root.RooRealVar(channelName+"_"+name+"_Alpha",
                               channelName+"_"+name+"_Alpha",
                               params['Alpha'])
      n = root.RooRealVar(channelName+"_"+name+"_n",
                               channelName+"_"+name+"_n",
                               params['n'])
      mix = root.RooRealVar(channelName+"_"+name+"_mix",
                               channelName+"_"+name+"_mix",
                               params['mix'])
      cb = root.RooCBShape(name+"_CB",name+"_CB",mMuMu,mean,width,alpha,n)
      gaus = root.RooGaussian(name+"_Gaus",name+"_Gaus",mMuMu,mean,width2)
      pdfMmumu = root.RooAddPdf(convertSigName(name),name,cb,gaus,mix)
      workspaceImportFn(pdfMmumu)
      rooParamList += [mean,width,alpha,n]

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())

    return debug

#makePDFSig = makePDFSigCBPlusGaus
makePDFSig = makePDFSigDG
makePDFBak = makePDFBakOld
#makePDFBak = makePDFBakMOverSq
#makePDFBak = makePDFBakExpMOverSq
#makePDFBak = makePDFBakExpLog

###################################################################################

class Analysis:
  def __init__(self,directory,signalNames,backgroundNames,dataNames,analysis,lumi,controlRegionVeryLow,controlRegionLow,controlRegionHigh,toyData=False,sigInject=0.0,sigInjectMass=125.0,energyStr="8TeV",cutString=""):
    if analysis[0] in "0123456789-.*+/":
       print("Error: Analysis: analysis name '{0}' begins with a forbidden character".format(analysis))
       sys.exit(1)
    getCutHist = getattr(self,"getCutHist")
    doSigInject = getattr(self,"doSigInject")
    self.treename = "outtree"
    self.directory = directory
    self.weightName = "puWeight"
    self.fileList = []
    self.treeList = []
    self.origHistList = []
    self.binSize = 0.5
    
    higgsPeakMean = args.higgsMass - 0.3
    self.higgsMass = 125.
    if args.higgsMass > 0.0:
      self.higgsMass = args.higgsMass

    self.energyStr = energyStr
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
    self.debug += "#Nominal Higgs Mass: "+str(self.higgsMass) +"\n"
    self.debug += "#Peak Centered at: "+str(higgsPeakMean) +"\n"

    self.workspace = root.RooWorkspace(analysis+energyStr)
    self.workspaceName = analysis+energyStr
    wImport = getattr(self.workspace,"import")
    self.sigInjectWorkspaces = []

    maxMass = controlRegionHigh[1]
    minMass = controlRegionLow[0]
    self.minMass = minMass
    self.maxMass = maxMass
    mMuMu = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)
    mMuMu.setRange("low",controlRegionLow[0],controlRegionLow[1])
    mMuMu.setRange("high",controlRegionHigh[0],controlRegionHigh[1])
    mMuMu.setRange("signal",controlRegionLow[1],controlRegionHigh[0])
    mMuMu.setRange("signalfit",SIGNALFIT[0],SIGNALFIT[1])
    self.mMuMu = mMuMu

    self.rooPUWeight = root.RooRealVar("puWeight","PU Weight",0.0,200.)
    self.observablesForRooDataset = root.RooArgSet(
                            self.mMuMu,
                            self.rooPUWeight,
                        )

    # Hack to Make makePDFBakOld work
    self.minMassZ = 88.
    self.maxMassZ = 94.
    mMuMuZ = root.RooRealVar("mMuMu","mMuMu",self.minMassZ,self.maxMassZ)
    self.mMuMuZ = mMuMuZ
    self.observablesForRooDatasetZ = root.RooArgSet(
                            mMuMuZ,
                            self.rooPUWeight,
                        )

    # Dealing with cut string
    KDString = '(sigMEPdf*{0}/(bakMEPdf*{1}+sigMEPdf*{0}))'.format(
                   MENormDict[energyStr]["sigMEPdf"],MENormDict[energyStr]["bakMEPdf"]
                )
    self.fullCutString = treeCut(analysis,cutString,eventWeights=False,muonRequirements=True,KDString=KDString)
    self.debug += "#  Full Cut String:\n"
    self.debug += "#  {0}\n".format(self.fullCutString)

    self.sigHistsRaw = []
    if not USEGPANNA:
      for name in signalNames:
        tmpH = self.getRooDataSample(name,self.observablesForRooDataset)
        self.debug += "#  N Events in {0}: {1}\n".format(name,tmpH.sumEntries())
        self.sigHistsRaw.append(tmpH)

    self.bakHistsRaw = []
    self.bakHistsRawZ = []
    for name in backgroundNames:
      tmpH = self.getRooDataSample(name,self.observablesForRooDataset)
      self.debug += "#  N Events in {0}: {1}\n".format(name,tmpH.sumEntries())
      self.bakHistsRaw.append(tmpH)
      tmpH = self.getRooDataSample(name,self.observablesForRooDatasetZ,True)
      self.bakHistsRawZ.append(tmpH)

    self.datHists = []
    self.datHistsZ = []
    for name in dataNames:
      tmpH = self.getRooDataSample(name,self.observablesForRooDataset)
      self.debug += "#  N Events in {0}: {1}\n".format(name,tmpH.sumEntries())
      self.datHists.append(tmpH)
      tmpH = self.getRooDataSample(name,self.observablesForRooDatasetZ,True)
      self.datHistsZ.append(tmpH)

    effMap = {}
    xsecMap = {}
    lowBin = 0
    massBounds = [controlRegionLow[0],controlRegionHigh[1]]
    self.massBounds = massBounds

    self.xsecBakTotal = 0.0
    self.xsecBakList = []
    self.effBakList = []
    self.bakHists = []
    self.bakHistTotal = None
    for h,name in zip(self.bakHistsRaw,backgroundNames):
      sys.exit(1)
      counts = h.sumEntries()
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

    #self.bakHistTotal.Scale(lumi)
    #self.bakHistTotalReal = self.bakHistTotal.Clone("data_obs")

    self.dataCountsTotal = None
    self.datHistTotal = None
    for h,name in zip(self.datHists,dataNames):
      counts = h.sumEntries()
      if self.dataCountsTotal == None:
        self.dataCountsTotal = counts
      else:
        self.dataCountsTotal += counts
      if self.datHistTotal == None:
        self.datHistTotal = h.Clone("data_obs")
      else:
        self.datHistTotal.add(h)
    assert(self.dataCountsTotal == self.datHistTotal.sumEntries())

    self.dataCountsTotalZ = None
    self.datHistTotalZ = None
    for h,name in zip(self.datHistsZ,dataNames):
      counts = h.sumEntries()
      if self.dataCountsTotalZ == None:
        self.dataCountsTotalZ = counts
      else:
        self.dataCountsTotalZ += counts
      if self.datHistTotalZ == None:
        self.datHistTotalZ = h.Clone("data_obs_Z")
      else:
        self.datHistTotalZ.add(h)
    assert(self.dataCountsTotalZ == self.datHistTotalZ.sumEntries())

    histToUseForBak = self.datHistTotal
    histToUseForBakZ = self.datHistTotalZ
    if histToUseForBak is None:
        print("Background MC rescaling RooDataSet to lumi doesn't work...\n  you must use real data...Exiting")
        sys.exit(1)
        doSigInject(self.bakHistTotal,sigInject,sigInjectMass)
        histToUseForBak = self.bakHistTotal
    else:
        self.dataCountsTotal += doSigInject(self.datHistTotal,sigInject,sigInjectMass)

    tmpBakParams, self.bakNormTup, tmpDebug = makePDFBak(
                                    analysis+energyStr,histToUseForBak,
                                    mMuMu, minMass,maxMass,wImport,
                                    mMuMuZ,histToUseForBakZ
                                       )
    self.params.extend(tmpBakParams)
    self.countsBakTotal = self.bakNormTup[0]*self.bakNormTup[1]
    self.debug += tmpDebug
    
    self.sigParamListList = []
    if USEGPANNA:
      for name in signalNames:
        self.debug +=  makePDFSigNew(analysis+energyStr,name,mMuMu,self.higgsMass,wImport)
    elif True: # Use ggH for non-VBF channels, and VBF shape for VBF
      sigNameToUse = None
      sigHistCounts = []
      for name, hist in zip(signalNames,self.sigHistsRaw):
         sigHistCounts += [hist.sumEntries()]
      maxCounts = max(sigHistCounts)
      maxIndex = sigHistCounts.index(maxCounts)
      self.debug += "# Signal Histogram for Shape: {0} Counts: {1:.1f}\n".format(signalNames[maxIndex],maxCounts)
      sigHistRawToUse = self.sigHistsRaw[maxIndex]
      for name in signalNames:
          sigParams, sigDebug, tmpDS = makePDFSig(name,sigHistRawToUse,mMuMu,minMass,maxMass,wImport,analysis+energyStr)
          self.sigParamListList.append(sigParams)
          self.debug += sigDebug
    else:
      for name, hist in zip(signalNames,self.sigHistsRaw):
          sigParams, sigDebug, tmpDS = makePDFSig(name,hist,mMuMu,minMass,maxMass,wImport,analysis+energyStr,forceMean=higgsPeakMean)
          self.sigParamListList.append(sigParams)
          self.debug += sigDebug
          
    self.xsecSigTotal = 0.0
    self.xsecSigList = []
    self.effSigList = []
    self.sigHists = []
    if USEGPANNA:
      for name in signalNames:
        higgsMassStr = "{0:.0f}".format(self.higgsMass)
        if self.higgsMass*10 % 1 > 0.0:
          higgsMassStr = "{0:.1f}".format(self.higgsMass)
        nameMatch = re.match(r"(.+)Hmumu[\d.]+_[\d]+TeV",name)
        assert(nameMatch)
        prodMode = nameMatch.group(1)
        eff, effErr = effReader(self.energyStr,prodMode,analysis,higgsMassStr)
        tmpName = name.replace("125",higgsMassStr)
        xs = eff*xsec[tmpName]
        self.xsecSigTotal += xs
        self.xsecSigList.append(xs)
        self.effSigList.append(eff)
    else:
      for h,name in zip(self.sigHistsRaw,signalNames):
        counts = h.sumEntries()
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
      assert(False)
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
      #toyDataHist = toyDataset.binnedClone("data_obs","Toy Data")
      self.dataCountsTotal = int(toyDataset.sumEntries())
      wImport(toyDataset)
    else:
      #realDataHist = root.RooDataHist("data_obs","Real Observed Data",root.RooArgList(mMuMu),self.datHistTotal)
      #wImport(realDataHist)
      #realDataHistNotVeryLow = realDataHist.reduce(root.RooFit.CutRange("low,signal,high"))
      #wImport(realDataHistNotVeryLow)
      self.datHistTotal.SetTitle("Real Observed Data")
      wImport(self.datHistTotal)

    self.debug += "#  Signal Efficiencies: \n"
    for sigName,sigEff in zip(signalNames,self.effSigList):
      self.debug += "#    {0}: {1:.2g}\n".format(sigName,sigEff)

  def getCutHist(self,inHist,outHist,cut,cutAbove=False):
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
      if keepBelow:
        for iY in range(0,cutBin):
          mySum += inHist.GetBinContent(iX,iY)
      else:
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

  def getRooDataSample(self,name,observables,aroundZ=False):
    tmpFLoc = self.directory+name+".root"
    tmpF = root.TFile(tmpFLoc)
    tmpTree = tmpF.Get(self.treename)
    tmpTree.SetCacheSize(10000000);
    tmpTree.AddBranchToCache("*");
    #self.treeList.append(tmpTree)
    #self.fileList.append(tmpF)
    minMass = self.minMass
    maxMass = self.maxMass
    mMuMu = self.mMuMu
    if aroundZ:
      minMass = self.minMassZ
      maxMass = self.maxMassZ
      mMuMu = self.mMuMuZ
    if False:
      #tmpSelTree = tmpTree.CopyTree(self.fullCutString)
      tmpSelTree = tmpTree
      print("Error: getRooDataSample RooDataSet Selection is broken!!!")
      tmpH = root.RooDataSet(name,name,tmpSelTree,
            observables, "",
            self.weightName
            )
      return tmpH
    else:
      histName = "hist{0:f}".format(time.time()).replace('.','')
      nBins = int((maxMass-minMass)/self.binSize)
      tmpHist = root.TH1F(histName,"",nBins,minMass,maxMass)
      drawString = "dimuonMass >> {0}".format(histName)
      tmpTree.Draw(drawString,self.fullCutString)
      self.origHistList.append(tmpHist)
      tmpH = root.RooDataHist(name,name,root.RooArgList(mMuMu),tmpHist)
      return tmpH

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
  def __init__(self,directory,analysisNames,signalNames,backgroundNames,dataNames,outfilename,lumi,nuisanceMap=None,controlRegionLow=[110.,115],controlRegionHigh=[135,150],controlRegionVeryLow=[80.,110.],sigInject=0.0,sigInjectMass=125.0,toyData=False,energyStr="8TeV",cutString=""):

    ########################
    ## Setup

    if len(analysisNames)==0:
      print("Error: DataCardMaker: len(analysisNames)==0")
      return

    channels = []
    self.is2D = False

    if type(energyStr) == list:
      for analysis in analysisNames:
        cutStringTmp = cutString
        if type(analysis)==list:
          if len(analysis)>1:
            if len(cutStringTmp)>0:
              cutStringTmp += " && "
            cutStringTmp += analysis[1]
          analysis = analysis[0]
        for es,sn,bn,dn,lu in zip(energyStr,signalNames,backgroundNames,dataNames,lumi):
          lu *= 1000.0
          tmp = Analysis(directory,sn,bn,dn,analysis,lu,controlRegionVeryLow,controlRegionLow,controlRegionHigh,toyData=toyData,sigInject=sigInject,sigInjectMass=sigInjectMass,energyStr=es,cutString=cutStringTmp)
          channels.append(tmp)
      self.channels = channels
    else:
      lumi *= 1000.0
      for analysis in analysisNames:
        cutStringTmp = cutString
        if type(analysis)==list:
          if len(analysis)>1:
            if len(cutStringTmp)>0:
              cutStringTmp += " && "
            cutStringTmp += analysis[1]
          analysis = analysis[0]
        tmp = Analysis(directory,signalNames,backgroundNames,dataNames,analysis,lumi,controlRegionVeryLow,controlRegionLow,controlRegionHigh,toyData=toyData,sigInject=sigInject,sigInjectMass=sigInjectMass,energyStr=energyStr,cutString=cutStringTmp)
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
        rootDebugString += "# "+channel.workspace.data("data_obs").GetTitle()+"\n"
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
          proc1FormatList.append(convertSigName(sigName))
  
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

      # Parameter Uncertainties
      for channel,channelName in zip(self.channels,self.channelNames):
        for nu in channel.params:
          if (not FREEBAKPARAMS) or ("voit" in nu.name):
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
    outfile.write("#\n#\n")
    outfile.write(rootDebugString)
    outfile.write("# Signal Injected: {0:.1f} X SM, mH = {1:.1f}\n".format(sigInject,sigInjectMass))

#    for channel,channelName in zip(self.channels,self.channelNames):
#        outfile.write("#\n")
#        outfile.write("#Efficiency: channel {0}: \n".format(channelName))
#        for sigName in channel.sigNames:
#          outfile.write("#  {0:>20}: {1:4f}\n".format(sigName,channel.getSigEff(sigName)))
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

  directory = "input/V00-01-10/"
  outDir = "statsCards/"
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
  analyses = [["VBFBDTCut"]]
  analyses += [["VBFPresel"]]
  puidTightString = " && jetLead_FullPUIDFlag >= 7 && jetSub_FullPUIDFlag >= 7"
  oldJetCutString = "jetLead_pt > 30. && jetSub_pt > 30."
  funJetCutString = "((jetLead_pt > 25.) || (jetLead_pt > 20. && abs(jetLead_eta)<2.4))"
  funJetCutString += " && ((jetSub_pt > 25.) || (jetSub_pt > 20. && abs(jetSub_eta)<2.4))"
  funJetCutString += puidTightString

  kindaFunJetCutString = oldJetCutString
  kindaFunJetCutString += puidTightString
  analyses += [["VBFCutBasedLoose",oldJetCutString]]
  analyses += [["VBFCutBasedTight",oldJetCutString]]
  #analyses += [["IncPreselPtG10"]+ x for x in categoriesInc]
  analyses = []
  #analyses += [["Inclusive",""]]
  #analyses += [["Jet0All","nJets == 0"]]
  #analyses += [["Jet1All","nJets == 1"]]
  #analyses += [["Jet2All","nJets >= 2"]]
  #analyses += [["Jet2AllPtMissL100.","nJets >= 2 && ptMiss<100."]]

  #analyses += [["Jet2Pass","nJets >= 2 && deltaEtaJets>3.5 && dijetMass>650. && ptMiss < 100."]]
  #analyses += [["Jet2Fail","nJets >= 2 && !(deltaEtaJets>3.5 && dijetMass>650.) && ptMiss < 100."]]
  #analyses += [["Jet1Pass","nJets == 1 && dimuonPt > 55."]]
  #analyses += [["Jet1Fail","nJets == 1 && !(dimuonPt > 55.)"]]
  #analyses += [["Jet0Pass","nJets == 0 && dimuonPt > 10."]]
  #analyses += [["Jet0Fail","nJets == 0 && !(dimuonPt > 10.)"]]

  #analyses += [["Jet0BB","nJets == 0"]]
  #analyses += [["Jet0BO","nJets == 0"]]
  #analyses += [["Jet0BE","nJets == 0"]]
  #analyses += [["Jet0OO","nJets == 0"]]
  #analyses += [["Jet0FF","nJets == 0"]]
  #analyses += [["Jet1BB","nJets == 1"]]
  #analyses += [["Jet1BO","nJets == 1"]]
  #analyses += [["Jet1BE","nJets == 1"]]
  #analyses += [["Jet1OO","nJets == 1"]]
  #analyses += [["Jet1FF","nJets == 1"]]



  #analyses += [["VBFJustinHighJetPtPtMissL25DEtaJetsG3p7MJJG400",kindaFunJetCutString+" && deltaEtaJets > 3.7 && dijetMass>400. && ptMiss < 25."]]
  #analyses += [["VBFJustinLowJetPtPtMissL25DEtaJetsG3p6MJJG450",funJetCutString+" && deltaEtaJets > 3.6 && dijetMass>450. && ptMiss < 25."]]
  combinations = []

  combinations.append((
    [
      ["Jet2Pass","nJets >= 2 && bdtVBF>0.05 && ptMiss < 100."],
      ["Jet2Fail","nJets >= 2 && !(bdtVBF>0.05) && ptMiss < 100."]
    ],"Jet2BDTSplitPtMiss100"
  ))
  combinations.append((
    [
      ["Jet2Pass","nJets >= 2 && bdtVBF>0.05 && ptMiss < 50."],
      ["Jet2Fail","nJets >= 2 && !(bdtVBF>0.05) && ptMiss < 50."]
    ],"Jet2BDTSplitPtMiss50"
  ))
  combinations.append((
    [
      ["Jet2Pass","nJets >= 2 && bdtVBF>0.05 && ptMiss < 40."],
      ["Jet2Fail","nJets >= 2 && !(bdtVBF>0.05) && ptMiss < 40."]
    ],"Jet2BDTSplitPtMiss40"
  ))
  combinations.append((
    [
      ["Jet2Pass","nJets >= 2 && bdtVBF>0.05 && ptMiss < 30."],
      ["Jet2Fail","nJets >= 2 && !(bdtVBF>0.05) && ptMiss < 30."]
    ],"Jet2BDTSplitPtMiss30"
  ))
  combinations.append((
    [
      ["Jet2Pass","nJets >= 2 && bdtVBF>0.05 && ptMiss < 25."],
      ["Jet2Fail","nJets >= 2 && !(bdtVBF>0.05) && ptMiss < 25."]
    ],"Jet2BDTSplitPtMiss25"
  ))

  #combinations.append((
  #  [
  #    ["Jet2Pass","nJets >= 2 && deltaEtaJets>3.5 && dijetMass>650. && ptMiss < 100."],
  #    ["Jet2Fail","nJets >= 2 && !(deltaEtaJets>3.5 && dijetMass>650.) && ptMiss < 100."]
  #  ],"Jet2SplitPtMiss100"
  #))
  #combinations.append((
  #  [
  #    ["Jet1Pass","nJets == 1 && dimuonPt > 55."],
  #    ["Jet1Fail","nJets == 1 && !(dimuonPt > 55.)"]
  #  ],"Jet1Split"
  #))
  #combinations.append((
  #  [
  #    ["Jet0Pass","nJets == 0 && dimuonPt > 10."],
  #    ["Jet0Fail","nJets == 0 && !(dimuonPt > 10.)"]
  #  ],"Jet0Split"
  #))
  #combinations.append((
  #  [
  #    ["Jet2Pass","nJets >= 2 && deltaEtaJets>3.5 && dijetMass>650. && ptMiss < 100."],
  #    ["Jet2Fail","nJets >= 2 && !(deltaEtaJets>3.5 && dijetMass>650.) && ptMiss < 100."],
  #    ["Jet1Pass","nJets == 1 && dimuonPt > 55."],
  #    ["Jet1Fail","nJets == 1 && !(dimuonPt > 55.)"],
  #    ["Jet0Pass","nJets == 0 && dimuonPt > 10."],
  #    ["Jet0Fail","nJets == 0 && !(dimuonPt > 10.)"],
  #  ],"JetNSplit"
  #))
  #combinations.append((
  #  [
  #    ["Jet2All","nJets >= 2 && ptMiss < 100."],
  #    ["Jet1All","nJets == 1"],
  #    ["Jet0All","nJets == 0"]
  #  ],"JetNAll"
  #))

  #combinations.append((
  #  [
  #    ["Jet1PassBB","nJets == 1 && dimuonPt > 55."],
  #    ["Jet1PassBO","nJets == 1 && dimuonPt > 55."],
  #    ["Jet1PassBE","nJets == 1 && dimuonPt > 55."],
  #    ["Jet1PassOO","nJets == 1 && dimuonPt > 55."],
  #    ["Jet1PassFF","nJets == 1 && dimuonPt > 55."],
  #    ["Jet1FailBB","nJets == 1 && !(dimuonPt > 55.)"],
  #    ["Jet1FailBO","nJets == 1 && !(dimuonPt > 55.)"],
  #    ["Jet1FailBE","nJets == 1 && !(dimuonPt > 55.)"],
  #    ["Jet1FailOO","nJets == 1 && !(dimuonPt > 55.)"],
  #    ["Jet1FailFF","nJets == 1 && !(dimuonPt > 55.)"],
  #  ],"Jet1SplitCat"
  #))
  #combinations.append((
  #  [
  #    ["Jet0PassBB","nJets == 0 && dimuonPt > 10."],
  #    ["Jet0PassBO","nJets == 0 && dimuonPt > 10."],
  #    ["Jet0PassBE","nJets == 0 && dimuonPt > 10."],
  #    ["Jet0PassOO","nJets == 0 && dimuonPt > 10."],
  #    ["Jet0PassFF","nJets == 0 && dimuonPt > 10."],
  #    ["Jet0FailBB","nJets == 0 && !(dimuonPt > 10.)"],
  #    ["Jet0FailBO","nJets == 0 && !(dimuonPt > 10.)"],
  #    ["Jet0FailBE","nJets == 0 && !(dimuonPt > 10.)"],
  #    ["Jet0FailOO","nJets == 0 && !(dimuonPt > 10.)"],
  #    ["Jet0FailFF","nJets == 0 && !(dimuonPt > 10.)"],
  #  ],"Jet0SplitCat"
  #))
  #combinations.append((
  #  [
  #    ["Jet2Pass","nJets >= 2 && deltaEtaJets>3.5 && dijetMass>650. && ptMiss < 100."],
  #    ["Jet2Fail","nJets >= 2 && !(deltaEtaJets>3.5 && dijetMass>650.) && ptMiss < 100."],
  #    ["Jet1PassBB","nJets == 1 && dimuonPt > 55."],
  #    ["Jet1PassBO","nJets == 1 && dimuonPt > 55."],
  #    ["Jet1PassBE","nJets == 1 && dimuonPt > 55."],
  #    ["Jet1PassOO","nJets == 1 && dimuonPt > 55."],
  #    ["Jet1PassFF","nJets == 1 && dimuonPt > 55."],
  #    ["Jet1FailBB","nJets == 1 && !(dimuonPt > 55.)"],
  #    ["Jet1FailBO","nJets == 1 && !(dimuonPt > 55.)"],
  #    ["Jet1FailBE","nJets == 1 && !(dimuonPt > 55.)"],
  #    ["Jet1FailOO","nJets == 1 && !(dimuonPt > 55.)"],
  #    ["Jet1FailFF","nJets == 1 && !(dimuonPt > 55.)"],
  #    ["Jet0PassBB","nJets == 0 && dimuonPt > 10."],
  #    ["Jet0PassBO","nJets == 0 && dimuonPt > 10."],
  #    ["Jet0PassBE","nJets == 0 && dimuonPt > 10."],
  #    ["Jet0PassOO","nJets == 0 && dimuonPt > 10."],
  #    ["Jet0PassFF","nJets == 0 && dimuonPt > 10."],
  #    ["Jet0FailBB","nJets == 0 && !(dimuonPt > 10.)"],
  #    ["Jet0FailBO","nJets == 0 && !(dimuonPt > 10.)"],
  #    ["Jet0FailBE","nJets == 0 && !(dimuonPt > 10.)"],
  #    ["Jet0FailOO","nJets == 0 && !(dimuonPt > 10.)"],
  #    ["Jet0FailFF","nJets == 0 && !(dimuonPt > 10.)"],
  #  ],"JetNSplitCat"
  #))
  #combinations.append((
  #      [["IncPreselPtG10"+x] for x in categoriesInc],"IncPreselCat"
  #))
  #combinations.append((
  #      [["VBFBDTCut",oldJetCutString]]+[["IncPreselPtG10"+x] for x in categoriesInc],"BDTCutCatVBFBDTOnly"
  #))
#  combinations.append((
#        ["VBFPresel"]+["IncPreselPtG10"],"BDTCutVBFBDTOnly"
#  ))

  # Multi-dimensional Optimization of Cuts
  # First two arguments are just like combinations
  # 3rd arg is a dictionary with the keys 
  # that are variables to be cut on with "L" or "G" appended
  # for cutting Less-than, or Greater-then.
  # Finally the last argument specifies whether to "Split"
  # events rather than cutting them.  If set to True,
  # A category passing the cuts will be made and combined
  # with a category failing the events.  If False, the
  # fail events are discarded.
  combinationsCutOpt = []
  #combinationsCutOpt.append((
  #  [["Yay","nJets>=2 && ptMiss < 100."]],"Jets2BDTSplitOptPtMissL100",{
  #      'bdtVBFG':[21,-0.5,0.5],
  #      },True
  #))
  #combinationsCutOpt.append((
  #  [["Yay","nJets>=2 && ptMiss < 50."]],"Jets2BDTSplitOptPtMissL50",{
  #      'bdtVBFG':[21,-0.5,0.5],
  #      },True
  #))
  #combinationsCutOpt.append((
  #  [["Yay","nJets>=2 && ptMiss < 40."]],"Jets2BDTSplitOptPtMissL40",{
  #      'bdtVBFG':[21,-0.5,0.5],
  #      },True
  #))
  #combinationsCutOpt.append((
  #  [["Yay","nJets>=2 && ptMiss < 30."]],"Jets2BDTSplitOptPtMissL30",{
  #      'bdtVBFG':[21,-0.5,0.5],
  #      },True
  #))
  #combinationsCutOpt.append((
  #  [["Yay","nJets>=2 && ptMiss < 25."]],"Jets2BDTSplitOptPtMissL25",{
  #      'bdtVBFG':[21,-0.5,0.5],
  #      },True
  #))
  #combinationsCutOpt.append((
  #  [["Yay","nJets>=2 && ptMiss < 100."]],"Jets2SplitOptPtMissL100",{
  #      'deltaEtaJetsG':[7,2.0,5.0],
  #      'dijetMassG':[11,300.,800.],
  #      },True
  #))
  combinationsCutOpt.append((
    [["Yay","nJets>=2 && ptMiss < 50."]],"Jets2SplitOptPtMissL50",{
        'deltaEtaJetsG':[7,2.0,5.0],
        'dijetMassG':[11,300.,800.],
        },True
  ))
  combinationsCutOpt.append((
    [["Yay","nJets>=2 && ptMiss < 30."]],"Jets2SplitOptPtMissL30",{
        'deltaEtaJetsG':[7,2.0,5.0],
        'dijetMassG':[11,300.,800.],
        },True
  ))
  #combinationsCutOpt.append((
  #  [["Yay","nJets>=2 && ptMiss < 40."]],"Jets2SplitOptPtMissL40",{
  #      'deltaEtaJetsG':[7,2.0,5.0],
  #      'dijetMassG':[11,300.,800.],
  #      },True
  #))
  #combinationsCutOpt.append((
  #  [["Yay","nJets==1"]],"Jets1SplitOpt",{
  #      'dimuonPtG':[16,0.,150.],
  #      },True
  #))
  #combinationsCutOpt.append((
  #  [["Yay","nJets==0"]],"Jets0SplitOpt",{
  #      'dimuonPtG':[13,0.,60.],
  #      },True
  #))
  #combinationsCutOpt.append((
  #  [["Yay","nJets==1"]],"Jets1CutOpt",{
  #      'dimuonPtG':[9,0.,40.],
  #      },False
  #))
  #combinationsCutOpt.append((
  #  [["Yay","nJets==0"]],"Jets0CutOpt",{
  #      'dimuonPtG':[9,0.,40.],
  #      },False
  #))
  #combinationsCutOpt.append((
  #  [["Yay","nJets==0"]],"Jets0OptSplitTest",{
  #      'dimuonPtG':[5,0.,100.],
  #      },True
  #))
  #combinationsCutOpt.append((
  #  [["Yay","nJets==0"]],"Jets0OptCutTest",{
  #      'dimuonPtG':[5,0.,100.],
  #      },False
  #))

  histPostFix="/mDiMu"
  signalNames=["ggHmumu125","vbfHmumu125","wHmumu125","zHmumu125"]
  signalNames=["ggHmumu125","vbfHmumu125"]
  #backgroundNames= ["DYJetsToLL","ttbar"]
  backgroundNames= []
  dataDict = {}
  dataDict["8TeV"] = [
    "SingleMuRun2012Av1",
    "SingleMuRun2012Av1Recover",
    "SingleMuRun2012Bv1",
    "SingleMuRun2012Cv1",
    "SingleMuRun2012Cv2",
    "SingleMuRun2012D",
  ]
  #dataDict["8TeV"] = []
  dataDict["7TeV"] = [
    "SingleMuRun2011Av1",
    "SingleMuRun2011Bv1"
  ]
  #dataDict["7TeV"] = []
  dataDict["14TeV"] = []
  lumiList = [lumiDict["8TeV"],20,25,30]
  lumiList = [lumiDict["8TeV"]]

  controlRegionVeryLow=[60,110]
  controlRegionLow=[110,120]
  controlRegionHigh=[130,160]
  if args.higgsMass > 0.0:
    controlRegionLow=[110,args.higgsMass-5.0]
    controlRegionHigh=[args.higgsMass+5.0,160]

  shape=True
  toyData=args.toyData

  if args.combinationsOnly:
    analyses = []
  threads = []
  if not args.cutOpt:
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
          cutsToDo = ""
          if len(ana) == 2:
            cutsToDo = ana[1]
          anaName = ana[0]
          tmp = ThreadedCardMaker(
            #__init__ args:
            directory,[anaName],
            appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),dataDict[p],
            controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,
            controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,sigInjectMass=args.signalInjectMass,
            energyStr=p, cutString=cutsToDo,
            #write args:
            outfilename=outDir+anaName+"_"+p+"_"+str(name)+".txt",lumi=i
            )
          threads.append(tmp)
        for comb in combinations:
         threads.append(
          ThreadedCardMaker(
            #__init__ args:
            directory,
            comb[0],
            appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),dataDict[p],
            controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,
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
          cutsToDo = ""
          if len(ana) == 2:
            cutsToDo = ana[1]
          anaName = ana[0]
          tmp = ThreadedCardMaker(
            #__init__ args:
            directory,[anaName],
            [appendPeriod(signalNames,p) for p in periods],
            [appendPeriod(backgroundNames,p) for p in periods],
            [dataDict[p] for p in periods],
            controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,
            controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,sigInjectMass=args.signalInjectMass,
            energyStr=periods,cutString=cutsToDo,
            #write args:
            outfilename=outDir+anaName+"_"+filenamePeriod+"_"+filenameLumi+".txt",lumi=[float(lumiDict[p]) for p in periods]
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
            controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,
            controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,sigInjectMass=args.signalInjectMass,
            energyStr=periods,
            #write args:
            outfilename=outDir+comb[1]+"_"+filenamePeriod+"_"+filenameLumi+".txt",lumi=[float(lumiDict[p]) for p in periods]
          )
         )
  
  else: #Do cutOpt Stuff
    print "combinationsCutOpt: {0}".format(combinationsCutOpt)
    for p in periods:
      i = lumiDict[p]
      for comb in combinationsCutOpt:
        print "\ncomb: {0}\n".format(comb)
        cutDict = comb[2]
        print "cutDict: {0}".format(cutDict)
        cutBaseString = "1"
        nameBaseString = comb[1]
        cutDictKeyList = cutDict.keys()
        nCutVars = len(cutDictKeyList)
        cutPoints = []
        cutPointsFail = []
        doSplit=comb[3]
        for key,cutIndex in zip(cutDictKeyList,range(nCutVars)):
          cutIndex = str(cutIndex)
          cut = key[:-1]
          cutType = key[-1]
          cutOp = ''
          cutOpFail = ''
          if cutType == "G":
            cutOp = '>'
            cutOpFail = '<='
          elif cutType == "L":
            cutOp = '<'
            cutOpFail = '>='
          else:
            print("Cut Type Not Recognized...exiting.")
            sys.exit(1)
          cutBaseString += " && "+cut+" "+cutOp+" {"+cutIndex+"}"
          nameBaseString += "_"+key+"{"+cutIndex+"}"
          cutData = cutDict[key]
          print "cutData: {0}".format(cutData)
          cutStep = (cutData[2]-cutData[1])/(cutData[0]-1)
          tmpCuts = [cutData[1]+iC*cutStep for iC in range(cutData[0])]
          cutPoints += [tmpCuts]
        print "cutBaseString: "+cutBaseString
        print "nameBaseString: "+nameBaseString
        print "cutPoints: {0}".format(cutPoints)
        assert(len(cutPoints)==nCutVars)
        
        notDone=True
        indexList = [0 for iC in range(nCutVars)]
        lenList = [len(iC) for iC in cutPoints]
        while notDone:
            cutList = [cL[iC] for iC, cL in zip(indexList,cutPoints)]
            cutString = cutBaseString.format(*cutList)
            cutStringFail = "!( "+cutString+" )"
            nameString = nameBaseString.format(*cutList)
            nameString = nameString.replace('.','p')
            print "cutString: "+cutString
            if doSplit:
              print "cutStringFail: "+cutStringFail
            analysisNames = []
            for ana in comb[0]:
              anaName = ana[0]
              anaCutsOrig = '1'
              if len(ana)>1:
                if ana[1] != "":
                  anaCutsOrig = ana[1]
              anaCuts = anaCutsOrig+" && "+cutString
              analysisNames.append([anaName,anaCuts])
              if doSplit:
                anaCuts = anaCutsOrig+" && "+cutStringFail
                analysisNames.append([anaName+"Fail",anaCuts])
            
            print("analysisNames: {0}".format(analysisNames))
            threads.append(
             ThreadedCardMaker(
               #__init__ args:
               directory,
               analysisNames,
               appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),dataDict[p],
               controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,
               controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,sigInjectMass=args.signalInjectMass,
               energyStr=p,
               #write args:
               outfilename=outDir+nameString+"_"+p+"_"+str(i)+".txt",lumi=i,
               cutString=""
             )
            )
            for ind,iLen,iInd in reversed(zip(indexList,lenList,range(nCutVars))):
              if ind < iLen-1:
                indexList[iInd] += 1
                break
              else:
                indexList[iInd] = 0
            if sum(indexList) == 0:
              notDone = False

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

  simpleScripts = True

  if simpleScripts:
    shutil.copy("etc/notlxbatch_simple.sh",outDir+"notlxbatch.sh")
    shutil.copy("etc/lxbatch_simple.sh",outDir+"lxbatch.sh")
  else:
    shutil.copy("etc/notlxbatch.sh",outDir+"notlxbatch.sh")
    shutil.copy("etc/lxbatch.sh",outDir+"lxbatch.sh")

  shutil.copy("etc/run.sh",outDir+"run.sh")
  shutil.copy("etc/uftrig01.sh",outDir+"uftrig01.sh")
      
  shutil.copy("etc/nuisanceDiff.py",outDir+"diffNuisances.py")
  shutil.copy("etc/fitNormsToText_mlfit.py",outDir+"fitNormsToText_mlfit.py")
  shutil.copy("etc/myNuisancePrinter.py",outDir+"myNuisancePrinter.py")
  shutil.copy("etc/hpcTemplate.sh",outDir+"hpcTemplate.sh")
  shutil.copy("etc/runHPC.sh",outDir+"runHPC.sh")
  shutil.copy("etc/getStatus2.sh",outDir+"getStatus2.sh")
  shutil.copy("etc/gof.sh",outDir+"gof.sh")
  shutil.copy("etc/gofHPC_Template.sh",outDir+"gofHPC_Template.sh")
  shutil.copy("etc/runHPC_GOF.sh",outDir+"runHPC_GOF.sh")
  shutil.copy("etc/runHPC_Compat.sh",outDir+"runHPC_Compat.sh")
  shutil.copy("etc/compatHPC_Template.sh",outDir+"compatHPC_Template.sh")
