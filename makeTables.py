#!/usr/bin/env python

import math
from math import sqrt
import ROOT as root
from helpers import *
import numpy
import glob
import re

from xsec import *
import makeCards
import makeShapePlots

from ROOT import gSystem
gSystem.Load('libRooFit')

root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

cats = [
    "Jets01PassPtG10BB",
    "Jets01PassPtG10BO",
    "Jets01PassPtG10BE",
    "Jets01PassPtG10OO",
    "Jets01PassPtG10OE",
    "Jets01PassPtG10EE",
    "Jets01FailPtG10BB",
    "Jets01FailPtG10BO",
    "Jets01FailPtG10BE",
    "Jets01FailPtG10OO",
    "Jets01FailPtG10OE",
    "Jets01FailPtG10EE",
    "Jet2CutsVBFPass",
    "Jet2CutsGFPass",
    "Jet2CutsFailVBFGF",
  ]

def getFWHM(catName,energyStr):

  x = root.RooRealVar("dimuonMass","dimuonMass",110.,160.)
  w = root.RooWorkspace("ws"+catName+energyStr)
  makeCards.makePDFSigNew(catName+energyStr,"sigPdf",x,125.0,getattr(w,"import"))
  #w.Print()
  pdf = w.pdf("sigPdf_hmm"+energyStr)
  return calcFWHM(pdf,x,110,160,0.02)

def getNSig(catName,energyStr):
  higgsMass = 125.
  for prodMode in ["gg","vbf","wh","zh"]:
    effReader = makeCards.effReaderFromFile('fitresults',
                                  prodMode,
                                  energyStr,
                                  catName)
    efficiencies = effReader.getEff()
    eff = efficiencies[ '{0:.1f}'.format(higgsMass) ]
  
    #print ' ####### higgsMassStr = %s' % higgsMassStr 
    tmpName = prodMode+"Hmumu"+"125"+"_"+energyStr
    #print ' #######  tmpName = %s' % tmpName
    xs = xsec[tmpName]
    lumi = lumiDict[energyStr]*1000.
    #print eff,xs,lumi
    return eff*xs*lumi

class DataCardGetter(makeShapePlots.ShapePlotter):

  def __init__(self,filename,fitDir,energyStr,fwhm):
  #def __init__(self,filename,outDir,fitDir,titleMap,signalInject=20.,binWidthOverride=0,energyStr=None):
    self.filename = filename # root part of data card
    self.textFileName = os.path.splitext(filename)[0]+".txt"
    self.fitFileName = fitDir+"/"+os.path.split(self.textFileName)[1]+".root"
    self.processNameMap, self.params, self.normErrMap = getattr(self,"readCard")(self.textFileName)

    self.energyStr = ""
    if energyStr:
      self.energyStr = energyStr
    else:
      tmpMatch = re.search(r"([\w]*)_(.+)_([.0-9]+)\.root",filename)
      if tmpMatch:
        self.energyStr = tmpMatch.group(2)
          
    self.data = {}
    self.f = root.TFile(filename)
    for channelKey in self.f.GetListOfKeys():
      if channelKey.GetClassName() != "RooWorkspace":
        continue
      channelNameOrig = channelKey.GetName()
      if not energyStr in channelNameOrig:
        continue
      channelName = re.sub("[\d]+TeV","",channelNameOrig)
      channelTitle = channelName
      self.channelName = channelName
      self.channelNameOrig = channelNameOrig
      channelWS = channelKey.ReadObj()
      mMuMu = channelWS.var("dimuonMass_CMShmm")
      bakPDF = channelWS.pdf("bak_"+channelNameOrig+"_CMShmm")
      sigPDF = channelWS.pdf("sigPdf_hmm"+energyStr+"_"+channelNameOrig+"_CMShmm")
      data_obs = channelWS.data("data_obs_"+channelNameOrig+"_CMShmm")
      
      rooDataTitle = data_obs.GetTitle()
      binWidth = 1
      if "Jet2" in channelName or "VBF" in channelName:
          binWidth *= 2.5
      elif "BO" in channelName:
          binWidth *= 1
      elif "BE" in channelName:
          binWidth *= 2.5
      elif "OO" in channelName:
          binWidth *= 2.5
      elif "OE" in channelName:
          binWidth *= 2.5
      elif "EE" in channelName:
          binWidth *= 2.5
      elif "FF" in channelName:
          binWidth *= 2.5
      elif "CC" in channelName:
          binWidth *= 2.5
      elif "BB" in channelName:
          binWidth *= 1

      binning = mMuMu.getBinning()
      xlow = binning.lowBound()
      xhigh = binning.highBound()
      #mMuMu.setRange("shapePlot",xlow,xhigh)
      mMuMu.setBins(int((xhigh-xlow)/binWidth))
      mMuMu.SetTitle("m_{#mu#mu} [GeV]")

      # Get Fit Result
      fr = None
      try:
        rf = root.TFile(self.fitFileName)
        if rf.IsZombie() or not rf.IsOpen():
          raise IOError("ROOT File not open or Zombie")
        self.fitFile = rf
        fr = rf.Get("fit_s")
        #fr = rf.Get("fit_b")
      except Exception as e:
        print("Warning, Couldn't find ML fit file: {0}\n{1}".format(self.fitFileName,e))
        # Backup Fit
        fr = bakPDF.fitTo(data_obs,root.RooFit.Save(True),root.RooFit.PrintLevel(-1))

      # Signal Stuff
      nSignal = 0.
      for key in self.processNameMap[channelNameOrig]:
        if not "bak" in key:
          nSignal += self.processNameMap[channelNameOrig][key]
      self.nSig = nSignal

      argList = fr.floatParsFinal()
      bestFitMu = argList.find("r").getVal()
      nSignal *= bestFitMu

      nTotalObs = data_obs.sumEntries()
      nBkgTotal = nTotalObs - nSignal
      #print nTotalObs,nBkgTotal,nSignal

      #Set the PDF pars value from the FitResults
      setPDFfromFR(fr,bakPDF,data_obs)

      # Do the work
      obsVarSet = root.RooArgSet(mMuMu)
      fwhmRangeName = "myIntRange_{0}".format(channelNameOrig)
      mMuMu.setRange(fwhmRangeName,125-0.5*fwhm,125.+0.5*fwhm) 
      pdfFrac = bakPDF.createIntegral(obsVarSet,obsVarSet,fwhmRangeName)
      self.nBkg = pdfFrac.getVal() * nBkgTotal

def getnBkgInFWHM(catName,energyStr,fwhm):
  fn = "/raid/raid7/jhugon/higgsDataCards/20140501/statsCards/" + catName + "_" + energyStr + "_125.0.root"
  fitDir = "/raid/raid7/jhugon/higgsDataCards/20140501/statsInput/"
  dcg = DataCardGetter(fn,fitDir,energyStr,fwhm)
  return dcg.nBkg, dcg.nSig

def getSignalCompFractions(catName,energyStr):
  fn = "/raid/raid7/jhugon/higgsDataCards/20140501/statsCards/" + catName + "_" + energyStr + "_125.0.root"
  fitDir = "/raid/raid7/jhugon/higgsDataCards/20140501/statsInput/"
  dcg = DataCardGetter(fn,fitDir,energyStr,2.5)
  gfNum = 0.
  vbfNum = 0.
  whNum = 0.
  zhNum = 0.
  for key in dcg.processNameMap[dcg.channelNameOrig]:
    if "ggH" in key:
      gfNum = dcg.processNameMap[dcg.channelNameOrig][key]
    elif "qqH" in key:
      vbfNum = dcg.processNameMap[dcg.channelNameOrig][key]
    elif "WH" in key:
      whNum = dcg.processNameMap[dcg.channelNameOrig][key]
    elif "ZH" in key:
      zhNum = dcg.processNameMap[dcg.channelNameOrig][key]
  totalNum = gfNum + vbfNum + whNum + zhNum
  gfFrac = gfNum/totalNum
  vbfFrac = vbfNum/totalNum
  whFrac = whNum/totalNum
  zhFrac = zhNum/totalNum
  return gfFrac,vbfFrac,whFrac,zhFrac

def getnDataInFWHM(catName,energyStr,fwhm):
  glbString = "/cms/data/store/user/jhugon/hmumu/stage2/SingleMuRun2011*.root"
  if energyStr == "8TeV":
    glbString = "/cms/data/store/user/jhugon/hmumu/stage2/SingleMuRun2012*.root"
  tree = root.TChain("outtree"+catName)
  for fn in glob.glob(glbString):
    tree.AddFile(fn)
  entryListName = "entryList"+catName+energyStr 
  tree.Draw(">> "+entryListName,"{0} > {1} && {0} < {2}".format("dimuonMass",125-0.5*fwhm,125.+0.5*fwhm),"entryList")
  entryList = root.gROOT.FindObject(entryListName)
  result = entryList.GetN()
  return result

def getEffs(catName,energyStr):
  gf = -1.
  vbf = -1.
  wh = -1.
  zh = -1. 
  fn = "/raid/raid7/jhugon/higgsDataCards/20140501/statsCards/" + catName + "_" + energyStr + "_125.0.txt"
  f = open(fn)
  foundSigEffs = False
  for l in f:
    if foundSigEffs:
      gfMatch = re.search(r"^#    ggHmumu125_[78]TeV: ([.0-9]+)",l)
      vbfMatch = re.search(r"^#    vbfHmumu125_[78]TeV: ([.0-9]+)",l)
      whMatch = re.search(r"^#    whHmumu125_[78]TeV: ([.0-9]+)",l)
      zhMatch = re.search(r"^#    zhHmumu125_[78]TeV: ([.0-9]+)",l)
      if gfMatch:
        gf = float(gfMatch.group(1))
      if vbfMatch:
        vbf = float(vbfMatch.group(1))
      if whMatch:
        wh = float(whMatch.group(1))
      if zhMatch:
        zh = float(zhMatch.group(1))
    elif re.search(r"^#  Signal Efficiencies:",l):
      foundSigEffs = True
  f.close()
  return gf,vbf,wh,zh


for energyStr in ["7TeV","8TeV"]:
  for cat in cats:
    fwhm     = getFWHM(cat,energyStr)
    nBkgInFWHM, nSig = getnBkgInFWHM(cat,energyStr,fwhm)
    nDataInFWHM = getnDataInFWHM(cat,energyStr,fwhm)
    print "{0:25} {1:3.1f} & {2:5.2f} & {3:9.1f} & {4:9.0f}".format(cat+energyStr,fwhm,nSig,nBkgInFWHM,nDataInFWHM)

print "\nNow gf, vbf, wh, zh efficiency * acceptance\n"
    
for energyStr in ["7TeV","8TeV"]:
  for cat in cats:
    gf,vbf,wh,zh     = getEffs(cat,energyStr)
    #print "{0:25} {1:5.1f}\% & {2:5.1f}\% & {3:5.1f}\% & {4:5.1f}\%".format(cat+energyStr,gf*100.,vbf*100.,wh*100.,zh*100.)
    print "{0:25} {1:5.1f}\% & {2:5.1f}\% & {3:5.1f}\%".format(cat+energyStr,gf*100.,vbf*100.,(wh+zh)/2.*100.)
   
print "\nNow fractional signal composition\n"
    
for energyStr in ["7TeV","8TeV"]:
  for cat in cats:
    gf,vbf,wh,zh     = getSignalCompFractions(cat,energyStr)
    #print "{0:25} {1:5.1f}\% & {2:5.1f}\% & {3:5.1f}\% & {4:5.1f}\%".format(cat+energyStr,gf*100.,vbf*100.,wh*100.,zh*100.)
    print "{0:25} {1:5.1f}\% & {2:5.1f}\% & {3:5.1f}\% & {4:5.1f}\%".format(cat+energyStr,gf*100.,vbf*100.,wh*100.,zh*100.)
    
