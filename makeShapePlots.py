#!/usr/bin/env python


import optparse
parser = optparse.OptionParser(description="Makes Shape Diagnostic Plots from Datacards")
parser.add_option("--signalInject", help="Sets a caption saying that signal was injected with strength",type=float,default=20.0)
parser.add_option("-b","--binWidthOverride", help="Overrides the default bin widths and sets all binning to this widht [GeV]",type=float,default=0.0)
#parser.add_option("--plotSignalStrength", help="Plots a signal bump with this strength",type=float,default=5.0)
#parser.add_option("--plotSignalBottom", help="Plots a signal bump on the bottom (bool)",action="store_true",default=True)
#parser.add_option("--signalInjectMass", help="Mass For Injected Signal",type=float,default=125.0)
#parser.add_option("-r","--rebinOverride", help="Rebin All plots with this rebinning, overriding all internal configuration",type=int,default=0)
args, fakeargs = parser.parse_args()

from helpers import *
from xsec import *
import math
import os.path
import glob
import random

from ROOT import gSystem
gSystem.Load('libRooFit')

root.gErrorIgnoreLevel = root.kWarning
#root.gROOT.SetBatch(True)
root.gStyle.SetOptStat(0)

#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

class ShapePlotter:
  def __init__(self,filename,outDir,fitDir,titleMap,signalInject=20.,binWidthOverride=0):
    self.signalInject=signalInject
    self.rmpList = []
    self.titleMap = titleMap
    self.filename = filename
    self.textFileName = os.path.splitext(filename)[0]+".txt"
    self.fitFileName = fitDir+"/"+os.path.split(self.textFileName)[1]+".root"
    self.processNameMap, self.params, self.normErrMap = getattr(self,"readCard")(self.textFileName)

    self.lumi = -1
    self.lumiStr = ""
    self.energyStr = ""
    tmpMatch = re.search(r"([\w]*)_(.+)_([.0-9]+)\.root",filename)
    if tmpMatch:
      self.energyStr = tmpMatch.group(2)
      self.lumi = float(lumiDict[self.energyStr])
        
      #self.lumi = float(tmpMatch.group(3))
      self.lumiStr = "L = {0:.1f} fb^{{-1}}".format(self.lumi)
      
    self.data = {}
    self.f = root.TFile(filename)
    for channelKey in self.f.GetListOfKeys():
      if channelKey.GetClassName() != "RooWorkspace":
        continue
      channelNameOrig = channelKey.GetName()
      channelName = re.sub("[\d]+TeV","",channelNameOrig)
      channelTitle = channelName
      if titleMap.has_key(channelName):
        channelTitle = titleMap[channelName]
      channelWS = channelKey.ReadObj()
      mMuMu = channelWS.var("dimuonMass")
      bakPDF = channelWS.pdf("bak")
      sigPDF = channelWS.pdf("ggH")
      data_obs = channelWS.data("data_obs")
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
      if binWidthOverride > 0:
        binWidth = binWidthOverride

      binning = mMuMu.getBinning()
      xlow = binning.lowBound()
      xhigh = binning.highBound()
      #mMuMu.setRange("shapePlot",xlow,xhigh)
      mMuMu.setBins(int((xhigh-xlow)/binWidth))
      mMuMu.SetTitle("M(#mu#mu) [GeV/c^{2}]")

      saveName = outDir+os.path.splitext(os.path.split(self.filename)[1])[0]+'_'+channelName
      saveName = re.sub(r"([\d]+)\.[\d]+",r"\1",saveName)

      # Get Fit Result
      fr = None
      try:
        rf = root.TFile(self.fitFileName)
        if rf.IsZombie() or not rf.IsOpen():
          raise IOError("ROOT File not open or Zombie")
        self.fitFile = rf
        #fr = rf.Get("fit_s")
        fr = rf.Get("fit_b")
      except Exception as e:
        print("Warning, Couldn't find ML fit file: {0}\n{1}".format(self.fitFileName,e))
        # Backup Fit
        fr = bakPDF.fitTo(data_obs,root.RooFit.Save(True),root.RooFit.PrintLevel(-1))

      # Signal Stuff
      nSignal = 0.
      for key in self.processNameMap[channelNameOrig]:
        if key != "bak":
          nSignal += self.processNameMap[channelNameOrig][key]
      nSignal *= signalInject
      signalLegEntry = "SM Higgs#times{0:.0f}".format(signalInject)

      #Set the PDF pars value from the FitResults
      setPDFfromFR(fr,bakPDF,data_obs)

      #Plot Time
      rmp = RooModelPlotter(mMuMu,bakPDF,data_obs,fr,
                            channelTitle,self.energyStr,self.lumi,
                            nSignal=nSignal,signalPdf=sigPDF,
                            signalLegEntry=signalLegEntry,
                            caption1="Analysis A"
                            )
      #rmp.draw(saveName)
      rmp.drawWithParams(saveName+"_params",["mixParam","bwWidth","bwmZ","expParam"])

      #Pull Distribution Time
      saveNameSplit = os.path.split(saveName)
      saveNamePulls = saveNameSplit[0]+"/"+"pulls_"+saveNameSplit[1]
      rmp.drawPulls(saveNamePulls)

      self.rmpList.append(rmp)

  def readCard(self,fn):
    f = open(fn)
    foundBin = False
    binList = []
    processList = []
    rateList = []
    paramMap = {}
    normErrMap = {}
    for line in f:
      if re.search("^bin",line):
        if foundBin:
          m =  re.findall("[\s]+[\w]+",line)
          binList.extend([i for i in m])
        else:
          foundBin = True
      if re.search("^process[\s]+[a-zA-Z]+",line):
          m =  re.findall("[\s]+[\w]+",line)
          processList.extend([i for i in m])
      if re.search("^rate[\s]+[-+eE.0-9]+",line):
          m =  re.findall("[\s]+[-+eE.0-9]+",line)
          rateList.extend([float(i) for i in m])
      paramMatch = re.search(r"([a-zA-Z0-9_]+)[\s]+param[\s]+([-.+eE0-9]+)[\s]+([-.+eE0-9]+)",line)
      if paramMatch:
        gs = paramMatch.groups()
        paramMap[gs[0]] = [gs[1],gs[2]]
      normMatch = re.search(r"bkN([a-zA-Z0-9_]+)[\s]+gmN[\s]+([-.+eE0-9]+)[\s-]+([.0-9]+)[\s-]+",line)
      if normMatch:
        gs = normMatch.groups()
        normErrMap[gs[0]] = 1.0/sqrt(float(gs[1]))
    binList = [x.replace(r" ","") for x in binList]
    processList = [x.replace(r" ","") for x in processList]
    result = {}
    for i in binList:
      if not result.has_key(i):
        result[i] = {}
    for b,p,r in zip(binList,processList,rateList):
      result[b][p] = r
    return result, paramMap, normErrMap

titleMap = {
  "AllCat":"All Categories Comb.",
  "IncCat":"Non-VBF Categories Comb.",
  "VBFCat":"VBF Categories Comb.",

  "IncPresel":"Non-VBF Preselection",
  "VBFPresel":"VBF Preselection",

  "Pt0to30":"p_{T}^{#mu#mu} #in [0,30]",
  "Pt30to50":"p_{T}^{#mu#mu} #in [30,50]",
  "Pt50to125":"p_{T}^{#mu#mu} #in [50,125]",
  "Pt125to250":"p_{T}^{#mu#mu} #in [125,250]",
  "Pt250":"p_{T}^{#mu#mu}>250",

  "VBFLoose":"VBFL",
  "VBFMedium":"VBFM",
  "VBFTight":"VBFT",
  "VBFVeryTight":"VBFVT",

  "BDTCut":"BDT Cut Combination",
  "IncBDTCut":"Non-VBF BDT Cut",
  "VBFBDTCut":"VBF BDT Cut",

  "BDTCutCat":"BDT Cut Cat. Combination",
  "IncBDTCutCat":"Non-VBF BDT Cut",
  "VBFBDTCutCat":"VBF BDT Cut",

  "IncPreselCat":"Non-VBF Cat. Preselection",
  "VBFPreselCat":"VBF Cat. Preselection",

  "IncBDTCutBB":"Non-VBF BDT Cut BB",
  "IncBDTCutBO":"Non-VBF BDT Cut BO",
  "IncBDTCutBE":"Non-VBF BDT Cut BE",
  "IncBDTCutOO":"Non-VBF BDT Cut OO",
  "IncBDTCutOE":"Non-VBF BDT Cut OE",
  "IncBDTCutEE":"Non-VBF BDT Cut EE",
  "IncBDTCutNotBB":"Non-VBF BDT Cut !BB",
  "VBFBDTCutBB":"VBF BDT Cut BB",
  "VBFBDTCutNotBB":"VBF BDT Cut !BB",
  "IncPreselBB":"Non-VBF Preselection BB",
  "IncPreselBO":"Non-VBF Preselection BO",
  "IncPreselBE":"Non-VBF Preselection BE",
  "IncPreselOO":"Non-VBF Preselection OO",
  "IncPreselOE":"Non-VBF Preselection OE",
  "IncPreselEE":"Non-VBF Preselection EE",
  "IncPreselNotBB":"Non-VBF Preselection !BB",
  "VBFPreselBB":"VBF Preselection BB",
  "VBFPreselNotBB":"VBF Preselection !BB",

  "IncPreselPtG10BB":"Non-VBF BB",
  "IncPreselPtG10BO":"Non-VBF BO",
  "IncPreselPtG10BE":"Non-VBF BE",
  "IncPreselPtG10OO":"Non-VBF OO",
  "IncPreselPtG10OE":"Non-VBF OE",
  "IncPreselPtG10EE":"Non-VBF EE",
  "IncPreselPtG10NotBB":"Non-VBF !BB",

  "IncPreselPtG":"Non-VBF Not Combined",

  "Jets01PassPtG10BB": "0,1-Jet Tight BB",
  "Jets01PassPtG10BO": "0,1-Jet Tight BO",
  "Jets01PassPtG10BE": "0,1-Jet Tight BE",
  "Jets01PassPtG10OO": "0,1-Jet Tight OO",
  "Jets01PassPtG10OE": "0,1-Jet Tight OE",
  "Jets01PassPtG10EE": "0,1-Jet Tight EE",
  "Jets01PassCatAll" : "0,1-Jet Tight Combination",
                        
  "Jets01FailPtG10BB": "0,1-Jet Loose BB",
  "Jets01FailPtG10BO": "0,1-Jet Loose BO",
  "Jets01FailPtG10BE": "0,1-Jet Loose BE",
  "Jets01FailPtG10OO": "0,1-Jet Loose OO",
  "Jets01FailPtG10OE": "0,1-Jet Loose OE",
  "Jets01FailPtG10EE": "0,1-Jet Loose EE",
  "Jets01FailCatAll" : "0,1-Jet Loose Combination",
                        
  "Jets01SplitCatAll": "0,1-Jet Combination",


  "Jet2CutsVBFPass":"2-Jet VBF Tight",
  "Jet2CutsGFPass":"2-Jet GF Tight",
  "Jet2CutsFailVBFGF":"2-Jet Loose",

  "Jet2SplitCutsGFSplit" : "2-Jet Combination",
  "CombSplitAll" : "H#rightarrow#mu#mu Combination",
}
        
if __name__ == "__main__":
  root.gROOT.SetBatch(True)

  dataDir = "statsCards/"
  outDir = "shapes/"
  fitDir = "statsInput/"

  plotRange= [110.,160]
  #plotRange= []
  normRange = [110.,120.,130.,160]

  rebin=1

  shapePlotterList = []
  for fn in glob.glob(dataDir+"*.root"):
    if re.search("P[\d.]+TeV",fn):
        continue

    if ("125" not in fn):
    #if ("148" not in fn):
    #if ("145" not in fn):
    #if ("147" not in fn):
        continue

    skip = False
    #if ("Jets01FailPtG10BB" in fn):
    #  skip = False
    #if ("Jets01PassPtG10BB" in fn):
    #  skip = False
    #if ("Jets01PassPtG10BO" in fn):
    #  skip = False
    #if ("Jet2CutsVBFPass"in fn):
    #  skip = False
    #if ("Jet2CutsGFPass" in fn):
    #  skip = False
    if ("CombSplitAll" in fn):
      skip = False

    if (skip):
      continue
    
    print fn
    s = ShapePlotter(fn,outDir,fitDir,titleMap,signalInject=args.signalInject,binWidthOverride=args.binWidthOverride)
    shapePlotterList.append(s)
