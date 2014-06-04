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

PRELIMINARYSTRING="CMS"

class ShapePlotter:
  def __init__(self,filename,outDir,fitDir,titleMap,signalInject=20.,binWidthOverride=0,energyStr=None):
    self.signalInject=signalInject
    self.rmpList = []
    self.titleMap = titleMap
    self.filename = filename
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
          
    self.lumi = float(lumiDict[self.energyStr])
    self.lumiStr = "L = {0:.1f} fb^{{-1}}".format(self.lumi)
      
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
      if titleMap.has_key(channelName):
        channelTitle = titleMap[channelName]
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
      if binWidthOverride > 0:
        binWidth = binWidthOverride

      binning = mMuMu.getBinning()
      xlow = binning.lowBound()
      xhigh = binning.highBound()
      #mMuMu.setRange("shapePlot",xlow,xhigh)
      mMuMu.setBins(int((xhigh-xlow)/binWidth))
      mMuMu.SetTitle("m_{#mu#mu} [GeV]")

      saveName = outDir+os.path.splitext(os.path.split(self.filename)[1])[0]+'_'+channelName
      saveName = re.sub(r"_[0-9P]+TeV_","_"+self.energyStr+"_",saveName)
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
        if not "bak" in key:
          nSignal += self.processNameMap[channelNameOrig][key]
      nSignal *= signalInject
      legEntrySignal = "SM Higgs#times{0:.0f}".format(signalInject)

      #Set the PDF pars value from the FitResults
      setPDFfromFR(fr,bakPDF,data_obs)

      #Plot Time
      rmp = RooModelPlotter(mMuMu,bakPDF,data_obs,fr,
                            channelTitle,self.energyStr.replace("TeV"," TeV"),self.lumi,
                            nSignal=nSignal,signalPdf=sigPDF,
                            legEntrySignal=legEntrySignal,
                            preliminaryString=PRELIMINARYSTRING
                            )
      rmp.draw(saveName)
      #rmp.drawWithParams(saveName+"_params",["mixParam","bwWidth","bwmZ","expParam"])

      #Pull Distribution Time
      rmp.drawPulls(saveName+"_pulls")

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
  for energyStr in ["7TeV","8TeV"]:
    for fn in glob.glob(dataDir+"*.root"):
      if re.search("P[\d.]+TeV",fn):
          continue

      if ("125" not in fn):
          continue

      skip = True
      if ("Jets01FailPtG10" in fn):
        skip = False
      if ("Jets01PassPtG10" in fn):
        skip = False
      if ("Jet2CutsFailVBFGF"in fn):
        skip = False
      if ("Jet2CutsVBFPass"in fn):
        skip = False
      if ("Jet2CutsGFPass" in fn):
        skip = False
      if ("All" in fn):
        skip = True

      #if "7P8TeV" in fn:
      #  skip = True

      if (skip):
        continue
      
      print fn
      s = ShapePlotter(fn,outDir,fitDir,TITLEMAP,signalInject=args.signalInject,binWidthOverride=args.binWidthOverride,energyStr=energyStr)
      shapePlotterList.append(s)
