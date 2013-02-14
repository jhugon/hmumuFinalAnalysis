#!/usr/bin/env python

from xsec import *
from helpers import *
import ROOT as root
import os
import sys

dataDir = "input/preApproveSample110to150/"
outDir = "output/"

RUNPERIOD="8TeV"
LUMI=lumiDict[RUNPERIOD]

scaleHiggsBy = 100.

LOGY=False
integralPlot=False
MCErrors=True
#PULLTYPE="adrian1"
PULLTYPE="pullMC"
allHiggsTogether = True
ylimitsRatio = [-4,4]
if PULLTYPE=="adrian1":
  ylimitsRatio = [-1.5,1.5]
elif PULLTYPE=="ratio":
  ylimitsRatio = [0.5,1.5]

#anotateText = "80 GeV < m_{#mu#mu} < 160 GeV; p_{T,#mu#mu}<20 GeV"
#anotateText = "110 GeV < m_{#mu#mu} < 160 GeV; p_{T,#mu#mu}<20 GeV"
anotateText = "110 GeV < m_{#mu#mu} < 150 GeV"
#anotateText = "80 GeV < m_{#mu#mu} < 160 GeV"
#anotateText = "VBF Preselection"

urLegendPos = [0.70,0.67,0.9,0.9]
ulLegendPos = [0.20,0.67,0.4,0.9]
ucLegendPos = [0.46,0.67,0.64,0.9]
lcLegendPos = [0.46,0.35,0.64,0.63]
llLegendPos = [0.20,0.35,0.4,0.63]
ccLegendPos = [0.46,0.47,0.64,0.7]
stdLegendPos = urLegendPos

#histDirs = ["","4GeVWindow/","PtDiMu100/","VBFPresel/","IncPresel/","NotBlindWindow/"]
histDirs = ["NotBlindWindow/"]
histDirs = ["VBFPreselDiMuPtL20/","IncPreselDiMuPtL20/"]
histDirs = ["VBFPresel/","IncPresel/"]
histDirs = ["IncPresel/"]
#histDirs = ["IncBDTCutBB/","VBFBDTCut/"]
#histDirs = ["VBFBDTCut/"]
histDirs = ["VBFPresel/"]

root.gErrorIgnoreLevel = root.kWarning

histNames = {}
if RUNPERIOD=="7TeV":
  print "Using 7TeV Settings"
  if len(histDirs) == 1 and histDirs[0] == "IncPresel/":
    print "Doing IncPresel"
    histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV/c^{2}]","xlimits":[110.0,149.99],"rebin":2,"ylimits":[0.1,5e5]}
    #histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV/c^{2}]","xlimits":[80.0,150.0],"rebin":2}
  
    #histNames["ptDiMu"] = {"xlabel":"p_{T,#mu#mu} [GeV/c]","xlimits":[0.0,20.0],"rebin":1,'ylimits':[0.,3.25e3]}
    histNames["ptDiMu"] = {"xlabel":"p_{T,#mu#mu} [GeV/c]","xlimits":[0.0,200.0],"rebin":2,"ylimits":[0.1,1e5]}

    histNames["yDiMu"] = {"xlabel":"y_{#mu#mu}","xlimits":[-2.2,2.2],"rebin":2,"ylimits":[0.1,3e6]}
    histNames["cosThetaStar"] = {"xlabel":"cos(#theta^{*})","xlimits":[-1.0,1.0],"rebin":2,"ylimits":[0.1,1e7]}
    histNames["BDTHistMuonOnly"] = {"xlabel":"BDT (Non-VBF Category)","xlimits":[-0.55,0.2],"rebin":2,"ylimits":[1e-2,1e8],'vertLines':{"8TeV":-0.55,"7TeV":-0.42}}
  elif len(histDirs) == 1 and histDirs[0] == "VBFPresel/":
    histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV/c^{2}]","xlimits":[110.0,149.99],"rebin":4,"ylimits":[0.,25]}
    #histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV/c^{2}]","xlimits":[80.0,150.0],"rebin":2}
  
    #histNames["ptDiMu"] = {"xlabel":"p_{T,#mu#mu} [GeV/c]","xlimits":[0.0,20.0],"rebin":1}
    histNames["ptDiMu"] = {"xlabel":"p_{T,#mu#mu} [GeV/c]","xlimits":[0.0,200.0],"rebin":5,"ylimits":[0.,25]}

    histNames["yDiMu"] = {"xlabel":"y_{#mu#mu}","xlimits":[-2.2,2.2],"rebin":4,"ylimits":[0.,30]}
    histNames["cosThetaStar"] = {"xlabel":"cos(#theta^{*})","xlimits":[-1.0,1.0],"rebin":5,"ylimits":[0.,50]}
  
    histNames["mDiJet"] = {"xlabel":"m_{jj} [GeV/c^{2}]","xlimits":[300.0,1400.0],"rebin":20,'ylimits':[0.1,50]}
    histNames["ptDiJet"] = {"xlabel":"p_{T,jj} [GeV/c]","xlimits":[0.0,400.0],"rebin":5,'ylimits':[0.1,35]}
    histNames["deltaEtaJets"] = {"xlabel":"#Delta#eta_{jj}","xlimits":[3.0,7.5],"rebin":2,"ylimits":[0.,50]}
  
  
    histNames["yDiJet"] = {"xlabel":"y_{jj}","xlimits":[-3.0,3.0],"rebin":4,"ylimits":[0.,40]}
  
    histNames["BDTHistVBF"] = {"xlabel":"BDT (VBF Category)","xlimits":[-0.4,0.25],"rebin":4,"ylimits":[0.0,55],'vertLines':{"8TeV":-0.04,"7TeV":-0.03}}
    histNames["ptmiss"] = {"xlabel":"p_{T}^{Miss} [GeV/c]","xlimits":[0,200],"leg":stdLegendPos,'rebin':2,'ylimits':[0.,35]}
  else:
    anotateText = "p_{T}(#mu#mu)>10 GeV"
    histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV/c^{2}]","xlimits":[110.0,149.99],"rebin":2,'ylimits':[0.1,5e5]}
    #anotateText = ""
    #histNames["nVtx"] = {"xlabel":"N_{vtx}","xlimits":[0,25],"leg":stdLegendPos,'ylimits':[0.1,3e3]}
elif RUNPERIOD=="8TeV":
  print "Using 8TeV Settings"
  if len(histDirs) == 1 and histDirs[0] == "IncPresel/":
    print "Doing IncPresel"
    histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV/c^{2}]","xlimits":[110.0,149.99],"rebin":2,'ylimits':[0.1,1e7]}
    #histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV/c^{2}]","xlimits":[80.0,150.0],"rebin":2}
    #histNames["ptDiMu"] = {"xlabel":"p_{T,#mu#mu} [GeV/c]","xlimits":[0.0,20.0],"rebin":1,'ylimits':[0.0,11000]}
    histNames["ptDiMu"] = {"xlabel":"p_{T,#mu#mu} [GeV/c]","xlimits":[0.0,200.0],"rebin":2,'ylimits':[0.1,1e5]}
    histNames["yDiMu"] = {"xlabel":"y_{#mu#mu}","xlimits":[-2.2,2.2],"rebin":2,"ylimits":[0.1,5e7]}
    
    histNames["cosThetaStar"] = {"xlabel":"cos(#theta^{*})","xlimits":[-1.0,1.0],"rebin":2,"ylimits":[0.1,5e7]}
    histNames["BDTHistMuonOnly"] = {"xlabel":"BDT (Non-VBF Category)","xlimits":[-1.0,0.2],"rebin":2,"ylimits":[1e-1,1e5],'vertLines':{"8TeV":-0.55,"7TeV":-0.42}}
  elif len(histDirs) == 1 and histDirs[0] == "VBFPresel/":
    histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV/c^{2}]","xlimits":[110.0,149.99],"rebin":4,"ylimits":[0.0,150]}
    #histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV/c^{2}]","xlimits":[80.0,150.0],"rebin":2}
    histNames["ptDiMu"] = {"xlabel":"p_{T,#mu#mu} [GeV/c]","xlimits":[0.0,200.0],"rebin":5,"ylimits":[0.0,160]}
    histNames["yDiMu"] = {"xlabel":"y_{#mu#mu}","xlimits":[-2.2,2.2],"rebin":5,"ylimits":[0.0,200]}
    
    histNames["cosThetaStar"] = {"xlabel":"cos(#theta^{*})","xlimits":[-1.0,1.0],"rebin":5,"ylimits":[0.,220]}
    histNames["mDiJet"] = {"xlabel":"m_{jj} [GeV/c^{2}]","xlimits":[300.0,1400.0],"rebin":20,"ylimits":[0.,250]}
    histNames["ptDiJet"] = {"xlabel":"p_{T,jj} [GeV/c]","xlimits":[0.0,400.0],"rebin":5,"ylimits":[0.1,200]}
    histNames["deltaEtaJets"] = {"xlabel":"#Delta#eta_{jj}","xlimits":[3.0,7.5],"rebin":2,"ylimits":[0.1,300]}
    histNames["yDiJet"] = {"xlabel":"y_{jj}","xlimits":[-3.0,3.0],"rebin":4,"ylimits":[0.0,300]}
    histNames["ptmiss"] = {"xlabel":"p_{T}^{Miss} [GeV/c]","xlimits":[0,200],"leg":stdLegendPos,"rebin":2,"ylimits":[0.1,160]}
    histNames["BDTHistVBF"] = {"xlabel":"BDT (VBF Category)","xlimits":[-0.4,0.25],"rebin":4,"ylimits":[0.,350],'vertLines':{"8TeV":-0.04,"7TeV":-0.03}}
  else:
    histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV/c^{2}]","xlimits":[110.0,149.99],"rebin":2,'ylimits':[0.1,1e5]}
    anotateText = "Non-VBF, BB, p_{T}(#mu#mu)>10 GeV"
    #histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV/c^{2}]","xlimits":[110.0,149.99],"rebin":2,'ylimits':[0.1,6e4]}
    #histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV/c^{2}]","xlimits":[110.0,149.99],"rebin":5,'ylimits':[0.1,6e2]}
    #histNames["nVtx"] = {"xlabel":"N_{vtx}","xlimits":[0,40],"leg":stdLegendPos,'ylimits':[0.0,8e3]}
else:
  print "Using Other Settings"
  histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV/c^{2}]","xlimits":[110.0,149.99],"rebin":2}
  #histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV/c^{2}]","xlimits":[80.0,150.0],"rebin":2}

  
  """
  #histNames["ptDiMu"] = {"xlabel":"p_{T,#mu#mu} [GeV/c]","xlimits":[0.0,20.0],"rebin":1}
  histNames["ptDiMu"] = {"xlabel":"p_{T,#mu#mu} [GeV/c]","xlimits":[0.0,200.0],"rebin":2}
  
  histNames["mDiJet"] = {"xlabel":"m_{jj} [GeV/c^{2}]","xlimits":[300.0,1400.0],"rebin":5}
  histNames["ptDiJet"] = {"xlabel":"p_{T,jj} [GeV/c]","xlimits":[0.0,400.0],"rebin":2}
  #histNames["deltaEtaJets"] = {"xlabel":"#Delta#eta_{jj}","xlimits":[0.0,8.5]}
  histNames["deltaEtaJets"] = {"xlabel":"#Delta#eta_{jj}","xlimits":[3.0,7.5],"ylimits":[0.1,1e3]}
  histNames["deltaRJets"] = {"xlabel":"#DeltaR_{jj}","xlimits":[0.0,8.5]}
  histNames["deltaPhiJets"] = {"xlabel":"#Delta#phi_{jj}","xlimits":[0.0,3.2],"rebin":2,"leg":ulLegendPos}
  histNames["deltaEtaMuons"] = {"xlabel":"#Delta#eta_{#mu#mu}","xlimits":[0.0,4.5]}
  histNames["deltaRMuons"] = {"xlabel":"#DeltaR_{#mu#mu}","xlimits":[0.0,6.0]}
  histNames["deltaPhiMuons"] = {"xlabel":"#Delta#phi_{#mu#mu}","xlimits":[0.0,3.2],"rebin":2,"leg":ulLegendPos}
  
  histNames["ptMu1"] = {"xlabel":"Leading Muon p_{T} [GeV/c]","xlimits":[20.0,200.0]}
  histNames["ptMu2"] = {"xlabel":"Sub-Leading Muon p_{T} [GeV/c]","xlimits":[20.0,200.0]}
  
  histNames["ptJet1"] = {"xlabel":"Leading Jet p_{T} [GeV/c]","xlimits":[30.0,250.0]}
  histNames["ptJet2"] = {"xlabel":"Sub-Leading Jet p_{T} [GeV/c]","xlimits":[0.0,250.0]}
  histNames["ht"] = {"xlabel":"Jet H_{T} [GeV], for Jets p_{T}>30 GeV","xlimits":[0.0,1000.0]}
  histNames["htInRapidityGap"] = {"xlabel":"Jet H_{T} in Jet Rapidity Gap [GeV], for Jets p_{T}>30 GeV","xlimits":[0.0,1000.0]}
  histNames["nJets"] = {"xlabel":"N_{jets}","xlimits":[0,10]}
  histNames["nJetsInRapidityGap"] = {"xlabel":"N_{jets} in Jet Rapidity Gap [GeV]","xlimits":[0,10]}
  
  histNames["etaMu1"] = {"xlabel":"Leading Muon #eta","xlimits":[-2.4,2.4],"rebin":2,"leg":ccLegendPos}
  histNames["etaMu2"] = {"xlabel":"Sub-Leading Muon #eta","xlimits":[-2.4,2.4],"rebin":2,"leg":ccLegendPos}
  
  histNames["etaJet1"] = {"xlabel":"Leading Jet #eta","xlimits":[-5.0,5.0],"rebin":2}
  histNames["etaJet2"] = {"xlabel":"Sub-Leading Jet #eta","xlimits":[-5.0,5.0],"rebin":2}
  
  histNames["yDiJet"] = {"xlabel":"y_{jj}","xlimits":[-5.0,5.0],"rebin":4,"ylimits":[0.1,1e3]}
  #histNames["yDiMu"] = {"xlabel":"y_{#mu#mu}","xlimits":[-2.2,2.2],"rebin":2,"leg":lcLegendPos}
  histNames["yDiMu"] = {"xlabel":"y_{#mu#mu}","xlimits":[-2.2,2.2],"rebin":2,"ylimits":[0.1,9e6]}
  
  #histNames["cosThetaStar"] = {"xlabel":"cos(#theta^{*})","xlimits":[-1.0,1.0],"rebin":2,"leg":lcLegendPos}
  histNames["cosThetaStar"] = {"xlabel":"cos(#theta^{*})","xlimits":[-1.0,1.0],"rebin":2,"ylimits":[0.1,1e7]}
  histNames["cosThetaStarCS"] = {"xlabel":"cos(#theta^{*}_{CS})","xlimits":[-1.0,1.0],"rebin":2,"leg":lcLegendPos}
  
  histNames["relIsoMu1"] = {"xlabel":"Leading Muon Relative PF Isolation","xlimits":[0,0.3],"rebin":2}
  histNames["relIsoMu2"] = {"xlabel":"Sub-Leading Muon Relative PF Isolation","xlimits":[0,0.3],"rebin":2}
  
  #histNames["BDTHistMuonOnly"] = {"xlabel":"BDT (Non-VBF Category)","xlimits":[-0.8,0.0],"rebin":2,"ylimits":[0.1,1e5],'vertLines':{"8TeV":-0.55,"7TeV":-0.48}}
  #histNames["BDTHistVBF"] = {"xlabel":"BDT (VBF Category)","xlimits":[-0.55,0.30],"rebin":2,"ylimits":[0.1,5e2],'vertLines':{"8TeV":-0.04,"7TeV":-0.03}}
  
  #histNames["BDTHistMuonOnly"] = {"xlabel":"BDT (Non-VBF Category)","xlimits":[-0.8,0.0],"rebin":2,'vertLines':{"8TeV":-0.55,"7TeV":-0.48}}
  #histNames["BDTHistVBF"] = {"xlabel":"BDT (VBF Category)","xlimits":[-0.55,0.30],"rebin":2,'vertLines':{"8TeV":-0.04,"7TeV":-0.03}}
  
  histNames["BDTHistMuonOnly"] = {"xlabel":"BDT (Non-VBF Category)","xlimits":[-0.8,0.0],"rebin":2,"ylimits":[1e-2,1e8],'vertLines':{"8TeV":-0.55,"7TeV":-0.42}}
  histNames["BDTHistVBF"] = {"xlabel":"BDT (VBF Category)","xlimits":[-0.55,0.30],"rebin":2,"ylimits":[1e-2,5e4],'vertLines':{"8TeV":-0.04,"7TeV":-0.03}}
  
  # 7TeV
  histNames["BDTHistMuonOnly"] = {"xlabel":"BDT (Non-VBF Category)","xlimits":[-0.55,0.2],"rebin":2,"ylimits":[1e-2,1e8],'vertLines':{"8TeV":-0.55,"7TeV":-0.42}}
  histNames["BDTHistVBF"] = {"xlabel":"BDT (VBF Category)","xlimits":[-0.3,0.4],"rebin":2,"ylimits":[1e-2,5e2],'vertLines':{"8TeV":-0.04,"7TeV":-0.03}}
  
  histNames["puJetIDSimpleDiscJet1"] = {"xlabel":"PU Jet ID Simple Discriminator--Leading Jet","xlimits":[-1,1],"rebin":2,"leg":ulLegendPos,"ylimits":[0.1,500.]}
  histNames["puJetIDSimpleDiscJet2"] = {"xlabel":"PU Jet ID Simple Discriminator--Sub-Leading Jet","xlimits":[-1,1],"rebin":2,"leg":ulLegendPos,"ylimits":[0.1,5e3]}
  histNames["puJetIDSimpleDiscJet3"] = {"xlabel":"PU Jet ID Simple Discriminator--3rd Leading Jet","xlimits":[-1,1],"rebin":2,"leg":ulLegendPos,"ylimits":[0.1,500.]}
  
  histNames["puJetIDSimpleJet1"] = {"xlabel":"PU Jet Simple Loose ID--Leading Jet","xlimits":[],"rebin":1,"leg":llLegendPos}
  histNames["puJetIDSimpleJet2"] = {"xlabel":"PU Jet Simple Loose ID--Sub-Leading Jet","xlimits":[],"rebin":1,"leg":llLegendPos}
  histNames["puJetIDSimpleJet3"] = {"xlabel":"PU Jet Simple Loose ID--3rd Leading Jet","xlimits":[],"rebin":1,"leg":llLegendPos}
  
  histNames["nVtx"] = {"xlabel":"N_{vtx}","xlimits":[0,40],"leg":stdLegendPos}
  histNames["met"] = {"xlabel":"E_{T}^{Miss}","xlimits":[0,300],"leg":stdLegendPos}
  histNames["ptmiss"] = {"xlabel":"p_{T}^{Miss} [GeV/c]","xlimits":[0,200],"leg":stdLegendPos}
  histNames["deltaPhiHJ1"] = {"xlabel":"#Delta#phi(j_{1},H)","xlimits":[0,3.2],'rebin':2,"leg":stdLegendPos,"ylimits":[5e-3,5e4]}
  """

for key in histNames:
  if re.match("pt.*",key):
    histNames[key]["units"] = " GeV/c"
  elif re.match("mDi.*",key):
    histNames[key]["units"] = " GeV/c^{2}"

tlatex = root.TLatex()
tlatex.SetNDC()
tlatex.SetTextSize(0.07)
tlatex.SetTextAlign(12)

#######################################
root.gROOT.SetBatch(True)
setStyle()
#######################################

scaleFactors = {}
#print "scale factors:"
for i in nEventsMap:
  if nEventsMap[i] ==0.0:
    scaleFactors[i] = 0.0
  else:
    p = getPeriod(i)
    scaleFactors[i] = xsec[i]*1000.0*LUMI/nEventsMap[i]*efficiencyMap[p]*mcPlotScaleFactorMap[p]
  #print "%s = %.2e" %(i,scaleFactors[i])

#######################################

class Dataset:
  def __init__(self,filename,legendEntry,color,scaleFactor,isData=False,isSignal=False):
    self.filename = filename
    self.legendEntry = legendEntry
    self.color = color
    self.scaleFactor = scaleFactor

    if filename != "":
      self.rootFile = root.TFile(filename)
    self.hists = {}
    self.datasetName = os.path.basename(filename)
    self.datasetName = self.datasetName.replace(".root","")
    self.isData=isData
    self.isSignal=isSignal

  def isZombie(self):
    return self.rootFile.IsZombie()

  def loadHistos(self,names,prefix=""):
    for name in names:
      #print("In datasetName: {0}, loading histogram: {1}".format(self.datasetName,name))
      tmpHistInfo = histNames[name]
      xlimits = tmpHistInfo["xlimits"]
      tmp = self.rootFile.Get(prefix+name)
      if type(tmp) != root.TH1F:
        print("Warning: In datasetName: {0}, loading histogram: {1}: Object type is not TH1F!!".format(self.datasetName,prefix+name))
        continue
      tmp.UseCurrentStyle()
      if histNames[name].has_key("rebin"):
        tmp.Rebin(histNames[name]["rebin"])
      if self.isSignal:
        tmp.SetLineColor(self.color)
      elif self.isData:
        tmp.SetMarkerColor(self.color)
        tmp.SetLineColor(self.color)
      else:
        tmp.SetFillColor(self.color)
      tmp.Scale(self.scaleFactor)
      self.hists[prefix+name] = tmp
  def eatOtherDatasets(self,listOfOthers):
    first = True
    for ds in listOfOthers:
      if first:
        first = False
        self.hists = ds.hists
        for key in self.hists:
          tmp = self.hists[key]
          if self.isSignal:
            tmp.SetLineColor(self.color)
          elif self.isData:
            tmp.SetMarkerColor(self.color)
            tmp.SetLineColor(self.color)
          else:
            tmp.SetFillColor(self.color)
      else:
        for key in self.hists:
          self.hists[key].Add(ds.hists[key])

#######################################

bkgDatasetList = []
for i in backgroundList:
  i += "_"+RUNPERIOD
  if i in scaleFactors:
    if scaleFactors[i]>0.0:
      filename = dataDir+i+".root"
      if not os.path.exists(filename):
          continue
      tmp = Dataset(filename,getLegendEntry(i),getColor(i),scaleFactors[i])
      if tmp.isZombie():
        print ("Warning: file for dataset {0} is Zombie!!".format(i))
        continue
      #print("Loading Dataset: {0}".format(i))
      for hDir in histDirs:
        tmp.loadHistos(histNames,prefix=hDir)
      bkgDatasetList.append(tmp)

sigDatasetList = []
for i in signalList:
  i += "_"+RUNPERIOD
  if i in scaleFactors:
    if scaleFactors[i]>0.0:
      filename = dataDir+i+".root"
      if not os.path.exists(filename):
          continue
      tmp = Dataset(filename,getLegendEntry(i),getColor(i),scaleFactors[i]*scaleHiggsBy,isSignal=True)
      if tmp.isZombie():
        print ("Warning: file for dataset {0} is Zombie!!".format(i))
        continue
      #print("Loading Dataset: {0}".format(i))
      for hDir in histDirs:
        tmp.loadHistos(histNames,prefix=hDir)
      sigDatasetList.append(tmp)
oldSigDatasetList = sigDatasetList
if allHiggsTogether:
  allHiggsDataset = Dataset("","H #rightarrow #mu#mu",root.kRed,1.0,isSignal=True)
  allHiggsDataset.eatOtherDatasets(sigDatasetList)
  sigDatasetList = [allHiggsDataset]

realDatasetList = []
for i in dataDict[RUNPERIOD]:
      filename = dataDir+i+".root"
      if not os.path.exists(filename):
        print("Error: Data file not found {}, exiting".format(filename))
        sys.exit(1)
      tmp = Dataset(filename,getLegendEntry(RUNPERIOD),1,1.0,isData=True)
      if tmp.isZombie():
        print ("Error: file for dataset {0} is Zombie!!".format(i))
        sys.exit(1)
      #print("Loading Dataset: {0}".format(i))
      for hDir in histDirs:
        tmp.loadHistos(histNames,prefix=hDir)
      realDatasetList.append(tmp)

#######################################

canvas = root.TCanvas("canvas")
leg = root.TLegend(*stdLegendPos)
leg.SetLineColor(0)
leg.SetFillColor(0)
uniqueLegendEntries = set()
for ds in bkgDatasetList:
  for hname in ds.hists:
    if ds.legendEntry not in uniqueLegendEntries:
      leg.AddEntry(ds.hists[hname],ds.legendEntry,"f")
      uniqueLegendEntries.add(ds.legendEntry)
    break
for ds in sigDatasetList:
  for hname in ds.hists:
    if ds.legendEntry not in uniqueLegendEntries:
      leg.AddEntry(ds.hists[hname],ds.legendEntry,"l")
      uniqueLegendEntries.add(ds.legendEntry)
    break
for ds in realDatasetList:
  for hname in ds.hists:
    if ds.legendEntry not in uniqueLegendEntries:
      leg.AddEntry(ds.hists[hname],ds.legendEntry,"pe")
      uniqueLegendEntries.add(ds.legendEntry)
    break
  break

#######################################
dataMCRatioStr = "############################################\n"

for histName in bkgDatasetList[0].hists:
  #print("Making Histo: %s" % histName)
  canvas.Clear()
  bkgHistList = []
  sigHistList = []
  for ds in bkgDatasetList:
    tmpHist = ds.hists[histName]
    bkgHistList.append(tmpHist)
  bkgHistList.reverse()
  for ds in sigDatasetList:
    tmpHist = ds.hists[histName]
    sigHistList.append(tmpHist)
  #bkgHistList.reverse()

  #dataHist = bkgDatasetList[0].hists[histName].Clone()
  #dataHist.Reset()
  dataHist = realDatasetList[0].hists[histName].Clone()
  dataHist.Reset()
  for realDS in realDatasetList:
    dataHist.Add(realDS.hists[histName])

  histBaseName = re.sub(r".*/","",histName)
  xtitle = histNames[histBaseName]["xlabel"]
  ylimits = []
  if histNames[histBaseName].has_key("ylimits"):
     ylimits = histNames[histBaseName]["ylimits"]
  stack = DataMCStack(bkgHistList, dataHist, canvas, xtitle,lumi=LUMI,logy=LOGY,xlimits=histNames[histBaseName]["xlimits"],signalsNoStack=sigHistList,integralPlot=integralPlot,energyStr=RUNPERIOD,ylimits=ylimits,ylimitsRatio=ylimitsRatio,pullType=PULLTYPE,doMCErrors=MCErrors)

  legLeftPos = stdLegendPos[0]
  legRightPos = stdLegendPos[2]
  scaleHiggsPos = "std"
  if histNames[histBaseName].has_key("leg"):
    setLegPos(leg,(histNames[histBaseName]["leg"]))
    legLeftPos = histNames[histBaseName]["leg"][0]
    legRightPos = histNames[histBaseName]["leg"][2]
    tmplegPos = histNames[histBaseName]["leg"]
    if tmplegPos == lcLegendPos:
      scaleHiggsPos = "lc"
    elif tmplegPos == llLegendPos:
      scaleHiggsPos = "ll"
    elif tmplegPos == ulLegendPos:
      scaleHiggsPos = "ul"
  else:
    setLegPos(leg,stdLegendPos)
  leg.Draw("same")
  #print("{},{},{},{}".format(leg.GetX1NDC(),leg.GetY1NDC(),leg.GetX2NDC(),leg.GetY2NDC()))

  #if histBaseName == "deltaPhiHJ1":
  if False:
    print "justin justin justin"
    print histNames[histBaseName]

  if scaleHiggsPos == "lc":
    if scaleHiggsBy != 1.0:
      tlatex.SetTextSize(0.07)
      tlatex.SetTextAlign(22)
      tlatex.DrawLatex(0.55,0.7,"Higgs #times {0:.0f}".format(scaleHiggsBy))

    tlatex.SetTextSize(0.03)
    tlatex.SetTextAlign(22)
    tlatex.DrawLatex(0.55,0.75,anotateText)
  elif scaleHiggsPos == "ll" or scaleHiggsPos == "ul":
    if scaleHiggsBy != 1.0:
      tlatex.SetTextSize(0.07)
      tlatex.SetTextAlign(22)
      tlatex.DrawLatex(0.55,0.8,"Higgs #times {0:.0f}".format(scaleHiggsBy))

    tlatex.SetTextSize(0.03)
    tlatex.SetTextAlign(22)
    tlatex.DrawLatex(0.55,1.0-gStyle.GetPadTopMargin()-0.02,anotateText)
  else:
    if scaleHiggsBy != 1.0:
      tlatex.SetTextSize(0.07)
      tlatex.SetTextAlign(32)
      tlatex.DrawLatex(legLeftPos-0.02,0.8,"Higgs #times {0:.0f}".format(scaleHiggsBy))

    tlatex.SetTextSize(0.03)
    tlatex.SetTextSize(0.04)
    tlatex.SetTextAlign(33)
    tlatex.DrawLatex(legLeftPos-0.02,1.0-gStyle.GetPadTopMargin()-0.02,anotateText)

  vertLine = None
  if histNames[histBaseName].has_key("vertLines"):
    #print("In vertLines for hist: {0}".format(histBaseName))
    padX1 = stack.pad1.GetX1()
    padX2 = stack.pad1.GetX2()
    padY1 = stack.pad1.GetY1()
    padY2 = stack.pad1.GetY2()
    stack.pad1.cd()
    vertLineX = histNames[histBaseName]["vertLines"][RUNPERIOD]
    vertLine = root.TLine()
    vertLine.SetLineColor(root.kRed+1)
    vertLine.SetLineWidth(3)
    pad1Width = stack.pad1.XtoPixel(stack.pad1.GetX2())-stack.pad1.XtoPixel(stack.pad1.GetX1())
    pad1Height = stack.pad1.YtoPixel(stack.pad1.GetY1())
    #vertLine.DrawLineNDC(0.1,0.1,stack.pad1.XtoPixel(0.0)/pad1Width,0.9)
    #vertLine.DrawLine(0.1,0.1,-0.4,10.)
    ybot = stack.stack.GetMinimum()
    ytop = stack.stack.GetMaximum()
    if stack.dataHist.GetMaximum()>ytop:
        ytop = stack.dataHist.GetMaximum()
    if LOGY:
      ytop = stack.stack.GetMaximum()*2
    #if histNames[histBaseName].has_key("ylimits"):
    #   ylimits = histNames[histBaseName]["ylimits"]
    #   ybot = ylimits[0]
    #   ytop = ylimits[1]
    vertLine.DrawLine(vertLineX,ybot,vertLineX,ytop)
    stack.pad1.RedrawAxis()

  if stack.mcSumHist.Integral() == 0:
    print("Warning: MC Sum Hist Was 0 for {0}, not saving image".format(histName))
    continue

  if histNames[histBaseName].has_key("units"):
    tmpAx = stack.stack.GetYaxis()
    tmpAx.SetTitle(tmpAx.GetTitle()+histNames[histBaseName]["units"])

  saveName = histName.replace("(","")
  saveName = saveName.replace(")","")
  saveName = saveName.replace("[","")
  saveName = saveName.replace("]","")
  saveName = saveName.replace("/","_")
  if integralPlot:
    saveName += "_IntPlot"
  saveAs(canvas,outDir+saveName+"_"+RUNPERIOD)

  match = re.match(r"(.*)_ptDiMu",saveName)
  if match:
    catName = match.group(1)
    dataMCRatioStr += "{0:<10} Data/MC Ratio: {1:.3f}\n".format(catName,float(stack.nDataEvents)/stack.nMCEvents)
  elif saveName == "ptDiMu":
    catName = "All"
    dataMCRatioStr += "{0:<10} Data/MC Ratio: {1:.3f}\n".format(catName,float(stack.nDataEvents)/stack.nMCEvents)

print(dataMCRatioStr)
