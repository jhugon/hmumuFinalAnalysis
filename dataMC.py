#!/usr/bin/python

from xsec import *
from helpers import *
import ROOT as root
import os

dataDir = "input/"
outDir = "output/"

LOGY=False

histNames = {}
histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV]","xlimits":[100.0,150.0],"rebin":4}
histNames["mDiMuVBFM"] = {"xlabel":"m_{#mu#mu}, After VBF Medium Selection [GeV]","xlimits":[100.0,150.0],"rebin":4}
histNames["mDiMuVBFT"] = {"xlabel":"m_{#mu#mu}, After VBF Tight Selection [GeV]","xlimits":[100.0,150.0],"rebin":4}
histNames["mDiMuVBFL"] = {"xlabel":"m_{#mu#mu}, After VBF Loose Selection [GeV]","xlimits":[100.0,150.0],"rebin":4}
histNames["mDiMuVBFVL"] = {"xlabel":"m_{#mu#mu}, After VBF Very Loose Selection [GeV]","xlimits":[100.0,150.0],"rebin":4}
histNames["mDiMuPtL30"] = {"xlabel":"m_{#mu#mu}, After p_{T}^{#mu#mu}<30 GeV Selection [GeV]","xlimits":[100.0,150.0],"rebin":4}
histNames["mDiMuPt30to50"] = {"xlabel":"m_{#mu#mu}, After 30<p_{T}^{#mu#mu}<50 GeV Selection [GeV]","xlimits":[100.0,150.0],"rebin":4}
histNames["mDiMuPt50to75"] = {"xlabel":"m_{#mu#mu}, After 50<p_{T}^{#mu#mu}<75 GeV Selection [GeV]","xlimits":[100.0,150.0],"rebin":4}
histNames["mDiMuPt75to125"] = {"xlabel":"m_{#mu#mu}, After 75<p_{T}^{#mu#mu}<125 GeV Selection [GeV]","xlimits":[100.0,150.0],"rebin":4}
histNames["mDiMuPt125"] = {"xlabel":"m_{#mu#mu}, After p_{T}^{#mu#mu}>125 GeV Selection [GeV]","xlimits":[100.0,150.0],"rebin":4}

histNames["ptDiMu"] = {"xlabel":"p_{T,#mu#mu} [GeV]","xlimits":[0.0,400.0]}
histNames["ptDiMuVBFM"] = {"xlabel":"p_{T,#mu#mu}, After VBF Selection [GeV]","xlimits":[0.0,400.0]}
histNames["ptDiMuVBFL"] = {"xlabel":"p_{T,#mu#mu}, After VBF-Loose Selection [GeV]","xlimits":[0.0,400.0]}

histNames["mDiJet"] = {"xlabel":"m_{jj} [GeV]","xlimits":[0.0,1200.0]}
histNames["deltaEtaJets"] = {"xlabel":"#Delta#eta_{jj}","xlimits":[0.0,10.0]}
histNames["deltaRJets"] = {"xlabel":"#DeltaR_{jj}","xlimits":[0.0,10.0]}
histNames["deltaPhiJets"] = {"xlabel":"#Delta#phi_{jj}","xlimits":[0.0,3.2]}
histNames["deltaEtaMuons"] = {"xlabel":"#Delta#eta_{#mumu}","xlimits":[0.0,5.0]}
histNames["deltaRMuons"] = {"xlabel":"#DeltaR_{#mu#mu}","xlimits":[0.0,6.0]}
histNames["deltaPhiMuons"] = {"xlabel":"#Delta#phi_{#mu#mu}","xlimits":[0.0,3.2]}

histNames["ptMu1"] = {"xlabel":"Leading Muon p_{T} [GeV]","xlimits":[0.0,400.0]}
histNames["ptMu2"] = {"xlabel":"Sub-Leading Muon p_{T} [GeV]","xlimits":[0.0,400.0]}

histNames["ptJet1"] = {"xlabel":"Leading Jet p_{T} [GeV]","xlimits":[0.0,400.0]}
histNames["ptJet2"] = {"xlabel":"Sub-Leading Jet p_{T} [GeV]","xlimits":[0.0,400.0]}

histNames["etaMu1"] = {"xlabel":"Leading Muon #eta","xlimits":[-2.4,2.4]}
histNames["etaMu2"] = {"xlabel":"Sub-Leading Muon #eta","xlimits":[-2.4,2.4]}

histNames["etaJet1"] = {"xlabel":"Leading Jet #eta","xlimits":[-5.0,5.0]}
histNames["etaJet2"] = {"xlabel":"Sub-Leading Jet #eta","xlimits":[-5.0,5.0]}

histNames["yDiMu"] = {"xlabel":"y_{#mu#mu} [GeV]","xlimits":[-3.0,3.0]}
histNames["yDiMuVBFM"] = {"xlabel":"y_{#mu#mu}, After VBF Selection [GeV]","xlimits":[-3.0,3.0]}
histNames["yDiMuVBFL"] = {"xlabel":"y_{#mu#mu}, After VBF-Loose Selection [GeV]","xlimits":[-3.0,3.0]}
histNames["yDiMuPt30"] = {"xlabel":"y_{#mu#mu}, After p_{T}^{#mu#mu}>30 GeV Selection [GeV]","xlimits":[-3.0,3.0]}
histNames["yDiMuPt50"] = {"xlabel":"y_{#mu#mu}, After p_{T}^{#mu#mu}>50 GeV Selection [GeV]","xlimits":[-3.0,3.0]}
histNames["yDiMuPt75"] = {"xlabel":"y_{#mu#mu}, After p_{T}^{#mu#mu}>75 GeV Selection [GeV]","xlimits":[-3.0,3.0]}

histNames["cosThetaStar"] = {"xlabel":"cos(#theta^{*})","xlimits":[-1.0,1.0]}
histNames["cosThetaStarVBFM"] = {"xlabel":"cos(#theta^{*}) After VBF Selection","xlimits":[-1.0,1.0]}
histNames["cosThetaStarVBFL"] = {"xlabel":"cos(#theta^{*}) After VBF-Loose Selection","xlimits":[-1.0,1.0]}
histNames["cosThetaStarVBFT"] = {"xlabel":"cos(#theta^{*}) After VBF-Tight Selection","xlimits":[-1.0,1.0]}
histNames["cosThetaStarPt30"] = {"xlabel":"cos(#theta^{*}) After p_{T}(#mu#mu)>30 GeV Selection","xlimits":[-1.0,1.0]}
histNames["cosThetaStarPt50"] = {"xlabel":"cos(#theta^{*}) After p_{T}(#mu#mu)>50 GeV Selection","xlimits":[-1.0,1.0]}
histNames["cosThetaStarPt75"] = {"xlabel":"cos(#theta^{*}) After p_{T}(#mu#mu)>75 GeV Selection","xlimits":[-1.0,1.0]}

histNames["likelihoodHistMuonOnly"] = {"xlabel":"Likelihood (Not-VBF Category)","xlimits":[-1,0.5],"rebin":20}
histNames["LDHistMuonOnly"] = {"xlabel":"LD (Not-VBF Category)","xlimits":[-0.2,1],"rebin":20}
histNames["BDTHistMuonOnly"] = {"xlabel":"BDT (Not-VBF Category)","xlimits":[-1,0.25],"rebin":20}
#histNames["BDTHistMuonOnly"] = {"xlabel":"BDT (Not-VBF Category)","xlimits":[-0.2,0.15],"rebin":10}

histNames["likelihoodHistVBF"] = {"xlabel":"Likelihood (VBF Category)","xlimits":[-0.5,1.0],"rebin":20}
histNames["LDHistVBF"] = {"xlabel":"LD (VBF Category)","xlimits":[0.0,0.75],"rebin":20}
histNames["BDTHistVBF"] = {"xlabel":"BDT (VBF Category)","xlimits":[-0.5,0.25],"rebin":20}
#histNames["BDTHistVBF"] = {"xlabel":"BDT (VBF Category)","xlimits":[-0.2,0.15],"rebin":10}

histNames["puJetIDSimpleDiscJet1"] = {"xlabel":"PU Jet ID Simple Discriminator--Leading Jet","xlimits":[-1,1],"rebin":1}
histNames["puJetIDSimpleDiscJet2"] = {"xlabel":"PU Jet ID Simple Discriminator--Sub-Leading Jet","xlimits":[-1,1],"rebin":1}
histNames["puJetIDSimpleDiscJet3"] = {"xlabel":"PU Jet ID Simple Discriminator--3rd Leading Jet","xlimits":[-1,1],"rebin":1}

histNames["puJetIDSimpleJet1"] = {"xlabel":"PU Jet Simple Loose ID--Leading Jet","xlimits":[],"rebin":1}
histNames["puJetIDSimpleJet2"] = {"xlabel":"PU Jet Simple Loose ID--Sub-Leading Jet","xlimits":[],"rebin":1}
histNames["puJetIDSimpleJet3"] = {"xlabel":"PU Jet Simple Loose ID--3rd Leading Jet","xlimits":[],"rebin":1}

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
    scaleFactors[i] = xsec[i]*1000.0*LUMI/nEventsMap[i]
  #print "%s = %.2e" %(i,scaleFactors[i])

#######################################

class Dataset:
  def __init__(self,filename,legendEntry,color,scaleFactor,isData=False,isSignal=False):
    self.filename = filename
    self.legendEntry = legendEntry
    self.color = color
    self.scaleFactor = scaleFactor

    self.rootFile = root.TFile(filename)
    self.hists = {}
    self.datasetName = os.path.basename(filename)
    self.datasetName = self.datasetName.replace(".root","")
    self.isData=isData
    self.isSignal=isSignal

  def isZombie(self):
    return self.rootFile.IsZombie()

  def loadHistos(self,names):
    for name in names:
      #print("In datasetName: {0}, loading histogram: {1}".format(self.datasetName,name))
      tmpHistInfo = histNames[name]
      xlimits = tmpHistInfo["xlimits"]
      tmp = self.rootFile.Get(name)
      if type(tmp) != root.TH1F:
        print("Warning: In datasetName: {0}, loading histogram: {1}: Object type is not TH1F!!".format(self.datasetName,name))
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
      self.hists[name] = tmp

#######################################

bkgDatasetList = []
for i in backgroundList:
  if i in scaleFactors:
    if scaleFactors[i]>0.0:
      filename = dataDir+i+".root"
      if not os.path.exists(filename):
          continue
      tmp = Dataset(filename,legendEntries[i],colors[i],scaleFactors[i])
      if tmp.isZombie():
        print ("Warning: file for dataset {0} is Zombie!!".format(i))
        continue
      print("Loading Dataset: {0}".format(i))
      tmp.loadHistos(histNames)
      bkgDatasetList.append(tmp)

sigDatasetList = []
for i in signalList:
  if i in scaleFactors:
    if scaleFactors[i]>0.0:
      filename = dataDir+i+".root"
      if not os.path.exists(filename):
          continue
      tmp = Dataset(filename,legendEntries[i],colors[i],scaleFactors[i],isSignal=True)
      if tmp.isZombie():
        print ("Warning: file for dataset {0} is Zombie!!".format(i))
        continue
      print("Loading Dataset: {0}".format(i))
      tmp.loadHistos(histNames)
      sigDatasetList.append(tmp)

#######################################

urLegendPos = [0.70,0.65,0.88,0.88]
ulLegendPos = [0.25,0.65,0.38,0.88]
stdLegendPos = urLegendPos
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
      leg.AddEntry(ds.hists[hname],ds.legendEntry,"f")
      uniqueLegendEntries.add(ds.legendEntry)
    break

#######################################

for histName in bkgDatasetList[0].hists:
  print("Making Histo: %s" % histName)
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

  dataHist = bkgDatasetList[0].hists[histName].Clone()  ##temporary
  dataHist.Reset()

  xtitle = histNames[histName]["xlabel"]
  stack = DataMCStack(bkgHistList, dataHist, canvas, xtitle,lumi=LUMI,logy=LOGY,xlimits=histNames[histName]["xlimits"],signalsNoStack=sigHistList)
  leg.Draw("same")

  if scaleHiggsBy != 1.0:
    tlatex.DrawLatex(0.33,0.8,"Higgs #times {0:.0f}".format(scaleHiggsBy))

  saveName = histName.replace("(","")
  saveName = saveName.replace(")","")
  saveName = saveName.replace("[","")
  saveName = saveName.replace("]","")
  saveAs(canvas,outDir+saveName)
