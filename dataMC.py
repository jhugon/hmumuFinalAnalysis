#!/usr/bin/python

from xsec import *
from helpers import *
import ROOT as root
import os

dataDir = "data/"
outDir = "output/"

LOGY=False

histNames = {}
histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV]","xlimits":[100.0,150.0]}
histNames["mDiMuVBFSelected"] = {"xlabel":"m_{#mu#mu}, After VBF Selection [GeV]","xlimits":[100.0,150.0]}
histNames["mDiMuVBFLooseSelected"] = {"xlabel":"m_{#mu#mu}, After VBF-Loose Selection [GeV]","xlimits":[100.0,150.0]}
histNames["mDiMuVBFTightSelected"] = {"xlabel":"m_{#mu#mu}, After VBF-Tight Selection [GeV]","xlimits":[100.0,150.0]}
histNames["mDiMuZPt30Selected"] = {"xlabel":"m_{#mu#mu}, After p_{T}^{#mu#mu}>30 GeV Selection [GeV]","xlimits":[100.0,150.0]}
histNames["mDiMuZPt50Selected"] = {"xlabel":"m_{#mu#mu}, After p_{T}^{#mu#mu}>50 GeV Selection [GeV]","xlimits":[100.0,150.0]}
histNames["mDiMuZPt75Selected"] = {"xlabel":"m_{#mu#mu}, After p_{T}^{#mu#mu}>75 GeV Selection [GeV]","xlimits":[100.0,150.0]}

histNames["ptDiMu"] = {"xlabel":"p_{T,#mu#mu} [GeV]","xlimits":[0.0,400.0]}
histNames["ptDiMuVBFSelected"] = {"xlabel":"p_{T,#mu#mu}, After VBF Selection [GeV]","xlimits":[0.0,400.0]}
histNames["ptDiMuVBFLooseSelected"] = {"xlabel":"p_{T,#mu#mu}, After VBF-Loose Selection [GeV]","xlimits":[0.0,400.0]}

histNames["mDiJet"] = {"xlabel":"m_{jj} [GeV]","xlimits":[0.0,1200.0]}
histNames["deltaEtaJets"] = {"xlabel":"#Delta#eta_{jj} [GeV]","xlimits":[0.0,10.0]}

histNames["ptMu1"] = {"xlabel":"Leading Muon p_{T} [GeV]","xlimits":[0.0,400.0]}
histNames["ptMu2"] = {"xlabel":"Sub-Leading Muon p_{T} [GeV]","xlimits":[0.0,400.0]}

histNames["ptJet1"] = {"xlabel":"Leading Jet p_{T} [GeV]","xlimits":[0.0,400.0]}
histNames["ptJet2"] = {"xlabel":"Sub-Leading Jet p_{T} [GeV]","xlimits":[0.0,400.0]}

histNames["etaMu1"] = {"xlabel":"Leading Muon #eta","xlimits":[-2.4,2.4]}
histNames["etaMu2"] = {"xlabel":"Sub-Leading Muon #eta","xlimits":[-2.4,2.4]}

histNames["etaJet1"] = {"xlabel":"Leading Jet #eta","xlimits":[-5.0,5.0]}
histNames["etaJet2"] = {"xlabel":"Sub-Leading Jet #eta","xlimits":[-5.0,5.0]}

histNames["yDiMu"] = {"xlabel":"y_{#mu#mu} [GeV]","xlimits":[-3.0,3.0]}
histNames["yDiMuVBFSelected"] = {"xlabel":"y_{#mu#mu}, After VBF Selection [GeV]","xlimits":[-3.0,3.0]}
histNames["yDiMuVBFLooseSelected"] = {"xlabel":"y_{#mu#mu}, After VBF-Loose Selection [GeV]","xlimits":[-3.0,3.0]}
histNames["yDiMuZPt30Selected"] = {"xlabel":"y_{#mu#mu}, After p_{T}^{#mu#mu}>30 GeV Selection [GeV]","xlimits":[-3.0,3.0]}
histNames["yDiMuZPt50Selected"] = {"xlabel":"y_{#mu#mu}, After p_{T}^{#mu#mu}>50 GeV Selection [GeV]","xlimits":[-3.0,3.0]}
histNames["yDiMuZPt75Selected"] = {"xlabel":"y_{#mu#mu}, After p_{T}^{#mu#mu}>75 GeV Selection [GeV]","xlimits":[-3.0,3.0]}

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
  def __init__(self,filename,legendEntry,color,scaleFactor):
    self.filename = filename
    self.legendEntry = legendEntry
    self.color = color
    self.scaleFactor = scaleFactor

    self.rootFile = root.TFile(filename)
    self.hists = {}
    self.datasetName = os.path.basename(filename)
    self.datasetName = self.datasetName.replace(".root","")

  def isZombie(self):
    return self.rootFile.IsZombie()

  def loadHistos(self,names):
    for name in names:
      print("In datasetName: {0}, loading histogram: {1}".format(self.datasetName,name))
      tmpHistInfo = histNames[name]
      xlimits = tmpHistInfo["xlimits"]
      tmp = self.rootFile.Get(name)
      print tmp
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

#######################################

for histName in bkgDatasetList[0].hists:
  print("Making Histo: %s" % histName)
  bkgHistList = []
  for ds in bkgDatasetList:
    tmpHist = ds.hists[histName]
    bkgHistList.append(tmpHist)
  bkgHistList.reverse()

  dataHist = bkgDatasetList[0].hists[histName].Clone()  ##temporary
  dataHist.Reset()

  xtitle = histNames[histName]["xlabel"]
  stack = DataMCStack(bkgHistList, dataHist, canvas, xtitle,lumi=LUMI,logy=LOGY,xlimits=histNames[histName]["xlimits"])
  leg.Draw("same")

  if scaleHiggsBy != 1.0:
    tlatex.DrawLatex(0.33,0.8,"Higgs #times {0:.0f}".format(scaleHiggsBy))

  saveName = histName.replace("(","")
  saveName = saveName.replace(")","")
  saveName = saveName.replace("[","")
  saveName = saveName.replace("]","")
  saveAs(canvas,outDir+saveName)
