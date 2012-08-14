#!/usr/bin/python

from xsec import *
from helpers import *
import ROOT as root
import os

dataDir = "input/4GeVWindow/"
#dataDir = "input/open/"
outDir = "output/"

LOGY=False
reverse=False

histNames = {}
histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV]","xlimits":[100.0,150.0]}
histNames["mDiMuVBFSelected"] = {"xlabel":"m_{#mu#mu}, After VBF Selection [GeV]","xlimits":[100.0,150.0]}
histNames["mDiMuVBFTightSelected"] = {"xlabel":"m_{#mu#mu}, After VBF-Tight Selection [GeV]","xlimits":[100.0,150.0]}
histNames["mDiMuZPt30Selected"] = {"xlabel":"m_{#mu#mu}, After p_{T}^{#mu#mu}>30 GeV Selection [GeV]","xlimits":[100.0,150.0]}
histNames["mDiMuZPt50Selected"] = {"xlabel":"m_{#mu#mu}, After p_{T}^{#mu#mu}>50 GeV Selection [GeV]","xlimits":[100.0,150.0]}
histNames["mDiMuZPt75Selected"] = {"xlabel":"m_{#mu#mu}, After p_{T}^{#mu#mu}>75 GeV Selection [GeV]","xlimits":[100.0,150.0]}

histNames["ptDiMu"] = {"xlabel":"p_{T,#mu#mu} [GeV]","xlimits":[0.0,200.0]}
histNames["ptDiMuVBFSelected"] = {"xlabel":"p_{T,#mu#mu}, After VBF Selection [GeV]","xlimits":[0.0,200.0],"rebin":10}
histNames["ptDiMuVBFLooseSelected"] = {"xlabel":"p_{T,#mu#mu}, After VBF-Loose Selection [GeV]","xlimits":[0.0,200.0],"rebin":5}

histNames["mDiJet"] = {"xlabel":"m_{jj} [GeV]","xlimits":[0.0,1200.0],"rebin":5}
histNames["deltaEtaJets"] = {"xlabel":"#Delta#eta_{jj}","xlimits":[0.0,10.0]}

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

histNames["cosThetaStar"] = {"xlabel":"cos(#theta^{*})","xlimits":[-1.0,1.0],"rebin":2}
histNames["cosThetaStarVBFSelected"] = {"xlabel":"cos(#theta^{*}) After VBF Selection","xlimits":[-1.0,1.0],"rebin":5}
histNames["cosThetaStarVBFLooseSelected"] = {"xlabel":"cos(#theta^{*}) After VBF-Loose Selection","xlimits":[-1.0,1.0],"rebin":5}
histNames["cosThetaStarVBFTightSelected"] = {"xlabel":"cos(#theta^{*}) After VBF-Tight Selection","xlimits":[-1.0,1.0],"rebin":5}
histNames["cosThetaStarZPt30Selected"] = {"xlabel":"cos(#theta^{*}) After p_{T}(#mu#mu)>30 GeV Selection","xlimits":[-1.0,1.0],"rebin":2}
histNames["cosThetaStarZPt50Selected"] = {"xlabel":"cos(#theta^{*}) After p_{T}(#mu#mu)>50 GeV Selection","xlimits":[-1.0,1.0],"rebin":2}
histNames["cosThetaStarZPt75Selected"] = {"xlabel":"cos(#theta^{*}) After p_{T}(#mu#mu)>75 GeV Selection","xlimits":[-1.0,1.0],"rebin":2}

histNames["mDiMuEta11"] = {"xlabel":"m_{#mu#mu} Both Muons |#eta|<1.0 [GeV]","xlimits":[100.0,150.0]}
histNames["mDiMuEta12"] = {"xlabel":"m_{#mu#mu} One Muons |#eta|<1.0, One Muon |#eta|>1.0 [GeV]","xlimits":[100.0,150.0]}
histNames["mDiMuEta22"] = {"xlabel":"m_{#mu#mu} Both Muons |#eta|>1.0[GeV]","xlimits":[100.0,150.0]}

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

if reverse:
  bkgDatasetList.reverse()

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
      leg.AddEntry(ds.hists[hname],ds.legendEntry,"l")
      uniqueLegendEntries.add(ds.legendEntry)
    break

#######################################

for histName in bkgDatasetList[0].hists:
  print("Making Histo: %s" % histName)
  bkgHistList = []
  for ds in bkgDatasetList:
    tmpHist = ds.hists[histName]
    tmpHist.GetXaxis().SetRangeUser(*histNames[histName]["xlimits"])
    if tmpHist.Integral() != 0.0:
      tmpHist.Scale(1.0/tmpHist.Integral())
    if histNames[histName].has_key("rebin"):
        tmpHist.Rebin(histNames[histName]["rebin"])
    bkgHistList.append(tmpHist)
  bkgHistList.reverse()

  xtitle = histNames[histName]["xlabel"]
  if LOGY:
    canvas.SetLogy(1)
  else:
    canvas.SetLogy(0)
  firstHist = True
  ymin = 9999999.0
  ymax = -9999999.0
  for hist in bkgHistList:
    hist.SetTitle("")
    hist.GetYaxis().SetTitle("Normalized Events")
    hist.GetXaxis().SetTitle(xtitle)
    hist.GetXaxis().SetRangeUser(*histNames[histName]["xlimits"])
    hist.SetFillStyle(0)
    hist.SetLineStyle(1)
    hist.SetLineWidth(2)
    hist.SetLineColor(hist.GetFillColor())
    if ymax < hist.GetMaximum():
      ymax = hist.GetMaximum()
    if ymin < hist.GetMinimum():
      ymin = hist.GetMinimum()
    if firstHist:
      hist.Draw()
      firstHist = False
    else:
      hist.Draw("same")
  firstHist = True
  print ymax
  for hist in bkgHistList:
    hist.SetTitle("")
    hist.GetYaxis().SetRangeUser(ymin*0.95,ymax*0.6)
    #hist.GetYaxis().SetRangeUser(ymin*0.95,ymax*0.8)
    if firstHist:
      hist.Draw()
      firstHist = False
    else:
      hist.Draw("same")
  leg.Draw("same")

  saveName = histName.replace("(","")
  saveName = saveName.replace(")","")
  saveName = saveName.replace("[","")
  saveName = saveName.replace("]","")
  saveAs(canvas,outDir+"dist_"+saveName)
