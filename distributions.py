#!/usr/bin/env python

from xsec import *
from helpers import *
import ROOT as root
import os

dataDir = "input/freezeSample110to150/"
caption1= "110 GeV < m_{#mu#mu} < 150 GeV"
dataDir = "input/preApproveSample/"
caption1= ""
outDir = "output/"

RUNPERIOD="8TeV"
LUMI=lumiDict[RUNPERIOD]

LOGY=False
reverse=False
allHiggsTogether = True
drawString="E"

histDirs = ["IncPreselPtG10/"]
histDirs = ["IncPresel/"]

backgroundList = [
"DYJetsToLL"#,
#"ttbar"
]
colors["DYJetsToLL"] = root.kBlue

urLegendPos = [0.68,0.65,0.88,0.88]
lcLegendPos = [0.45,0.25,0.65,0.48]
lrLegendPos = [0.68,0.25,0.88,0.48]
ulLegendPos = [0.25,0.65,0.45,0.88]
stdLegendPos = urLegendPos

histNames = {}
histNames["mDiMu"] = {"xlabel":"m_{#mu#mu} [GeV/c^{2}]","xlimits":[100.0,150.0],"rebin":2}

histNames["ptDiMu"] = {"xlabel":"p_{T,#mu#mu} [GeV/c]","xlimits":[0.0,50.0],"rebin":1}
histNames["yDiMu"] = {"xlabel":"y_{#mu#mu}","xlimits":[-3.0,3.0],"leg":lcLegendPos,"rebin":2}
histNames["cosThetaStar"] = {"xlabel":"cos(#theta^{*})","xlimits":[-1.0,1.0],"leg":lcLegendPos,"rebin":2}

"""
histNames["mDiJet"] = {"xlabel":"m_{jj} [GeV]","xlimits":[0.0,1200.0],"rebin":5}
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

histNames["cosThetaStarCS"] = {"xlabel":"cos(#theta^{*}_{CS})","xlimits":[-1.0,1.0],"rebin":2}

histNames["likelihoodHistMuonOnly"] = {"xlabel":"Likelihood (Not-VBF Category)","xlimits":[-1,0.5],"rebin":20}
histNames["BDTHistMuonOnly"] = {"xlabel":"BDT (Non-VBF Category)","xlimits":[-1,0.25],"rebin":20}

histNames["likelihoodHistVBF"] = {"xlabel":"Likelihood (VBF Category)","xlimits":[-0.5,1.0],"rebin":20}
histNames["BDTHistVBF"] = {"xlabel":"BDT (VBF Category)","xlimits":[-0.5,0.25],"rebin":20}

histNames["puJetIDSimpleDiscJet1"] = {"xlabel":"PU Jet ID Simple Discriminator--Leading Jet","xlimits":[-1,1],"rebin":1}
histNames["puJetIDSimpleDiscJet2"] = {"xlabel":"PU Jet ID Simple Discriminator--Sub-Leading Jet","xlimits":[-1,1],"rebin":1}
histNames["puJetIDSimpleDiscJet3"] = {"xlabel":"PU Jet ID Simple Discriminator--3rd Leading Jet","xlimits":[-1,1],"rebin":1}

histNames["puJetIDSimpleJet1"] = {"xlabel":"PU Jet Simple Loose ID--Leading Jet","xlimits":[],"rebin":1}
histNames["puJetIDSimpleJet2"] = {"xlabel":"PU Jet Simple Loose ID--Sub-Leading Jet","xlimits":[],"rebin":1}
histNames["puJetIDSimpleJet3"] = {"xlabel":"PU Jet Simple Loose ID--3rd Leading Jet","xlimits":[],"rebin":1}
"""

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

    if filename != "":
      self.rootFile = root.TFile(filename)
    self.hists = {}
    self.datasetName = os.path.basename(filename)
    self.datasetName = self.datasetName.replace(".root","")

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
      tmp.SetFillColor(self.color)
      tmp.SetFillStyle(1)
      tmp.SetLineColor(self.color)
      tmp.SetMarkerColor(self.color)
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
          tmp.SetFillColor(self.color)
          tmp.SetFillStyle(1)
          tmp.SetMarkerColor(self.color)
          tmp.SetLineColor(self.color)
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
      tmp = Dataset(filename,getLegendEntry(i),getColor(i),scaleFactors[i])
      if tmp.isZombie():
        print ("Warning: file for dataset {0} is Zombie!!".format(i))
        continue
      #print("Loading Dataset: {0}".format(i))
      for hDir in histDirs:
        tmp.loadHistos(histNames,prefix=hDir)
      sigDatasetList.append(tmp)
oldSigDatasetList = sigDatasetList
if allHiggsTogether:
  allHiggsDataset = Dataset("","H #rightarrow #mu#mu",root.kRed,1.0)
  allHiggsDataset.eatOtherDatasets(sigDatasetList)
  sigDatasetList = [allHiggsDataset]

bkgDatasetList += sigDatasetList

if reverse:
  bkgDatasetList.reverse()

#######################################

canvas = root.TCanvas("canvas")
#leg = root.TLegend(*stdLegendPos)
leg = root.TLegend(*lcLegendPos)
leg.SetLineColor(0)
leg.SetFillColor(0)
uniqueLegendEntries = set()
for ds in bkgDatasetList:
  for hname in ds.hists:
    if ds.legendEntry not in uniqueLegendEntries:
      leg.AddEntry(ds.hists[hname],ds.legendEntry,"lpe")
      uniqueLegendEntries.add(ds.legendEntry)
    break

#######################################

for histName in bkgDatasetList[0].hists:
  canvas.Clear()
  histBaseName = re.sub(r".*/","",histName)
  print("Making Histo: %s" % histName)
  bkgHistList = []
  for ds in bkgDatasetList:
    tmpHist = ds.hists[histName]
    if histNames[histBaseName].has_key("xlimits"):
      if histNames[histBaseName]["xlimits"] != []:
        tmpHist.GetXaxis().SetRangeUser(*histNames[histBaseName]["xlimits"])
    if histNames[histBaseName].has_key("rebin"):
        tmpHist.Rebin(histNames[histBaseName]["rebin"])
    if tmpHist.Integral() != 0.0:
      tmpHist.Scale(1.0/tmpHist.Integral())
    bkgHistList.append(tmpHist)
  bkgHistList.reverse()

  xtitle = histNames[histBaseName]["xlabel"]
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
    if histNames[histBaseName].has_key("xlimits"):
      if histNames[histBaseName]["xlimits"] != []:
        hist.GetXaxis().SetRangeUser(*histNames[histBaseName]["xlimits"])
    hist.SetFillStyle(0)
    hist.SetLineStyle(1)
    hist.SetLineWidth(2)
    hist.SetLineColor(hist.GetFillColor())
    if ymax < hist.GetMaximum():
      ymax = hist.GetMaximum()
    if ymin < hist.GetMinimum():
      ymin = hist.GetMinimum()
    if firstHist:
      hist.Draw(drawString)
      firstHist = False
    else:
      hist.Draw(drawString+" same")
  firstHist = True
  print ymax
  for hist in bkgHistList:
    hist.SetTitle("")
    hist.GetYaxis().SetRangeUser(ymin*0.95,ymax*0.6)
    #hist.GetYaxis().SetRangeUser(ymin*0.95,ymax*0.8)
    if firstHist:
      hist.Draw(drawString)
      firstHist = False
    else:
      hist.Draw(drawString+" same")
  canvas.RedrawAxis()
  if histNames[histBaseName].has_key("leg"):
    setLegPos(leg,(histNames[histBaseName]["leg"]))
  else:
    setLegPos(leg,(stdLegendPos))
  leg.Draw()

  tlatex = root.TLatex()
  tlatex.SetNDC()
  tlatex.SetTextFont(root.gStyle.GetLabelFont())
  tlatex.SetTextSize(0.04)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,"#sqrt{{s}}={0}".format(RUNPERIOD))
  tlatex.SetTextAlign(13)
  tlatex.DrawLatex(gStyle.GetPadLeftMargin()+0.03,1.0-gStyle.GetPadTopMargin()-0.04,caption1)

  saveName = histName.replace("(","")
  saveName = saveName.replace(")","")
  saveName = saveName.replace("[","")
  saveName = saveName.replace("]","")
  saveName = saveName.replace("/","_")
  saveAs(canvas,outDir+"dist_"+saveName+"_"+RUNPERIOD)
