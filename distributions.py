#!/usr/bin/env python

from xsec import *
from helpers import *
import ROOT as root
import os

dataDir = "input/freezeSample110to150/"
caption1= "110 GeV < m_{#mu#mu} < 150 GeV"
dataDir = "input/jets20/"
caption1= ""
outDir = "output/"

RUNPERIOD="8TeV"
LUMI=lumiDict[RUNPERIOD]

LOGY=False
reverse=False
allHiggsTogether = False
drawString="E"

#histDirs = ["IncPreselPtG10/"]
histDirs = ["VBFPresel/"]

backgroundList = [
#"DYJetsToLL"#,
#"ttbar"
]

signalList = [
"ggHmumu125",
"vbfHmumu125",
"ggHmumu125ChangeEvents",
"vbfHmumu125ChangeEvents"
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


if True:
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

for i in histNames:
  histNames[i]["rebin"] = int(0.5*histNames[i]["rebin"])

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
  print i
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
leg = root.TLegend(*stdLegendPos)
#leg = root.TLegend(*lcLegendPos)
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
