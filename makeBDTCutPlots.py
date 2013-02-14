#!/usr/bin/env python

from xsec import *
from helpers import *
import ROOT as root
import os
import sys
import random

dataDir = "input/preApproveSample//"
outDir = "output/"

RUNPERIOD="7TeV"
LUMI=lumiDict[RUNPERIOD]

LOGY=True
integralPlot=False
ylimitsRatio = [0.5,1.5]

mRange = [120.0,130.0]

urLegendPos = [0.70,0.67,0.9,0.9]
ulLegendPos = [0.20,0.67,0.4,0.9]
ucLegendPos = [0.46,0.67,0.64,0.9]
lcLegendPos = [0.46,0.35,0.64,0.63]
llLegendPos = [0.20,0.35,0.4,0.63]
ccLegendPos = [0.46,0.47,0.64,0.7]
stdLegendPos = urLegendPos
drawStr = "hist"
drawStr = "E1"

histDirs = ["VBFPresel/"]

root.gErrorIgnoreLevel = root.kWarning

histNames = {}
if RUNPERIOD == "8TeV":
  #histNames["BDTHistMuonOnlyVMass"] = {"xlabel":"BDT Cut (Non-VBF Category)","xlimits":[-1.,0.1],
  #                                      'ylimits':[1e-6,10.0],
  #                                      'ylimitsSqrt':[1e-3,10],
  #                              'vertLines':{"8TeV":-0.55,"7TeV":-0.42 },
  #                                      "rebin":1}
  histNames["BDTHistVBFVMass"] = {"xlabel":"BDT Cut (VBF Category)","xlimits":[-0.4,0.25],
                                        'ylimits':[1e-3,1.0],
                                        'ylimitsSqrt':[1e-3,1.0],
                                        'ylimitsSqr':[1e-5,1.0],
                                'vertLines':{"8TeV":-0.04,"7TeV":-0.03},
                                        "rebin":4}
if RUNPERIOD == "7TeV":
#  histNames["BDTHistMuonOnlyVMass"] = {"xlabel":"BDT Cut (Non-VBF Category)","xlimits":[-0.55,0.2],
#                                        'ylimits':[1e-6,10.0],
#                                        'ylimitsSqrt':[1e-3,10],
#                                'vertLines':{"8TeV":-0.55,"7TeV":-0.42 },
#                                        "rebin":1}
  histNames["BDTHistVBFVMass"] = {"xlabel":"BDT Cut (VBF Category)","xlimits":[-0.4,0.25],
                                        'ylimits':[1e-6,100.0],
                                        'ylimitsSqrt':[1e-2,1.0],
                                        'ylimitsSqr':[1e-5,1.0],
                                'vertLines':{"8TeV":-0.04,"7TeV":-0.03},
                                        "rebin":4}

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

def drawLatex():
    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(root.gStyle.GetLabelFont())
    #tlatex.SetTextSize(0.05)
    tlatex.SetTextSize(0.04)
    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(0.98-gStyle.GetPadRightMargin(),0.88,"{0:.0f} GeV < m_{{#mu#mu}} < {1:.0f} GeV".format(*mRange))
    tlatex.DrawLatex(0.98-gStyle.GetPadRightMargin(),0.83,"L={0:.1f} fb^{{-1}}".format(LUMI))
    tlatex.DrawLatex(0.98-gStyle.GetPadRightMargin(),0.78,"#sqrt{s} = "+RUNPERIOD)

    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,"")

class Dataset:
  def __init__(self,filename,legendEntry,color,scaleFactor,isData=False,isSignal=False):
    self.filename = filename
    self.legendEntry = legendEntry
    self.color = color
    self.color = root.kBlue
    self.scaleFactor = scaleFactor

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
      tmp.SetName(tmp.GetName()+str(random.randint(0,10000)))
      if type(tmp) != root.TH2F:
        print("Warning: In datasetName: {0}, loading histogram: {1}: Object type is not TH2F!!".format(self.datasetName,prefix+name))
        continue
      tmp = hist2to1CollapseY(tmp,mRange)
      if type(tmp) != root.TH1F:
        print("Warning: In datasetName: {0}, loading histogram: {1}: Object type is not TH1F!!".format(self.datasetName,prefix+name))
        continue
      tmp.UseCurrentStyle()
      if histNames[name].has_key("rebin"):
        tmp.Rebin(histNames[name]["rebin"])
      tmp.SetLineColor(self.color)
      tmp.SetFillColor(self.color)
      tmp.Scale(self.scaleFactor)
      self.hists[prefix+name] = tmp

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
  if "wH" in i or "zH" in i:
    continue
  i += "_"+RUNPERIOD
  if i in scaleFactors:
    if scaleFactors[i]>0.0:
      filename = dataDir+i+".root"
      if not os.path.exists(filename):
          continue
      tmp = Dataset(filename,getLegendEntry(i),getColor(i),scaleFactors[i],isSignal=True)
      if tmp.isZombie():
        print ("Warning: file for dataset {0} is Zombie!!".format(i))
        continue
      #print("Loading Dataset: {0}".format(i))
      for hDir in histDirs:
        tmp.loadHistos(histNames,prefix=hDir)
      sigDatasetList.append(tmp)

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
canvas.SetLogy(LOGY)
leg = root.TLegend(*ulLegendPos)
leg.SetLineColor(0)
leg.SetFillColor(0)
uniqueLegendEntries = set()
#for ds in bkgDatasetList:
#  for hname in ds.hists:
#    if ds.legendEntry not in uniqueLegendEntries:
#      leg.AddEntry(ds.hists[hname],ds.legendEntry,"f")
#      uniqueLegendEntries.add(ds.legendEntry)
#    break
for ds in sigDatasetList:
  for hname in ds.hists:
    if ds.legendEntry not in uniqueLegendEntries:
      leg.AddEntry(ds.hists[hname],ds.legendEntry,"l")
      uniqueLegendEntries.add(ds.legendEntry)
    break
#for ds in realDatasetList:
#  for hname in ds.hists:
#    if ds.legendEntry not in uniqueLegendEntries:
#      leg.AddEntry(ds.hists[hname],ds.legendEntry,"pe")
#      uniqueLegendEntries.add(ds.legendEntry)
#    break
#  break

#######################################

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

  histName = os.path.split(histName)[1]

  bkgSumHist = bkgHistList[0].Clone()
  bkgSumHist.Reset()
  for i in bkgHistList:
    bkgSumHist.Add(i)

  if bkgSumHist.GetEntries() == 0:
    continue

  bkgIntHist = getIntegralHist(bkgSumHist)
  sigIntHistList = []
  for i in sigHistList:
    sigIntHistList.append(getIntegralHist(i))

  # Makes all together
  sumSigIntHist = None
  for i in sigIntHistList:
    if sumSigIntHist == None:
        sumSigIntHist = i.Clone("SigIntHistAll")
    else:
        sumSigIntHist.Add(i)
  sigIntHistList = [sumSigIntHist]

  sobHistList = []
  for hist in sigIntHistList:
    hist2 = hist.Clone("sob_"+hist.GetName())
    hist2.SetMarkerStyle(0)
    hist2.SetTitle('')
    hist2.GetXaxis().SetTitle(histNames[histName]['xlabel'])
    hist2.GetXaxis().SetRangeUser(*histNames[histName]['xlimits'])
    hist2.GetYaxis().SetTitle("S/B")
    hist2.GetYaxis().SetRangeUser(*histNames[histName]['ylimits'])
    hist2.Divide(bkgIntHist)
    sobHistList.append(hist2)

  drawn = False
  for sig in sobHistList:
    if drawn:
      sig.Draw(drawStr+" same")
    else:
      sig.Draw(drawStr)
      drawn = True

  drawLatex()
  setLegPos(leg,(ulLegendPos))
  #leg.Draw()

  vertLine = None
  if histNames[histName].has_key("vertLines"):
    vertLineX = histNames[histName]["vertLines"][RUNPERIOD]
    vertLine = root.TLine()
    vertLine.SetLineColor(root.kRed+1)
    vertLine.SetLineWidth(3)
    vertLine.DrawLine(vertLineX,sobHistList[0].GetMinimum(),vertLineX,sobHistList[0].GetMaximum())

  saveName = histName.replace("(","")
  saveName = saveName.replace(")","")
  saveName = saveName.replace("[","")
  saveName = saveName.replace("]","")
  saveName = saveName.replace("/","_")
  saveAs(canvas,outDir+"sob_"+saveName+"_"+RUNPERIOD)

  canvas.Clear()

  #sqrtTH1(bkgSumHist)
  sosqrtbHistList = []
  for hist in sigIntHistList:
    hist2 = hist.Clone("sosqrtb_"+hist.GetTitle())
    tmpBakHist = bkgIntHist.Clone("tmpBakInt"+str(random.randint(0,1000)))
    tmpBakHist.Add(hist2)
    sqrtTH1(tmpBakHist)
    hist2.Divide(tmpBakHist)
    hist2.SetMarkerStyle(0)
    hist2.SetTitle('')
    hist2.GetXaxis().SetTitle(histNames[histName]['xlabel'])
    hist2.GetYaxis().SetTitle("S/#sqrt{S+B}")
    hist2.GetXaxis().SetRangeUser(*histNames[histName]['xlimits'])
    hist2.GetYaxis().SetRangeUser(*histNames[histName]['ylimitsSqrt'])
    sosqrtbHistList.append(hist2)
  
  drawn = False
  for sig in sosqrtbHistList:
    if drawn:
      sig.Draw(drawStr+" same")
    else:
      sig.Draw(drawStr)
      drawn = True

  drawLatex()
  setLegPos(leg,(ulLegendPos))
  #leg.Draw()

  vertLine = None
  if histNames[histName].has_key("vertLines"):
    vertLineX = histNames[histName]["vertLines"][RUNPERIOD]
    vertLine = root.TLine()
    vertLine.SetLineColor(root.kRed+1)
    vertLine.SetLineWidth(3)
    vertLine.DrawLine(vertLineX,sosqrtbHistList[0].GetMinimum(),vertLineX,sosqrtbHistList[0].GetMaximum())

  saveAs(canvas,outDir+"sosqrtb_"+saveName+"_"+RUNPERIOD)

  # Square
  squareHistList = []
  sqrtsquareHistList = []
  for hist in sigIntHistList:
    hist2 = hist.Clone("ssqrob_"+hist.GetTitle())
    tmpBakHist = bkgIntHist.Clone("tmpBakInt"+str(random.randint(0,1000)))
    tmpBakHist.Add(hist2)
    hist2.Multiply(hist2)
    hist2.Divide(tmpBakHist)
    hist2.SetMarkerStyle(0)
    hist2.SetTitle('')
    hist2.GetXaxis().SetTitle(histNames[histName]['xlabel'])
    hist2.GetYaxis().SetTitle("#frac{S^{2}}{S+B}")
    hist2.GetXaxis().SetRangeUser(*histNames[histName]['xlimits'])
    hist2.GetYaxis().SetRangeUser(*histNames[histName]['ylimitsSqr'])
    squareHistList.append(hist2)
    sqrtsquareHist = hist2.Clone("sqrtsquare_"+hist.GetTitle())
    sqrtTH1(sqrtsquareHist)
    sqrtsquareHist.GetYaxis().SetTitle("#sqrt{#frac{S^{2}}{S+B}}")
    sqrtsquareHist.GetYaxis().SetRangeUser(*histNames[histName]['ylimitsSqrt'])
    sqrtsquareHistList.append(sqrtsquareHist)
  
  drawn = False
  for sig in squareHistList:
    if drawn:
      sig.Draw(drawStr+" same")
    else:
      sig.Draw(drawStr)
      drawn = True

  drawLatex()
  setLegPos(leg,(ulLegendPos))
  #leg.Draw()

  vertLine = None
  if histNames[histName].has_key("vertLines"):
    vertLineX = histNames[histName]["vertLines"][RUNPERIOD]
    vertLine = root.TLine()
    vertLine.SetLineColor(root.kRed+1)
    vertLine.SetLineWidth(3)
    vertLine.DrawLine(vertLineX,squareHistList[0].GetMinimum(),vertLineX,squareHistList[0].GetMaximum())

  saveAs(canvas,outDir+"ssqrob_"+saveName+"_"+RUNPERIOD)

  drawn = False
  for sig in sqrtsquareHistList:
    if drawn:
      sig.Draw(drawStr+" same")
    else:
      sig.Draw(drawStr)
      drawn = True
  saveAs(canvas,outDir+"sqrtsquare_"+saveName+"_"+RUNPERIOD)
