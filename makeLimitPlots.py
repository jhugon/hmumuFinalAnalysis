#!/usr/bin/python

from helpers import *
import ROOT as root
import glob
import re

dirName = "statsInput/"
caption2 = "#pm2 GeV Search Window"

outDir = "statsOutput/"

root.gROOT.SetBatch(True)
setStyle()
canvas = root.TCanvas()
canvas.SetLogx(1)
canvas.SetLogy(1)
#######################################

def getData(fileString,matchString=r"_([\d]+).txt.out"):
  def sortfun(s):
    match = re.search(matchString,s)
    result = 1e12
    if match:
      result = float(match.group(1))
    return result

  result = []
  fNames =  glob.glob(fileString)
  fNames.sort(key=sortfun)
  #print fileString
  #print fNames
  for fname in fNames: 
    tmpF = open(fname)
    match = re.search(matchString,fname)
    obs = -10.0
    median = -10.0
    low2sig = -10.0
    low1sig = -10.0
    high1sig = -10.0
    high2sig = -10.0
    xNum = -10.0
    if match:
      xNum = match.group(1)
    for line in tmpF:
      obsMatch = re.search(r"Observed[\s]Limit:[^.\d]*< ([.\deE]+)",line)
      low2sigMatch = re.search(r"Expected.*2\.5.:[^.\d]*< ([.\deE]+)",line)
      low1sigMatch = re.search(r"Expected.*16.0[^.\d]*< ([.\deE]+)",line)
      medianMatch = re.search(r"Expected.*50\.0.*< ([.\deE]+)",line)
      high1sigMatch = re.search(r"Expected.*84.0.*< ([.\deE]+)",line)
      high2sigMatch = re.search(r"Expected.*97.5.*< ([.\deE]+)",line)
      if obsMatch:
        obs = obsMatch.group(1)
      if low2sigMatch:
        low2sig = low2sigMatch.group(1)
      if low1sigMatch:
        low1sig = low1sigMatch.group(1)
      if medianMatch:
        median = medianMatch.group(1)
      if high1sigMatch:
        high1sig = high1sigMatch.group(1)
      if high2sigMatch:
        high2sig = high2sigMatch.group(1)
    thisPoint = [xNum,obs,low2sig,low1sig,median,high1sig,high2sig]
    if thisPoint.count(-10.0)>0:
        continue
    #print thisPoint
    result.append(thisPoint)
  return result

class RelativePlot:
  def __init__(self,dataPoints, canvas, legend, caption, ylabel="Expected 95% CL Limit #sigma/#sigma_{SM}", xlabel="Integrated Luminosity [fb^{-1}]",caption2=""):
    expGraph = root.TGraph()
    oneSigGraph = root.TGraphAsymmErrors()
    oneSigGraph.SetFillColor(root.kGreen)
    oneSigGraph.SetLineStyle(0)
    twoSigGraph = root.TGraphAsymmErrors()
    twoSigGraph.SetFillColor(root.kYellow)
    twoSigGraph.SetLineStyle(0)
    oneGraph = root.TGraph()
    oneGraph.SetLineColor(root.kRed)
    oneGraph.SetLineStyle(2)
    self.expGraph = expGraph
    self.oneSigGraph = oneSigGraph
    self.twoSigGraph = twoSigGraph
    self.oneGraph = oneGraph
    iPoint = 0
    ymax = 1e-20
    ymin = 1e20
    for point in dataPoints:
      xNum = float(point[0])
      obs = float(point[1])
      median = float(point[4])
      low2sig = median - float(point[2])
      low1sig = median - float(point[3])
      high1sig = float(point[5]) - median
      high2sig = float(point[6]) - median
      expGraph.SetPoint(iPoint,xNum,median) 
      oneSigGraph.SetPoint(iPoint,xNum,median) 
      twoSigGraph.SetPoint(iPoint,xNum,median) 
      oneSigGraph.SetPointError(iPoint,0.0,0.0,low1sig,high1sig) 
      twoSigGraph.SetPointError(iPoint,0.0,0.0,low2sig,high2sig) 
      oneGraph.SetPoint(iPoint,xNum,1.0)
      iPoint += 1

    twoSigGraph.Draw("a3")
    twoSigGraph.GetXaxis().SetTitle(xlabel)
    twoSigGraph.GetYaxis().SetTitle(ylabel)
    oneSigGraph.Draw("3")
    expGraph.Draw("l")
    oneGraph.Draw("l")

    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(root.gStyle.GetLabelFont())
    #tlatex.SetTextSize(0.05)
    tlatex.SetTextSize(0.04)
    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,"CMS Preliminary")
    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,caption)


titleMap = {
  "combined":"Combined H#rightarrow#mu#mu",
  "combinedVBFOnly":"Combined H#rightarrow#mu#mu VBF Channels",
  "combinedMuOnly":"Combined H#rightarrow#mu#mu Non-VBF Channels",
  "BDTCombination":"Combined BDT H#rightarrow#mu#mu",
  "BDTVBF":"H#rightarrow#mu#mu VBF BDT",
  "BDTMuonOnly":"H#rightarrow#mu#mu Non-VBF BDT",
  "VBFL":"H#rightarrow#mu#mu, VBFL",
  "VBFM":"H#rightarrow#mu#mu, VBFM",
  "VBFT":"H#rightarrow#mu#mu, VBFT",
  "PtL30":"H#rightarrow#mu#mu, p_{T}(#mu#mu) < 30 GeV",
  "Pt30to50":"H#rightarrow#mu#mu, p_{T}(#mu#mu) #in [30,50] GeV",
  "Pt50to75":"H#rightarrow#mu#mu, p_{T}(#mu#mu) #in [50,75] GeV",
  "Pt75to125":"H#rightarrow#mu#mu, p_{T}(#mu#mu) #in [75,125] GeV",
  "Pt125":"H#rightarrow#mu#mu, p_{T}(#mu#mu) > 125 GeV"
}

allfiles = glob.glob(dirName+"*.txt.out")
## Limit v. Lumi
plots = set()
for fn in allfiles:
  match = re.search(r".*/(.*)_[\d]+.txt.out",fn)
  badPlot = re.search(r"PM",fn)
  badPlot2 = re.search(r"BDT.*BDT",fn)
  if match and not (badPlot or badPlot2):
    plots.add(match.group(1))

legend = root.TLegend(0.58,0.70,0.9,0.9)
legend.SetFillColor(0)
legend.SetLineColor(0)
for plotName in plots:
  data = getData(dirName+plotName+"_*.txt.out")
  if len(data)==0:
    continue
  incPlot = RelativePlot(data,canvas,legend,titleMap[plotName],caption2=caption2)
  saveAs(canvas,outDir+plotName)

#Limit v. Window
plots = set()
for fn in allfiles:
  match = re.search(r".*/PM(.*)PM",fn)
  if match and not badPlot:
    plots.add(match.group(1))

canvas.SetLogx(0)
for plotName in plots:
  data = getData(dirName+"PM"+plotName+"PM*.txt.out",matchString=r"PM([.\d]+).*.txt.out")
  if len(data)==0:
    continue
  incPlot = RelativePlot(data,canvas,legend,titleMap[plotName]+" L=20fb^{-1}",caption2=caption2,xlabel="Search Window: m_{#mu#mu} #in 125 GeV #pm X [GeV]")
  saveAs(canvas,outDir+"PM"+plotName)

#Limit v. BDT Cut
plots = set()
for fn in allfiles:
  match = re.search(r".*/BDT(.*)BDT",fn)
  if match and not badPlot:
    plots.add(match.group(1))

titleMap["Mu"] = "H#rightarrow#mu#mu Non-VBF BDT"
titleMap["VBF"] = "H#rightarrow#mu#mu VBF BDT"

for plotName in plots:
  data = getData(dirName+"BDT"+plotName+"BDT*.txt.out",matchString=r"BDT([.\d-]+).*.txt.out")
  if len(data)==0:
    continue
  incPlot = RelativePlot(data,canvas,legend,titleMap[plotName]+" L=20fb^{-1}",caption2=caption2,xlabel="BDT Cut")
  saveAs(canvas,outDir+"BDTCut"+plotName)
