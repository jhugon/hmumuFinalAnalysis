#!/usr/bin/python

from helpers import *
import ROOT as root
import glob
import re

import matplotlib.pyplot as mpl
import numpy

dirName = "statsInput/"
caption2 = "#pm2 GeV Search Window"

outDir = "statsOutput/"

root.gROOT.SetBatch(True)
setStyle()
canvas = root.TCanvas()
canvas.SetLogx(1)
canvas.SetLogy(1)
#######################################

def getData(fileString,matchString=r"_([\d]+).txt.out",dontMatchStrings=[],doSort=True):
  def sortfun(s):
    match = re.search(matchString,s)
    result = 1e12
    if match:
      result = float(match.group(1))
    return result

  result = []
  fNames =  glob.glob(fileString)
  if doSort:
    fNames.sort(key=sortfun)
  #print fileString
  #print fNames
  for fname in fNames: 
    dontMatch = False
    for dont in dontMatchStrings:
      if re.search(dont,fname):
        dontMatch = True
    if dontMatch:
        continue
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
    if thisPoint.count("-10.0")>0:
        continue
    #print thisPoint
    result.append(thisPoint)
  return result

class RelativePlot:
  def __init__(self,dataPoints, canvas, legend, caption, ylabel="Expected 95% CL Limit #sigma/#sigma_{SM}", xlabel="Integrated Luminosity [fb^{-1}]",caption2="",ylimits=[],xlimits=[],vertLines=[20.0]):
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
      if len(xlimits)==2:
        if xNum<xlimits[0]:
            continue
        if xNum>xlimits[1]:
            continue
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

    self.vertLines = []
    self.vertLabel = []
    for xPos in vertLines:
      tmp = root.TGraph()
      tmp.SetPoint(0,xPos,ymin)
      tmp.SetPoint(1,xPos,ymax)
      tmp.SetLineColor(root.kRed)
      self.vertLines.append(tmp)
      self.vertLabel.append(xPos)

    label = root.TLatex()
    #label.SetNDC()
    label.SetTextFont(root.gStyle.GetLabelFont("X"))
    label.SetTextSize(root.gStyle.GetLabelSize("X"))
    label.SetTextColor(root.kRed)
    label.SetTextAlign(22)
    self.label=label

    twoSigGraph.Draw("a3")
    twoSigGraph.GetXaxis().SetTitle(xlabel)
    twoSigGraph.GetYaxis().SetTitle(ylabel)
    if len(ylimits)==2:
        twoSigGraph.GetYaxis().SetRangeUser(*ylimits)
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

    for g in self.vertLines:
      g.Draw("l")


class ComparePlot:
  def __init__(self,data,ylabel="Expected 95% CL Limit $\sigma/\sigma_{SM}$",titleMap={}):
    fig = mpl.figure()
    self.fig = fig
    ax1 = fig.add_subplot(111)
    self.ax1 = ax1
    ax1bounds = ax1.get_position().bounds
    ax1.set_position([0.25,0.1,0.7,0.85])

    data.sort(key=lambda x: x[0].lower())

    medians = [float(point[4]) for point in data]
    high1sigs = [float(point[5])-float(point[4]) for point in data]
    low1sigs = [float(point[4])-float(point[3]) for point in data]

    xPos = numpy.arange(len(medians))
    xLabels = [point[0] for point in data]
    xLabels = [re.sub(r".*/","",s) for s in xLabels]
    xLabels = [titleMap[s] if titleMap.has_key(s) else s for s in xLabels]

    ax1.set_yticks(xPos+0.25)
    ax1.set_yticklabels(tuple(xLabels),size="small")
    ax1.set_xlabel(ylabel)
    bars = ax1.barh(xPos,medians, 0.5, xerr=[low1sigs,high1sigs],ecolor="k")
    self.bars = bars
    writeInValues = getattr(self,"writeInValues")
    writeInValues(bars)
  def save(self,saveName):
    self.fig.savefig(saveName+".png")
    self.fig.savefig(saveName+".pdf")
  def writeInValues(self,rects):
    size="medium"
    if len(rects)<5:
      size="large"
    elif len(rects)>12:
      size="xx-small"
    elif len(rects)>9:
      size="x-small"
    for rect in rects:
      width = rect.get_width()
      if (width < 5): # The bars aren't wide enough to print the ranking inside
        xloc = width + 1 # Shift the text to the right side of the right edge
        clr = 'black' # Black against white background
        align = 'left'
      else:
        xloc = 0.98*width # Shift the text to the left side of the right edge
        clr = 'white' # White on magenta
        align = 'right'
      yloc = rect.get_y()+rect.get_height()/2.0 
      valueStr = "{0:.1f}".format(width)
      self.ax1.text(xloc, yloc, valueStr, horizontalalignment=align,
            verticalalignment='center', color=clr, weight='bold',size=size)

    

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

comparisonMap = {
  "combined":"C&C Comb.",
  "combinedVBFOnly":"C&C VBF Comb.",
  "combinedMuOnly":"C&C !VBF Comb.",
  "BDTCombination":"BDT Comb.",
  "BDTVBF":"BDT VBF",
  "BDTMuonOnly":"BDT !VBF",
  "VBFL":"VBFL",
  "VBFM":"VBFM",
  "VBFT":"VBFT",
  "PtL30":"$p_{T}^{\mu\mu} < 30$ GeV",
  "Pt30to50":"$p_{T}^{\mu\mu} \in [30,50]$ GeV",
  "Pt50to75":"$p_{T}^{\mu\mu} \in [50,75]$ GeV",
  "Pt75to125":"$p_{T}^{\mu\mu} \in [75,125]$ GeV",
  "Pt125":"$p_{T}^{\mu\mu} > 125$ GeV"
}

ylimits=[1.0,500.0]

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
  incPlot = RelativePlot(data,canvas,legend,titleMap[plotName],caption2=caption2,ylimits=ylimits)
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
  incPlot = RelativePlot(data,canvas,legend,titleMap[plotName]+" L=20fb^{-1}",caption2=caption2,xlabel="BDT Cut",xlimits=[-1000,0.12])
  saveAs(canvas,outDir+"BDTCut"+plotName)

## Compare all types of limits
compareData = getData(dirName+"*_20.txt.out",matchString=r"(.*)_[\d]+.txt.out",dontMatchStrings=[r"BDT.+BDT",r"PM"],doSort=False)
#print compareData
comparePlot = ComparePlot(compareData,titleMap=comparisonMap)
comparePlot.fig.text(0.9,0.2,"$L=20$ fb$^{-1}$",horizontalalignment="right",size="x-large")
comparePlot.save(outDir+"compare")

## Compare Best Performing Limits
compareData = getData(dirName+"*_20.txt.out",matchString=r"(.*)_[\d]+.txt.out",dontMatchStrings=[r"BDT.+BDT",r"PM","VBFT","VBFL","VBFM","Pt"],doSort=False)
comparePlot = ComparePlot(compareData,titleMap=comparisonMap)
comparePlot.fig.text(0.9,0.2,"$L=20$ fb$^{-1}$",horizontalalignment="right",size="x-large")
comparePlot.save(outDir+"compareGood")

## Compare BDT Limits
compareData = getData(dirName+"*BDT*_20.txt.out",matchString=r"(.*)_[\d]+.txt.out",dontMatchStrings=[r"BDT.+BDT",r"PM"],doSort=False)
comparePlot = ComparePlot(compareData,titleMap=comparisonMap)
comparePlot.fig.text(0.9,0.2,"$L=20$ fb$^{-1}$",horizontalalignment="right",size="x-large")
comparePlot.save(outDir+"compareBDT")

## Compare Not-BDT Limits
compareData = getData(dirName+"*_20.txt.out",matchString=r"(.*)_[\d]+.txt.out",dontMatchStrings=[r"BDT.+BDT",r"PM",r"BDT"],doSort=False)
comparePlot = ComparePlot(compareData,titleMap=comparisonMap)
comparePlot.fig.text(0.9,0.2,"$L=20$ fb$^{-1}$",horizontalalignment="right",size="x-large")
comparePlot.save(outDir+"compareNonBDT")
