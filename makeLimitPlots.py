#!/usr/bin/python

from helpers import *
import ROOT as root
import glob
import re

import matplotlib.pyplot as mpl
import numpy

dirName = "statsInput/"
caption2 = "1D Shape Analysis"

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
    maxWidth = 0.0
    for rect in rects:
      width = rect.get_width()
      if width>maxWidth:
        maxWidth = width
    for rect in rects:
      width = rect.get_width()
      if (width/maxWidth < 0.1): # The bars aren't wide enough to print the ranking inside
        xloc = width + 0.1/maxWidth # Shift the text to the right side of the right edge
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
  "BDTComb":"Combined BDT H#rightarrow#mu#mu",
  "BDTComb1d":"Combined BDT 1D H#rightarrow#mu#mu",
  "LHComb":"Combined BDT H#rightarrow#mu#mu",
  "BDTHistMuonOnly":"H#rightarrow#mu#mu: BDT Inclusive",
  "BDTHistVBF":"H#rightarrow#mu#mu: BDT VBF",
  "likelihoodHistMuonOnly":"H#rightarrow#mu#mu: Likelihood Inclusive",
  "likelihoodHistVBF":"H#rightarrow#mu#mu: Likelihood VBF",
  "BDTHistMuonOnlyVMass":"H#rightarrow#mu#mu: BDT v. m_{#mu#mu} Inclusive",
  "BDTHistVBFVMass":"H#rightarrow#mu#mu: BDT v. m_{#mu#mu} VBF",
  "likelihoodHistMuonOnlyVMass":"H#rightarrow#mu#mu: Likelihood v. m_{#mu#mu} Inclusive",
  "likelihoodHistVBFVMass":"H#rightarrow#mu#mu: Likelihood v. m_{#mu#mu} VBF",
  "mDiMu":"H#rightarrow#mu#mu: m_{#mu#mu}",

  "THBDTHistMuonOnlyVMass":"TH2 H#rightarrow#mu#mu: BDT v. m_{#mu#mu} Inclusive",
  "THBDTHistVBFVMass":"TH2 H#rightarrow#mu#mu: BDT v. m_{#mu#mu} VBF",
  "THlikelihoodHistMuonOnlyVMass":"TH2 H#rightarrow#mu#mu: Likelihood v. m_{#mu#mu} Inclusive",
  "THlikelihoodHistVBFVMass":"TH2 H#rightarrow#mu#mu: Likelihood v. m_{#mu#mu} VBF"
}

comparisonMap = {
  "combined":"Combination",
  "BDTComb":"BDT Comb.",
  "BDTComb1d":"BDT 1D Comb.",
  "LHComb":"LH Comb.",
  "BDTHistMuonOnly":"BDT Inc.",
  "BDTHistVBF":"BDT VBF",
  "likelihoodHistMuonOnly":"LH Inc.",
  "likelihoodHistVBF":"LH VBF",

  "BDTHistMuonOnlyVMass":"BDT v. $m_{\mu\mu}$ Inc.",
  "BDTHistVBFVMass":"BDT v. $m_{\mu\mu}$ VBF",
  "likelihoodHistMuonOnlyVMass":"LH v. $m_{\mu\mu}$ Inc.",
  "likelihoodHistVBFVMass":"LH v. $m_{\mu\mu}$ VBF",
  "mDiMu":"$m_{\mu\mu}$ Shape",

  "THBDTHistMuonOnlyVMass":"TH2 BDT v. $m_{\mu\mu}$ Inc.",
  "THBDTHistVBFVMass":"TH2 BDT v. $m_{\mu\mu}$ VBF",
  "THlikelihoodHistMuonOnlyVMass":"TH2 LH v. $m_{\mu\mu}$ Inc.",
  "THlikelihoodHistVBFVMass":"TH2 LH v. $m_{\mu\mu}$ VBF"
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

## Compare all types of limits
compareData = getData(dirName+"*_20.txt.out",matchString=r"(.*)_[\d]+.txt.out",dontMatchStrings=[r"BDT.+BDT",r"PM"],doSort=False)
#print compareData
comparePlot = ComparePlot(compareData,titleMap=comparisonMap)
comparePlot.fig.text(0.9,0.2,"$L=20$ fb$^{-1}$",horizontalalignment="right",size="x-large")
comparePlot.save(outDir+"compare")
