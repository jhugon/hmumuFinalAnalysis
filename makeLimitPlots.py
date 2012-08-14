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

def getData(fileString):
  def tryint(s):
    try:
        return int(s)
    except:
        return s

  def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

  result = []
  fNames =  glob.glob(fileString)
  fNames.sort(key=alphanum_key)
  #print fileString
  #print fNames
  for fname in fNames: 
    tmpF = open(fname)
    match = re.search(r"_([\d]+).txt.out",fname)
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
    tlatex.SetTextSize(0.05)
    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,"CMS Preliminary")
    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,caption)


allfiles = glob.glob(dirName+"*.txt.out")
plots = set()
for fn in allfiles:
  match = re.search(r".*/(.*)_[\d]+.txt.out",fn)
  badPlot = re.search(r"combine125PM",fn)
  if match and not badPlot:
    plots.add(match.group(1))

legend = root.TLegend(0.58,0.70,0.9,0.9)
legend.SetFillColor(0)
legend.SetLineColor(0)
for plotName in plots:
  data = getData(dirName+plotName+"_*.txt.out")
  if len(data)==0:
    continue
  incPlot = RelativePlot(data,canvas,legend,"H->#mu#mu",caption2=caption2)
  saveAs(canvas,outDir+plotName)
    

