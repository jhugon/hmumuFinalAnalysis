#!/usr/bin/env python

from singleHelpers import *
from helpers import *
import ROOT as root
import glob
import re
import os.path
import cPickle
from copy import deepcopy
from array import array
import sys

from xsec import *

#######################################

def getData(fileString,matchString=r"_N([\d]+)\.txt\.out",dontMatchStrings=[],doSort=True,xMax=1.0e20):
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
      if float(xNum) > xMax:
        continue
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
    if thisPoint.count(-10.0)>0:
        continue
    #print thisPoint
    result.append(thisPoint)
      
  return result

if __name__ == "__main__":
  root.gROOT.SetBatch(True)

  data = getData("statsCards/*.out")
  canvas = root.TCanvas()
  axisHist = root.TH2F("axisHist","",1,0,100,1,0,30)
  setHistTitles(axisHist,"N_{bias}","Expected Limit #sigma/#sigma_{SM}")
  graph = root.TGraph()
  graphPred = root.TGraph()
  for i,p in enumerate(data):
    nb = float(p[0])
    limitMu = float(p[4])
    limitN = limitMu*9.
    graph.SetPoint(i,nb,limitMu)
    graphPred.SetPoint(i,nb,sqrt(12.9**2+(1.65*nb/9.)**2))

  #fr = graph.Fit("pol2","MEFS","",0.,100.)
  
  graphPred.SetLineColor(root.kBlue)

  leg = root.TLegend(0.2,0.15,0.8,0.45)
  leg.SetFillColor(0)
  leg.SetLineColor(0)
  leg.AddEntry(graph,"Computed Limits","p")
  leg.AddEntry(graphPred,"My Estimate","l")


  axisHist.Draw()
  leg.Draw()
  graph.Draw("P")
  graphPred.Draw("L")
  tlatex = drawStandardCaptions(canvas,"0,1-Jet Tight BB","#sqrt{s} = 8 TeV, L = 19.7 fb^{-1}","0,1-Jet Tight BB",preliminaryString="CMS Internal")
  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.88,"")
  saveAs(canvas,"~/public_html/hmumu/toyTest/Stupid")
