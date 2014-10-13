#!/usr/bin/env python

import optparse
parser = optparse.OptionParser(description="Makes Limit Plots from output text from combine tool.")
parser.add_option("--xs", help="Limit on XS*BR",action="store_true",default=True)
args, fakeargs = parser.parse_args()

import singleHelpers
from helpers import *
import ROOT as root
import glob
import re
import os.path
from copy import deepcopy
from array import array
import sys
import csv

from xsec import *

##############################################
# Only use Matplotlib if not in CMSSW_6*
cmsswVersion = ""
mplGood = True
if os.environ.has_key("CMSSW_VERSION"):
  cmsswVersion = os.environ["CMSSW_VERSION"]
#def ComparePlotTable(*arg):
#  #Dummy so that there are no errors
#  pass
if "CMSSW_6" in cmsswVersion:
  mplGood = False
if mplGood:
  from etc.evilMatPlotlibFunctions import *
##############################################


def getData():

   heeF = open("etc/limits_xsbr_hee.csv")
   data = csv.reader(heeF)
   result = []
   for line in data:
     if re.match(r"[0-9.-]+",line[0]):
       line = [float(i) for i in line]
       x = line[0]
       obs = line[6]
       exp = line[5]
       p1sig = exp+line[3]
       m1sig = exp-line[4]
       p2sig = exp+line[1]
       m2sig = exp-line[2]
       result.append([x,obs,m2sig,m1sig,exp,p1sig,p2sig])
   for i in result:
    print i
   return result

class RelativePlot:
  def __init__(self,dataPoints, canvas, legend, caption, ylabel="95% CL limit on #sigma/#sigma_{SM} (H #rightarrow #mu#mu)", xlabel="Integrated Luminosity [fb^{-1}]",caption2="",caption3="",ylimits=[],xlimits=[],vertLines=[],showObs=False,energyStr="8TeV",caption4=""):
    expGraph = root.TGraph()
    expGraph.SetLineStyle(2)
    expGraph.SetMarkerStyle(20)
    expGraph.SetMarkerSize(0.9)
    oneSigGraph = root.TGraphAsymmErrors()
    oneSigGraph.SetFillColor(root.kGreen)
    oneSigGraph.SetLineColor(root.kGreen)
    oneSigGraph.SetLineStyle(0)
    twoSigGraph = root.TGraphAsymmErrors()
    twoSigGraph.SetFillColor(root.kYellow)
    twoSigGraph.SetLineColor(root.kYellow)
    twoSigGraph.SetLineStyle(0)
    oneGraph = root.TGraph()
    oneGraph.SetLineColor(root.kRed)
    oneGraph.SetLineStyle(3)
    obsGraph = root.TGraph()
    obsGraph.SetLineColor(1)
    obsGraph.SetLineStyle(1)
    obsGraph.SetLineWidth(3)
    obsGraph.SetMarkerStyle(20)
    obsGraph.SetMarkerSize(1.1)
    obsGraph.SetMarkerColor(1)
    self.expGraph = expGraph
    self.oneSigGraph = oneSigGraph
    self.twoSigGraph = twoSigGraph
    self.oneGraph = oneGraph
    self.obsGraph = obsGraph
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
      obsGraph.SetPoint(iPoint,xNum,obs)
      iPoint += 1

      if float(point[6]) > ymax:
        ymax = float(point[6])
      if obs > ymax:
        ymax = obs

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
    if args.xs:
      twoSigGraph.GetYaxis().SetTitleSize(0.9*twoSigGraph.GetYaxis().GetTitleSize())
      twoSigGraph.GetYaxis().SetTitleOffset(1.15*twoSigGraph.GetYaxis().GetTitleOffset())
    if len(ylimits)==2:
        twoSigGraph.GetYaxis().SetRangeUser(*ylimits)
    else:
        twoSigGraph.GetYaxis().SetRangeUser(0.0,ymax*1.1)
        if args.xs:
          twoSigGraph.GetYaxis().SetRangeUser(0.0,ymax*2.)
    if len(xlimits)==2:
        twoSigGraph.GetXaxis().SetRangeUser(*xlimits)
    oneSigGraph.Draw("3")
    expGraph.Draw("l")
    expGraph.Draw("p")
    if showObs:
      obsGraph.Draw("l")
      obsGraph.Draw("p")

    legPos = [gStyle.GetPadLeftMargin()+0.03,0.55,0.68,0.93-gStyle.GetPadTopMargin()]
    if args.xs:
      legPos = [0.45,0.55,0.97-gStyle.GetPadRightMargin(),0.93-gStyle.GetPadTopMargin()]
    if caption=="" or not caption:
      legPos = [0.179598,0.592262,0.679598,0.892857]
    if args.xs and (caption=="" or not caption):
      legPos = [0.451149,0.592262,0.920977,0.892857]
    leg = root.TLegend(*legPos)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    leg.AddEntry(obsGraph,"Observed limit","lp")
    leg.AddEntry(expGraph,"Median expected limit","lp")
    leg.AddEntry(oneSigGraph,"68% CL expected limit","f")
    leg.AddEntry(twoSigGraph,"95% CL expected limit","f")
    self.legPos = legPos
    self.leg = leg
    leg.Draw()

    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(62)
    tlatex.SetTextSize(0.06)
    tlatex.SetTextAlign(11)
    if args.xs:
      tlatex.DrawLatex(0.18,0.84,PRELIMINARYSTRING)
    else:
      tlatex.DrawLatex(0.8,0.84,PRELIMINARYSTRING)
    tlatex.SetTextFont(42)
    tlatex.SetTextSize(0.04)
    tlatex.DrawLatex(0.15,0.94,"H #rightarrow e^{+}e^{-}")
    tlatex.SetTextAlign(31)
    tlatex.SetTextSize(0.04)
    tlatex.DrawLatex(0.95,0.94,caption2)
    tlatex.SetTextSize(0.04)
    tlatex.SetTextAlign(11)
    tlatex.DrawLatex(0.19,0.86,caption)

    self.vertLine = root.TLine()
    self.vertLine.SetLineColor(root.kRed)
    self.vertLine.SetLineWidth(3)
    self.arrows = []
    for xPos in vertLines:
      self.vertLine.DrawLine(xPos,ylimits[0],xPos,ylimits[1])
      xAxis = twoSigGraph.GetXaxis()
      arrowLength = (xAxis.GetXmax() - xAxis.GetXmin())/10.
      arrowY = (ylimits[1]-ylimits[0])*0.8
      arrowHeadSize = 0.025
      arrow = root.TArrow(xPos,arrowY,xPos+arrowLength,arrowY,
                                arrowHeadSize,"|>")
      arrow.SetLineWidth(3)
      arrow.SetLineColor(root.kRed)
      arrow.SetFillColor(root.kRed)
      #arrow.SetAngle(40)
      arrow.Draw()
      self.arrows += [arrow]

    canvas.RedrawAxis()

if __name__ == "__main__":

  dirName = "statsInput/"
  
  outDir = "statsOutput/"
  
  root.gErrorIgnoreLevel = root.kWarning
  root.gROOT.SetBatch(True)
  setStyle()

  canvas = root.TCanvas()

  PRELIMINARYSTRING="CMS"

  for period in ["8TeV"]:
    energyStr = "8TeV"
    energyStrWrite = energyStr
    if energyStr == "7P8TeV":
      energyStrWrite = "7 & 8 TeV"
    else:
      energyStrWrite = energyStr.replace("TeV"," TeV")
    caption2 = "#sqrt{s} = "+energyStrWrite
    legend = root.TLegend(0.58,0.70,0.9,0.9)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    for plotName in [""]:
      data = getData()
      xlimits = []
      ylimits = [0.,0.14]
      vertLines = []
      if len(data)<=1:
        continue
      xlabel="m_{H} [GeV]"
      ylabel="95% CL upper limit on #sigma #times B(H #rightarrow e^{+}e^{-}) [pb]"
      caption2 = "{0:.1f} fb^{{-1}} (8 TeV)".format(19.7)
      caption3 = ""
      caption4 = ""
      title = ""
      incPlot = RelativePlot(data,canvas,legend,title,caption2=caption2,ylimits=ylimits,energyStr=energyStrWrite,xlabel=xlabel,caption3=caption3,showObs=True,xlimits=xlimits,vertLines = vertLines,ylabel=ylabel,caption4=caption4)
      saveAs(canvas,outDir+"xsbr_EE_"+energyStr)
