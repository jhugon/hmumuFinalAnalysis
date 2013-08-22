#!/usr/bin/env python

from helpers import *
import ROOT as root
import glob
import re
import os.path
from copy import deepcopy
from array import array
import sys

from xsec import *

#######################################

class RelativePlot:
  def __init__(self,dataPoints, canvas, legend, caption, ylabel="95% CL Limit on #sigma/#sigma_{SM} (H#rightarrow#mu#mu)", xlabel="Integrated Luminosity [fb^{-1}]",caption2="",caption3="",ylimits=[],xlimits=[],vertLines=[],showObs=False,energyStr="8TeV",caption4=""):
    expGraph = root.TGraph()
    expGraph.SetLineStyle(2)
    expGraph.SetMarkerStyle(20)
    expGraph.SetMarkerSize(0.9)
    oneSigGraph = root.TGraphAsymmErrors()
    oneSigGraph.SetFillColor(root.kGreen)
    oneSigGraph.SetLineStyle(0)
    twoSigGraph = root.TGraphAsymmErrors()
    twoSigGraph.SetFillColor(root.kYellow)
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
    if len(ylimits)==2:
        twoSigGraph.GetYaxis().SetRangeUser(*ylimits)
    else:
        twoSigGraph.GetYaxis().SetRangeUser(0.0,ymax*1.1)
    if len(xlimits)==2:
        twoSigGraph.GetXaxis().SetRangeUser(*xlimits)
    oneSigGraph.Draw("3")
    expGraph.Draw("l")
    expGraph.Draw("p")
    if showObs:
      obsGraph.Draw("l")
      obsGraph.Draw("p")

    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(root.gStyle.GetLabelFont())
    #tlatex.SetTextSize(0.05)
    tlatex.SetTextSize(0.04)
    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.88,caption2)
    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.82,caption3)
    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.88,caption4)

    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,caption)

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
  import csv
  root.gErrorIgnoreLevel = root.kWarning
  root.gROOT.SetBatch(True)
  setStyle()
  outDir = "statsOutput/"

  ###########################################################################

  eeDataFile = open("table_limit_hee.csv")
  mmDataFile = open("xsbr.txt")
  eeDataReader = csv.reader(eeDataFile)
  massList = []
  mmMassList = []
  eeExpList = []
  eeObsList = []
  mmExpList = []
  mmObsList = []
  for row in eeDataReader:
    if row[0]=="Mass":
      continue
    massList.append(float(row[0]))
    eeObsList.append(float(row[1]))
    eeExpList.append(float(row[2]))
  for row in mmDataFile:
    if row[0] == '#':
      continue
    row = row[:-1]
    row = row.split(" ")
    for i in range(row.count('')):
      row.remove('')
    if row[7] != "8TeV":
      continue
    if float(row[0])<120. or float(row[0])>150.:
        continue
    mmMassList.append(float(row[0]))
    mmExpList.append(float(row[4]))
    mmObsList.append(float(row[1]))
  eeDataFile.close()
  mmDataFile.close()

  mmExpGraph = root.TGraph()
  mmObsGraph = root.TGraph()
  eeExpGraph = root.TGraph()
  eeObsGraph = root.TGraph()
  eeOvermmExpGraph = root.TGraph()
  for i,mass,mmMass,eeExp,eeObs,mmExp,mmObs in zip(range(len(massList)),massList,mmMassList,eeExpList,eeObsList,mmExpList,mmObsList):
    assert(mass==mmMass)
   #print("mass: {0:3.0f} eeExp: {1:5.3f} mmExp: {2:5.3f} eeObs: {3:5.3f} mmObs: {4:5.3f}".format(mass,eeExp,mmExp,eeObs,mmObs))
    if mass==125.:
      print("mass: {0:3.0f} eeExp: {1:5.3f} mmExp: {2:5.3f} ratio: {3:5.3f} ".format(mass,eeExp,mmExp,eeExp/mmExp))
    mmExpGraph.SetPoint(i,mass,mmExp)
    mmObsGraph.SetPoint(i,mass,mmObs)
    eeExpGraph.SetPoint(i,mass,eeExp)
    eeObsGraph.SetPoint(i,mass,eeObs)
    eeOvermmExpGraph.SetPoint(i,mass,eeExp/mmExp)

  mmExpGraph.SetLineStyle(2)
  eeExpGraph.SetLineStyle(2)
  eeExpGraph.SetLineColor(root.kRed+1)
  eeObsGraph.SetLineColor(root.kRed+1)
  eeExpGraph.SetMarkerColor(root.kRed+1)
  eeObsGraph.SetMarkerColor(root.kRed+1)
  mmExpGraph.SetMarkerStyle(20)
  mmObsGraph.SetMarkerStyle(20)
  eeExpGraph.SetMarkerStyle(20)
  eeObsGraph.SetMarkerStyle(20)
  mmExpGraph.SetMarkerSize(0.8)
  mmObsGraph.SetMarkerSize(0.8)
  eeExpGraph.SetMarkerSize(0.8)
  eeObsGraph.SetMarkerSize(0.8)

  eeOvermmExpGraph.SetMarkerStyle(20)
  eeOvermmExpGraph.SetMarkerSize(0.8)

  ###########################################################################
  
  gStyle.SetPadLeftMargin(gStyle.GetPadLeftMargin()*1.1)
  gStyle.SetTitleOffset(gStyle.GetTitleOffset("Y")*1.22,"Y")
  canvas = root.TCanvas()
  xlabel="m_{H} [GeV/c^{2}]"
  ylabel="95% CL Limit on #sigma #times BR [pb]"
  caption1 = ""
  caption2 = "#sqrt{s} = 8 TeV"
  caption3 = ""
  caption4 = ""
  #caption1 = "#sqrt{s} = 8 TeV"
  #caption2 = "H#rightarrowe^{+}e^{-} L = 19.6 fb^{-1}"
  #caption3 = "H#rightarrow#mu^{+}#mu^{-} L = 19.8 fb^{-1}"
  #caption4 = ""
  #legend = root.TLegend(0.58,0.70,0.9,0.9)
  legend = root.TLegend(0.40,0.75,0.9,0.9)
  legend.SetFillColor(0)
  legend.SetLineColor(0)
  legend.AddEntry(mmObsGraph,"H#rightarrowe^{+}e^{-} L = 19.6 fb^{-1}","lp")
  legend.AddEntry(eeObsGraph,"H#rightarrow#mu^{+}#mu^{-} L = 19.8 fb^{-1}","lp")
  #mg = root.TMultiGraph()
  #mg.Add(mmExpGraph)
  #mg.Add(eeExpGraph)
  #mg.Add(mmObsGraph)
  #mg.Add(eeObsGraph)
  #mg.Draw("APL")
  #setHistTitles(mg,xlabel,ylabel)
  #incPlot = RelativePlot(data,canvas,legend,title,caption2=caption2,ylimits=ylimits,energyStr=energyStrWrite,xlabel=xlabel,caption3=caption3,showObs=args.higgsMass,xlimits=xlimits,vertLines = vertLines,ylabel=ylabel,caption4=caption4)
  axisHist = root.TH2F("axisHistMG","",1,115,155,1,0,0.07)
  setHistTitles(axisHist,xlabel,ylabel)
  axisHist.Draw()
  mmExpGraph.Draw("L")
  mmObsGraph.Draw("LP")
  eeExpGraph.Draw("L")
  eeObsGraph.Draw("LP")
  legend.Draw()

  tlatex = root.TLatex()
  tlatex.SetNDC()
  tlatex.SetTextFont(root.gStyle.GetLabelFont())
  #tlatex.SetTextSize(0.05)
  tlatex.SetTextSize(0.04)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(0.03+gStyle.GetPadLeftMargin(),0.88,caption2)
  tlatex.DrawLatex(0.03+gStyle.GetPadLeftMargin(),0.82,caption3)
  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.88,caption4)

  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,caption1)

  tlatex.SetTextSize(0.045)
  tlatex.SetTextAlign(21)
  tlatex.DrawLatex((1.0-gStyle.GetPadRightMargin()+gStyle.GetPadLeftMargin())/2.,gStyle.GetPadBottomMargin()+0.05,"Dotted Lines: Expected Limit")
  saveAs(canvas,outDir+"eeVmm_xsbr_8TeV")

  #########################################################3

  gStyle.SetTitleOffset(gStyle.GetTitleOffset("Y")*0.85,"Y")
  canvas.Clear()
  axisHist = root.TH2F("axisHistRatio","",1,115,155,1,1.,1.7)
  #setHistTitles(axisHist,xlabel,"H#rightarrowe^{+}e^{-}/H#rightarrow#mu^{+}#mu^{-} Ratio of Upper Limit on #sigma #times BR")
  setHistTitles(axisHist,xlabel,"Ratio of #sigma #times BR Upper Limits: H#rightarrowee / H#rightarrow#mu#mu")
  axisHist.Draw()
  eeOvermmExpGraph.Draw("LP")

  caption1 = ""
  caption2 = "H#rightarrowe^{+}e^{-} L = 19.6 fb^{-1}"
  caption3 = "H#rightarrow#mu^{+}#mu^{-} L = 19.8 fb^{-1}"
  caption4 = "#sqrt{s} = 8 TeV"
  
  tlatex = root.TLatex()
  tlatex.SetNDC()
  tlatex.SetTextFont(root.gStyle.GetLabelFont())
  #tlatex.SetTextSize(0.05)
  tlatex.SetTextSize(0.04)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(0.03+gStyle.GetPadLeftMargin(),0.88,caption2)
  tlatex.DrawLatex(0.03+gStyle.GetPadLeftMargin(),0.82,caption3)
  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.88,caption4)

  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,caption1)
  saveAs(canvas,outDir+"eeOvermm_xsbr_8TeV")
