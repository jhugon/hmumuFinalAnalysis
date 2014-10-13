#!/usr/bin/env python

import singleHelpers
from helpers import *
import ROOT as root
import glob
import re
import os.path
import cPickle

from xsec import *

from singleUseScripts.biasPklToMu import getSMSigCounts

#######################################

PRELIMINARYSTRING="CMS"

#~48 Charactars Max

titleMap = {
  "AllCat":"All Categories Comb.",
  "IncCat":"Non-VBF Categories Comb.",
  "VBFCat":"VBF Categories Comb.",

  "IncPresel":"Non-VBF Preselection",
  "VBFPresel":"VBF Preselection",

  "Pt0to30":"p_{T}^{#mu#mu} #in [0,30]",
  "Pt30to50":"p_{T}^{#mu#mu} #in [30,50]",
  "Pt50to125":"p_{T}^{#mu#mu} #in [50,125]",
  "Pt125to250":"p_{T}^{#mu#mu} #in [125,250]",
  "Pt250":"p_{T}^{#mu#mu}>250",

  "VBFLoose":"VBFL",
  "VBFMedium":"VBFM",
  "VBFTight":"VBFT",
  "VBFVeryTight":"VBFVT",

  "BDTCut":"BDT Cut Combination",
  "IncBDTCut":"Non-VBF BDT Cut",
  "VBFBDTCut":"VBF BDT Cut",

  "BDTCutCat":"BDT Cut Cat. Combination",
  "IncBDTCutCat":"Non-VBF BDT Cut",
  "VBFBDTCutCat":"VBF BDT Cut",

  "IncPreselCat":"Non-VBF Cat. Preselection",
  "VBFPreselCat":"VBF Cat. Preselection",

  "IncBDTCutBB":"Non-VBF BDT Cut BB",
  "IncBDTCutBO":"Non-VBF BDT Cut BO",
  "IncBDTCutBE":"Non-VBF BDT Cut BE",
  "IncBDTCutOO":"Non-VBF BDT Cut OO",
  "IncBDTCutOE":"Non-VBF BDT Cut OE",
  "IncBDTCutEE":"Non-VBF BDT Cut EE",
  "IncBDTCutNotBB":"Non-VBF BDT Cut !BB",
  "VBFBDTCutBB":"VBF BDT Cut BB",
  "VBFBDTCutNotBB":"VBF BDT Cut !BB",
  "IncPreselBB":"Non-VBF Preselection BB",
  "IncPreselBO":"Non-VBF Preselection BO",
  "IncPreselBE":"Non-VBF Preselection BE",
  "IncPreselOO":"Non-VBF Preselection OO",
  "IncPreselOE":"Non-VBF Preselection OE",
  "IncPreselEE":"Non-VBF Preselection EE",
  "IncPreselNotBB":"Non-VBF Preselection !BB",
  "VBFPreselBB":"VBF Preselection BB",
  "VBFPreselNotBB":"VBF Preselection !BB",

  "IncPreselPtG10BB":"Non-VBF BB",
  "IncPreselPtG10BO":"Non-VBF BO",
  "IncPreselPtG10BE":"Non-VBF BE",
  "IncPreselPtG10OO":"Non-VBF OO",
  "IncPreselPtG10OE":"Non-VBF OE",
  "IncPreselPtG10EE":"Non-VBF EE",
  "IncPreselPtG10NotBB":"Non-VBF !BB",

  "IncPreselPtG":"Non-VBF Not Combined",

  "Jets01PassPtG10BB": "0,1-Jet Tight BB",
  "Jets01PassPtG10BO": "0,1-Jet Tight BO",
  "Jets01PassPtG10BE": "0,1-Jet Tight BE",
  "Jets01PassPtG10OO": "0,1-Jet Tight OO",
  "Jets01PassPtG10OE": "0,1-Jet Tight OE",
  "Jets01PassPtG10EE": "0,1-Jet Tight EE",
  "Jets01PassCatAll" : "0,1-Jet Tight Combination",
                        
  "Jets01FailPtG10BB": "0,1-Jet Loose BB",
  "Jets01FailPtG10BO": "0,1-Jet Loose BO",
  "Jets01FailPtG10BE": "0,1-Jet Loose BE",
  "Jets01FailPtG10OO": "0,1-Jet Loose OO",
  "Jets01FailPtG10OE": "0,1-Jet Loose OE",
  "Jets01FailPtG10EE": "0,1-Jet Loose EE",
  "Jets01FailCatAll" : "0,1-Jet Loose Combination",
                        
  "Jets01SplitCatAll": "0,1-Jet Combination",


  "Jet2CutsVBFPass":"2-Jet VBF Tight",
  "Jet2CutsGFPass":"2-Jet GF Tight",
  "Jet2CutsFailVBFGF":"2-Jet Loose",

  "Jet2SplitCutsGFSplit" : "2-Jet Combination",
  "CombSplitAll" : "Combination",
}

#######################################


def getDataMu(fileString,matchString=r"_([-\d.]+)\.txt\.mu",dontMatchStrings=[],doSort=True,xMax=1e20,xMin=-1e20):
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
    low1sig = -10.0
    high1sig = -10.0
    xNum = -10.0
    if match:
      xNum = match.group(1)
      if float(xNum) > xMax:
        continue
      if float(xNum) < xMin:
        continue
    for line in tmpF:
      #obsMatch = re.search(r"Best fit r:[\s]+([.\deE-]+)[\s]+([.\deE-]+)/+([.\deE-]+)",line)
      #obsMatch = re.search(r"Best fit r:[\s]+([.\deE-]+)[\s]+\-([.\deE-]+)/\+([.\deE-]+)",line)
      obsMatch = re.search(r"r[\s]*:[\s]*([.\deE+-]+)[\s]+\-([.\deE-]+)/\+([.\deE-]+)",line)
      if obsMatch:
        obs = obsMatch.group(1)
        low1sig = obsMatch.group(2)
        high1sig = obsMatch.group(3)
        #print("original line: {0}\n  Extracted: {1} -{2} +{3}".format(obsMatch.group(0),obs,low1sig,high1sig))
    thisPoint = [xNum,obs,low1sig,high1sig]
    if thisPoint.count("-10.0")>0:
        continue
    if thisPoint.count(-10.0)>0:
        continue
    #print thisPoint
    result.append(thisPoint)
  return result


class RelativePlot:
  def __init__(self,dataPoints, canvas, legend, caption, ylabel="Error on #sigma/#sigma_{SM}", xlabel="Integrated Luminosity [fb^{-1}]",caption2="",caption3="",ylimits=[],xlimits=[],vertLines=[],showObs=False,energyStr="8TeV",doMuExtraPlot=False):
    errGraph = root.TGraph()
    #errGraph.SetLineStyle(2)
    errGraph.SetLineColor(root.kRed)
    errGraph.SetMarkerStyle(20)
    errGraph.SetMarkerSize(1.1)
    errGraph.SetMarkerColor(root.kRed)
    self.errGraph = errGraph
    obsGraph = root.TGraph()
    obsGraph.SetLineColor(root.kRed)
    obsGraph.SetLineStyle(2)
    obsGraph.SetMarkerStyle(20)
    obsGraph.SetMarkerSize(1.1)
    obsGraph.SetMarkerColor(root.kRed)
    self.obsGraph = obsGraph
    muExtraGraphErr = root.TGraphAsymmErrors()
    muExtraGraph = root.TGraph()
    muExtraGraph.SetLineColor(root.kBlack)
    muExtraGraph.SetMarkerStyle(20)
    muExtraGraph.SetMarkerSize(1.0)
    muExtraGraph.SetMarkerColor(root.kBlack)
    muExtraGraphErr.SetFillColor(root.kGreen)
    self.muExtraGraph = muExtraGraph
    self.muExtraGraphErr = muExtraGraphErr
    iPoint = 0
    ymax = 1e-20
    ymin = 1e20
    for point in dataPoints:
      value = None
      obs = 0.0
      xNum = float(point[0])
      if len(xlimits)==2:
        if xNum<xlimits[0]:
            continue
        if xNum>xlimits[1]:
            continue
      if len(point)==4: #mu
        #thisPoint = [xNum,obs,low1sig,high1sig]
        high1sig = float(point[3])
        low1sig = float(point[2])
        if high1sig > low1sig:
            value = high1sig
        else:
            value = low1sig
      else: #sig len=3
        #thisPoint = [xNum,obs,exp]
        value = float(point[2])
        obs = float(point[1])
        if value == -20.0:
            value = obs
      errGraph.SetPoint(iPoint,xNum,value)
      obsGraph.SetPoint(iPoint,xNum,obs)
      if doMuExtraPlot:
        muExtraGraph.SetPoint(iPoint,xNum,float(point[1]))
        muExtraGraphErr.SetPoint(iPoint,xNum,float(point[1]))
        muExtraGraphErr.SetPointError(iPoint,0.,0.,float(point[2]),float(point[3]))
        #canvas.SetGrid(0,1)
      iPoint += 1

    self.vertLines = []
    self.vertLabel = []
    for xPos in vertLines:
      tmp = root.TGraph()
      tmp.SetPoint(0,xPos,ymin)
      tmp.SetPoint(1,xPos,ymax)
      tmp.SetLineColor(root.kRed)
      tmp.SetLineStyle(3)
      self.vertLines.append(tmp)
      self.vertLabel.append(xPos)

    label = root.TLatex()
    #label.SetNDC()
    label.SetTextFont(root.gStyle.GetLabelFont("X"))
    label.SetTextSize(root.gStyle.GetLabelSize("X"))
    label.SetTextColor(root.kRed)
    label.SetTextAlign(22)
    self.label=label

    if doMuExtraPlot:
      muExtraGraphErr.GetXaxis().SetTitle(xlabel)
      muExtraGraphErr.GetYaxis().SetTitle(ylabel)
      if len(ylimits)==2:
        muExtraGraphErr.GetYaxis().SetRangeUser(*ylimits)
      muExtraGraphErr.Draw("a3")
      muExtraGraph.Draw("l")
      muExtraGraph.Draw("p")
      #muExtraGraphErr.Draw("a4")
      #muExtraGraph.Draw("c")
    else:
      errGraph.GetXaxis().SetTitle(xlabel)
      errGraph.GetYaxis().SetTitle(ylabel)
      if showObs:
        if len(ylimits)==2:
          errGraph.GetYaxis().SetRangeUser(*ylimits)
        else:
          maxVal = 0.0
          tmpY = root.Double()
          tmpX = root.Double()
          for i in range(obsGraph.GetN()):
            obsGraph.GetPoint(i,tmpX,tmpY)
            if maxVal < float(tmpY):
              maxVal = float(tmpY)
          for i in range(errGraph.GetN()):
            errGraph.GetPoint(i,tmpX,tmpY)
            if maxVal < float(tmpY):
              maxVal = float(tmpY)
          errGraph.SetMaximum(maxVal)
        errGraph.Draw("al")
        obsGraph.Draw("l")
        obsGraph.Draw("p")
      else:
        if len(ylimits)==2:
          errGraph.GetYaxis().SetRangeUser(*ylimits)
        errGraph.Draw("al")
        errGraph.Draw("p")

    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(62)
    tlatex.SetTextSize(0.05)
    tlatex.SetTextAlign(11)
    tlatex.DrawLatex(0.15,0.94,PRELIMINARYSTRING)
    tlatex.SetTextFont(42)
    tlatex.SetTextSize(0.04)
    tlatex.DrawLatex(0.27,0.94,"H #rightarrow #mu^{+}#mu^{-}")
    tlatex.SetTextAlign(31)
    tlatex.SetTextSize(0.035)
    tlatex.DrawLatex(0.95,0.94,caption2)

    tlatex.SetTextSize(0.04)
    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(gStyle.GetPadLeftMargin()+0.03,0.87,caption)
    tlatex.DrawLatex(gStyle.GetPadLeftMargin()+0.03,0.82,caption3)

    for g in self.vertLines:
      g.Draw("l")

    canvas.RedrawAxis()

class PValuePlotTogether:
  def __init__(self,dataDict, canvas, caption="Standard Model H #rightarrow #mu^{+}#mu^{-}", ylabel="p-Value", xlabel="m_{H} [GeV/c^{2}]",caption2="",caption3="",ylimits=[],xlimits=[],energyStr="8TeV"):
    graphs = []
    ymax = 1.0
    ymin = 1e20
    ymin = 2e-3 # To see 3sigma line
    xmin = 1e20
    xmax = -1e20
    sortedChannels = sorted(dataDict.keys())
    sortedChannels.reverse()
    for channel in sortedChannels:
      dataPoints = dataDict[channel]
      graph = root.TGraph()
      graph.SetName("Pvalue_"+energyStr+"_"+channel)
      #graph.SetLineStyle(2)
      graph.SetLineColor(colorMap[channel])
      graph.SetMarkerStyle(20)
      graph.SetMarkerSize(1.1)
      graph.SetMarkerColor(colorMap[channel])
      iPoint = 0
      for point in dataPoints:
        value = None
        obs = 0.0
        xNum = float(point[0])
        if len(xlimits)==2:
          if xNum<xlimits[0]:
              continue
          if xNum>xlimits[1]:
              continue
        if xNum < xmin:
            xmin = xNum
        if xNum > xmax:
            xmax = xNum
        #thisPoint = [xNum,obs,exp]
        obs = float(point[1])
        graph.SetPoint(iPoint,xNum,obs)
        iPoint += 1
        if obs < ymin:
          ymin = obs
      graphs += [graph]
    ymin *=0.5
    for graph in graphs:
      if len(ylimits)==2:
        graph.GetYaxis().SetRangeUser(*ylimits)
      else:
        graph.GetYaxis().SetRangeUser(ymin*0.5,ymax)

    self.hLine = root.TLine()
    self.hLine.SetLineColor(root.kBlack)
    self.hLine.SetLineWidth(3)
    self.hLine.SetLineStyle(3)
    sigmaVals = [0.8413,0.9772,0.9987]
    sigmaVals = [1.0-x for x in sigmaVals]
    hLabelX = None

    label = root.TLatex()
    label.SetTextFont(root.gStyle.GetLabelFont("X"))
    label.SetTextSize(root.gStyle.GetLabelSize("X"))
    label.SetTextAlign(11)
    label.SetTextColor(1)
    self.label=label

    drawn = False
    for graph in graphs:
      if drawn == False:
        graph.Draw("al")
        drawn = True
        graph.GetXaxis().SetTitle(xlabel)
        graph.GetYaxis().SetTitle(ylabel)
        graph.Draw("al")
        xmin = graph.GetXaxis().GetXmin()
        xmax = graph.GetXaxis().GetXmax()
        hLabelX = (xmax-xmin)/0.02+xmin
        for yPos in sigmaVals:
         if yPos < ymin:
           break
         self.hLine.DrawLine(xmin,yPos,xmax,yPos)
        graph.Draw("l")
      else:
        graph.Draw("l")
      graph.Draw("p")

    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(62)
    tlatex.SetTextSize(0.05)
    tlatex.SetTextAlign(11)
    tlatex.DrawLatex(0.15,0.94,PRELIMINARYSTRING)
    tlatex.SetTextFont(42)
    tlatex.SetTextSize(0.04)
    tlatex.DrawLatex(0.27,0.94,"H #rightarrow #mu^{+}#mu^{-}")
    tlatex.SetTextFont(root.gStyle.GetLabelFont())
    tlatex.SetTextSize(0.035)
    tlatex.SetTextAlign(31)
    tlatex.DrawLatex(0.95,0.94,caption)

    tlatex.SetTextSize(0.04)
    tlatex.SetTextAlign(12)
    caption23 = caption2
    if caption3 != "":
      caption23 += ", "+caption3
    tlatex.DrawLatex(gStyle.GetPadLeftMargin()+0.03,0.88,caption23)

    #tlatex.SetTextAlign(12)
    #tlatex.SetTextSize(0.05)
    #tlatex.DrawLatex(gStyle.GetPadLeftMargin()+0.03,gStyle.GetPadBottomMargin()+0.23,"p-Value of #sigma/#sigma_{SM}=1")
    #tlatex.DrawLatex(gStyle.GetPadLeftMargin()+0.03,gStyle.GetPadBottomMargin()+0.17,"Hypothesis")

    #legend = root.TLegend(0.58,0.70,0.9,0.9) # UR
    baselineY = gStyle.GetPadBottomMargin()+0.02
    marginR = 1.0-gStyle.GetPadRightMargin()+0.01
    marginL = gStyle.GetPadLeftMargin()+0.01
    #legend = root.TLegend(marginR-0.28,baselineY,marginR,baselineY+0.4)
    #legend = root.TLegend(marginL+0.1,baselineY,marginL+0.1+0.45,baselineY+0.37)
    legend = root.TLegend(marginL+0.1,baselineY,marginL+0.1+0.45,baselineY+0.32)
    if len(graphs) < 5:
      #legend = root.TLegend(marginR-0.32,baselineY,marginR,baselineY+0.25)
      legend = root.TLegend(marginL+0.1,baselineY,marginL+0.1+0.35,baselineY+0.27)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    self.legend = legend
    for channel,graph in zip(reversed(sortedChannels),reversed(graphs)):
      entry = titleMap[channel]
      if channel == "BDTCutCatVBFBDTOnly":
        entry = "VBF & Non-VBF"
      legend.AddEntry(graph,entry,"lp")

    legend.Draw()

    for yPos,iSigma in zip(sigmaVals,range(len(sigmaVals))):
      if yPos < ymin:
        break
      label.DrawLatex(xmin+0.05*(xmax-xmin),yPos*1.05,"{0}#sigma".format(iSigma+1))

    canvas.RedrawAxis()
    self.graphs = graphs

if __name__ == "__main__":

  dirName = "statsCards/"
  #dirName = "/data/uftrig01b/digiovan/baselinePP/exppluspolin_fixEff_DiffMeans_unbinned_loosePUID_fixBkgVBFFit_syst_22Jan2013ReReco_assProd/statsInputUltimate/"
  #dirName = "/data/uftrig01b/kropiv/HiggsMuMu/CMSSW_6_1_1/Results/hmumuFinalAnalysis/statsInput/"
  #dirName = "/data/uftrig01b/digiovan/baselinePP/m110to160_pixelLumi/hmumuFinalAnalysis/statsInput/"
  #dirName = "/raid/raid7/jhugon/higgsDataCards/20140501/statsInput/"
  
  root.gROOT.SetBatch(True)
  
  for period in ["7TeV","8TeV","14TeV","7P8TeV"]:
    fnToGlob = dirName+"*_"+period+"_*.txt.mu2"
    allfiles = glob.glob(fnToGlob)

    energyStr = ""
    plots = set()
    for fn in allfiles:
      match = re.search(r".*/(.+)_(.+)_[.\d]+.txt.mu2",fn)
      badPlot = re.search(r"Silly",fn)
      badPlot2 = re.search(r"Silly",fn)
      if match and not (badPlot or badPlot2):
        plots.add(match.group(1))
        energyStr = match.group(2)
    plots = list(plots)
    plots = sortCatNames(plots)
    #print plots
  
    print "="*40
    print period
    for plotName in plots:
      data = None
      doMuExtraPlot=True
      showObs=False
      
      data = getDataMu(dirName+plotName+"_"+period+"_125.0.txt.mu*")
      data = [float(i) for i in data[0]]
      data.append(data[1]-data[3])
      data.append(data[1]+data[2])
      #print "{0:26}".format(titleMap[plotName])+": {1:6.1f} +{2:<6.1f} -{3:<6.1f}".format(*data)
      #print "{0:26}".format(titleMap[plotName])+": {1:6.1f} +{2:<6.1f} -{3:<6.1f}  68% region: [{4:6.1f},{5:6.1f}]".format(*data)
      #print "  {0:25}".format('"'+plotName+'"')+": [{4:6.1f},{5:6.1f}], # {1:6.1f}".format(*data)

      newData = [data[1],max(data[2],data[3])]
      smCounts = getSMSigCounts(plotName,125)
      newData = [ i*smCounts for i in newData ]
      print "{0:26}".format(titleMap[plotName])+": {0:6.1f} +/- {1:6.1f}".format(*newData)
