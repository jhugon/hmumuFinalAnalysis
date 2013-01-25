#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description="Makes Limit Plots from output text from combine tool.")
parser.add_argument("--bdtCut", help="Makes plots v. BDT Cut Instead of Luminosity",action="store_true",default=False)
args = parser.parse_args()


from helpers import *
import ROOT as root
import glob
import re

import matplotlib.pyplot as mpl
import numpy

from xsec import *

#######################################

#~48 Charactars Max
titleMap = {
  "AllCat":"H->#mu#mu Catagories Combination",
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

  "BDTCut":"H->#mu#mu BDT Combination",
  "IncBDTCut":"Non-VBF BDT ",
  "VBFBDTCut":"VBF BDT ",

  "BDTCutCat":"BDT Res. Cat. Combination",
  "IncBDTCutCat":"Non-VBF BDT Resolution Categories",
  "VBFBDTCutCat":"VBF BDT Resolution Categories",

  "PreselCat":"Res. Cat. Preselection Combination",
  "IncPreselCat":"Non-VBF Resolution Cat. Preselection",
  "VBFPreselCat":"VBF Cat. Resolution Preselection",

  "IncBDTCutBB":"Non-VBF BDT BB",
  "IncBDTCutBO":"Non-VBF BDT BO",
  "IncBDTCutBE":"Non-VBF BDT BE",
  "IncBDTCutOO":"Non-VBF BDT OO",
  "IncBDTCutOE":"Non-VBF BDT OE",
  "IncBDTCutEE":"Non-VBF BDT EE",
  "IncBDTCutNotBB":"Non-VBF BDT !BB",
  "VBFBDTCutBB":"VBF BDT BB",
  "VBFBDTCutNotBB":"VBF BDT !BB",
  "IncPreselBB":"Non-VBF Preselection BB",
  "IncPreselBO":"Non-VBF Preselection BO",
  "IncPreselBE":"Non-VBF Preselection BE",
  "IncPreselOO":"Non-VBF Preselection OO",
  "IncPreselOE":"Non-VBF Preselection OE",
  "IncPreselEE":"Non-VBF Preselection EE",
  "IncPreselNotBB":"Non-VBF Preselection !BB",
  "VBFPreselBB":"VBF Preselection BB",
  "VBFPreselNotBB":"VBF Preselection !BB"
}

comparisonMap = {
  "AllCat":"All Cat. Comb.",
  "IncCat":"!VBF Cat. Comb.",
  "VBFCat":"VBF Cat. Comb.",

  "Presel":"Presel. Comb.",
  "IncPresel":"!VBF Presel.",
  "VBFPresel":"VBF Presel.",

  "Pt0to30":"$p_{T}^{\mu\mu} \in [0,30]$",
  "Pt30to50":"$p_{T}^{\mu\mu} \in [30,50]$",
  "Pt50to125":"$p_{T}^{\mu\mu} \in [50,125]$",
  "Pt125to250":"$p_{T}^{\mu\mu} \in [125,250]$",
  "Pt250":"$p_{T}^{\mu\mu}>250$",

  "VBFLoose":"VBFL",
  "VBFMedium":"VBFM",
  "VBFTight":"VBFT",
  "VBFVeryTight":"VBFVT",

  "BDTCut":"BDT Comb.",
  "IncBDTCut":"!VBF BDT",
  "VBFBDTCut":"VBF BDT",

  "BDTCutCat":"BDT Res. Comb.",
  "IncBDTCutCat":"!VBF BDT Res.",
  "VBFBDTCutCat":"VBF BDT Res.",

  "BDTCutCatIncOnly":"BDT Comb. Inc Res.",

  #"BDTCutCat":"Combination",
  #"IncBDTCutCat":"Non-VBF",
  #"VBFBDTCutCat":"VBF",

  "PreselCat":"Presel. Res. Comb.",
  "IncPreselCat":"!VBF Res. Presel.",
  "VBFPreselCat":"VBF Res. Presel.",

  "IncBDTCutBB":"!VBF BDT BB",
  "IncBDTCutBO":"!VBF BDT BO",
  "IncBDTCutBE":"!VBF BDT BE",
  "IncBDTCutOO":"!VBF BDT OO",
  "IncBDTCutOE":"!VBF BDT OE",
  "IncBDTCutEE":"!VBF BDT EE",
  "IncBDTCutNotBB":"!VBF BDT !BB",
  "VBFBDTCutBB":"VBF BDT BB",
  "VBFBDTCutNotBB":"VBF BDT !BB",
  "IncPreselBB":"!VBF Presel. BB",
  "IncPreselBO":"!VBF Presel. BO",
  "IncPreselBE":"!VBF Presel. BE",
  "IncPreselOO":"!VBF Presel. OO",
  "IncPreselOE":"!VBF Presel. OE",
  "IncPreselEE":"!VBF Presel. EE",
  "IncPreselNotBB":"!VBF Presel. !BB",
  "VBFPreselBB":"VBF Presel. BB",
  "VBFPreselNotBB":"VBF Presel. !BB"
}

#######################################

def getData(fileString,matchString=r"_([-\d.]+)\.txt\.out",dontMatchStrings=[],doSort=True):
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
    if thisPoint.count(-10.0)>0:
        continue
    #print thisPoint
    result.append(thisPoint)
  return result

class RelativePlot:
  def __init__(self,dataPoints, canvas, legend, caption, ylabel="Expected 95% CL Limit #sigma/#sigma_{SM}", xlabel="Integrated Luminosity [fb^{-1}]",caption2="",caption3="",ylimits=[],xlimits=[],vertLines=[],showObs=False,energyStr="8TeV"):
    expGraph = root.TGraph()
    expGraph.SetLineStyle(2)
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

    twoSigGraph.Draw("a3")
    twoSigGraph.GetXaxis().SetTitle(xlabel)
    twoSigGraph.GetYaxis().SetTitle(ylabel)
    if len(ylimits)==2:
        twoSigGraph.GetYaxis().SetRangeUser(*ylimits)
    oneSigGraph.Draw("3")
    expGraph.Draw("l")
    oneGraph.Draw("l")
    if showObs:
      obsGraph.Draw("l")

    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(root.gStyle.GetLabelFont())
    #tlatex.SetTextSize(0.05)
    tlatex.SetTextSize(0.04)
    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.88,caption2)
    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.82,caption3)

    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,caption)

    for g in self.vertLines:
      g.Draw("l")

    canvas.RedrawAxis()

class ComparePlot:
  def __init__(self,data,ylabel="Expected 95% CL Limit $\sigma/\sigma_{SM}$",titleMap={},showObs=False,xlimits=[]):
    fig = mpl.figure()
    self.fig = fig
    ax1 = fig.add_subplot(111)
    self.ax1 = ax1
    ax1bounds = ax1.get_position().bounds
    ax1.set_position([0.25,0.1,0.7,0.85])

    data.sort(key=lambda x: x[0].lower())

    obs = [float(point[1]) for point in data]
    medians = [float(point[4]) for point in data]
    high1sigs = [float(point[5])-float(point[4]) for point in data]
    low1sigs = [float(point[4])-float(point[3]) for point in data]
    high2sigs = [float(point[6])-float(point[4]) for point in data]
    low2sigs = [float(point[4])-float(point[2]) for point in data]

    xPos = numpy.arange(len(medians))
    xLabels = [point[0] for point in data]
    xLabels = [re.sub(r".*/","",s) for s in xLabels]
    xLabels = [titleMap[s] if titleMap.has_key(s) else s for s in xLabels]

    ax1.set_yticks(xPos+0.25)
    ax1.set_yticklabels(tuple(xLabels),size="small")
    ax1.set_xlabel(ylabel)
    if len(xlimits)==2:
      ax1.set_xlim(*xlimits)
    bars = ax1.barh(xPos,medians, 0.5, xerr=[low2sigs,high2sigs],ecolor="k")
    self.bars = bars
    xPosObs = [x+0.25 for x in xPos]
    # Now for the 1sigma error bars
    self.oneSig = ax1.errorbar(medians,xPosObs,xerr=[low1sigs,high1sigs],color="g",linestyle="None")
    if showObs:
      #self.obs = ax1.plot(obs,xPosObs,marker="|",color="r",markersize=42,linestyle="None",markeredgewidth=2)
      self.obs = ax1.plot(obs,xPosObs,marker="o",color="r",markersize=10,linestyle="None",markeredgecolor='r')
    writeInValues = getattr(self,"writeInValues")
    writeInValues(bars,high2sigs)
    if showObs:
      writeInValues(bars,high2sigs,obs=obs)
  def save(self,saveName):
    self.fig.savefig(saveName+".png")
    self.fig.savefig(saveName+".pdf")
  def writeInValues(self,rects,high2sigs,obs=None):
    axisWidth = self.ax1.get_xlim()[1]-self.ax1.get_xlim()[0]
    sizelen = len(rects)
    if obs != None:
      sizelen -=5
    size="medium"
    if sizelen<5:
      size="large"
    elif sizelen>12:
      size="xx-small"
    elif sizelen>9:
      size="x-small"
    maxWidth = axisWidth
    for rect,sigUp,i in zip(rects,high2sigs,range(len(rects))):
      width = rect.get_width()
      if (width/maxWidth < 0.1): # The bars aren't wide enough to print the ranking inside
        xloc = width + sigUp + maxWidth*0.01 # Shift the text to the right side of the right edge
        clr = 'black' # Black against white background
        align = 'left'
        if obs != None:
          clr = 'red' # Black against white background
          align = 'right'
          xloc += maxWidth*0.04
      else:
        xloc = width - maxWidth*0.01 # Shift the text to the left side of the right edge
        clr = 'white' # White on magenta
        align = 'right'
        if obs != None:
          clr = 'red' # Black against white background
          align = 'left'
          xloc += maxWidth*0.04
      yloc = rect.get_y()+rect.get_height()/2.0 
      valueStr = "{0:.1f}".format(width)
      if obs != None:
        valueStr = "{0:.1f}".format(obs[i])
        if i == len(rects)-1:
          yloc -= 0.15
        else:
          yloc += 0.15
      self.ax1.text(xloc, yloc, valueStr, horizontalalignment=align,
            verticalalignment='center', color=clr, weight='bold',size=size)

class ComparePlotTable:
  def __init__(self,data,ylabel="Expected 95% CL Limit $\sigma/\sigma_{SM}$",titleMap={},xlimits=[]):
    fig = mpl.figure(figsize=(16,8))
    self.fig = fig
    ax1 = fig.add_subplot(111)
    self.ax1 = ax1
    ax1bounds = ax1.get_position().bounds
    ax1bounds = list(ax1bounds)
    self.plotxLimit = 0.6
    ax1bounds[2] *= self.plotxLimit
    ax1.set_position(ax1bounds)

    data.sort(key=lambda x: x[0].lower())

    obs = [float(point[1]) for point in data]
    medians = [float(point[4]) for point in data]
    high1sigs = [float(point[5])-float(point[4]) for point in data]
    low1sigs = [float(point[4])-float(point[3]) for point in data]
    high2sigs = [float(point[6])-float(point[4]) for point in data]
    low2sigs = [float(point[4])-float(point[2]) for point in data]

    xPos = numpy.arange(len(medians))
    self.xPos = xPos
    xLabels = [point[0] for point in data]
    xLabels = [re.sub(r".*/","",s) for s in xLabels]
    xLabels = [titleMap[s] if titleMap.has_key(s) else s for s in xLabels]

    ax1.set_yticks(xPos+0.25)
    ax1.set_yticklabels(tuple(xLabels),size="small")
    ax1.set_xlabel(ylabel)
    if len(xlimits)==2:
      ax1.set_xlim(*xlimits)
    bars = ax1.barh(xPos,medians, 0.5, xerr=[low2sigs,high2sigs],ecolor="k")
    self.bars = bars
    xPosObs = [x+0.25 for x in xPos]
    # Now for the 1sigma error bars
    self.oneSig = ax1.errorbar(medians,xPosObs,xerr=[low1sigs,high1sigs],color="g",linestyle="None")
    self.obs = ax1.plot(obs,xPosObs,marker="o",color="r",markersize=25,linestyle="None")
  def writeTable(self):
    for i in self.xPos:
      self.fig.text(self.plotxLimit,1.0/len(self.xPos)*i,"Red Lines: Observed Limit",horizontalalignment="right",size="medium",color="r")
  def save(self,saveName):
    self.fig.savefig(saveName+".png")
    self.fig.savefig(saveName+".pdf")

if __name__ == "__main__":

  dirName = "statsInput/"
  
  outDir = "statsOutput/"
  
  root.gErrorIgnoreLevel = root.kWarning
  root.gROOT.SetBatch(True)
  setStyle()
  canvas = root.TCanvas()
  if not args.bdtCut:
    canvas.SetLogx(1)
    canvas.SetLogy(1)
  
  mpl.rcParams["font.family"] = "sans-serif"
  #print mpl.rcParams["backend"]

  ylimits=[0.1,100.0]
  if args.bdtCut:
    ylimits=[1.,30.0]
    ylimits=[1.,70.0]

  lumisToUse={"7TeV":lumiDict["7TeV"],"8TeV":lumiDict["8TeV"],"7P8TeV":lumiDict["8TeV"]+lumiDict["7TeV"]}
  
  for period in ["7TeV","8TeV","14TeV","7P8TeV"]:
    fnToGlob = dirName+"*_"+period+"_*.txt.out"
    allfiles = glob.glob(fnToGlob)
    
    ## Limit v. Lumi
    energyStr = ""
    plots = set()
    for fn in allfiles:
      match = re.search(r".*/(.+)_(.+)_[-.\d]+\.txt\.out",fn)
      badPlot = re.search(r"Silly",fn)
      badPlot2 = re.search(r"Silly",fn)
      if match and not (badPlot or badPlot2):
        plots.add(match.group(1))
        energyStr = match.group(2)
  
    caption2 = "#sqrt{s}="+energyStr
    legend = root.TLegend(0.58,0.70,0.9,0.9)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    for plotName in plots:
      data = getData(dirName+plotName+"_"+energyStr+"_*.txt.out")
      if len(data)<=1:
        continue
      xlabel="Integrated Luminosity [fb^{-1}]"
      caption3 = ""
      if args.bdtCut:
        xlabel="BDT Discriminant Cut"
        match = re.match(r"([^0-9.]*)([0-9.]*)",plotName)
        assert(match)
        caption3 = "L = {0:.1f} fb^{{-1}}".format(float(match.group(2)))
        plotName = match.group(1)
      #elif period == "14TeV":
      #  title = "Standard Model H#rightarrow#mu#mu"
      title = titleMap[plotName]
      incPlot = RelativePlot(data,canvas,legend,title,caption2=caption2,ylimits=ylimits,energyStr=energyStr,xlabel=xlabel,caption3=caption3)
      saveAs(canvas,outDir+plotName+"_"+energyStr)

    if args.bdtCut:
        continue

    ## Compare all types of limits
    if period == "14TeV":
        continue

    veto = []
    #veto = [r"CNC",r"PM","BB","BO","BE","OO","OE","EE","NotBB"]
    #veto = [r"CNC",r"PM","Presel","BB","BO","BE","OO","OE","EE","NotBB"]
    #veto = [r"CNC",r"PM","BDT","BB","BO","BE","OO","OE","EE","NotBB"]
    #veto = [r"CNC",r"PM","Presel"]
    #veto = [r"CNC",r"PM","BDT"]

    #veto = ["Cat","Presel"]
    #veto = ["Cat","Comb","Presel"]
    #veto = ["Comb","Presel"]
    veto = ["Presel"]
    veto = ["Presel","BB","BO","BE","OO","OE","EE"]

    mustBe = r"(.+)_(.+)_[.\d]+.txt.out"
    #mustBe = r"(.+Cat)_(.+)_[.\d]+.txt.out"
    
    desiredLumiStr=str(lumisToUse[period])
    fnGlobStr = dirName+"*_"+energyStr+"_"+desiredLumiStr+".txt.out"
    compareData = getData(fnGlobStr,matchString=mustBe,dontMatchStrings=veto,doSort=False)
    energyStrWrite = None
    compareXlims = [0.,25.]
    if energyStr == "7TeV":
      compareXlims = [0.,65.]
    if energyStr == "7P8TeV":
      energyStrWrite = "7 & 8 TeV"
    else:
      energyStrWrite = energyStr.replace("TeV"," TeV")
    #print compareData
    if len(compareData)==0:
        print("No Data to Compare for {0}!!".format(period))
        continue
    lumiStrWrite = "{0:.1f}".format(float(desiredLumiStr))
    comparePlot = ComparePlot(compareData,titleMap=comparisonMap,showObs=True,xlimits=compareXlims)
    comparePlot.fig.text(0.9,0.2,"$\mathcal{L}="+lumiStrWrite+"$ fb$^{-1}$",horizontalalignment="right",size="x-large")
    comparePlot.fig.text(0.9,0.27,"$\sqrt{s}=$"+energyStrWrite,horizontalalignment="right",size="x-large")
    #comparePlot.fig.text(0.9,0.13,"Red Lines: Observed Limit",horizontalalignment="right",size="medium",color="r")
    comparePlot.fig.text(0.9,0.425,"Red: Observed Limit",horizontalalignment="right",size="large",color="r")
    comparePlot.save(outDir+"compareObs"+"_"+energyStr)
    
    #comparePlot = ComparePlot(compareData,titleMap=comparisonMap,showObs=False,xlimits=compareXlims)
    #comparePlot.fig.text(0.9,0.2,"$\mathcal{L}="+desiredLumiStr+"$ fb$^{-1}$",horizontalalignment="right",size="x-large")
    #comparePlot.fig.text(0.9,0.27,"$\sqrt{s}=$"+energyStr,horizontalalignment="right",size="x-large")
    #comparePlot.fig.text(0.9,0.50,"Calculated from MC",horizontalalignment="right",size="x-large")
    #comparePlot.save(outDir+"compare"+"_"+energyStr)

    #comparePlot = ComparePlotTable(compareData,titleMap=comparisonMap,xlimits=compareXlims)
    #comparePlot.fig.text(0.9,0.2,"$\mathcal{L}="+lumiStrWrite+"$ fb$^{-1}$",horizontalalignment="right",size="x-large")
    #comparePlot.fig.text(0.9,0.27,"$\sqrt{s}=$"+energyStrWrite,horizontalalignment="right",size="x-large")
    #comparePlot.fig.text(0.9,0.13,"Red Lines: Observed Limit",horizontalalignment="right",size="medium",color="r")
    #comparePlot.fig.text(0.9,0.425,"Red Lines: Observed Limit",horizontalalignment="right",size="medium",color="r")
    #comparePlot.save(outDir+"compareObs"+"_"+energyStr)
