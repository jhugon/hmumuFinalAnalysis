#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description="Makes Limit Plots from output text from combine tool.")
parser.add_argument("--bdtCut", help="Makes plots v. BDT Cut Instead of Luminosity",action="store_true",default=False)
parser.add_argument("-m","--higgsMass", help="Makes plots v. Higgs Mass",action="store_true",default=False)
parser.add_argument("--signalInject", help="Sets a caption saying that signal was injected with strength",type=float,default=0.0)
parser.add_argument("--signalInjectMass", help="Mass For Injected Signal",type=float,default=125.0)
args = parser.parse_args()


from helpers import *
import ROOT as root
import glob
import re
import os.path

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
  "IncPreselCat":"Non-VBF",
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
  "VBFPreselNotBB":"VBF Preselection !BB",

  "IncPtCut":"Non-VBF",

  "IncPreselPtG10BB":"Non-VBF BB",
  "IncPreselPtG10BO":"Non-VBF BO",
  "IncPreselPtG10BE":"Non-VBF BE",
  "IncPreselPtG10OO":"Non-VBF OO",
  "IncPreselPtG10OE":"Non-VBF OE",
  "IncPreselPtG10EE":"Non-VBF EE",
  "IncPreselPtG10":"Non-VBF",
  "BDTCutCatVBFBDTOnly": "VBF & Non-VBF Combination"
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
  "VBFPreselNotBB":"VBF Presel. !BB",

  "BDTCutCatVBFBDTOnly":'Comb. VBFBDT IncPreRes',
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
    obsGraph.SetLineWidth(3)
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
  def __init__(self,data,ylabel="95% CL Limit $\sigma/\sigma_{SM}$",titleMap={},showObs=False,xlimits=[]):
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
  def __init__(self,data,ylabel="95% CL Limit $\sigma/\sigma_{SM}$",titleMap={},xlimits=[],brazil=True,obsCircles=True,vertLine1=True,toHighlight=[],anotation1='',anotation2=''):
    self.toHighlight = toHighlight
    #data.sort(key=lambda x: x[0].lower())
    data.sort(key=lambda x: float(x[4]))
    fig = mpl.figure(figsize=(16,8))
    self.axRightBound = 0.60
    self.axLeftBound = 0.132
    self.tabRightBound = 0.05
    self.obsColWidth = 0.25
    self.linewidth = 2.5
    mpl.rcParams["lines.linewidth"] = self.linewidth
    mpl.rcParams["axes.linewidth"] = self.linewidth
    mpl.rcParams["patch.linewidth"] = self.linewidth
    mpl.rcParams["xtick.major.width"] = self.linewidth
    mpl.rcParams["xtick.minor.width"] = self.linewidth
    mpl.rcParams["ytick.major.width"] = self.linewidth
    mpl.rcParams["ytick.minor.width"] = self.linewidth
    mpl.rcParams['font.size'] = 22.0
    self.legendFontSize = 24.0
    self.tickLabelFontSize = 20.0
    self.legXMin = 0.74
    self.legYMax = 0.2
    self.legYMin = 0.13
    #mpl.rcParams["font.weight"]= 700
    #mpl.rcParams["mathtext.fontset"]= 'custom'
    #mpl.rcParams["mathtext.default"]= 'bf'
    self.data = data
    self.fig = fig
    ax1 = fig.add_subplot(111)
    self.ax1 = ax1
    ax1bounds = ax1.get_position().bounds
    ax1bounds = list(ax1bounds)
    ax1bounds[0] = self.axLeftBound
    ax1bounds[2] = self.axRightBound -self.axLeftBound
    ax1.set_position(ax1bounds)
    self.maxX = ax1.transAxes.inverted().transform(
            fig.transFigure.transform((1.0-self.tabRightBound,0.))
            )[0]
    self.obsColWidthAx = ax1.transAxes.inverted().transform(
            fig.transFigure.transform((self.obsColWidth,0.))
            )[0] - ax1.transAxes.inverted().transform(
            fig.transFigure.transform((0.0,0.))
            )[0]

    obs = [float(point[1]) for point in data]
    medians = [float(point[4]) for point in data]
    high1sigs = [float(point[5])-float(point[4]) for point in data]
    low1sigs = [float(point[4])-float(point[3]) for point in data]
    high2sigs = [float(point[6])-float(point[4]) for point in data]
    low2sigs = [float(point[4])-float(point[2]) for point in data]
    self.obs = obs
    self.medians = medians
    self.high1sigs = high1sigs
    self.high2sigs = high2sigs
    self.low1sigs = low1sigs
    self.low2sigs = low2sigs
    maxDataPoint = max([float(point[6]) for point in data])
    xmax = maxDataPoint*1.1

    xPos = numpy.arange(len(medians))
    self.xPos = xPos
    xLabels = [point[0] for point in data]
    xLabels = [re.sub(r".*/","",s) for s in xLabels]
    xLabels = [titleMap[s] if titleMap.has_key(s) else s for s in xLabels]
    self.xLabels = xLabels

    ax1.set_yticks(xPos+0.25)
    ax1.set_xlabel(ylabel)
    ax1.set_xlim(0.0,xmax)
    ax1.set_ylim(-0.25,len(xLabels)-0.25)
    xPosObs = [x+0.25 for x in xPos]
    self.xPosObs = xPosObs
    getattr(self,'highlight')()
    if brazil:
      getattr(self,'writeBrazil')()
    else:
      getattr(self,'writeBars')()
    if obsCircles:
      getattr(self,'writeObsCircles')()
    getattr(self,'writeTable')()
    getattr(self,'writeLegend')()
    if vertLine1:
      ax1.axvline(1.0,color='k',ls='--')
    ax1.text(0.0,1.01,PRELIMINARYSTRING,ha="left",va="baseline",size="x-large",transform=ax1.transAxes)

    ax1.text(0.99,0.03,anotation2,ha="right",va="baseline",size=24,transform=ax1.transAxes)
    if len(medians)>5:
      anotation1 =  anotation1.replace("&",'')
      anotation1 =  anotation1.replace("and",'')
      anotation1 =  anotation1.replace("And",'')
      ax1.text(1.0,self.legYMax+0.05,anotation1,
                    va="baseline",ha='right',size=26.,color='k',
                    transform=ax1.transAxes)
    else:
      anotation1 =  anotation1.replace("\n",'')
      ax1.text(0.5,0.99,anotation1,
                    va="top",ha='center',size=20.,color='k',
                    transform=ax1.transAxes)
    #ax1.set_yticklabels(tuple(xLabels),size="medium")
    ax1.set_yticklabels(tuple(['' for i in xLabels]))
    getattr(self,'writeTickLabels')()
  def writeObsCircles(self):
    self.obsPoints = self.ax1.plot(self.obs,self.xPosObs,marker="o",color="r",markersize=15,linestyle="None",markeredgecolor='r')
  def writeObsLines(self):
    self.obsPoints = self.ax1.plot(self.obs,self.xPosObs,marker="|",color="r",markersize=49,markeredgewidth=4.,linestyle="None")
  def writeLegend(self):
    ymax = self.legYMax
    ymin = self.legYMin
    xmin = self.legXMin
    self.ax1.plot([xmin],[ymax],marker='o',color='r',markersize=10,linestyle="None",markeredgecolor='r',transform=self.ax1.transAxes)
    self.ax1.plot([xmin],[(ymax-ymax)/2.+ymin],marker='|',color='k',markersize=30,markeredgewidth=self.linewidth,linestyle="None",transform=self.ax1.transAxes)
    self.ax1.text(xmin+0.02,ymax,'Observed',ha='left',va='center',size=self.legendFontSize,transform=self.ax1.transAxes,
                                                color='r')
    self.ax1.text(xmin+0.02,(ymax-ymax)/2.+ymin,'Expected',ha='left',va='center',size=self.legendFontSize,transform=self.ax1.transAxes)
    
  def writeTable(self):
    dispToFig = self.fig.transFigure.inverted()
    axToDisp = self.ax1.transAxes
    hLineList = []
    nHLines = len(self.data)
    maxX = self.maxX
    for i in range(nHLines+1):
      yCord = float(i)/nHLines
      l = mpl.Line2D([1,maxX],[yCord,yCord],color='k')
      self.ax1.add_artist(l)
      l.set_transform(self.ax1.transAxes)
      l.set_clip_on(False)
      hLineList.append(l)
    vLineList = []
    vlineXCords = []
    columnXCords = []
    nVLines = 4
    columnWidth = (maxX-1.0)/(nVLines)
    for i in range(nVLines):
      vlineXCords.append(1.0+(i+1.0)*columnWidth)
      columnXCords.append(1.0+(i+0.5)*columnWidth)
    for xCord in vlineXCords:
      l = mpl.Line2D([xCord,xCord],[0,1],color='k')
      self.ax1.add_artist(l)
      l.set_transform(self.ax1.transAxes)
      l.set_clip_on(False)
      vLineList.append(l)
    reversedData = reversed(self.data)
    for i in range(nHLines):
      yCord = float(i+0.5)/nHLines
      median = float(self.data[i][4])
      for xCord,j in zip(columnXCords,range(nVLines)):
        s = None
        color = 'k'
        if j == 0:
          color = 'r'
          s = "{0:.1f}".format(float(self.data[i][j+1]))
        elif j == 1:
          s = "{0:.1f}".format(median)
        elif j == 2:
          s = "+{0:.1f}\n-{1:.1f}".format(float(self.data[i][5])-median,median-float(self.data[i][3]))
        elif j == 3:
          s = "+{0:.1f}\n-{1:.1f}".format(float(self.data[i][6])-median,median-float(self.data[i][2]))
        self.ax1.text(xCord,yCord,s,va="center",ha='center',size=25.,
            transform=axToDisp,
            color=color)
    colLabels = ["Observed","Expected",r"$1\sigma$",r"$2\sigma$"]
    colLabelColors = ['r','k','k','k','k','k']
    for xCord,lab,color in zip(columnXCords,colLabels,colLabelColors):
      yCord = 1.01
      self.ax1.text(xCord,yCord,lab,va="baseline",ha='center',size=20.,
            transform=axToDisp,
            color=color)
  def writeBars(self):
    self.bars = self.ax1.barh(self.xPos,self.medians, 0.5, 
                xerr=[self.low2sigs,self.high2sigs],ecolor="k")
    # Now for the 1sigma error bars
    self.oneSig = self.ax1.errorbar(self.medians,self.xPosObs,
        xerr=[self.low1sigs,self.high1sigs],color="g",linestyle="None")
  def writeBrazil(self):
    medianLines = []
    boxes1 = []
    boxes2 = []
    barHeight = 0.5
    for point,yLow in zip(self.data,range(len(self.data))):
      xMedian = float(point[4])
      xHigh1Sig = float(point[5])
      xHigh2Sig = float(point[6])
      xLow1Sig = float(point[3])
      xLow2Sig = float(point[2])

      box2Sig = mpl.Rectangle([xLow2Sig,yLow],xHigh2Sig-xLow2Sig,barHeight,color='#FFFF00')
      self.ax1.add_artist(box2Sig)
      boxes1.append(box2Sig)

      box1Sig = mpl.Rectangle([xLow1Sig,yLow],xHigh1Sig-xLow1Sig,barHeight,color='#00FF00')
      self.ax1.add_artist(box1Sig)
      boxes1.append(box1Sig)

      l = mpl.Line2D([xMedian,xMedian],[yLow,yLow+barHeight],color='k',lw=3.)
      self.ax1.add_artist(l)
      medianLines.append(l)
  def writeTickLabels(self):
    #self.axLeftBound = 0.127
    for lab,i in zip(self.xLabels,range(len(self.xLabels))):
        yInFig = self.fig.transFigure.inverted().transform(
            self.ax1.transData.transform((1.0,i+0.25))
            )[1]
        size = self.tickLabelFontSize
        ha = 'right'
        xPos = self.axLeftBound-0.01
        if len(lab)> 10:
          size *= 0.95
          ha = 'center'
          xPos = self.axLeftBound/2.0
        self.fig.text(xPos,yInFig,lab,va='center',ha=ha,size=size)
  def highlight(self):
    nCats = len(self.data)
    boxes = []
    axTrans = self.ax1.transAxes
    color = 'lightSkyBlue'
    color = 'PowderBlue'
    color = 'lightBlue'
    for i in range(len(self.data)):
      if os.path.split(self.data[i][0])[1] in self.toHighlight:
        tmpBox =  mpl.Rectangle([0.0,i/float(nCats)],self.maxX,1./nCats,fill=True,fc=color,ec=color,transform=axTrans,clip_on=False)
        self.ax1.add_artist(tmpBox)
        boxes.append(tmpBox)

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
  if (not args.bdtCut) and (not args.higgsMass):
    canvas.SetLogx(1)
    canvas.SetLogy(1)
  
  mpl.rcParams["font.family"] = "sans-serif"
  #print mpl.rcParams["backend"]

  ylimits=[0.1,100.0]

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
  
    energyStrWrite = energyStr
    if energyStr == "7P8TeV":
      energyStrWrite = "7 & 8 TeV"
    else:
      energyStrWrite = energyStr.replace("TeV"," TeV")
    if args.signalInject > 0.0:
      energyStrWrite += "   Signal Injected: {0:.1f}#times SM".format(args.signalInject)
      energyStrWrite += " m_{H} = "+"{0:.1f}".format(args.signalInjectMass)+" GeV/c^{2}"
    caption2 = "#sqrt{s} = "+energyStrWrite
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
        if plotName == "IncPtCut":
          xlabel="p_{T}(#mu#mu) Cut [GeV/c]"
          if energyStr == "8TeV":
            ylimits = [0.,16.]
          elif energyStr == "7TeV":
            ylimits = [0.,32.]
      elif args.higgsMass:
        if energyStr == "8TeV":
            ylimits = [0.,15.]
        elif energyStr == "7TeV":
            ylimits = [0.,40.]
        xlabel="m_{H} [GeV/c^{2}]"
        caption3 = "L = {0:.1f} fb^{{-1}}".format(float(lumisToUse[energyStr]))
      #elif period == "14TeV":
      #  title = "Standard Model H#rightarrow#mu#mu"
      title = titleMap[plotName]
      incPlot = RelativePlot(data,canvas,legend,title,caption2=caption2,ylimits=ylimits,energyStr=energyStrWrite,xlabel=xlabel,caption3=caption3,showObs=args.higgsMass)
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
    if energyStr == "7P8TeV":
      energyStrWrite = "7 & 8 TeV"
    else:
      energyStrWrite = energyStr.replace("TeV"," TeV")
    #print compareData
    if len(compareData)==0:
        print("No Data to Compare for {0}!!".format(period))
        continue
    lumiStrWrite = "{0:.1f}".format(float(desiredLumiStr))

    anotation1 = "$\sqrt{s}$ = "+energyStrWrite
    anotation1 += ", $\mathcal{L}$ = "+lumiStrWrite+" fb$^{-1}$"
    if energyStr == "7P8TeV":
      anotation1 = "$\sqrt{s}$ = 7 TeV, $\mathcal{L}$ = "
      anotation1 += "{0:.1f}".format(lumiDict["7TeV"])
      anotation1 += " fb$^{-1}$"
      anotation1 += " &\n "
      anotation1 += "$\sqrt{s}$ = 8 TeV, $\mathcal{L}$ = "
      anotation1 += "{0:.1f}".format(lumiDict["8TeV"])
      anotation1 += " fb$^{-1}$"
    anotation2 = '$m_H=125$ GeV/$c^2$'

    ## Inclusive Categories
    
    veto = ["VBF","IncBDT"]
    mustBe="(IncPresel.*)_(.+)_[.\d]+.txt.out"
    tmpMap = { 
  "IncPreselPtG10BB":"Non-VBF BB",
  "IncPreselPtG10BO":"Non-VBF BO",
  "IncPreselPtG10BE":"Non-VBF BE",
  "IncPreselPtG10OO":"Non-VBF OO",
  "IncPreselPtG10OE":"Non-VBF OE",
  "IncPreselPtG10EE":"Non-VBF EE",
  "IncPreselPtG10":"Non-VBF \nNo Categories",
  "IncPreselCat":"Non-VBF\nCategory\nCombination",
        }
    compareData = getData(fnGlobStr,matchString=mustBe,dontMatchStrings=veto,doSort=False)
    if len(compareData)>0:
      comparePlot = ComparePlotTable(compareData,titleMap=tmpMap,vertLine1=False,anotation1=anotation1,anotation2=anotation2)
      comparePlot.save(outDir+"compareNonVBFCats"+"_"+energyStr)

    ## VBF v. Inclusive Categories
    
    veto = ["EE",'BB',"BO","BE","OE","OO","BDTCutVBF"]
    mustBe="(.*)_(.+)_[.\d]+.txt.out"
    veto2 = ["BDTCut","IncPresel","Presel","PreselCat","VBFPresel","IncPreselPtG10","BDTCutCat"]
    tmpMap = { 
        "VBFBDTCut":'VBF',
        "IncPreselCat":'Non-VBF',
        "BDTCutCatVBFBDTOnly":'VBF & Non-VBF\nCombination',
        }
    compareData = getData(fnGlobStr,matchString=mustBe,dontMatchStrings=veto,doSort=False)
    veto3 = []
    for i in range(len(compareData)):
        for j in veto2:
          if os.path.split(compareData[i][0])[1] == j:
            veto3.append(i)
    for i in reversed(veto3):
        compareData.pop(i)
    if len(compareData)>0:
      comparePlot = ComparePlotTable(compareData,titleMap=tmpMap,vertLine1=False,anotation1=anotation1,anotation2=anotation2)
      comparePlot.save(outDir+"compareFinal"+"_"+energyStr)

    ## All Categories

    mustBe = r"(.+)_(.+)_[.\d]+.txt.out"
    veto = ["EE",'BB',"BO","BE","OE","OO"]
    compareData = getData(fnGlobStr,matchString=mustBe,dontMatchStrings=veto,doSort=False)
    comparePlot = ComparePlotTable(compareData,titleMap={},vertLine1=False,anotation1=anotation1,anotation2=anotation2)
    comparePlot.save(outDir+"all"+"_"+energyStr)

    ## Inc BDT v. Presel

    mustBe = r"(Inc.+Cat)_(.+)_[.\d]+.txt.out"
    veto = ["EE",'BB',"BO","BE","OE","OO"]
    veto2 = ["BDTCut"]
    compareData = getData(fnGlobStr,matchString=mustBe,dontMatchStrings=veto,doSort=False)
    veto3 = []
    for i in range(len(compareData)):
        for j in veto2:
          if os.path.split(compareData[i][0])[1] == j:
            veto3.append(i)
    for i in reversed(veto3):
        compareData.pop(i)
    if len(compareData)>0:
      comparePlot = ComparePlotTable(compareData,titleMap={},vertLine1=False,anotation1=anotation1,anotation2=anotation2)
      comparePlot.save(outDir+"IncBDTVPresel"+"_"+energyStr)
