#!/usr/bin/env python

from helpers import *
import ROOT as root
import glob
import re
import os.path

import matplotlib.pyplot as mpl
import numpy

from xsec import *

class ComparePlotTable:
  def __init__(self,data,ylabel="95% CL Limit $\sigma/\sigma_{SM}$",titleMap={},xlimits=[],brazil=True,obsCircles=True,vertLine1=True,toHighlight=[],anotation1='',anotation2='',showObs=True,signalInject=0.0):
    self.showObs = showObs
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
    if obsCircles and showObs:
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
    if signalInject > 0.0:
      ax1.text(xmax/2.0,len(medians)-1.25, r"Signal Injected: {0:.1f}$\times$ SM".format(signalInject),
                    va="center",ha='center',size=30.,color='r')
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
          if self.showObs:
            s = "{0:.1f}".format(float(self.data[i][j+1]))
          else:
            s = "{0}".format("XX")
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
