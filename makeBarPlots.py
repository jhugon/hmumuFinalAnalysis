#!/usr/bin/env python

from helpers import *
import ROOT as root
import glob
import re

import matplotlib.pyplot as mpl
import numpy

from xsec import *

#######################################

class BarPlot:
  def __init__(self,data,xlabel="Expected 95% CL Limit $\sigma/\sigma_{SM}$",titleMap={},showObs=False,xlim=[0,20]):
    fig = mpl.figure()
    self.fig = fig
    ax1 = fig.add_subplot(111)
    self.ax1 = ax1
    ax1bounds = ax1.get_position().bounds
    ax1.set_position([0.25,0.1,0.7,0.85])

    #data.sort(key=lambda x: x[0].lower())
    data.reverse()

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
    ax1.set_xlabel(xlabel)
    ax1.set_xlim(xlim)
    bars = ax1.barh(xPos,medians, 0.5, xerr=[low2sigs,high2sigs],ecolor="k")
    self.bars = bars
    xPosObs = [x+0.25 for x in xPos]
    # Now for the 1sigma error bars
    self.oneSig = ax1.errorbar(medians,xPosObs,xerr=[low1sigs,high1sigs],color="g",linestyle="None")
    if showObs:
      self.obs = ax1.plot(obs,xPosObs,marker="|",color="r",markersize=10,linestyle="None")
    writeInValues = getattr(self,"writeInValues")
    writeInValues(bars,high2sigs)
  def save(self,saveName):
    self.fig.savefig(saveName+".png")
    self.fig.savefig(saveName+".pdf")
  def writeInValues(self,rects,high2sigs):
    axisWidth = self.ax1.get_xlim()[1]-self.ax1.get_xlim()[0]
    size="medium"
    if len(rects)<5:
      size="large"
    elif len(rects)>12:
      size="xx-small"
    elif len(rects)>9:
      size="x-small"
    maxWidth = axisWidth
    for rect,sigUp in zip(rects,high2sigs):
      width = rect.get_width()
      if (width/maxWidth < 0.1): # The bars aren't wide enough to print the ranking inside
        xloc = width + sigUp + maxWidth*0.01 # Shift the text to the right side of the right edge
        clr = 'black' # Black against white background
        align = 'left'
      else:
        xloc = width - maxWidth*0.01 # Shift the text to the left side of the right edge
        clr = 'white' # White on magenta
        align = 'right'
      yloc = rect.get_y()+rect.get_height()/2.0 
      valueStr = "{0:.1f}".format(width)
      self.ax1.text(xloc, yloc, valueStr, horizontalalignment=align,
            verticalalignment='center', color=clr, weight='bold',size=size)

if __name__ == "__main__":

  # Wrong Dataformat: title, obs, median, 1sig high, 1sig low, 2sig high, 2sig low

  data = [
    #("MuScle PFIso",4.2, 2.2, 3.0,4.33,6.2,8.8), #Cat ISMEAR=2
    #("MuScle TrkIso",4.98, 2.84,3.84,5.4,7.8,10.9) #Cat ISMEAR=2

    #BDT Comb no cat
    #("PFIsoLooseBDT",6.399,2.67,3.689,5.359,7.944,11.37), #No Smearing
    #("PFIsoTightBDT",4.352,3.011,4.126,5.953,8.78,12.46), #no smear
    #("TrkIsoLooseBDT",11.02,3.55,4.833,6.9,10.1,14.13), #no smear
    #("TrkIsoTightBDT",12.63,4.04,5.5,7.8,11.3,15.8), #no smear

    #Presel Comb no cat
    ("PF Loose Iso",5.31,4.2,5.7,8.16,11.8,16.5), #No Smearing
    ("PF Tight Iso",5.6,4.14,5.6,8.,11.5,16.19), #no smear
    ("Trk Loose Iso",6.0,4.6,6.18,8.66,12.35,17.2), #no smear
    ("Trk Tight Iso",7.23,4.62,6.15,8.66,12.3,17.2), #no smear
  ]

  data2 = [
    ("Normal MuScle",3.21, 2.5,3.4,4.83,7.,9.7), #Cat
#    ("+Shape Res Unc.", 3.9,2.4,3.22,4.516,6.37,8.55), #Cat
    ("+3% $\ln(N)$ Unc.", 3.2, 2.49, 3.39,4.859,7.,9.8 ) #Cat
  ]

  energyStr = "8 TeV"
  desiredLumiStr = "12"


  myPlot = BarPlot(data)
  myPlot.fig.text(0.9,0.2,"$\mathcal{L}="+desiredLumiStr+"$ fb$^{-1}$",horizontalalignment="right",size="x-large")
  myPlot.fig.text(0.9,0.27,"$\sqrt{s}=$"+energyStr,horizontalalignment="right",size="x-large")
  #myPlot.fig.text(0.9,0.54,"Preliminary!!",horizontalalignment="right",size="x-large")
  #myPlot.fig.text(0.9,0.47,"Using Untested Smearing!",horizontalalignment="right",size="x-large")
  myPlot.fig.text(0.9,0.55,"No BDTs",horizontalalignment="right",size="x-large")
  myPlot.fig.text(0.9,0.48,"No Resolution Categories",horizontalalignment="right",size="x-large")
  myPlot.save("trkIsoCompare")

  myPlot = BarPlot(data2)
  myPlot.fig.text(0.9,0.2,"$\mathcal{L}="+desiredLumiStr+"$ fb$^{-1}$",horizontalalignment="right",size="x-large")
  myPlot.fig.text(0.9,0.27,"$\sqrt{s}=$"+energyStr,horizontalalignment="right",size="x-large")
  myPlot.save("resErrLimits")


