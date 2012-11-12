#!/usr/bin/env python

from helpers import *
import math
import os.path
import glob
import numpy
import matplotlib.pyplot as mpl

root.gErrorIgnoreLevel = root.kWarning
root.gROOT.SetBatch(True)
root.gStyle.SetOptStat(0)

class ResStudy:
  def __init__(self,fileNames,titles):
    getData = getattr(self,"getData")
    countsList = []
    resList = []
    quantList = []
    categoriesList = []
    for fn in fileNames:
        prefix = "IncPresel"
        categories = ["BB","BO","BE","OO","OE","EE"]
        if re.search(r"VBF",fn) or re.search(r"vbf",fn) :
          prefix = "VBFPresel"
          categories = ["BB","NotBB"]
        tmpCounts, tmpRes, tmpQuants = getData(fn,categories,prefix)
        countsList.append(tmpCounts)
        resList.append(tmpRes)
        quantList.append(tmpQuants)
        categoriesList.append(categories)
    self.titles = titles
    self.countsList = countsList
    self.resList = resList
    self.quantList = quantList
    self.categoriesList = categoriesList

  def plot(self):
    plotPies = getattr(self,"plotPies")
    plotRes = getattr(self,"plotRes")
    plotPies()
    plotRes()

  def plotRes(self):
    colors = ["b","r","g","k","p"]
    fig = mpl.figure()
    ax = fig.add_subplot(111)
    allCatList = []
    for cat in self.categoriesList:
      allCatList.extend(cat)
    catSoFar = 0
    iPlot = 0
    amountUp = mpl.rcParams["figure.subplot.bottom"]
    dataCoords = []
    rectList = []
    for title, resMap, categories in zip(self.titles,self.resList,self.categoriesList):
      color = colors[iPlot]
      xPos = numpy.arange(catSoFar,catSoFar+len(categories))+0.25
      rects = ax.bar(xPos,[resMap[i] for i in categories],0.5,color=color)
      catSoFar += len(categories)
      dataCoords.append((xPos[len(xPos)-1]+xPos[0]+0.5)/2.)
      rectList.append(rects)
      iPlot += 1

    tsize="small"
    if len(allCatList)<5:
      tsize="medium"
    elif len(allCatList)>12:
      tsize="xx-small"
    elif len(allCatList)>9:
      tsize="x-small"
    maxHeight = ax.get_ylim()[1]-ax.get_ylim()[0]
    for rects in rectList:
      for rect in rects:
        yloc = rect.get_height() - maxHeight*0.01
        clr = 'white'
        xloc = rect.get_x()+rect.get_width()/2.0 
        valueStr = "{0:.1f}".format(rect.get_height())
        ax.text(xloc, yloc, valueStr, horizontalalignment='center',
              verticalalignment="top", color=clr, weight='bold',size=tsize)


    tickPos = numpy.arange(len(allCatList))+0.25
    ax.set_xticks(tickPos+0.25)
    ax.set_xticklabels(tuple(allCatList))
    ax.set_title(r"$H \rightarrow \mu\mu$ Width")
    ax.set_ylabel(r"$\sigma$  [GeV]")

    iPlot = 0
    for x, title, in zip(dataCoords, self.titles):
      color = colors[iPlot]
      iPlot += 1
      dataCord = (x,1.)
      displayCord = ax.transData.transform(dataCord)
      figCord = fig.transFigure.inverted().transform(displayCord)
      text = r"$gg \rightarrow H \rightarrow \mu\mu$"
      if title == "VBF" or title == "vbf":
        text = r"VBF $H \rightarrow \mu\mu$"
      fig.text(figCord[0],amountUp*0.3,text,color=color,
            ha="center",va="center")

    fig.savefig("resPlot.png")

  def plotPies(self):
    #fig = mpl.figure(figsize=(8,16))
    fig = mpl.figure(figsize=(8,12))
    subplotBase = len(self.titles)*100+11
    iPlot = 0
    for title, countsMap, categories in zip(self.titles,self.countsList,self.categoriesList):
      ax = fig.add_subplot(subplotBase+iPlot)
      ax.pie([countsMap[i] for i in categories],labels=tuple(categories),
            shadow=False,autopct="%1.0f%%")
      if title == "VBF" or title == "vbf":
        ax.set_title(r"VBF $H \rightarrow \mu\mu$ Fractions")
      else:
        ax.set_title(r"$gg \rightarrow H \rightarrow \mu\mu$ Fractions")
      iPlot += 1
    fig.savefig("resFractions.png")

  def getData(self,infilename,categories, prefix=""):
    keyList = categories + ["All"]
    histMap = {}
    countsMap = {}
    colors = [root.kRed,root.kBlue,root.kGreen,root.kOrange,root.kPink,root.kCyan,root.kBlue+1,root.kBlack]
    infile = root.TFile(infilename)
    histBase = "mDiMu"
    
    histMap["All"] = infile.Get(prefix+"/"+histBase)
    for c in categories:
      histMap[c] = infile.Get(prefix+c+"/"+histBase)
    
    for k in histMap:
      countsMap[k] = histMap[k].Integral()
    
    fractionsString = ""
    fractionsString += "{0:<8}{1:}\n".format("","Fraction Of Signal")
    countsAll = countsMap["All"]
    for k in keyList:
      #countsMap[k] = (countsMap[k]) / countsAll
      countsMap[k] /= countsAll
      fractionsString += "{0:<8}{1:>6.1%}\n".format(k,countsMap[k])
    
    for c in histMap:
      histMap[c].Scale(1.0/histMap[c].Integral())
    
    meanRMSString = ""
    meanRMSString += "{0:<8}{1:<10}{2:<6}\n".format("","Mean","RMS")
    for c in keyList:
      meanRMSString += "{0:<8}{1:<10.2f}{2:<6.2f}\n".format(c,histMap[c].GetMean(),histMap[c].GetRMS())
    
    resMap={}
    quantMap={}
    quantilesString = ""
    quantilesString += "{0:<8}{1:<10}{2:^24}{3:^24}\n".format("","Median","1 Sigma Quantile Width","2 Sigma Quantile Width")
    for c in keyList:
      quantiles = getMedianAndQuantileInterval(histMap[c],0.159)
      err = (quantiles[2]-quantiles[0])/2.
      quantMap[c] = quantiles
      resMap[c] = err
      quantiles2 = getMedianAndQuantileInterval(histMap[c],0.023)
      err2 = (quantiles2[2]-quantiles2[0])/2.
      quantilesString += "{0:<8}{1:<10.2f}{2:^24.2f}{3:^24.2f}\n".format(c,quantiles[1],err,err2)
    
    """
    for k,c in zip(keyList,colors):
      histMap[k].SetFillStyle(0)
      histMap[k].UseCurrentStyle()
      histMap[k].SetLineColor(c)
    
    drawn = False
    canvas = root.TCanvas()
    keyList.reverse()
    for k in keyList:
      if drawn:
        histMap[k].Draw("hist same")
      else:
        histMap[k].SetTitle("")
        histMap[k].GetXaxis().SetRangeUser(100,150)
        histMap[k].GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
        histMap[k].GetYaxis().SetTitle("Normalized Events")
        histMap[k].Draw("hist")
        drawn = True
    
    urLegendPos = [0.70,0.67,0.9,0.9]
    leg = root.TLegend(*urLegendPos)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    for k in keyList:
      leg.AddEntry(histMap[k],k,"l")
    
    leg.Draw()
    saveAs(canvas,"resolutionCategories")
    """

    return countsMap, resMap, quantMap

if __name__ == "__main__":
  infilename = "input/ggHmumu125_8TeV.root"

  infiles = []
  titles = []
  infiles.append("input/smearing/ggHmumu125_8TeV.root")
  titles.append("gg")
  infiles.append("input/smearing/vbfHmumu125_8TeV.root")
  titles.append("VBF")

  rs = ResStudy(infiles,titles)
  rs.plot()
  print rs.quantList
