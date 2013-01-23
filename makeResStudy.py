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

from makeCards import makePDFSig, SIGNALFIT

urLegendPos = [0.65,0.67,0.9,0.9]
minMass = 110.
maxMass = 160.
vminMass = 60.
vmaxMass = 160.
signalRange = [120.,130.]
plotRange = [110,140.]
rooPlotRange = root.RooFit.Range(*plotRange)

class ResStudy:
  def __init__(self,fileName,title):
    f = root.TFile(fileName)
    self.f = f
    self.colors = [1,root.kRed+1,root.kBlue+1,root.kGreen+1,root.kCyan,root.kPink]
    self.title = title
    self.saveTitle = re.sub(r"rightarrow","",title)
    self.saveTitle = re.sub(r"#","",self.saveTitle)
    self.saveTitle = re.sub(r"{\W}+","",self.saveTitle)
    self.saveTitle = re.sub(r" +","",self.saveTitle)
    print self.title
    print self.saveTitle
    # Old Code
    getData = getattr(self,"getData")
    countsList = []
    resList = []
    quantList = []
    categoriesList = []
    prefix = "IncPresel"
    categories = ["BB","BO","BE","OO","OE","EE"]
    if re.search(r"VBF",fileName) or re.search(r"vbf",fileName) :
      prefix = "VBFPresel"
      categories = ["BB","NotBB"]
    self.categories = categories
    if True:
      tmpCounts, tmpRes, tmpQuants, hists = getData(f,categories,prefix)
      countsList.append(tmpCounts)
      resList.append(tmpRes)
      quantList.append(tmpQuants)
      categoriesList.append(categories)
      self.histsCatDict = hists
    self.titles = [title]
    self.countsList = countsList
    self.resList = resList
    self.quantList = quantList
    self.categoriesList = categoriesList

    # New Code
    self.workspace = root.RooWorkspace("w")
    wImport = getattr(self.workspace,"import")
    mMuMu = root.RooRealVar("mMuMu","m_{#mu#mu} [GeV]",vminMass,vmaxMass)
    mMuMu.setRange("signal",signalRange[0],signalRange[1])
    mMuMu.setRange("signalfit",SIGNALFIT[0],SIGNALFIT[1])
    mMuMu.setRange("all",minMass,maxMass)
    mMuMu.setRange("vall",vminMass,vmaxMass)
    mMuMu.setRange("draw",plotRange[0],plotRange[1])
    wImport(mMuMu)
    self.mMuMu = mMuMu
    for i in categories:
       tmpParamList, tmpDebug = makePDFSig(i,self.histsCatDict[i],mMuMu,minMass,maxMass,wImport,i)
       #print tmpDebug

    #self.workspace.Print()

    self.canvas = root.TCanvas("canvas")

  def plot(self,ymax=None):
    plotPies = getattr(self,"plotPies")
    plotRes = getattr(self,"plotRes")
    plotPDF = getattr(self,"plotPDF")
    plotData = getattr(self,"plotData")
    plotPies()
    plotRes(ymax)
    plotPDF()
    plotData()

  def plotRes(self,ymax=None):
    colors = ["b","r","g","k","p"]
    fig = mpl.figure()
    ax = fig.add_subplot(111)
    if ymax != None:
      ax.set_ylim(0.0,ymax)
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
    title2 = r"$gg \rightarrow H \rightarrow \mu\mu$ Resolution"
    if "VBF" in title or "vbf" in title:
      title2 = r"VBF $H \rightarrow \mu\mu$ Resolution"
    ax.set_title(title2)
    ax.set_ylabel(r"$1\sigma$ Quantile Resolution  [GeV]")

    iPlot = 0
    for x in dataCoords:
      color = colors[iPlot]
      iPlot += 1
      dataCord = (x,1.)
      displayCord = ax.transData.transform(dataCord)
      figCord = fig.transFigure.inverted().transform(displayCord)

    fig.savefig("resPlot"+self.saveTitle+".png")
    fig.savefig("resPlot"+self.saveTitle+".pdf")

  def plotPies(self):
    vertical = True
    fig = None
    subplotBase = None
    if len(self.titles) == 1:
      fig = mpl.figure(figsize=(8,8))
    elif vertical:
      fig = mpl.figure(figsize=(8,12))
      subplotBase = len(self.titles)*100+11
    else:
      pass
    iPlot = 0
    for title, countsMap, categories in zip(self.titles,self.countsList,self.categoriesList):
      ax = None
      if len(self.titles) == 1:
        ax = fig.add_subplot(111)
      elif vertical:
        ax = fig.add_subplot(subplotBase+iPlot)
      else:
        pass
      ax.pie([countsMap[i] for i in categories],labels=tuple(categories),
            shadow=False,autopct="%1.0f%%")
      if "VBF" in title or "vbf" in title:
        ax.set_title(r"VBF $H \rightarrow \mu\mu$ Fractions")
      else:
        ax.set_title(r"$gg \rightarrow H \rightarrow \mu\mu$ Fractions")
      iPlot += 1
    fig.savefig("resFractions"+self.saveTitle+".png")
    fig.savefig("resFractions"+self.saveTitle+".pdf")

  def getData(self,infile,categories, prefix=""):
    keyList = categories + ["All"]
    histMap = {}
    countsMap = {}
    colors = [root.kRed,root.kBlue,root.kGreen,root.kOrange,root.kPink,root.kCyan,root.kBlue+1,root.kBlack]
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
    quantilesString += "{0:<8}{1:<10}{2:^24}{3:^24}\n".format("","Median","1 Sigma Quantile Resolution","2 Sigma Quantile Resolution")
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

    return countsMap, resMap, quantMap, histMap

  def plotData(self):
    saveNameBase = "resSigHist_"+self.saveTitle
    datasets = []
    iColor = 0
    if True:
      self.canvas.Clear()
      frame = self.mMuMu.frame(rooPlotRange)
      frame.SetTitle("")
      frame.SetYTitle("Signal MC Events")
      leg = root.TLegend(*urLegendPos)
      leg.SetFillColor(0)
      leg.SetLineColor(0)
      for i in self.categories:
        tmpDataset = self.workspace.data(i+"_Template")
        rooLCol = root.RooFit.LineColor(self.colors[iColor])
        rooMCol = root.RooFit.MarkerColor(self.colors[iColor])
        rooNameStr = "Curve_"+i+"_Template"
        rooName = root.RooFit.Name(rooNameStr)
        tmpDataset.plotOn(frame,rooLCol,rooMCol,rooName)
        tmpDatasetH = frame.getHist(rooNameStr)
        leg.AddEntry(tmpDatasetH,i,"p")
        datasets.append(tmpDataset)
        iColor +=1
      frame.Draw()
      leg.Draw()
      tlatex = root.TLatex()
      tlatex.SetNDC()
      tlatex.SetTextFont(root.gStyle.GetLabelFont())
      tlatex.SetTextSize(0.05)
      tlatex.SetTextAlign(12)
      tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
      tlatex.SetTextAlign(32)
      tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,self.title)
      saveAs(self.canvas,saveNameBase)

  def plotPDF(self):
    saveNameBase = "resSigShape_"+self.saveTitle
    curves = []
    iColor = 0
    self.canvas.Clear()
    frame = self.mMuMu.frame(rooPlotRange)
    frame.SetTitle("")
    frame.SetYTitle("Signal MC Events")
    leg = root.TLegend(*urLegendPos)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    for i in self.categories:
      if True:
        tmpCurve = self.workspace.pdf(i)
        rooLCol = root.RooFit.LineColor(self.colors[iColor])
        rooNameStr = "Curve_"+i
        rooName = root.RooFit.Name(rooNameStr)
        rooRange = root.RooFit.Range("all")
        tmpCurve.plotOn(frame,rooLCol,rooName,rooRange)
        tmpCurveH = frame.getCurve(rooNameStr)
        leg.AddEntry(tmpCurveH,i,"l")
        curves.append(tmpCurve)
        iColor += 1
    frame.Draw()
    leg.Draw()
    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(root.gStyle.GetLabelFont())
    tlatex.SetTextSize(0.05)
    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,self.title)

    saveAs(self.canvas,saveNameBase)

if __name__ == "__main__":

  infiles = []
  titles = []
  #infiles.append("input/ggHmumu125_8TeV.root")
  #infiles.append("input/vbfHmumu125_8TeV.root")
  infiles.append("input/muscle/ggHmumu125_8TeV.root")
  infiles.append("input/muscle/vbfHmumu125_8TeV.root")
  #infiles.append("input/smearing/ggHmumu125_8TeV.root")
  #infiles.append("input/smearing/vbfHmumu125_8TeV.root")
  #infiles.append("input/rochester/ggHmumu125_8TeV.root")
  #infiles.append("input/rochester/vbfHmumu125_8TeV.root")
  titles.append("gg #rightarrow H #rightarrow #mu#mu")
  titles.append("VBF H #rightarrow #mu#mu")

  for f,t in zip(infiles,titles):
    rs = ResStudy(f,t)
    rs.plot(3.5)
