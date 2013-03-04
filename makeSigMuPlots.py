#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description="Makes cards for use in the CMS Combine tool.")
parser.add_argument("--signalInject", help="Inject Signal with Strength into data_obs",type=float,default=0.0)
parser.add_argument("--signalInjectMass", help="Mass For Injected Signal",type=float,default=125.0)
parser.add_argument("-m","--higgsMass", help="Makes plots v. Higgs Mass",action="store_true",default=False)
parser.add_argument("-p","--pValue", help="Makes p-value plots instead of significance",action="store_true",default=False)
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

  "IncPreselPtG10BB":"Non-VBF BB",
  "IncPreselPtG10BO":"Non-VBF BO",
  "IncPreselPtG10BE":"Non-VBF BE",
  "IncPreselPtG10OO":"Non-VBF OO",
  "IncPreselPtG10OE":"Non-VBF OE",
  "IncPreselPtG10EE":"Non-VBF EE",
  "IncPreselPtG10":"Non-VBF",
  "BDTCutCatVBFBDTOnly": "VBF & Non-VBF Combination"
}

colorMap = {
  "IncPresel":root.kRed,
  "VBFPresel":root.kBlue,
  "VBFBDTCut":root.kBlue,
  "IncPreselCat":root.kRed,

  "IncPreselPtG10BB":root.kCyan,
  "IncPreselPtG10BO":root.kOrange,
  "IncPreselPtG10BE":root.kGreen,
  "IncPreselPtG10OO":root.kOrange+3,
  "IncPreselPtG10OE":root.kMagenta,
  "IncPreselPtG10EE":root.kViolet+2,
  "IncPreselPtG10":root.kBlue,
  "BDTCutCatVBFBDTOnly": 1
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

# Do match and don't match w/o extensions.  \.expsig and \.sig are added automatically
def getDataSig(fileString,matchString=r"_([-\d.]+)\.txt",dontMatchStrings=[],doSort=True,getPValue=False):
  def sortfun(s):
    match = re.search(matchString,s)
    result = 1e12
    if match:
      result = float(match.group(1))
    return result

  #print fileString
  #print matchString
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
    match = re.search(matchString+r"\.sig",fname)
    obs = -10.0
    exp = -20.0
    xNum = -10.0
    if match:
      xNum = match.group(1)
    for line in tmpF:
      obsMatch = None
      if getPValue:
        obsMatch = re.search(r"p-value = ([.\deE]+)",line)
      else:
        obsMatch = re.search(r"^Significance:[\s]+([.\deE]+)",line)
      if obsMatch:
        obs = obsMatch.group(1)
    expFname = os.path.splitext(fname)[0]+".expsig"
    """
    try:
      tmpFExp = open(expFname)
      for line in tmpFExp:
        obsMatch = None
        if getPValue:
          obsMatch = re.search(r"p-value = ([.\deE]+)",line)
        else:
          obsMatch = re.search(r"^Significance:[\s]+([.\deE]+)",line)
        if obsMatch:
          exp = obsMatch.group(1)
    except Exception:
      print("Expected Significance Not Found: "+expFname)
    """
    thisPoint = [xNum,obs,exp]
    #print thisPoint
    if thisPoint.count("-10.0")>0:
        continue
    if thisPoint.count(-10.0)>0:
        continue
    #print thisPoint
    result.append(thisPoint)
  #print result
  return result

def getDataMu(fileString,matchString=r"_([-\d.]+)\.txt\.mu",dontMatchStrings=[],doSort=True):
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
    for line in tmpF:
      #obsMatch = re.search(r"Best fit r:[\s]+([.\deE-]+)[\s]+([.\deE-]+)/+([.\deE-]+)",line)
      obsMatch = re.search(r"Best fit r:[\s]+([.\deE-]+)[\s]+\-([.\deE-]+)/\+([.\deE-]+)",line)
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
    muExtraGraph.SetLineColor(root.kBlue+1)
    muExtraGraph.SetMarkerStyle(20)
    muExtraGraph.SetMarkerSize(1.1)
    muExtraGraph.SetMarkerColor(root.kBlue+1)
    muExtraGraphErr.SetFillColor(root.kCyan)
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
    tlatex.SetTextFont(root.gStyle.GetLabelFont())
    #tlatex.SetTextSize(0.05)
    tlatex.SetTextSize(0.04)
    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,caption)

    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin()-0.03,0.88,caption2)
    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin()-0.03,0.82,caption3)

    for g in self.vertLines:
      g.Draw("l")

    canvas.RedrawAxis()

class PValuePlotTogether:
  def __init__(self,dataDict, canvas, caption="Standard Model H#rightarrow#mu#mu", ylabel="p-Value", xlabel="m_{H} [GeV/c^{2}]",caption2="",caption3="",ylimits=[],xlimits=[],energyStr="8TeV"):
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
        print xmin, xmax
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
    tlatex.SetTextFont(root.gStyle.GetLabelFont())
    #tlatex.SetTextSize(0.05)
    tlatex.SetTextSize(0.04)
    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,caption)

    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin()-0.03,0.88,caption2)
    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin()-0.03,0.82,caption3)

    #legend = root.TLegend(0.58,0.70,0.9,0.9) # UR
    baselineY = gStyle.GetPadBottomMargin()+0.02
    marginR = 1.0-gStyle.GetPadRightMargin()+0.01
    legend = root.TLegend(marginR-0.28,baselineY,marginR,baselineY+0.4)
    if len(graphs) < 5:
      legend = root.TLegend(marginR-0.32,baselineY,marginR,baselineY+0.25)
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

class ComparePlot:
  def __init__(self,data,ylabel="Expected 95% CL Limit $\sigma/\sigma_{SM}$",titleMap={},showObs=False):
    fig = mpl.figure()
    self.fig = fig
    ax1 = fig.add_subplot(111)
    self.ax1 = ax1
    ax1bounds = ax1.get_position().bounds
    ax1.set_position([0.25,0.1,0.7,0.85])

    data.sort(key=lambda x: x[0].lower())

    obs = [float(point[1]) for point in data]
    exp = None
    high1sigs = None
    low1sigs = None
    if len(data[0])==4: #mu
      #thisPoint = [xNum,obs,low1sig,high1sig]
      high1sigs = [float(point[3]) for point in data]
      low1sigs = [float(point[2]) for point in data]
    else: #sig len=3
      #thisPoint = [xNum,obs,exp]
      exp = [float(point[2]) for point in data]
      high1sigs = [0.0 for point in data]

    xPos = numpy.arange(len(obs))
    xLabels = [point[0] for point in data]
    xLabels = [re.sub(r".*/","",s) for s in xLabels]
    xLabels = [titleMap[s] if titleMap.has_key(s) else s for s in xLabels]

    ax1.set_yticks(xPos+0.25)
    ax1.set_yticklabels(tuple(xLabels),size="small")
    ax1.set_xlabel(ylabel)
    #ax1.set_xlim([0,20])
    bars = None
    if len(data[0])==4: #mu
      bars = ax1.barh(xPos,obs, 0.5, xerr=[low1sigs,high1sigs],ecolor="k")
    else: #sig
      bars = ax1.barh(xPos,exp, 0.5)
    self.bars = bars
    xPosObs = [x+0.25 for x in xPos]
    if showObs and len(data[0]) != 4: #is sig
      self.obs = ax1.plot(obs,xPosObs,marker="|",color="r",markersize=10,linestyle="None")
    writeInValues = getattr(self,"writeInValues")
    writeInValues(bars,high1sigs)
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

  dirName = "statsInput/"
  
  outDir = "statsOutput/"
  
  root.gErrorIgnoreLevel = root.kWarning
  root.gROOT.SetBatch(True)
  setStyle()
  canvas = root.TCanvas()
  if not args.higgsMass:
    canvas.SetLogx(1)
  canvas.SetLogy(0)
  
  mpl.rcParams["font.family"] = "sans-serif"
  #print mpl.rcParams["backend"]

  ylimits=[]

  lumisToUse={"7TeV":lumiDict["7TeV"],"8TeV":lumiDict["8TeV"],"7P8TeV":lumiDict["8TeV"]+lumiDict["7TeV"]}
  
  for period in ["7TeV","8TeV","14TeV","7P8TeV"]:
    fnToGlob = dirName+"*_"+period+"_*.txt.out"
    allfiles = glob.glob(fnToGlob)

    ## Limit v. Lumi
    energyStr = ""
    plots = set()
    for fn in allfiles:
      match = re.search(r".*/(.+)_(.+)_[.\d]+.txt.out",fn)
      badPlot = re.search(r"Silly",fn)
      badPlot2 = re.search(r"Silly",fn)
      if match and not (badPlot or badPlot2):
        plots.add(match.group(1))
        energyStr = match.group(2)

    if energyStr == "7P8TeV":
      energyStr = "7 & 8 TeV"
    else:
      energyStr.replace("TeV"," TeV")
  
    caption2 = "#sqrt{s} = "+energyStr
    caption3 = ""
    if args.signalInject>0.0:
      caption3 = "Signal Injected {0:.1f}#timesSM".format(args.signalInject)
      caption3 += " m_{H} = "+"{0:.1f}".format(args.signalInjectMass)+" GeV/c^{2}"
    legend = root.TLegend(0.58,0.70,0.9,0.9)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    
    for ytitle,fnPref in zip(["Significance","Error on #sigma_{Measured}/#sigma_{SM}","#sigma_{Measured}/#sigma_{SM}"],["sig","mu","muToy"]):
      if args.higgsMass and fnPref == "mu":
        continue
      if args.pValue and fnPref == "sig":
        canvas.SetLogy(1)
        ytitle = "p-Value of Background Only Hypothesis"
      else:
        canvas.SetLogy(0)
      for plotName in plots:
        data = None
        doMuExtraPlot=False
        showObs=False
        if fnPref == "sig":
          data = getDataSig(dirName+plotName+"_"+period+"_*.txt*",getPValue=args.pValue)
          showObs=True
        else:
          data = getDataMu(dirName+plotName+"_"+period+"_*.txt*")
        if fnPref == "muToy":
          doMuExtraPlot = True
        #print("{0} {1} v. Lumi: {2}".format(period,fnPref, data))
        if len(data)<=1:
          continue
        title = "Standard Model H#rightarrow#mu#mu"
        xlabel="Integrated Luminosity [fb^{-1}]"
        if args.higgsMass:
          title = titleMap[plotName]
          xlabel="m_{H} [GeV/c^{2}]"
        incPlot = RelativePlot(data,canvas,legend,title,caption2=caption2,caption3=caption3,ylabel=ytitle,energyStr=energyStr,doMuExtraPlot=doMuExtraPlot,showObs=showObs,xlabel=xlabel)
        saveAs(canvas,outDir+fnPref+plotName+"_"+period)

    ## All p-values together plot
    canvas.SetLogy(1)
    pValueVetos = [
        [
          "VBFBDTCut",
          "BDTCutCatVBFBDTOnly"
        ],
        [
          "IncPreselPtG10BB",
          "IncPreselPtG10BE",
          "IncPreselPtG10BO",
          "IncPreselPtG10EE",
          "IncPreselPtG10OE",
          "IncPreselPtG10OO"
        ]
    ]
    for saveName,vetos in zip(["NonVBF","Final"],pValueVetos):
      if len(plots)==0 or not args.higgsMass:
        continue
      pValueDict = {}
      for plotName in plots:
        if plotName in vetos:
            continue

        data = getDataSig(dirName+plotName+"_"+period+"_*.txt*",getPValue=True)
        pValueDict[plotName] = data
      pValueAllPlot = PValuePlotTogether(pValueDict,canvas,caption2=caption2,caption3=caption3,energyStr=energyStr)
      saveAs(canvas,outDir+"pValues_"+saveName+period)
    canvas.SetLogy(0)
    
    ## Compare all types of limits
    print "Doing Compare..."
    if period == "14TeV":
        continue
    #veto = [r"CNC",r"PM","BB","BO","BE","OO","OE","EE","NotBB"]
    #veto = [r"CNC",r"PM","Presel","BB","BO","BE","OO","OE","EE","NotBB"]
    #veto = [r"CNC",r"PM","BDT","BB","BO","BE","OO","OE","EE","NotBB"]
    #veto = [r"CNC",r"PM","Presel"]
    #veto = [r"CNC",r"PM","BDT"]
    veto = []

    mustBe = r"(.+)_(.+)_[.\d]+.txt"
    #mustBe = r"(.+Cat)_(.+)_[.\d]+.txt"

    desiredLumiStr=str(lumisToUse[period])
    fnGlobStr = dirName+"*_"+period+"_"+desiredLumiStr+".txt*"
    compareDataSig = getDataSig(fnGlobStr,matchString=mustBe,dontMatchStrings=veto,doSort=False)
    compareDataMu = getDataMu(fnGlobStr,matchString=mustBe,dontMatchStrings=veto,doSort=False)
    #print compareDataSig
    #print compareDataMu

    for datCase,ytitle,fnPref in zip([compareDataSig,compareDataMu],["Expected Significance","$\sigma/\sigma_{SM}$"],["sig","mu"]):
      if len(datCase)==0:
        print("No Data to Compare for {0} {1}!!".format(period,fnPref))
        continue
      comparePlot = ComparePlot(datCase,titleMap=comparisonMap,showObs=False,ylabel=ytitle)
      comparePlot.ax1.set_xlim(*[0,20])
      comparePlot.fig.text(0.9,0.2,"$\mathcal{L}="+desiredLumiStr+"$ fb$^{-1}$",horizontalalignment="right",size="x-large")
      comparePlot.fig.text(0.9,0.27,"$\sqrt{s}=$"+energyStr,horizontalalignment="right",size="x-large")
      comparePlot.save(outDir+fnPref+"Compare"+"_"+energyStr)
