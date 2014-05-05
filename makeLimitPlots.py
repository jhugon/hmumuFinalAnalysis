#!/usr/bin/env python

import optparse
parser = optparse.OptionParser(description="Makes Limit Plots from output text from combine tool.")
parser.add_option("--bdtCut", help="Makes plots v. BDT Cut Instead of Luminosity",action="store_true",default=False)
parser.add_option("--cutOpt", help="Makes plots for cutOptimization",action="store_true",default=False)
parser.add_option("-m","--higgsMass", help="Makes plots v. Higgs Mass",action="store_true",default=False)
parser.add_option("--signalInject", help="Sets a caption saying that signal was injected with strength",type=float,default=0.0)
parser.add_option("--signalInjectMass", help="Mass For Injected Signal",type=float,default=125.0)
parser.add_option("-p","--printLimits", help="Just Print The Limits",action="store_true",default=False)
parser.add_option("--xs", help="Limit on XS*BR",action="store_true",default=False)
parser.add_option("--blind",help="Don't Show Observed Limit",action="store_true",default=False)
args, fakeargs = parser.parse_args()

from helpers import *
import ROOT as root
import glob
import re
import os.path
import cPickle
from copy import deepcopy
from array import array
import sys

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

PRELIMINARYSTRING="CMS"

#######################################

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
                        
  #"Jets01SplitCatAll": "H#rightarrow#mu#mu 0,1-Jet Combination",
  "Jets01SplitCatAll": "0,1-Jet Combination",


  "Jet2CutsVBFPass":"2-Jet VBF Tight",
  "Jet2CutsGFPass":"2-Jet GF Tight",
  "Jet2CutsFailVBFGF":"2-Jet Loose",

  #"Jet2SplitCutsGFSplit" : "H#rightarrow#mu#mu 2-Jet Combination",
  #"CombSplitAll" : "H#rightarrow#mu#mu Combination",
  "Jet2SplitCutsGFSplit" : "2-Jet Combination",
  "CombSplitAll" : "Standard Model H #rightarrow #mu^{+}#mu^{-}",
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

def getData(fileString,matchString=r"_([-\d.]+)\.txt\.out",dontMatchStrings=[],doSort=True,xMax=1.0e20,xMin=-1.0e20):
  def sortfun(s):
    match = re.search(matchString,s)
    result = 1e12
    if match:
      result = float(match.group(1))
    return result

  print fileString
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
      if float(xNum) > xMax:
        continue
      if float(xNum) < xMin:
        continue
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
    if args.xs and args.higgsMass:
      matchEn = re.match(r".+_([1478P]+TeV)_.+",fname)
      if not matchEn:
        print("Error:Trying to make XS*BR limits w.r.t. higgs Mass, can't match energy for filename: "+fname)
        print("Exiting.")
        sys.exit(1)
      energy = matchEn.group(1)
      energyOrig = energy
      if energy == "7P8TeV":
        energy = "8TeV"
      elif not (energy == "8TeV" or energy == "7TeV" or energy == "14TeV"):
        print("Error:Trying to make XS*BR limits w.r.t. higgs Mass, don't recognize energy '"+energy+"' for filename: "+fname)
        print("Exiting.")
        sys.exit(1)
      xsForPoint = xsec['ggHmumu'+xNum+'_'+energy]
      xsForPoint += xsec['vbfHmumu'+xNum+'_'+energy]
      xsForPoint += xsec['whHmumu'+xNum+'_'+energy]
      xsForPoint += xsec['zhHmumu'+xNum+'_'+energy]
      outStr = thisPoint[0]+" "
      for i in range(1,7):
        thisPoint[i] = str(float(thisPoint[i])*xsForPoint)
        outStr += "{0:10.5f} ".format(float(thisPoint[i]))
      outStr += energyOrig
      print outStr
    #print thisPoint
    result.append(thisPoint)
      
  return result

class RelativePlot:
  def __init__(self,dataPoints, canvas, legend, caption, ylabel="95% CL Limit on #sigma/#sigma_{SM} (H #rightarrow #mu^{+}#mu^{-})", xlabel="Integrated Luminosity [fb^{-1}]",caption2="",caption3="",ylimits=[],xlimits=[],vertLines=[],showObs=False,energyStr="8TeV",caption4=""):
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
          twoSigGraph.GetYaxis().SetRangeUser(0.0,ymax*1.65)
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
    leg = root.TLegend(*legPos)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    leg.AddEntry(obsGraph,"Observed Limit","lp")
    leg.AddEntry(expGraph,"Median Expected Limit","lp")
    leg.AddEntry(oneSigGraph,"#pm1 #sigma Expected Limit","f")
    leg.AddEntry(twoSigGraph,"#pm2 #sigma Expected Limit","f")
    self.legPos = legPos
    self.leg = leg
    leg.Draw()

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

class CutOptPlots:
  def __init__(self,globName):
    self.goodDrawn=False
    self.histList = []
    self.graphList = []
    self.globName = globName
    self.canvas = root.TCanvas("CutOptCanvas")
    self.tlatex = root.TLatex()
    self.tlatex.SetNDC()
    self.tlatex.SetTextFont(root.gStyle.GetLabelFont())
    self.tlatex.SetTextSize(0.04)
    #root.gStyle.SetPalette(53)
    #root.gStyle.SetPaintTextFormat(".2f")
    root.gStyle.SetPaintTextFormat(".1f")
    files = glob.glob(globName)
    self.energyStr = ""
    self.data = {}
    for f in sorted(files):
      data =  getData(f)
      if len(data)==0:
        continue
      data =  data[0]
      data = [float(x) for x in data]
      metamatch = re.search(r"_([0-9P]+TeV)_",f)
      if metamatch:
        self.energyStr = metamatch.group(1)
        #print("Energy Str "+self.energyStr)
      else:
        print("Error: couldn't match energy string: "+f)
      f = re.sub(".*/","",f)
      f = re.sub("_8TeV.*","",f)
      f = re.sub("_7TeV.*","",f)
      f = re.sub("_7P8TeV.*","",f)
      match = re.match(r"([a-zA-Z0-9]+)_.*",f)
      assert(match)
      testName = match.group(1)
      match = re.findall(r"_([a-zA-Z0-9]+[GLS])([0-9p-]+)",f)
      assert(match)
      limit = data[4]
      #limitErr = abs(limit-data[3])
      #limitErr = max(abs(limit-data[5]),limitErr)
      #print("Relative Limit Error: {0:.1%}".format(limitErr/limit))
      tmpDict= {}
      tmpDict['limit'] = limit
      for i in match:
        val = float((i[1]).replace('p','.'))
        tmpDict[i[0]] = val
      if not self.data.has_key(testName):
        self.data[testName] = []
      self.data[testName] += [tmpDict]
    for c in self.data:
      d = self.data[c]
      d.sort(key=lambda x: x['limit'])
    if self.energyStr == "7P8TeV":
      self.energyStrWrite = "#sqrt{s} = 7 & 8 TeV"
    else:
      self.energyStrWrite = "#sqrt{s} = "+re.sub("TeV"," TeV",self.energyStr)

  def printBest(self,N):
    for c in self.data:
      print("The Best {1} For Test: {0}".format(c,N))
      for i in range(min(len(self.data[c]),N)):
        d = self.data[c][i]
        printStr = "  {limit:.2f}"
        for key in d:
          if key == 'limit':
            continue
          printStr += " "+key +' {'+key+':.1f}'
        print(printStr.format(**d))

  def plot1D(self,dataName,xName,holdConstDict,rootHistParamList):
    self.goodDrawn = False
    try:
      yMax = -1e6
      yMin = 1e6
      self.canvas.cd()
      graph = root.TGraph()
      graph.SetMarkerColor(1)
      #graph.SetMarkerSize(2*graph.GetMarkerSize())
      graph.SetLineColor(1)
      graph.SetLineWidth(2)
      data = self.data[dataName]
      iPoint = 0
      data.sort(key=lambda x: x[xName])
      for d in data:
        doContinue = False
        for sliceVar in holdConstDict:
          if not d.has_key(sliceVar):
            break
          if d[sliceVar] != holdConstDict[sliceVar]:
            doContinue = True
            break
        if doContinue:
          break
        x = d[xName]
        limit = d['limit']
        yMax = max(limit,yMax)
        yMin = min(limit,yMin)
        #print("{0}: {1:4.0f} limit: {2:4.1f}".format(xName,x,limit))
        graph.SetPoint(iPoint,x,limit)
        iPoint += 1
      xTitle = xName
      yTitle = "Median Expected Limit on #sigma/#sigma_{SM}"
      if xTitle[-1] == 'G':
          xTitle = xTitle[:-1]+" > X"
      elif xTitle[-1] == 'L':
          xTitle = xTitle[:-1]+" < X"
      else:
          raise
      graph.Draw('AC')
      xMin = graph.GetXaxis().GetXmin()
      xMax = graph.GetXaxis().GetXmax()
      padding = yMax-yMin
      yMin -= padding
      yMax += padding
      yMin = math.floor(yMin)
      yMax = math.ceil(yMax)
      hist = root.TH2F(dataName+"dummy","",1,xMin,xMax,1,yMin,yMax)
      setHistTitles(graph,xTitle,yTitle)
      setHistTitles(hist,xTitle,yTitle)
      hist.Draw()
      graph.Draw('C')
      #graph.Draw('P')
      self.graphList += [graph]
      self.histList += [hist]
      self.tlatex.SetTextAlign(12)
      self.tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
      self.tlatex.SetTextAlign(13)
      self.tlatex.DrawLatex(gStyle.GetPadLeftMargin()+0.03,1.0-gStyle.GetPadTopMargin()-0.02,self.energyStrWrite)
      self.tlatex.DrawLatex(gStyle.GetPadLeftMargin()+0.03,1.0-gStyle.GetPadTopMargin()-0.06,"L = {0:.1f} fb^{{-1}}".format(float(lumiDict[self.energyStr])))
      self.goodDrawn = True
    except LookupError as e:
      print("Warning: could not draw {0}, because key not found".format(e))

  def plot2D(self,dataName,xName,yName,holdConstDict,rootHistParamList):
    self.goodDrawn = False
    try:
      self.canvas.cd()
      hist = root.TH2F(*rootHistParamList)
      #hist.SetMarkerColor(0)
      #hist.SetMarkerSize(2*hist.GetMarkerSize())
      data = self.data[dataName]
      #print data
      for d in data:
        doContinue = False
        for sliceVar in holdConstDict:
          if not d.has_key(sliceVar):
            break
          if d[sliceVar] != holdConstDict[sliceVar]:
            doContinue = True
            break
        if doContinue:
          break
        x = d[xName]+0.00001
        y = d[yName]+0.00001
        limit = d['limit']
        ix = hist.GetXaxis().FindBin(x)
        iy = hist.GetYaxis().FindBin(y)
        hist.SetBinContent(ix,iy,limit)
        #hist.Fill(x,y)
        #print x,y,limit, d['ptMissL']
      xTitle = xName
      yTitle = yName
      if xTitle[-1] == 'G':
          xTitle = xTitle[:-1]+" > X"
      elif xTitle[-1] == 'L':
          xTitle = xTitle[:-1]+" < X"
      else:
          raise
      if yTitle[-1] == 'G':
          yTitle = yTitle[:-1]+" > Y"
      elif yTitle[-1] == 'L':
          yTitle = yTitle[:-1]+" < Y"
      else:
          raise
      setHistTitles(hist,xTitle,yTitle)
      hist.Draw('coltext')
      self.histList += [hist]
      self.tlatex.SetTextAlign(12)
      self.tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
      self.tlatex.DrawLatex(gStyle.GetPadLeftMargin(),
            gStyle.GetPadBottomMargin()*0.3,
            self.energyStrWrite+", "+"L = {0:.1f} fb^{{-1}}".format(float(lumiDict[self.energyStr]))
            )
      self.goodDrawn = True
    except LookupError as e:
      print("Warning: could not draw {0}, because key not found".format(e))

  def annotatePlot(self,text):
    if self.goodDrawn:
      optPlots.tlatex.SetTextAlign(32)
      optPlots.tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,text)
  def save(self,fn):
    if self.goodDrawn:
      saveAs(self.canvas,fn)

if __name__ == "__main__":

  dirName = "statsInput/"
  #dirName = "/data/uftrig01b/digiovan/baselinePP/exppluspolin_fixEff_DiffMeans_unbinned_loosePUID_fixBkgVBFFit_syst_22Jan2013ReReco_assProd/statsInputUltimate/"
  #dirName = "/data/uftrig01b/kropiv/HiggsMuMu/CMSSW_6_1_1/Results/hmumuFinalAnalysis/statsInput/"
  #dirName = "/data/uftrig01b/digiovan/baselinePP/m110to160_pixelLumi/hmumuFinalAnalysis/statsInput/"
  
  outDir = "statsOutput/"
  
  root.gErrorIgnoreLevel = root.kWarning
  root.gROOT.SetBatch(True)
  setStyle()
  if args.cutOpt:
    for period in ["7TeV","8TeV"]:
      optPlots = CutOptPlots(dirName+"*"+period+"*.txt.out")
      optPlots.printBest(10)
      tlatex = root.TLatex()

      dataName = 'NonVBFCutOpt'
      xName = 'dimuonPtG'
      holdConstDict = {}
      rootHistParamList = [dataName,'',15,0.,75.]
      optPlots.plot1D(dataName,xName,holdConstDict,rootHistParamList)
      optPlots.annotatePlot("Non-VBF")
      optPlots.save(outDir+dataName+"_"+period)


      dataName = 'BDTCutOptPtMissL40'
      xName = 'bdtVBFG'
      holdConstDict = {}
      rootHistParamList = [dataName,'',20,-0.5,0.5]
      optPlots.plot1D(dataName,xName,holdConstDict,rootHistParamList)
      optPlots.annotatePlot("VBF BDT p_{T}^{Miss}<40 GeV")
      optPlots.save(outDir+dataName+"_"+period)

      dataName = 'VBFCutBasedOptPtMissL40'
      xName = 'dijetMassG'
      yName = 'deltaEtaJetsG'
      holdConstDict = {}
      rootHistParamList = [dataName,'',8,300.,700.,4,3.,5.]
      optPlots.plot2D(dataName,xName,yName,holdConstDict,rootHistParamList)
      optPlots.annotatePlot("VBF Cut Based, p_{T}^{Miss}<40 GeV/c")
      optPlots.save(outDir+dataName+"_"+period)

    sys.exit(0)

  canvas = root.TCanvas()
  if (not args.bdtCut) and (not args.higgsMass):
    canvas.SetLogx(1)
    canvas.SetLogy(1)
  
  #mpl.rcParams["font.family"] = "sans-serif"
  #print mpl.rcParams["backend"]

  ylimits=[0.1,100.0]
  ylimits=[0.1,60.0]

  lumisToUse={"7TeV":lumiDict["7TeV"],"8TeV":lumiDict["8TeV"],"7P8TeV":lumiDict["8TeV"]+lumiDict["7TeV"]}

  if args.printLimits:
    fnToGlob = dirName+"*_*TeV_*.txt.out"
    allfiles = glob.glob(fnToGlob)
    fLens = [len(re.sub(".*/","",f)) for f in allfiles]
    maxLen = max(fLens)+1
    maxLen = str(maxLen)
    print(("\n\n{0:"+maxLen+"}  {1:>4}  {2:>4}  {3:>4} {4:>4} {5:>4} {6:>4}").format("file","obs","exp","-1s","+1s","-2s","+2s"))
    for f in sorted(allfiles):
      data =  getData(f)
      if len(data)==0:
        continue
      data =  data[0]
      data = [float(x) for x in data]
      f = re.sub(".*/","",f)
      f += ":"
      print(("{0:"+maxLen+"}  {1:4.1f}  {2:4.2f}  {3:4.1f} {4:4.1f} {5:4.1f} {6:4.1f}").format(f,data[1],data[4],data[3],data[5],data[2],data[6]))
      #massStr = f[20:23]
      #exp = float(data[4])
      #obs = float(data[1])
      #p1sig = float(data[5])-exp
      #p2sig = float(data[6])-exp
      #m1sig = exp-float(data[3])
      #m2sig = exp-float(data[2])
      #print(("{0:4} & {1:4.1f} & {2:4.1f} & ${{}}^{{+{3:4.1f}}}_{{-{4:3.1f}}}$ & ${{}}^{{+{5:4.1f}}}_{{-{6:4.1f}}}$ \\\\\n\\hline").format(massStr,exp,obs,p1sig,m1sig,p2sig,m2sig))
    print
    sys.exit(0)
  
  for period in ["7TeV","8TeV","14TeV","7P8TeV"]:
    fnToGlob = dirName+"*_"+period+"_*.txt.out"
    if args.xs:
      fnToGlob = dirName+"CombSplitAll_"+period+"_*.txt.out"
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
      xMax = 1e20
      xMin = -1e20
      if args.higgsMass:
        xMax = 150.
        xMin = 120.
      data = getData(dirName+plotName+"_"+energyStr+"_*.txt.out",xMax=xMax,xMin=xMin)
      if plotName == "CombSplitAll" and period=="7P8TeV":
        pklFile = open("limitsCombSplitAll7P8TeV.pkl","w")
        cPickle.dump(data,pklFile)
        pklFile.close()
      xlimits = []
      vertLines = []
      if len(data)<=1:
        continue
      xlabel="Integrated Luminosity [fb^{-1}]"
      ylabel="95% CL_{s} Upper Limit on #sigma/#sigma_{SM} (H #rightarrow #mu^{+}#mu^{-})"
      if args.xs:
        ylabel="95% CL_{s} Upper Limit on #sigma #times BR (H #rightarrow #mu^{+}#mu^{-}) [pb]"
        if energyStr == "7P8TeV":
          ylabel="95% CL_{s} Limit on #sigma(8 TeV) #times BR (H #rightarrow #mu^{+}#mu^{-}) [pb]"
      caption3 = ""
      caption4 = ""
      if args.bdtCut:
        xlimits = [-0.9,1.0]
        xlabel="BDT Discriminant Cut"
        match = re.match(r"([^0-9.]*)([0-9.]*)",plotName)
        assert(match)
        caption3 = "L = {0:.1f} fb^{{-1}}".format(float(match.group(2)))
        plotName = match.group(1)
        if plotName == "IncPtCut":
          xlabel="p_{T}(#mu#mu) Cut [GeV/c]"
          xlimits = []
          vertLines += [10.0]
          if energyStr == "8TeV":
            ylimits = [0.,16.]
          elif energyStr == "7TeV":
            ylimits = [0.,32.]
        elif plotName == "VBFmDiJetCut":
          xlabel="m_{jj} Cut [GeV/c]"
          xlimits = []
          vertLines += [550.0]
          if energyStr == "8TeV":
            ylimits = [0.,50.]
          elif energyStr == "7TeV":
            ylimits = [0.,50.]
        elif plotName == "VBFdeltaEtaJetCut":
          xlabel="#Delta#eta(jj) Cut"
          xlimits = []
          vertLines += [3.0]
          if energyStr == "8TeV":
            ylimits = [0.,50.]
          elif energyStr == "7TeV":
            ylimits = [0.,50.]
        else:
          if energyStr == "8TeV":
            ylimits = [0.,35.]
            #vertLines += [-0.04]
            vertLines += [0.0]
          elif energyStr == "7TeV":
            ylimits = [0.,70.]
            vertLines += [0.0]
      elif args.higgsMass:
        if energyStr == "8TeV":
            caption2 = "#sqrt{{s}} = 8 TeV L = {0:.1f} fb^{{-1}}".format(float(lumiDict[energyStr]))
            caption3 = ""
        elif energyStr == "7TeV":
            caption2 = "#sqrt{{s}} = 7 TeV L = {0:.1f} fb^{{-1}}".format(float(lumiDict[energyStr]))
            caption3 = ""
        elif energyStr == "7P8TeV":
            caption2 = "#sqrt{{s}} = 7 TeV L = {0:.1f} fb^{{-1}}".format(float(lumiDict["7TeV"]))+", "+ "#sqrt{{s}} = 8 TeV L = {0:.1f} fb^{{-1}}".format(float(lumiDict["8TeV"]))
            caption3 = ""
        ylimits = []
        xlabel="m_{H} [GeV/c^{2}]"
      #elif period == "14TeV":
      #  title = "Standard Model H#rightarrow#mu#mu"
      title = titleMap[plotName]
      #title = "GF Only "+title
      #ylabel="95% CL Limit on #sigma/#sigma_{SM} (GF H#rightarrow#mu#mu)"
      #caption4 = "#sigma_{VBF}= #sigma_{WH}= #sigma_{ZH}= 0"
      #title = "VBF Only "+title
      #ylabel="95% CL Limit on #sigma/#sigma_{SM} (VBF H#rightarrow#mu#mu)"
      #caption4 = "#sigma_{GF}= #sigma_{WH}= #sigma_{ZH}= 0"
      showObs = args.higgsMass and not args.blind
      incPlot = RelativePlot(data,canvas,legend,title,caption2=caption2,ylimits=ylimits,energyStr=energyStrWrite,xlabel=xlabel,caption3=caption3,showObs=showObs,xlimits=xlimits,vertLines = vertLines,ylabel=ylabel,caption4=caption4)
      saveAs(canvas,outDir+plotName+"_"+energyStr)
      #saveAs(canvas,outDir+"gfOnly_"+plotName+"_"+energyStr)
      #saveAs(canvas,outDir+"vbfOnly_"+plotName+"_"+energyStr)
      if args.xs:
        saveAs(canvas,outDir+"xsbr_"+plotName+"_"+energyStr)

    if not mplGood:
        continue

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
      comparePlot = ComparePlotTable(compareData,titleMap=tmpMap,vertLine1=False,anotation1=anotation1,anotation2=anotation2,signalInject=args.signalInject)
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
      comparePlot = ComparePlotTable(compareData,titleMap=tmpMap,vertLine1=False,anotation1=anotation1,anotation2=anotation2,signalInject=args.signalInject)
      comparePlot.save(outDir+"compareFinal"+"_"+energyStr)

    ## All Categories

    mustBe = r"(.+)_(.+)_[.\d]+.txt.out"
    veto = ["EE",'BB',"BO","BE","OE","OO"]
    compareData = getData(fnGlobStr,matchString=mustBe,dontMatchStrings=veto,doSort=False)
    comparePlot = ComparePlotTable(compareData,titleMap={},vertLine1=False,anotation1=anotation1,anotation2=anotation2,signalInject=args.signalInject)
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
      comparePlot = ComparePlotTable(compareData,titleMap={},vertLine1=False,anotation1=anotation1,anotation2=anotation2,signalInject=args.signalInject)
      comparePlot.save(outDir+"IncBDTVPresel"+"_"+energyStr)
