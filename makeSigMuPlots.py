#!/usr/bin/env python

import optparse
parser = optparse.OptionParser(description="Makes cards for use in the CMS Combine tool.")
parser.add_option("--signalInject", help="Inject Signal with Strength into data_obs",type=float,default=0.0)
parser.add_option("--signalInjectMass", help="Mass For Injected Signal",type=float,default=125.0)
parser.add_option("-m","--higgsMass", help="Makes plots v. Higgs Mass",action="store_true",default=False)
parser.add_option("-p","--pValue", help="Makes p-value plots instead of significance",action="store_true",default=False)
parser.add_option("-e","--expected", help="Makes Expected sig/mu plots v. lumi",action="store_true",default=False)
args, fakeargs = parser.parse_args()

from helpers import *
import ROOT as root
import glob
import re
import os.path

from xsec import *

import numpy

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
                        
  "Jets01SplitCatAll": "0,1-Jet Combination",


  "Jet2CutsVBFPass":"2-Jet VBF Tight",
  "Jet2CutsGFPass":"2-Jet GF Tight",
  "Jet2CutsFailVBFGF":"2-Jet VBF Loose",

  "Jet2SplitCutsGFSplit" : "2-Jet Combination",
  "CombSplitAll" : "Combination",
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
  "BDTCutCatVBFBDTOnly": 1,

  "Jets01PassPtG10BB": root.kCyan,
  "Jets01PassPtG10BO": root.kOrange,
  "Jets01PassPtG10BE": root.kGreen,
  "Jets01PassPtG10OO": root.kOrange+3,
  "Jets01PassPtG10OE": root.kMagenta,
  "Jets01PassPtG10EE": root.kViolet+2,
  "Jets01PassCatAll" : root.kBlue+3,
                                      
  "Jets01FailPtG10BB": root.kCyan,
  "Jets01FailPtG10BO": root.kOrange,
  "Jets01FailPtG10BE": root.kGreen,
  "Jets01FailPtG10OO": root.kOrange+3,
  "Jets01FailPtG10OE": root.kMagenta,
  "Jets01FailPtG10EE": root.kViolet+2,
  "Jets01FailCatAll" : root.kBlue+3,
                                      
  "Jets01SplitCatAll": root.kBlue,

  "Jet2CutsVBFPass":root.kGreen,
  "Jet2CutsGFPass": root.kBlue+3,
  "Jet2CutsFailVBFGF":root.kOrange+3,

  "Jet2SplitCutsGFSplit" : root.kRed,
  "CombSplitAll" : 1,
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
    obs = -15.0
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
      if xNum == -10.:
        match = re.search(matchString+r"\.sigexp",expFname)
        if match:
          xNum = match.group(1)
    except Exception:
      #print("Expected Significance Not Found: "+expFname)
      pass
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

def getDataMu2(fileString,matchString=r".*_([-\d.]+)\.txt\.mu\.root",dontMatchStrings=[],doSort=True):
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
    xNum = -10.0
    match = re.match(matchString,fname)
    if match:
      xNum = match.group(1)
    else:
      continue
    #print "x: ",xNum
    tmpF = root.TFile(fname)
    avgMu = 0.
    muSigma = 0.
    avgErrMu = 0.
    errMuSigma = 0.
    tmpTree = tmpF.Get("limit")
    nEntries = tmpTree.GetEntries()
    if nEntries == 0:
        continue
    nToys = nEntries/4
    nToysReal = 0.
    muVals = []
    muErrVals = []
    for iToy in range(nToys):
      iEntry = iToy*4+1
      tmpTree.GetEntry(iEntry)
      if tmpTree.limit == -50.:
        continue
      iEntry = iToy*4+3
      tmpTree.GetEntry(iEntry)
      avgMu += tmpTree.limit
      avgErrMu += tmpTree.limitErr/tmpTree.limit
      muVals.append(tmpTree.limit)
      muErrVals.append(tmpTree.limitErr/tmpTree.limit)
      #print("Toy: {0}, Mu: {1:.2f} +/- {2:.2f}".format(iToy,tmpTree.limit,tmpTree.limitErr))
      nToysReal += 1.
    avgMu /= float(nToysReal)
    avgErrMu /= float(nToysReal)
    for iToy in range(nToys):
      iEntry = iToy*4+1
      tmpTree.GetEntry(iEntry)
      if tmpTree.limit == -50.:
        continue
      iEntry = iToy*4+3
      tmpTree.GetEntry(iEntry)
      muSigma += (tmpTree.limit-avgMu)**2
      errMuSigma += (tmpTree.limitErr/tmpTree.limit-avgErrMu)**2
    muSigma /= float(nToysReal-1)
    errMuSigma /= float(nToysReal-1)
    muSigma = sqrt(muSigma)
    errMuSigma = sqrt(errMuSigma)

    medMu = computeMedian(muVals)
    medErrMu = computeMedian(muErrVals)
    #print("Original Values for mu: {0} {1}".format(avgMu,avgErrMu))
    #print("Median Values for mu  : {0} {1}".format(medMu,medErrMu))

    #thisPoint = [xNum,avgMu,muSigma,avgErrMu,errMuSigma]
    thisPoint = [xNum,medMu,muSigma,medErrMu,errMuSigma]
    if thisPoint.count("-10.0")>0:
        continue
    if thisPoint.count(-10.0)>0:
        continue
    #print thisPoint
    result.append(thisPoint)
  return result

class RelativePlot:
  def __init__(self,dataPoints, canvas, legend, caption, ylabel="Error on #sigma/#sigma_{SM}", xlabel="Integrated Luminosity [fb^{-1}]",caption2="",caption3="",ylimits=[],xlimits=[],vertLines=[],showObs=False,energyStr="8TeV",doMuExtraPlot=False,horizLines=[],drawErrBands=True,yToFind=[],xToFind=[]):
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
    if not drawErrBands:
      muExtraGraph.SetLineColor(root.kRed)
      muExtraGraph.SetMarkerColor(root.kRed)
    muExtraGraphErr.SetFillColor(root.kGreen)
    self.muExtraGraph = muExtraGraph
    self.muExtraGraphErr = muExtraGraphErr
    iPoint = 0
    ymax = 1e-20
    ymin = 1e20
    xmax = 1e-20
    xmin = 1e20
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
        high1sig = float(point[3])*100.
        low1sig = float(point[2])*100.
        if high1sig > low1sig:
            value = high1sig
        else:
            value = low1sig
      elif len(point)==5: #mu2
        value = point[3]
        high1sig = point[4]
        low1sig = point[4]
        value *= 100.
        high1sig *= 100.
        low1sig *= 100.
        print("Error on mu for: {0:.2f} fb^-1: {1:.2f} +/- {2:.2f}".format(xNum,value,high1sig))
      else: #sig len=3
        #thisPoint = [xNum,obs,exp]
        value = float(point[2])
        obs = float(point[1])
        if value == -20.0:
            value = obs
      errGraph.SetPoint(iPoint,xNum,value)
      obsGraph.SetPoint(iPoint,xNum,obs)
      if doMuExtraPlot:
        if len(point)==4: #mu
          muExtraGraph.SetPoint(iPoint,xNum,float(point[1]))
          muExtraGraphErr.SetPoint(iPoint,xNum,float(point[1]))
          muExtraGraphErr.SetPointError(iPoint,0.,0.,float(point[2]),float(point[3]))
        elif len(point)==5: #mu2
          muExtraGraph.SetPoint(iPoint,xNum,value)
          muExtraGraphErr.SetPoint(iPoint,xNum,value)
          muExtraGraphErr.SetPointError(iPoint,0.,0.,low1sig,high1sig)
      if xNum < xmin:
        xmin = xNum
      if xNum > xmax:
        xmax = xNum
      if value < ymin:
        ymin = value 
      if value > ymax:
        ymax = value
      iPoint += 1

    xWidth = xmax-xmin
    yWidth = ymax-ymin
    xEdge = xWidth*0.1
    yEdge = yWidth*0.1
    xmin -= xEdge
    xmax += xEdge
    ymin -= yEdge
    ymax += yEdge
    if len(xlimits) == 2:
      xmin = xlimits[0]
      xmax = xlimits[1]
    if len(ylimits) == 2:
      ymin = ylimits[0]
      ymax = ylimits[1]
    self.axisHist = root.TH2F("axisHist","",1,xmin,xmax,1,ymin,ymax)
    setHistTitles(self.axisHist,xlabel,ylabel)

    self.axisHist.Draw()

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

    pairsToPlot = []
    for x in xToFind:
      y = tgraphFindY(errGraph,x)
      pairsToPlot.append((x,y))

    for y in yToFind:
      x = tgraphFindX(errGraph,y)
      pairsToPlot.append((x,y))

    self.line = root.TLine()
    self.line.SetLineWidth(2)
    self.line.SetLineStyle(2)
    self.lines = []
    for pair in pairsToPlot:
      x,y = pair
      self.lines.append(self.line.DrawLine(xmin,y,x,y))
      self.lines.append(self.line.DrawLine(x,ymin,x,y))

    label = root.TLatex()
    #label.SetNDC()
    label.SetTextFont(root.gStyle.GetLabelFont("X"))
    label.SetTextSize(root.gStyle.GetLabelSize("X"))
    label.SetTextColor(root.kRed)
    label.SetTextAlign(22)
    self.label=label

    if doMuExtraPlot:
      if drawErrBands:
        muExtraGraphErr.Draw("3")
        muExtraGraph.Draw("l")
        muExtraGraph.Draw("p")
      else:
        muExtraGraph.Draw("l")
        muExtraGraph.Draw("p")
    else:
      if showObs:
        errGraph.Draw("l")
        obsGraph.Draw("l")
        obsGraph.Draw("p")
      else:
        errGraph.Draw("l")
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

    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(gStyle.GetPadLeftMargin()+0.03,0.88,caption2)
    tlatex.DrawLatex(gStyle.GetPadLeftMargin()+0.03,0.82,caption3)

    for g in self.vertLines:
      g.Draw("l")

    if doMuExtraPlot and drawErrBands:
      for l in self.lines:
        l.Draw()

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

  dirName = "statsInput/"
  
  outDir = "statsOutput/"
  
  root.gErrorIgnoreLevel = root.kWarning
  root.gROOT.SetBatch(True)
  setStyle()
  canvas = root.TCanvas()
  if not args.higgsMass:
    #canvas.SetLogx(1)
    pass
  canvas.SetLogy(0)
  
  ylimits=[]

  lumisToUse={"7TeV":lumiDict["7TeV"],"8TeV":lumiDict["8TeV"],"7P8TeV":lumiDict["8TeV"]+lumiDict["7TeV"]}
  
  for period in ["7TeV","8TeV","14TeV","7P8TeV"]:
    fnToGlob = dirName+"*_"+period+"_*.txt.sig"
    allfiles = glob.glob(fnToGlob)

    ## Limit v. Lumi
    energyStr = ""
    plots = set()
    for fn in allfiles:
      match = re.search(r".*/(.+)_(.+)_[.\d]+.txt.sig",fn)
      badPlot = re.search(r"Silly",fn)
      badPlot2 = re.search(r"Silly",fn)
      if match and not (badPlot or badPlot2):
        plots.add(match.group(1))
        energyStr = match.group(2)

    caption2 = ""
    caption3 = ""
    if energyStr == '':
      pass
    elif energyStr == "7P8TeV":
      caption2 = "#sqrt{{s}} = 7 TeV L = {0:.1f} fb^{{-1}}".format(float(lumiDict["7TeV"]))
      caption3 = "#sqrt{{s}} = 8 TeV L = {0:.1f} fb^{{-1}}".format(float(lumiDict["8TeV"]))
    elif energyStr == "14TeV":
      caption2 = "#sqrt{{s}} = {0}".format(energyStr)
      caption3 = ""
    else:
      caption2 = "#sqrt{{s}} = {0} L = {1:.1f} fb^{{-1}}".format(energyStr,float(lumiDict[energyStr]))
      caption3 = ""

    if energyStr == "7P8TeV":
      energyStr = "7 & 8 TeV"
    else:
      energyStr.replace("TeV"," TeV")
  
    if args.signalInject>0.0:
      caption2 += "Signal Injected {0:.1f}#timesSM".format(args.signalInject)
      caption2 += " m_{H} = "+"{0:.1f}".format(args.signalInjectMass)+" GeV/c^{2}"
    legend = root.TLegend(0.58,0.70,0.9,0.9)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    
    #print "allfiles:"
    #print allfiles
    #print "plots:"
    #print plots
    
    canvas.SetLogy(0)
    for plotName in plots:
      data = None
      doMuExtraPlot=False
      showObs=False
      data = getDataMu(dirName+plotName+"_"+period+"_*.txt*")
      #print "muData:"
      #print data
      if len(data)<=1:
        continue
      title = "Standard Model H#rightarrow#mu#mu"
      xlabel="Integrated Luminosity [fb^{-1}]"
      ytitle = "Error on #sigma/#sigma_{SM} [%]"
      drawErrBands=True
      xlimits=[0.,3200.]
      #xlimits=[0.,1700.]
      ylimits=[0.,100.]
      if args.higgsMass:
        title = titleMap[plotName]
        doMuExtraPlot=True
        ytitle = "Best Fit #sigma/#sigma_{SM}"
        xlabel="m_{H} [GeV/c^{2}]"
        xlimits=[]
        ylimits=[]
      incPlot = RelativePlot(data,canvas,legend,title,caption2=caption2,caption3=caption3,ylabel=ytitle,energyStr=energyStr,doMuExtraPlot=doMuExtraPlot,showObs=showObs,xlabel=xlabel,drawErrBands=drawErrBands,xlimits=xlimits,ylimits=ylimits)
      saveAs(canvas,outDir+"mu"+plotName+"_"+period)

      if args.higgsMass:
        continue

      data = getDataMu2(dirName+plotName+"_"+period+"_*.txt*")
      #print "muData2:"
      #print data
      doMuExtraPlot=True
      drawErrBands=False
      incPlot = RelativePlot(data,canvas,legend,title,caption2=caption2,caption3=caption3,ylabel=ytitle,energyStr=energyStr,doMuExtraPlot=doMuExtraPlot,showObs=showObs,xlabel=xlabel,drawErrBands=drawErrBands,xlimits=xlimits,ylimits=ylimits,xToFind=[300.,3000.])
      saveAs(canvas,outDir+"mu2"+plotName+"_"+period)

      data = getDataSig(dirName+plotName+"_"+period+"_*.txt*")
      #print "SigData:"
      #print data
      if len(data)<=1:
        continue

      title = "Standard Model H#rightarrow#mu#mu"
      xlabel="Integrated Luminosity [fb^{-1}]"
      ytitle = "Expected Significance"
      xlimits=[0.,1700.]
      ylimits=[0.,6.0]
      xlimits=[10.,3200.]
      ylimits=[0.,8.0]
      incPlot = RelativePlot(data,canvas,legend,title,caption2=caption2,caption3=caption3,ylabel=ytitle,energyStr=energyStr,showObs=showObs,xlabel=xlabel,xlimits=xlimits,ylimits=ylimits,yToFind=[3.,5.])
      saveAs(canvas,outDir+"sig"+plotName+"_"+period)

    ## All p-values together plot
    canvas.SetLogy(1)
    pValueVetos = [
        [
            "Jets01PassPtG10BB",
            "Jets01PassPtG10BO",
            "Jets01PassPtG10BE",
            "Jets01PassPtG10OO",
            "Jets01PassPtG10OE",
            "Jets01PassPtG10EE",
            "Jets01PassCatAll" ,
                                                
            "Jets01FailPtG10BB",
            "Jets01FailPtG10BO",
            "Jets01FailPtG10BE",
            "Jets01FailPtG10OO",
            "Jets01FailPtG10OE",
            "Jets01FailPtG10EE",
            "Jets01FailCatAll" ,
                                                
            "Jets01SplitCatAll",
          
            "CombSplitAll" 
        ],
        [
            "Jets01FailPtG10BB",
            "Jets01FailPtG10BO",
            "Jets01FailPtG10BE",
            "Jets01FailPtG10OO",
            "Jets01FailPtG10OE",
            "Jets01FailPtG10EE",
            "Jets01FailCatAll" ,
                                                
            "Jets01SplitCatAll",
          
          
            "Jet2CutsVBFPass",
            "Jet2CutsGFPass",
            "Jet2CutsFailVBFGF",
          
            "Jet2SplitCutsGFSplit" ,
            "CombSplitAll" 
        ],
        [
            "Jets01PassPtG10BB",
            "Jets01PassPtG10BO",
            "Jets01PassPtG10BE",
            "Jets01PassPtG10OO",
            "Jets01PassPtG10OE",
            "Jets01PassPtG10EE",
            "Jets01PassCatAll" ,
                                                
            "Jets01SplitCatAll",
          
          
            "Jet2CutsVBFPass",
            "Jet2CutsGFPass",
            "Jet2CutsFailVBFGF",
          
            "Jet2SplitCutsGFSplit" ,
            "CombSplitAll" 
        ],
        [
            "Jets01PassPtG10BB",
            "Jets01PassPtG10BO",
            "Jets01PassPtG10BE",
            "Jets01PassPtG10OO",
            "Jets01PassPtG10OE",
            "Jets01PassPtG10EE",
            "Jets01PassCatAll" ,
                                                
            "Jets01FailPtG10BB",
            "Jets01FailPtG10BO",
            "Jets01FailPtG10BE",
            "Jets01FailPtG10OO",
            "Jets01FailPtG10OE",
            "Jets01FailPtG10EE",
            "Jets01FailCatAll" ,
                                                
            "Jet2CutsVBFPass",
            "Jet2CutsGFPass",
            "Jet2CutsFailVBFGF"
        ]
    ]
    for saveName,vetos in zip(["VBF","NonVBFTight","NonVBFLoose","Final"],pValueVetos):
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
