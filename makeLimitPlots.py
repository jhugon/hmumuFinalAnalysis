#!/usr/bin/env python

import optparse
parser = optparse.OptionParser(description="Makes Limit Plots from output text from combine tool.")
parser.add_option("--bdtCut", help="Makes plots v. BDT Cut Instead of Luminosity",action="store_true",default=False)
parser.add_option("-m","--higgsMass", help="Makes plots v. Higgs Mass",action="store_true",default=False)
parser.add_option("--signalInject", help="Sets a caption saying that signal was injected with strength",type=float,default=0.0)
parser.add_option("--signalInjectMass", help="Mass For Injected Signal",type=float,default=125.0)
args, fakeargs = parser.parse_args()

from helpers import *
import ROOT as root
import glob
import re
import os.path

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
  "VBFmDiJetCut":"VBF Cut-Based",
  "VBFdeltaEtaJetCut":"VBF Cut-Based",

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
    expGraph.SetMarkerStyle(20)
    expGraph.SetMarkerSize(0.9)
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
    if len(ylimits)==2:
        twoSigGraph.GetYaxis().SetRangeUser(*ylimits)
    else:
        twoSigGraph.GetYaxis().SetRangeUser(0.0,ymax*1.1)
    if len(xlimits)==2:
        twoSigGraph.GetXaxis().SetRangeUser(*xlimits)
    oneSigGraph.Draw("3")
    expGraph.Draw("l")
    expGraph.Draw("p")
    if showObs:
      obsGraph.Draw("l")
      obsGraph.Draw("P")

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
  
  #mpl.rcParams["font.family"] = "sans-serif"
  #print mpl.rcParams["backend"]

  ylimits=[0.1,100.0]
  ylimits=[0.1,60.0]

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
      xlimits = []
      vertLines = []
      if len(data)<=1:
        continue
      xlabel="Integrated Luminosity [fb^{-1}]"
      caption3 = ""
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
            ylimits = [0.,40.]
        elif energyStr == "7TeV":
            ylimits = [0.,70.]
        ylimits = []
        xlabel="m_{H} [GeV/c^{2}]"
        caption3 = "L = {0:.1f} fb^{{-1}}".format(float(lumisToUse[energyStr]))
      #elif period == "14TeV":
      #  title = "Standard Model H#rightarrow#mu#mu"
      title = titleMap[plotName]
      incPlot = RelativePlot(data,canvas,legend,title,caption2=caption2,ylimits=ylimits,energyStr=energyStrWrite,xlabel=xlabel,caption3=caption3,showObs=args.higgsMass,xlimits=xlimits,vertLines = vertLines)
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
