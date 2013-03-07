#!/usr/bin/env python

from helpers import *
import ROOT as root
import glob
import re
import os.path
import random

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

def getDataGOF(basename,outDir):
  def getPValue(filename):
    filename = re.sub(r"\.txt.*","",filename)
    outFilename = filename
    outFilename = os.path.basename(outFilename)
    endStr = [".txt.CCC.root",".txt.CCC-Toys-*.root"]
    observed = None
    nToys = None
    nNum = 0
    tmpFiles = []
    tmpTrees = []
    tmpHists = []
    for postfix in endStr:
      tree = root.TChain("limit")
      for fname in glob.glob(filename+postfix):
        tree.AddFile(fname)
      nEntries = tree.GetEntries()
      if nEntries > 1:
        nToys = nEntries
      for i in range(nEntries):
        tree.GetEntry(i)
        toy = tree.limit
        if nEntries == 1:
          observed = toy
        elif toy > observed:
          nNum += 1
      tmpFiles += [f]
      tmpTrees += [tree]
    result = float(nNum)/nToys
    print("{} obs: {:.2f} nToys: {} result: {:.2f}".format(filename,observed,nToys,result))
    ## Now Drawing Distribution
    canvas = root.TCanvas("tmpCanvas"+str(random.randint(0,10000)))
    histName = "tmpHist"+str(random.randint(0,10000))
    tmpHist = None
    first = True
    for tree in tmpTrees:
      if first:
        first = False
        continue
      if tmpHist == None:
        tree.Draw("limit >> "+histName)
        tmpHist = root.gDirectory.Get(histName)
        tmpHist.SetTitle(filename)
        binWidthStr = getBinWidthStr(tmpHist)
        setHistTitles(tmpHist,"Channel Compatability Test Statistic","Toys/"+binWidthStr)
        tmpHist.SetLineWidth(2)
      else:
        tree.Draw("limit >>+ "+histName)
    tmpHist.Draw()

    # Observed Arrow
    yAxis = tmpHist.GetYaxis()
    arrowLength = (1.0-2e-3)/10.
    arrowYmax = tmpHist.GetMaximum()
    arrowHeadSize = 0.025
    arrow = root.TLine(observed,0.0,observed,arrowYmax
                                )
    arrow.SetLineWidth(3)
    arrow.SetLineColor(root.kRed)
    #arrow.SetFillColor(root.kRed)
    #arrow.SetAngle(40)
    arrow.Draw()

    label = root.TLatex()
    label.SetNDC()
    label.SetTextFont(root.gStyle.GetLabelFont("X"))
    label.SetTextSize(root.gStyle.GetLabelSize("X"))
    label.SetTextAlign(23)
    label.DrawLatex(
            (1.0+root.gStyle.GetPadLeftMargin()-root.gStyle.GetPadRightMargin())/2.,
            1.0-root.gStyle.GetPadTopMargin()-0.02,
            "Observed: {0:.2f}, nToys: {1}, p-Value: {2:.2g}".format(
                        observed,nToys,result)
            )
    saveAs(canvas,outDir+"testStat_"+outFilename)
    return result

  result = {}

  result['sig'] = []
  files = glob.glob(basename+"*.CCC.root")
  for f in files:
    massMatch = re.match(r".+_.+TeV_([0-9.]+).txt.*",f)
    assert(massMatch)
    mass = massMatch.group(1)
    result["sig"] += [[float(mass),getPValue(f)]]
  return result

class PValuePlotTogether:
  def __init__(self,dataDict, canvas, caption="Standard Model H#rightarrow#mu#mu", ylabel="Channel Compatability p-Value", xlabel="m_{H} [GeV/c^{2}]",caption2="",caption3="",ylimits=[],xlimits=[],energyStr="8TeV"):
    graphs = []
    ymax = 1.0
    ymin = 1e20
    ymin = 2e-3 # To see 3sigma line
    xmin = 1e20
    xmax = -1e20
    sortedChannels = sorted(dataDict.keys())
    sortedChannels.reverse()
    for channel in sortedChannels:
      dataPoints = dataDict[channel]['sig']
      dataPoints.sort(key=lambda x:x[0])
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

    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin()-0.03,0.88,caption2)
    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin()-0.03,0.82,caption3)

    #tlatex.SetTextAlign(12)
    #tlatex.SetTextSize(0.05)
    #tlatex.DrawLatex(gStyle.GetPadLeftMargin()+0.03,gStyle.GetPadBottomMargin()+0.23,"p-Value of #sigma/#sigma_{SM}=1")
    #tlatex.DrawLatex(gStyle.GetPadLeftMargin()+0.03,gStyle.GetPadBottomMargin()+0.17,"Hypothesis")

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

class ChannelPlot:
  def __init__(self,filename,outname):
    self.canvas = root.TCanvas()
    getattr(self,"getData")(filename)
    getattr(self,"plot")(outname)

  def getData(self,filename):
    f = root.TFile(filename)
    fit_nom = f.Get("fit_nominal")
    fpf_nom = fit_nom.floatParsFinal()
    fit_alt = f.Get("fit_alternate")
    fpf_alt = fit_alt.floatParsFinal()

    r_nom = fpf_nom.find("r")
    data = {'r':[r_nom.getVal(),-r_nom.getErrorLo(),r_nom.getErrorHi()]}

    xMin = r_nom.getVal()-r_nom.getErrorLo()
    xMax = r_nom.getVal()+r_nom.getErrorHi()

    self.differentEnergies = False

    for i in range(fpf_alt.getSize()):
      param = fpf_alt.at(i)
      name = param.GetName()
      match = re.match(r".*ChannelCompatibilityCheck.*_r_(.*)([\d])TeV",name)
      energy = None
      if match:
        channel = match.group(1)+match.group(2)
        tmpEnergy = match.group(2)
        data[channel] = [param.getVal(),-param.getErrorLo(),param.getErrorHi(),match.group(1),match.group(2)]
        print channel
        print data[channel]
        high = data[channel][0]+data[channel][2]
        low = data[channel][0]-data[channel][1]
        if high > xMax:
            xMax = high
        if low < xMin:
            xMin = low


        if energy == None:
            energy = tmpEnergy
        else:
          if energy != tmpEnergy:
            self.differentEnergies = True
        
    self.data = data
    self.channels = self.data.keys()
    self.channels.pop(self.channels.index("r"))
    self.nChannels = len(self.channels)
    self.xMax = xMax
    self.xMin = xMin

  def plot(self,savename):
    self.canvas.cd()
    self.canvas.SetLeftMargin(0.3)
    frame = root.TH2F("frame","",1,self.xMin*0.5,self.xMax*1.1,self.nChannels+1,0,self.nChannels+1)
    if self.differentEnergies:
      frame.GetYaxis().SetLabelSize(frame.GetYaxis().GetLabelSize()*1.5)
    else:
      frame.GetYaxis().SetLabelSize(frame.GetYaxis().GetLabelSize()*1.75)
    frame.GetXaxis().SetLabelSize(frame.GetXaxis().GetLabelSize()*1.4)
    frame.GetXaxis().SetTitleSize(frame.GetXaxis().GetTitleSize()*1.2)
    frame.GetXaxis().SetTitleOffset(frame.GetXaxis().GetTitleOffset()*0.8)
    setHistTitles(frame,"Best Fit #sigma/#sigma_{SM}","")
    points = root.TGraphAsymmErrors()
    points.SetMarkerStyle(21)
    points.SetMarkerColor(1)
    points.SetMarkerSize(2)
    points.SetLineColor(root.kRed)
    points.SetLineWidth(5)
    
    for i in range(self.nChannels):
        datapoint = self.data[self.channels[i]]
        points.SetPoint(i,datapoint[0],i+0.5)
        points.SetPointError(i,datapoint[1],datapoint[2],0.,0.)
        extra = ""
        if self.differentEnergies:
            extra = " {0} TeV".format(datapoint[4])
        frame.GetYaxis().SetBinLabel(i+1,titleMap[datapoint[3]]+extra)

    globalFitVal = self.data['r'][0]
    globalFitP1Sig = self.data['r'][2]+globalFitVal
    globalFitM1Sig = -self.data['r'][1]+globalFitVal

    globalFitBand = root.TBox(globalFitM1Sig, 0., globalFitP1Sig, self.nChannels);
    globalFitBand.SetFillColor(root.kGreen);

    globalFitLine = root.TLine(globalFitVal, 0., globalFitVal, self.nChannels);
    globalFitLine.SetLineWidth(5);
    globalFitLine.SetLineColor(1);

    frame.Draw()
    globalFitBand.Draw()
    globalFitLine.Draw()
    points.Draw("PZ")

    self.canvas.RedrawAxis() # Updates Axis Lines

    saveAs(self.canvas,savename)

if __name__ == "__main__":

  dirName = "statsInput/"
  
  outDir = "statsOutput/"
  
  root.gErrorIgnoreLevel = root.kWarning
  root.gROOT.SetBatch(True)
  setStyle()
  canvas = root.TCanvas()
  
  mpl.rcParams["font.family"] = "sans-serif"
  #print mpl.rcParams["backend"]

  ylimits=[]

  lumisToUse={"7TeV":lumiDict["7TeV"],"8TeV":lumiDict["8TeV"],"7P8TeV":lumiDict["8TeV"]+lumiDict["7TeV"]}
  
  for period in ["7TeV","8TeV","14TeV","7P8TeV"]:
    fnToGlob = dirName+"*_"+period+"_*.txt.*root"
    allfiles = glob.glob(fnToGlob)

    """
    ## Limit v. Lumi
    energyStr = ""
    plots = set()
    for fn in allfiles:
      match = re.search(r".*/(.+)_(.+)_[.\d]+.txt.*root",fn)
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
      if len(plots)==0:
        continue
      pValueDict = {}
      for plotName in plots:
        if plotName in vetos:
            continue

        data = getDataGOF(dirName+plotName+"_"+period,outDir)
        pValueDict[plotName] = data
      pValueAllPlot = PValuePlotTogether(pValueDict,canvas,caption2=caption2,caption3=caption3,energyStr=energyStr)
      saveAs(canvas,outDir+"ccc_"+saveName+period)
    canvas.SetLogy(0)
    """

    compareFileGlob = "BDTCutCatVBFBDTOnly*"+period+"*125.0.txt.CCC.root"
    compareFiles = glob.glob(dirName+compareFileGlob)
    for f in compareFiles:
      ChannelPlot(f,outDir+"cccCompare_"+period)
    
