#!/usr/bin/env python

from xsec import *
from helpers import *
import ROOT as root
import os
import sys

#dataDir = "input/separateSamplesTrainOnlyVBFLarge110to150/"
dataDir = "input/lowPtCuts/"
#dataDir = "input/jets20f25/"
#dataDir = "input/separateSamplesTrainOnlyVBFLargeBDTG110to150/"
outDir = "output/"

RUNPERIOD="8TeV"
LUMI=lumiDict[RUNPERIOD]

scaleHiggsBy = 1000.

LOGY=False
integralPlot=False
MCErrors=True
#PULLTYPE="adrian1"
PULLTYPE="pullMC"
allHiggsTogether = False
ylimitsRatio = [-5,5]
if PULLTYPE=="adrian1":
  ylimitsRatio = [-1.5,1.5]
elif PULLTYPE=="ratio":
  ylimitsRatio = [0.5,1.5]

#anotateText = "80 GeV < m_{#mu#mu} < 160 GeV; p_{T,#mu#mu}<20 GeV"
#anotateText = "110 GeV < m_{#mu#mu} < 160 GeV; p_{T,#mu#mu}<20 GeV"
anotateText = "110 GeV < m_{#mu#mu} < 150 GeV"
#anotateText = "80 GeV < m_{#mu#mu} < 160 GeV"
#anotateText = "VBF Preselection"

urLegendPos = [0.70,0.67,0.9,0.9]
ulLegendPos = [0.20,0.67,0.4,0.9]
ucLegendPos = [0.46,0.67,0.64,0.9]
lcLegendPos = [0.46,0.35,0.64,0.63]
llLegendPos = [0.20,0.35,0.4,0.63]
ccLegendPos = [0.46,0.47,0.64,0.7]
stdLegendPos = urLegendPos

#histDirs = ["","4GeVWindow/","PtDiMu100/","VBFPresel/","IncPresel/","NotBlindWindow/"]
histDirs = ["NotBlindWindow/"]
histDirs = ["VBFPreselDiMuPtL20/","IncPreselDiMuPtL20/"]
#histDirs = ["VBFPresel/","IncPresel/"]
#histDirs = ["IncPresel/"]
#histDirs = ["IncPreselPtG10BB/","VBFBDTCut/"]
#histDirs = ["VBFBDTCut/"]
#histDirs = ["VBFBDTCut/"]
histDirs = ["","BB/","IncPreselBB/","VBFPresel/"]
histDirs = [""]

CUTS="dimuonMass < 150. && dimuonMass > 110."

## New Jet Pt Cuts Justin-Special
#CUTS += " && ((jetLead_pt > 25.) || (jetLead_pt > 20. && abs(jetLead_eta)<2.4))"
#CUTS += " && ((jetSub_pt > 25.) || (jetSub_pt > 20. && abs(jetSub_eta)<2.4))"
#CUTS += " && (abs(jetLead_eta) < 2.4 || jetLead_PUIDDisc > -0.5)"
#CUTS += " && (abs(jetSub_eta) < 2.4 || jetSub_PUIDDisc > -0.5)"

# New Jet Pt Cuts
#CUTS += " && ((jetLead_pt > 25.) || (jetLead_pt > 20. && abs(jetLead_eta)<2.4))"
#CUTS += " && ((jetSub_pt > 25.) || (jetSub_pt > 20. && abs(jetSub_eta)<2.4))"
#CUTS += " && ((jetLead_pt > 30.) || (jetLead_pt > 20. && abs(jetLead_eta)<2.4))"
#CUTS += " && ((jetSub_pt > 30.) || (jetSub_pt > 20. && abs(jetSub_eta)<2.4))"
#CUTS += " && ((jetLead_pt > 25.) || (jetLead_pt > 20. && abs(jetLead_eta)<2.6))"
#CUTS += " && ((jetSub_pt > 25.) || (jetSub_pt > 20. && abs(jetSub_eta)<2.6))"
#CUTS += " && jetLead_pt > 30."
#CUTS += " && jetSub_pt > 30."
CUTS += " && jetLead_pt > 25."
CUTS += " && jetSub_pt > 25."
# New Jet Pt Cuts PUID Tight
CUTS += " && (jetLead_FullPUIDFlag >= 7)"
CUTS += " && (jetSub_FullPUIDFlag >= 7)"
# New Jet Pt Cuts PUID Medium
#CUTS += " && (abs(jetLead_eta) < 2.4 || jetLead_FullPUIDFlag >= 6)"
#CUTS += " && (abs(jetSub_eta) < 2.4 || jetSub_FullPUIDFlag >= 6)"
# New Jet Pt Cuts PUID Loose
#CUTS += " && (abs(jetLead_eta) < 2.4 || jetLead_FullPUIDFlag >= 4)"
#CUTS += " && (abs(jetSub_eta) < 2.4 || jetSub_FullPUIDFlag >= 4)"

## Old Jet Pt Cuts
#CUTS += " && jetLead_pt > 30. && jetSub_pt > 30."

# VBF Cut Based 1
#CUTS += " && deltaEtaJets > 3.5 && dijetMass > 550. && ptMiss < 100."
# VBF Cut Based 2
#CUTS += " && deltaEtaJets > 3.5 && dijetMass > 500. && ptMiss < 50."

root.gErrorIgnoreLevel = root.kWarning
GLOBALCOUNTER=0

histNames = {}
if True:
    """
    histNames["dimuonMass"] = {"xlabel":"m_{#mu#mu} [GeV/c^{2}]","xlimits":[110.0,150.],"nbins":40}#,"ylimits":[0.1,5e5]}
    histNames["dimuonPt"] = {"xlabel":"p_{T,#mu#mu} [GeV/c]","xlimits":[0.0,200.0],"nbins":20}#,"ylimits":[0.1,1e5]}
    histNames["dimuonY"] = {"xlabel":"y_{#mu#mu}","xlimits":[-2.2,2.2],"nbins":10}#,"ylimits":[0.1,3e6]}
    histNames["cosThetaStar"] = {"xlabel":"cos(#theta^{*})","xlimits":[-2.2,2.2],"nbins":10}#,"ylimits":[0.1,3e6]}
    histNames["muonLead_pt"] = {"xlabel":"Leading Muon p_{T} [GeV/c]","xlimits":[25.,150.],"nbins":25}#,"ylimits":[0.1,3e6]}
    histNames["muonSub_pt"] = {"xlabel":"Sub-Leading Muon p_{T} [GeV/c]","xlimits":[25.,150.],"nbins":25}#,"ylimits":[0.1,3e6]}
    histNames["muonLead_eta"] = {"xlabel":"Leading Muon #eta","xlimits":[-2.1,2.1],"nbins":25}#,"ylimits":[0.1,3e6]}
    histNames["muonSub_eta"] = {"xlabel":"Sub-Leading Muon #eta","xlimits":[-2.1,2.1],"nbins":10}#,"ylimits":[0.1,3e6]}
    """

    histNames["nJets"] = {"xlabel":"N_{jets}","xlimits":[-0.5,5.5],"nbins":6}#,"ylimits":[0.1,3e6]}
    histNames["ptMiss"] = {"xlabel":"Missing p_{T} [GeV/c]","xlimits":[0.0,300.0],"nbins":12}#,"ylimits":[0.1,3e6]}
    histNames["deltaEtaJets"] = {"xlabel":"#Delta#eta(j_{1},j_{2})","xlimits":[0.0,7.0],"nbins":14}#,"ylimits":[0.1,3e6]}

    histNames["dijetMass"] = {"xlabel":"m_{jj} [GeV/c^{2}]","xlimits":[0.,1500.],"nbins":30}#,"ylimits":[0.1,5e5]}
    histNames["dijetPt"] = {"xlabel":"p_{T,jj} [GeV/c]","xlimits":[0.0,1000.0],"nbins":50}#,"ylimits":[0.1,1e5]}
    histNames["dijetY"] = {"xlabel":"y_{jj}","xlimits":[-5.0,5.0],"nbins":20}#,"ylimits":[0.1,3e6]}

    histNames["jetLead_pt"] = {"xlabel":"Leading Jet p_{T} [GeV/c]","xlimits":[20.,400.],"nbins":19}#,"ylimits":[0.1,3e6]}
    histNames["jetSub_pt"] = {"xlabel":"Sub-Leading Jet p_{T} [GeV/c]","xlimits":[20.,400.],"nbins":19}#,"ylimits":[0.1,3e6]}
    histNames["jetLead_eta"] = {"xlabel":"Leading Jet #eta","xlimits":[-5,5],"nbins":20}#,"ylimits":[0.1,3e6]}
    histNames["jetSub_eta"] = {"xlabel":"Sub-Leading Jet #eta","xlimits":[-5,5],"nbins":20}#,"ylimits":[0.1,3e6]}
    histNames["jetLead_abseta"] = {"xlabel":"Leading Jet |#eta|","xlimits":[0,5],"nbins":25}#,"ylimits":[0.1,3e6]}
    histNames["jetSub_abseta"] = {"xlabel":"Sub-Leading Jet |#eta|","xlimits":[0,5],"nbins":25}#,"ylimits":[0.1,3e6]}
    histNames["jetLead_PUIDDisc"] = {"xlabel":"Leading Jet PUID","xlimits":[-1,1],"nbins":20,"leg":ulLegendPos}#,"ylimits":[0.1,3e6]}
    histNames["jetSub_PUIDDisc"] = {"xlabel":"Sub-Leading Jet PUID","xlimits":[-1,1],"nbins":20,"leg":ulLegendPos}#,"ylimits":[0.1,3e6]}

    histNames["jetLead_PUIDDisc_Forward"] = {"xlabel":"Leading Jet PUID (|#eta|>2.4)","xlimits":[-1,1],"nbins":20,"leg":ulLegendPos}#,"ylimits":[0.1,3e6]}
    histNames["jetSub_PUIDDisc_Forward"] = {"xlabel":"Sub-Leading Jet PUID (|#eta|>2.4)","xlimits":[-1,1],"nbins":20,"leg":ulLegendPos}#,"ylimits":[0.1,3e6]}

    histNames["jetLead_PUIDDisc_Central"] = {"xlabel":"Leading Jet PUID (|#eta|<2.4)","xlimits":[-1,1],"nbins":20,"leg":ulLegendPos}#,"ylimits":[0.1,3e6]}
    histNames["jetSub_PUIDDisc_Central"] = {"xlabel":"Sub-Leading Jet PUID (|#eta|<2.4)","xlimits":[-1,1],"nbins":20,"leg":ulLegendPos}#,"ylimits":[0.1,3e6]}

    """
    histNames["KD"] = {"xlabel":"MEKD","xlimits":[0.0,1.0],"nbins":20}#,"ylimits":[0.1,3e6]}
    histNames["KDPdf"] = {"xlabel":"MEKD","xlimits":[0.0,1.0],"nbins":20}#,"ylimits":[0.1,3e6]}
    """

for key in histNames:
  if "BDTG" in dataDir and "BDT" in key:
    histNames[key]["xlimits"] = [-1,1]
    histNames[key]['vertLines']={"8TeV":0.2,"7TeV":0.2}
  if re.search("Pt",key) or re.match(".*pt.*",key):
    histNames[key]["units"] = " GeV/c"
  elif re.search("Mass",key):
    histNames[key]["units"] = " GeV/c^{2}"
  if LOGY and histNames[key].has_key("ylimits"):
    histNames[key]["ylimits"][1] *= 100
    if key == "yDiMu" or key == "yDiJet" or key == "cosThetaStar":
      histNames[key]["ylimits"][1] *= 10
    elif RUNPERIOD == "7TeV" and key != "mDiMu":
      histNames[key]["ylimits"][1] /= 10

tlatex = root.TLatex()
tlatex.SetNDC()
tlatex.SetTextSize(0.07)
tlatex.SetTextAlign(12)

#######################################
root.gROOT.SetBatch(True)
setStyle()
#######################################

scaleFactors = {}
#print "scale factors:"
for i in nEventsMap:
  if nEventsMap[i] ==0.0:
    scaleFactors[i] = 0.0
  else:
    p = getPeriod(i)
    scaleFactors[i] = xsec[i]*1000.0*LUMI/nEventsMap[i]*efficiencyMap[p]*mcPlotScaleFactorMap[p]
  #print "%s = %.2e" %(i,scaleFactors[i])

#######################################

class Dataset:
  def __init__(self,filename,legendEntry,color,scaleFactor,isData=False,isSignal=False):
    self.filename = filename
    self.legendEntry = legendEntry
    self.color = color
    self.scaleFactor = scaleFactor

    if filename != "":
      self.rootFile = root.TFile(filename)
      self.tree = self.rootFile.Get("outtree")
      self.tree.SetCacheSize(10000000);
      self.tree.AddBranchToCache("*");
    self.hists = {}
    self.datasetName = os.path.basename(filename)
    self.datasetName = self.datasetName.replace(".root","")
    self.isData=isData
    self.isSignal=isSignal

  def isZombie(self):
    return self.rootFile.IsZombie()

  def loadHistos(self,names,prefix=""):
    global GLOBALCOUNTER
    for name in names:
      #print("In datasetName: %, loading histogram: %" % (self.datasetName,name))
      tmpHistInfo = histNames[name]
      xlimits = tmpHistInfo["xlimits"]
      nbins = tmpHistInfo["nbins"]
      tmpHistName = name+str(GLOBALCOUNTER)
      GLOBALCOUNTER += 1
      varToDraw = name
      tmpCUTS = CUTS

      if "abs" in name:
        varToDraw = varToDraw.replace("abs","")
        varToDraw = "abs("+varToDraw+")"
      if "_Forward" in name:
        varToDraw = varToDraw.replace("_Forward","")
        if "jetLead" in varToDraw:
          tmpCUTS += " && abs(jetLead_eta) > 2.4"
        elif "jetSub" in varToDraw:
          tmpCUTS += " && abs(jetSub_eta) > 2.4"
      if "_Central" in name:
        varToDraw = varToDraw.replace("_Central","")
        if "jetLead" in varToDraw:
          tmpCUTS += " && abs(jetLead_eta) < 2.4"
        elif "jetSub" in varToDraw:
          tmpCUTS += " && abs(jetSub_eta) < 2.4"
      if name == "KD":
        sigNorm = MENormDict[RUNPERIOD]["sigME"]
        bakNorm = MENormDict[RUNPERIOD]["bakME"]
        varToDraw = "sigME*%f/(bakME*%f+sigME*%f)" % (sigNorm,bakNorm,sigNorm)
      if name == "KDPdf":
        sigNorm = MENormDict[RUNPERIOD]["sigMEPdf"]
        bakNorm = MENormDict[RUNPERIOD]["bakMEPdf"]
        varToDraw = "sigMEPdf*%f/(bakMEPdf*%f+sigMEPdf*%f)" %(sigNorm,bakNorm,sigNorm)
      drawStr = varToDraw+" >> "+tmpHistName+"("+str(nbins)+","+str(xlimits[0])+","+str(xlimits[1])+")"
      print drawStr
      cutStr = treeCut(prefix,tmpCUTS)
      print cutStr
      self.tree.Draw(drawStr,cutStr)
      tmp = root.gDirectory.Get(tmpHistName)
      if type(tmp) != root.TH1F:
        print("Warning: In datasetName: %s, loading histogram: %s: Object type is not TH1F!!" % (self.datasetName,prefix+name))
        continue
      tmp.UseCurrentStyle()
      if self.isSignal:
        tmp.SetLineColor(self.color)
      elif self.isData:
        tmp.SetMarkerColor(self.color)
        tmp.SetLineColor(self.color)
      else:
        tmp.SetFillColor(self.color)
      tmp.Scale(self.scaleFactor)
      self.hists[prefix+name] = tmp
  def eatOtherDatasets(self,listOfOthers):
    first = True
    for ds in listOfOthers:
      if first:
        first = False
        self.hists = ds.hists
        for key in self.hists:
          tmp = self.hists[key]
          if self.isSignal:
            tmp.SetLineColor(self.color)
            tmp.SetLineWidth(3)
          elif self.isData:
            tmp.SetMarkerColor(self.color)
            tmp.SetLineColor(self.color)
          else:
            tmp.SetFillColor(self.color)
      else:
        for key in self.hists:
          self.hists[key].Add(ds.hists[key])

#######################################

bkgDatasetList = []
for i in backgroundList:
  i += "_"+RUNPERIOD
  if i in scaleFactors:
    if scaleFactors[i]>0.0:
      filename = dataDir+i+".root"
      if not os.path.exists(filename):
          continue
      tmp = Dataset(filename,getLegendEntry(i),getColor(i),scaleFactors[i])
      if tmp.isZombie():
        print ("Warning: file for dataset %i is Zombie!!" % (i))
        continue
      #print("Loading Dataset: %s" % (i))
      for hDir in histDirs:
        tmp.loadHistos(histNames,prefix=hDir)
      bkgDatasetList.append(tmp)

sigDatasetList = []
for i in signalList:
  i += "_"+RUNPERIOD
  if i in scaleFactors:
    if scaleFactors[i]>0.0:
      filename = dataDir+i+".root"
      if not os.path.exists(filename):
          continue
      tmp = Dataset(filename,getLegendEntry(i),getColor(i),scaleFactors[i]*scaleHiggsBy,isSignal=True)
      if tmp.isZombie():
        print ("Warning: file for dataset %s is Zombie!!" %s (i))
        continue
      #print("Loading Dataset: %s" % (i))
      for hDir in histDirs:
        tmp.loadHistos(histNames,prefix=hDir)
      sigDatasetList.append(tmp)
oldSigDatasetList = sigDatasetList
if allHiggsTogether:
  #allHiggsDataset = Dataset("","H #rightarrow #mu#mu",root.kCyan-4,1.0,isSignal=True)
  #allHiggsDataset = Dataset("","H #rightarrow #mu#mu",root.kGray,1.0,isSignal=True)
  allHiggsDataset = Dataset("","H #rightarrow #mu#mu",1,1.0,isSignal=True)
  allHiggsDataset.eatOtherDatasets(sigDatasetList)
  sigDatasetList = [allHiggsDataset]

realDatasetList = []
for i in dataDict[RUNPERIOD]:
      filename = dataDir+i+".root"
      if not os.path.exists(filename):
        print("Error: Data file not found %s, exiting" % (filename))
        sys.exit(1)
      tmp = Dataset(filename,getLegendEntry(RUNPERIOD),1,1.0,isData=True)
      if tmp.isZombie():
        print ("Error: file for dataset %s is Zombie!!" % (i))
        sys.exit(1)
      #print("Loading Dataset: %s" % (i))
      for hDir in histDirs:
        tmp.loadHistos(histNames,prefix=hDir)
      realDatasetList.append(tmp)

#######################################

canvas = root.TCanvas("canvas")
leg = root.TLegend(*stdLegendPos)
leg.SetLineColor(0)
leg.SetFillColor(0)
uniqueLegendEntries = set()
for ds in bkgDatasetList:
  for hname in ds.hists:
    if ds.legendEntry not in uniqueLegendEntries:
      leg.AddEntry(ds.hists[hname],ds.legendEntry,"f")
      uniqueLegendEntries.add(ds.legendEntry)
    break
for ds in sigDatasetList:
  for hname in ds.hists:
    if ds.legendEntry not in uniqueLegendEntries:
      leg.AddEntry(ds.hists[hname],ds.legendEntry,"l")
      uniqueLegendEntries.add(ds.legendEntry)
    break
for ds in realDatasetList:
  for hname in ds.hists:
    if ds.legendEntry not in uniqueLegendEntries:
      leg.AddEntry(ds.hists[hname],ds.legendEntry,"pe")
      uniqueLegendEntries.add(ds.legendEntry)
    break
  break

#######################################
dataMCRatioStr = "############################################\n"

for histName in bkgDatasetList[0].hists:
  #print("Making Histo: %s" % histName)
  canvas.Clear()
  bkgHistList = []
  sigHistList = []
  for ds in bkgDatasetList:
    tmpHist = ds.hists[histName]
    bkgHistList.append(tmpHist)
  bkgHistList.reverse()
  for ds in sigDatasetList:
    tmpHist = ds.hists[histName]
    sigHistList.append(tmpHist)
  #bkgHistList.reverse()

  #dataHist = bkgDatasetList[0].hists[histName].Clone()
  #dataHist.Reset()
  dataHist = realDatasetList[0].hists[histName].Clone()
  dataHist.Reset()
  for realDS in realDatasetList:
    dataHist.Add(realDS.hists[histName])

  histBaseName = re.sub(r".*/","",histName)
  xtitle = histNames[histBaseName]["xlabel"]
  ylimits = []
  if histNames[histBaseName].has_key("ylimits"):
     ylimits = histNames[histBaseName]["ylimits"]
  legLeftPos = stdLegendPos[0]
  legRightPos = stdLegendPos[2]
  legBotPos = stdLegendPos[1]
  scaleHiggsPos = "std"
  if histNames[histBaseName].has_key("leg"):
    setLegPos(leg,(histNames[histBaseName]["leg"]))
    legLeftPos = histNames[histBaseName]["leg"][0]
    legRightPos = histNames[histBaseName]["leg"][2]
    legBotPos = histNames[histBaseName]["leg"][1]
    tmplegPos = histNames[histBaseName]["leg"]
    if tmplegPos == lcLegendPos:
      scaleHiggsPos = "lc"
    elif tmplegPos == llLegendPos:
      scaleHiggsPos = "ll"
    elif tmplegPos == ulLegendPos:
      scaleHiggsPos = "ul"
  else:
    setLegPos(leg,stdLegendPos)
  # For Auto y-limits
  yMaxVals = [legBotPos,0.75]
  yMaxXRanges = [[legLeftPos,legRightPos]]
  if scaleHiggsPos == "ul":
    yMaxXRanges += [[legRightPos,0.7]]
  else:
    yMaxXRanges += [[0.3,legLeftPos]]
  stack = DataMCStack(bkgHistList, dataHist, canvas, xtitle,lumi=LUMI,logy=LOGY,xlimits=histNames[histBaseName]["xlimits"],signalsNoStack=sigHistList,integralPlot=integralPlot,energyStr=RUNPERIOD,ylimits=ylimits,ylimitsRatio=ylimitsRatio,pullType=PULLTYPE,doMCErrors=MCErrors,yMaxVals=yMaxVals,yMaxXRanges=yMaxXRanges)


  leg.Draw("same")

  if scaleHiggsPos == "lc":
    if scaleHiggsBy != 1.0:
      tlatex.SetTextSize(0.07)
      tlatex.SetTextAlign(22)
      tlatex.DrawLatex(0.55,0.7,"Higgs #times %.0f" % (scaleHiggsBy))

    tlatex.SetTextSize(0.03)
    tlatex.SetTextAlign(22)
    tlatex.DrawLatex(0.55,0.75,anotateText)
  elif scaleHiggsPos == "ll" or scaleHiggsPos == "ul":
    if scaleHiggsBy != 1.0:
      tlatex.SetTextSize(0.07)
      tlatex.SetTextAlign(22)
      tlatex.DrawLatex(0.55,0.8,"Higgs #times %.0f" % (scaleHiggsBy))

    tlatex.SetTextSize(0.04)
    tlatex.SetTextAlign(23)
    tlatex.DrawLatex(0.55,1.0-gStyle.GetPadTopMargin()-0.02,anotateText)
  else:
    if scaleHiggsBy != 1.0:
      tlatex.SetTextSize(0.07)
      tlatex.SetTextAlign(32)
      tlatex.DrawLatex(legLeftPos-0.02,0.8,"Higgs #times %.0f" % (scaleHiggsBy))

    tlatex.SetTextSize(0.04)
    tlatex.SetTextAlign(33)
    tlatex.DrawLatex(legLeftPos-0.02,1.0-gStyle.GetPadTopMargin()-0.02,anotateText)

  vertLine = None
  arrow = None
  if histNames[histBaseName].has_key("vertLines"):
    #print("In vertLines for hist: %s" % (histBaseName))
    padX1 = stack.pad1.GetX1()
    padX2 = stack.pad1.GetX2()
    padY1 = stack.pad1.GetY1()
    padY2 = stack.pad1.GetY2()
    stack.pad1.cd()
    vertLineX = histNames[histBaseName]["vertLines"][RUNPERIOD]
    vertLine = root.TLine()
    vertLine.SetLineColor(root.kRed+1)
    vertLine.SetLineWidth(3)
    pad1Width = stack.pad1.XtoPixel(stack.pad1.GetX2())-stack.pad1.XtoPixel(stack.pad1.GetX1())
    pad1Height = stack.pad1.YtoPixel(stack.pad1.GetY1())
    #vertLine.DrawLineNDC(0.1,0.1,stack.pad1.XtoPixel(0.0)/pad1Width,0.9)
    #vertLine.DrawLine(0.1,0.1,-0.4,10.)
    ybot = stack.stack.GetMinimum()
    ytop = stack.stack.GetMaximum()
    if stack.dataHist.GetMaximum()>ytop:
        ytop = stack.dataHist.GetMaximum()
    if LOGY:
      ytop = stack.stack.GetMaximum()*2
    #if histNames[histBaseName].has_key("ylimits"):
    #   ylimits = histNames[histBaseName]["ylimits"]
    #   ybot = ylimits[0]
    #   ytop = ylimits[1]
    vertLine.DrawLine(vertLineX,ybot,vertLineX,ytop)

    xAxis = stack.stack.GetXaxis()
    arrowLength = (xAxis.GetXmax() - xAxis.GetXmin())/10.
    print(xAxis.GetXmax() , xAxis.GetXmin())
    arrowY = (ytop-ybot)*0.5
    arrowHeadSize = 0.025
    arrow = root.TArrow(vertLineX,arrowY,vertLineX+arrowLength,arrowY,0.025,"|>")
    arrow.SetLineWidth(3)
    #arrow.SetAngle(40)
    arrow.SetLineColor(root.kRed)
    arrow.SetFillColor(root.kRed)
    #arrow.Draw()

    stack.pad1.RedrawAxis()

  if stack.mcSumHist.Integral() == 0:
    print("Warning: MC Sum Hist Was 0 for %s, not saving image" % (histName))
    continue

  if histNames[histBaseName].has_key("units"):
    tmpAx = stack.stack.GetYaxis()
    tmpAx.SetTitle(tmpAx.GetTitle()+histNames[histBaseName]["units"])

  saveName = histName.replace("(","")
  saveName = saveName.replace(")","")
  saveName = saveName.replace("[","")
  saveName = saveName.replace("]","")
  saveName = saveName.replace("/","_")
  if integralPlot:
    saveName += "_IntPlot"
  saveAs(canvas,outDir+saveName+"_"+RUNPERIOD)

  match = re.match(r"(.*)_ptDiMu",saveName)
  if match:
    catName = match.group(1)
    dataMCRatioStr += "%-10s Data/MC Ratio: %.3f\n" % (catName,float(stack.nDataEvents)/stack.nMCEvents)
  elif saveName == "ptDiMu":
    catName = "All"
    dataMCRatioStr += "%-10s Data/MC Ratio: %.3f\n" % (catName,float(stack.nDataEvents)/stack.nMCEvents)

print(dataMCRatioStr)
