#!/usr/bin/env python

import optparse
parser = optparse.OptionParser(description="Makes data and MC comparison plots.")
parser.add_option("--scaleHiggsBy", help="Scale Factor For Higgs Samples",type=float,default=1000.0)
parser.add_option("--energy", help="Energy to Use, either 7TeV or 8TeV",default="8TeV")
parser.add_option("--disableMCErrors", help="Energy to Use, either 7TeV or 8TeV",action="store_true",default=False)
parser.add_option("--outPrefix", help="Output filename prefix",default="")
parser.add_option("--n01Jet", help="Preset to run 0,1 Jet Presel Plots",action="store_true",default=False)
parser.add_option("--n01JetMassOnly", help="Preset to run 0,1 Jet Presel Mass Plot",action="store_true",default=False)
parser.add_option("--n2Jet", help="Preset to run 2 Jet Presel Plots",action="store_true",default=False)
parser.add_option("--n2JetMassOnly", help="Preset to run 2 Jet Presel Mass Plot",action="store_true",default=False)
parser.add_option("--n2JetVBFTight", help="Preset to run 2 Jet VBF Tight Plots",action="store_true",default=False)
parser.add_option("--n2JetVBFLoose", help="Preset to run 2 Jet VBF Loose Plots",action="store_true",default=False)
parser.add_option("--n2JetGFTight", help="Preset to run 2 Jet GF Tight Plots",action="store_true",default=False)
args, fakeargs = parser.parse_args()

if args.n01Jet:
  args.scaleHiggsBy = 500.
  args.outPrefix = "nonVBFPresel"
  args.disableMCErrors = True
  print("Running 01Jet Presel")

if args.n01JetMassOnly:
  args.scaleHiggsBy = 500.
  args.outPrefix = "nonVBFPresel"
  args.disableMCErrors = True
  print("Running 01Jet Presel Mass Only")

if args.n2Jet:
  args.scaleHiggsBy = 1000.
  args.outPrefix = ""
  print("Running 2Jet Presel")

if args.n2JetMassOnly:
  args.scaleHiggsBy = 200.
  args.outPrefix = "presel"
  print("Running 2Jet Presel Mass Only")

if args.n2JetVBFTight:
  args.scaleHiggsBy = 50.
  args.outPrefix = "vbfTight"
  print("Running 2Jet VBF Tight")

if args.n2JetVBFLoose:
  args.scaleHiggsBy = 50.
  args.outPrefix = "loose"
  print("Running 2Jet Loose")

if args.n2JetGFTight:
  args.scaleHiggsBy = 50.
  args.outPrefix = "gfTight"
  print("Running 2Jet GF Tight")

from xsec import *
from helpers import *
import ROOT as root
import os
import sys

dataDir = getDataStage2Directory()
#dataDir = "/data/uftrig01b/jhugon/hmumu/analysisV00-01-10/forGPReRecoMuScleFit/"
outDir = "output/"

RUNPERIOD=args.energy
LUMI=lumiDict[RUNPERIOD]

scaleHiggsBy = args.scaleHiggsBy

SCALEMC2DATA=True
JETErrors=True
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

if args.disableMCErrors:
  JETErrors=False
  MCErrors=False

anotateText2 = ""
anotateText3 = "Analysis A"

anotateText = "M(#mu#mu) #in [110,160] GeV/c^{2}"
#if args.n2JetMassOnly or args.n2JetVBFTight or args.n2JetGFTight or args.n2JetVBFLoose or args.n01JetMassOnly:
#  anotateText = ""

urLegendPos = [0.70,0.67,0.9,0.9]
ulLegendPos = [0.20,0.67,0.4,0.9]
ucLegendPos = [0.46,0.67,0.64,0.9]
lcLegendPos = [0.46,0.35,0.64,0.63]
llLegendPos = [0.20,0.35,0.4,0.63]
ccLegendPos = [0.46,0.47,0.64,0.7]
stdLegendPos = urLegendPos

if args.outPrefix != "":
  args.outPrefix += "/"
histDirs = [args.outPrefix]

CUTS="dimuonMass < 160. && dimuonMass > 110."

if args.n01Jet or args.n01JetMassOnly:
  anotateText2 = "0,1-Jet Preselection"
  CUTS+=" && !(jetLead_pt>40. && jetSub_pt>30. && ptMiss<40.)"

if args.n2Jet or args.n2JetMassOnly:
  anotateText2 = "2-Jet Preselection"
  CUTS+=" && jetLead_pt>40. && jetSub_pt>30. && ptMiss<40."
  CUTS+=" && jetLead_pt>40. && jetSub_pt>30."

if args.n2JetVBFTight:
  anotateText2 = "2-Jet VBF Tight"
  CUTS+=" && jetLead_pt>40. && jetSub_pt>30. && ptMiss<40. && dijetMass > 650. && deltaEtaJets>3.5"

if args.n2JetVBFLoose:
  anotateText2 = "2-Jet Loose"
  CUTS+=" && jetLead_pt>40. && jetSub_pt>30. && ptMiss<40. && !(dijetMass > 650. && deltaEtaJets>3.5) && !(dijetMass>250. && dimuonPt>50.)"

if args.n2JetGFTight:
  anotateText2 = "2-Jet GF Tight"
  CUTS+=" && jetLead_pt>40. && jetSub_pt>30. && ptMiss<40. && !(dijetMass > 650. && deltaEtaJets>3.5) && dijetMass>250. && dimuonPt>50."

root.gErrorIgnoreLevel = root.kWarning
GLOBALCOUNTER=0

histNames = {}
if True:
    if args.n2JetVBFTight or args.n2JetGFTight or args.n2JetVBFLoose:
      histNames["dimuonMass"] = {"xlabel":"M(#mu#mu) [GeV/c^{2}]","xlimits":[110.0,160.],"nbins":20}#,"ylimits":[0.1,5e5]}
    elif not args.n2Jet:
      histNames["dimuonMass"] = {"xlabel":"M(#mu#mu) [GeV/c^{2}]","xlimits":[110.0,160.],"nbins":50}#,"ylimits":[0.1,5e5]}
    if not ( args.n2JetVBFTight or args.n2JetGFTight or args.n2JetVBFLoose or args.n2JetMassOnly or args.n01JetMassOnly):
      histNames["dimuonPt"] = {"xlabel":"p_{T}(#mu#mu) [GeV/c]","xlimits":[0.0,200.0],"nbins":20}#,"ylimits":[0.1,1e5]}
#    histNames["dimuonY"] = {"xlabel":"y_{#mu#mu}","xlimits":[-2.2,2.2],"nbins":22}#,"ylimits":[0.1,3e6]}
#    histNames["cosThetaStar"] = {"xlabel":"cos(#theta^{*})","xlimits":[-1,1],"nbins":20}#,"ylimits":[0.1,3e6]}
    #histNames["muonLead_pt"] = {"xlabel":"Leading Muon p_{T} [GeV/c]","xlimits":[25.,150.],"nbins":25}#,"ylimits":[0.1,3e6]}
    #histNames["muonSub_pt"] = {"xlabel":"Sub-Leading Muon p_{T} [GeV/c]","xlimits":[25.,150.],"nbins":25}#,"ylimits":[0.1,3e6]}
    #histNames["muonLead_eta"] = {"xlabel":"Leading Muon #eta","xlimits":[-2.1,2.1],"nbins":25}#,"ylimits":[0.1,3e6]}
    #histNames["muonSub_eta"] = {"xlabel":"Sub-Leading Muon #eta","xlimits":[-2.1,2.1],"nbins":10}#,"ylimits":[0.1,3e6]}

    #histNames["nJets"] = {"xlabel":"N_{jets}","xlimits":[-0.5,5.5],"nbins":6}#,"ylimits":[0.1,3e6]}
    #histNames["ptMiss"] = {"xlabel":"Missing p_{T} [GeV/c]","xlimits":[0.0,300.0],"nbins":12}#,"ylimits":[0.1,3e6]}
    if not (args.n01Jet or args.n2JetMassOnly or args.n2JetVBFTight or args.n2JetGFTight or args.n2JetVBFLoose or args.n01JetMassOnly):
      histNames["deltaEtaJets"] = {"xlabel":"|#Delta#eta(jj)|","xlimits":[0.0,7.0],"nbins":14}#,"ylimits":[0.1,3e6]}

      histNames["dijetMass"] = {"xlabel":"M(jj) [GeV/c^{2}]","xlimits":[0.,1000.],"nbins":20}#,"ylimits":[0.1,5e5]}
    #histNames["dijetPt"] = {"xlabel":"p_{T,jj} [GeV/c]","xlimits":[0.0,1000.0],"nbins":50}#,"ylimits":[0.1,1e5]}
    #histNames["dijetY"] = {"xlabel":"y_{jj}","xlimits":[-5.0,5.0],"nbins":20}#,"ylimits":[0.1,3e6]}

    #histNames["jetLead_pt"] = {"xlabel":"Leading Jet p_{T} [GeV/c]","xlimits":[20.,400.],"nbins":19}#,"ylimits":[0.1,3e6]}
    #histNames["jetSub_pt"] = {"xlabel":"Sub-Leading Jet p_{T} [GeV/c]","xlimits":[20.,400.],"nbins":19}#,"ylimits":[0.1,3e6]}
    #histNames["jetLead_pt_Central"] = {"xlabel":"Leading Jet p_{T} [GeV/c] (|#eta|<2.4)","xlimits":[20.,400.],"nbins":19}#,"ylimits":[0.1,3e6]}
    #histNames["jetSub_pt_Central"] = {"xlabel":"Sub-Leading Jet p_{T} [GeV/c] (|#eta|<2.4)","xlimits":[20.,400.],"nbins":19}#,"ylimits":[0.1,3e6]}
    #histNames["jetLead_pt_lowpt_Forward"] = {"xlabel":"Leading Jet p_{T} [GeV/c] (|#eta|>2.4)","xlimits":[20.,100.],"nbins":16}#,"ylimits":[0.1,3e6]}
    #histNames["jetSub_pt_lowpt_Forward"] = {"xlabel":"Sub-Leading Jet p_{T} [GeV/c] (|#eta|>2.4)","xlimits":[20.,100.],"nbins":16}#,"ylimits":[0.1,3e6]}
    #histNames["jetLead_pt_lowpt_Forward"] = {"xlabel":"Leading Jet p_{T} [GeV/c] (|#eta|>2.4)","xlimits":[20.,100.],"nbins":16}#,"ylimits":[0.1,3e6]}
    #histNames["jetSub_pt_lowpt_Forward"] = {"xlabel":"Sub-Leading Jet p_{T} [GeV/c] (|#eta|>2.4)","xlimits":[20.,100.],"nbins":16}#,"ylimits":[0.1,3e6]}
    #histNames["jetLead_pt_lowpt_Central"] = {"xlabel":"Leading Jet p_{T} [GeV/c] (|#eta|<2.4)","xlimits":[20.,100.],"nbins":16}#,"ylimits":[0.1,3e6]}
    #histNames["jetSub_pt_lowpt_Central"] = {"xlabel":"Sub-Leading Jet p_{T} [GeV/c] (|#eta|<2.4)","xlimits":[20.,100.],"nbins":16}#,"ylimits":[0.1,3e6]}
    #histNames["jetLead_eta"] = {"xlabel":"Leading Jet #eta","xlimits":[-5,5],"nbins":20}#,"ylimits":[0.1,3e6]}
    #histNames["jetSub_eta"] = {"xlabel":"Sub-Leading Jet #eta","xlimits":[-5,5],"nbins":20}#,"ylimits":[0.1,3e6]}
    #histNames["jetLead_abseta"] = {"xlabel":"Leading Jet |#eta|","xlimits":[0,5],"nbins":25}#,"ylimits":[0.1,3e6]}
    #histNames["jetSub_abseta"] = {"xlabel":"Sub-Leading Jet |#eta|","xlimits":[0,5],"nbins":25}#,"ylimits":[0.1,3e6]}
    #histNames["jetLead_PUIDDisc"] = {"xlabel":"Leading Jet PUID","xlimits":[-1,1],"nbins":20,"leg":ulLegendPos}#,"ylimits":[0.1,3e6]}
    #histNames["jetSub_PUIDDisc"] = {"xlabel":"Sub-Leading Jet PUID","xlimits":[-1,1],"nbins":20,"leg":ulLegendPos}#,"ylimits":[0.1,3e6]}

    #histNames["jetLead_PUIDDisc_Forward"] = {"xlabel":"Leading Jet PUID (|#eta|>2.4)","xlimits":[-1,1],"nbins":20,"leg":ulLegendPos}#,"ylimits":[0.1,3e6]}
    #histNames["jetSub_PUIDDisc_Forward"] = {"xlabel":"Sub-Leading Jet PUID (|#eta|>2.4)","xlimits":[-1,1],"nbins":20,"leg":ulLegendPos}#,"ylimits":[0.1,3e6]}

    #histNames["jetLead_PUIDDisc_Central"] = {"xlabel":"Leading Jet PUID (|#eta|<2.4)","xlimits":[-1,1],"nbins":20,"leg":ulLegendPos}#,"ylimits":[0.1,3e6]}
    #histNames["jetSub_PUIDDisc_Central"] = {"xlabel":"Sub-Leading Jet PUID (|#eta|<2.4)","xlimits":[-1,1],"nbins":20,"leg":ulLegendPos}#,"ylimits":[0.1,3e6]}

    #histNames["KD"] = {"xlabel":"MEKD","xlimits":[0.0,1.0],"nbins":20}#,"ylimits":[0.1,3e6]}
    #histNames["KDPdf"] = {"xlabel":"MEKD","xlimits":[0.0,1.0],"nbins":20}#,"ylimits":[0.1,3e6]}

    #histNames["bdtVBF"] = {"xlabel":"2-Jet BDT","xlimits":[-0.5,0.5],"nbins":40}#,"ylimits":[0.1,3e6]}

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

if "KD" in CUTS:
  sigNorm = MENormDict[RUNPERIOD]["sigME"]
  bakNorm = MENormDict[RUNPERIOD]["bakME"]
  varToPutIn = "sigME*%f/(bakME*%f+sigME*%f)" % (sigNorm,bakNorm,sigNorm)
  CUTS = CUTS.replace("KD",varToPutIn)

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
  def __init__(self,filename,legendEntry,color,scaleFactor,isData=False,isSignal=False,error=None):
    self.filename = filename
    self.legendEntry = legendEntry
    self.color = color
    self.scaleFactor = scaleFactor
    self.error=error

    print "filename: {0}".format(filename)
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

      if "_lowpt" in name:
        varToDraw = varToDraw.replace("_lowpt","")
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

      if self.error!= None:
       cutList = re.split('\W+', tmpCUTS)
       varList = re.split('\W+', varToDraw)
       cutSet = set()
       varSet = set()
       for i in cutList:
         if (not re.match(r'\d+',i)) and (not i==''):
            if not (i in cutSet):
              cutSet.add(i)
       for i in varList:
         if (not re.match(r'\d+',i)) and (not i==''):
            if not (i in varSet):
              varSet.add(i)
       cutList = list(cutSet)
       varList = list(varSet)
       newCutList = [i + "_"+self.error for i in cutList]
       newVarList = [i + "_"+self.error for i in varList]
       branchNameList = [i.GetName() for i in self.tree.GetListOfBranches()]
       for i in reversed(range(len(newCutList))):
         newVarTmp = newCutList[i]
         if not (newVarTmp in branchNameList):
            newCutList.pop(i)
            cutList.pop(i)
       for i in reversed(range(len(newVarList))):
         newVarTmp = newVarList[i]
         if not (newVarTmp in branchNameList):
            newVarList.pop(i)
            varList.pop(i)

       for i,iNew in zip(cutList,newCutList):
         tmpCUTS = tmpCUTS.replace(i,iNew)
       for i,iNew in zip(varList,newVarList):
         varToDraw = varToDraw.replace(i,iNew)

      drawStr = varToDraw+" >> "+tmpHistName+"("+str(nbins)+","+str(xlimits[0])+","+str(xlimits[1])+")"
      #print drawStr
      cutStr = treeCut(prefix,tmpCUTS)
      #print cutStr
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

bkgErrDatasetMap = {}
if JETErrors:
  errorList = ["JESUp","JESDown","JERUp","JERDown"]
  for e in errorList:
    bkgErrDatasetListTmp = []
    for i in backgroundList:
      i += "_"+RUNPERIOD
      if i in scaleFactors:
        if scaleFactors[i]>0.0:
          filename = dataDir+i+".root"
          if not os.path.exists(filename):
              continue
          tmp = Dataset(filename,getLegendEntry(i),getColor(i),scaleFactors[i],error=e)
          if tmp.isZombie():
            print ("Warning: file for dataset %i is Zombie!!" % (i))
            continue
          print("Loading Variation Dataset: %s %s" % (i,e))
          for hDir in histDirs:
            tmp.loadHistos(histNames,prefix=hDir)
          bkgErrDatasetListTmp.append(tmp)
    bkgErrDatasetMap[e] = bkgErrDatasetListTmp

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

  mcVariations=None
  if JETErrors:
    mcVariations={}
    for e in bkgErrDatasetMap:
      bkgErrList = []
      for ds in bkgErrDatasetMap[e]:
        bkgErrList.append(ds.hists[histName])
      mcVariations[e] = bkgErrList

  histBaseName = re.sub(r".*/","",histName)
  print("Making hist: "+histBaseName)
  categoryName = re.sub(r"/.*","",histName)+'/'
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
  stack = DataMCStack(bkgHistList, dataHist, canvas, xtitle,lumi=LUMI,logy=LOGY,xlimits=histNames[histBaseName]["xlimits"],signalsNoStack=sigHistList,integralPlot=integralPlot,energyStr=RUNPERIOD,ylimits=ylimits,ylimitsRatio=ylimitsRatio,pullType=PULLTYPE,doMCErrors=MCErrors,yMaxVals=yMaxVals,yMaxXRanges=yMaxXRanges,mcVariations=mcVariations,scaleMC2Data=SCALEMC2DATA)


  leg.Draw("same")

  if scaleHiggsPos == "lc":
    if scaleHiggsBy != 1.0:
      tlatex.SetTextSize(0.07)
      tlatex.SetTextAlign(22)
      tlatex.DrawLatex(0.55,0.7,"Higgs #times %.0f" % (scaleHiggsBy))

    tlatex.SetTextSize(0.03)
    tlatex.SetTextAlign(22)
    if histBaseName != "dimuonMass":
      tlatex.DrawLatex(0.55,0.75,anotateText)
    tlatex.DrawLatex(0.55,0.55,anotateText2)
    tlatex.DrawLatex(0.55,0.6,anotateText3)
  elif scaleHiggsPos == "ll" or scaleHiggsPos == "ul":
    if scaleHiggsBy != 1.0:
      tlatex.SetTextSize(0.07)
      tlatex.SetTextAlign(22)
      tlatex.DrawLatex(0.55,0.82,"Higgs #times %.0f" % (scaleHiggsBy))

    tlatex.SetTextSize(0.04)
    tlatex.SetTextAlign(23)
    if histBaseName != "dimuonMass":
      tlatex.DrawLatex(0.55,1.0-gStyle.GetPadTopMargin()-0.02,anotateText)
    tlatex.DrawLatex(0.55,0.72,anotateText2)
    tlatex.DrawLatex(0.55,0.77,anotateText3)
  else:
    if scaleHiggsBy != 1.0:
      tlatex.SetTextSize(0.07)
      tlatex.SetTextAlign(32)
      tlatex.DrawLatex(legLeftPos-0.02,0.82,"Higgs #times %.0f" % (scaleHiggsBy))

    tlatex.SetTextSize(0.04)
    tlatex.SetTextAlign(33)
    tlatex.SetTextSize(0.035)
    if histBaseName != "dimuonMass":
      tlatex.DrawLatex(legLeftPos-0.02,1.0-gStyle.GetPadTopMargin()-0.02,anotateText)
    tlatex.SetTextSize(0.04)
    tlatex.DrawLatex(legLeftPos-0.02,0.77,anotateText2)
    tlatex.DrawLatex(legLeftPos-0.02,0.72,anotateText3)
    #tlatex.SetTextAlign(23)
    #tlatex.DrawLatex(0.55,0.77,anotateText3)

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
    tmpAx = stack.histForAxis.GetYaxis()
    tmpAx.SetTitle(tmpAx.GetTitle()+histNames[histBaseName]["units"])

  saveName = histName.replace("(","")
  saveName = saveName.replace(")","")
  saveName = saveName.replace("[","")
  saveName = saveName.replace("]","")
  saveName = saveName.replace("/","_")
  if integralPlot:
    saveName += "_IntPlot"
  saveAs(canvas,outDir+saveName+"_"+RUNPERIOD)

  dataMCRatioStr += "%-10s Data/MC Ratio: %.3f\n" % (histName,float(stack.nDataEvents)/stack.nMCEvents)

print(dataMCRatioStr)
