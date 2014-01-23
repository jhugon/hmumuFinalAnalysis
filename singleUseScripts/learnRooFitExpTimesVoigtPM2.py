#! /usr/bin/env python

from ROOT import gSystem

import sys
import os
import re
import math
from ROOT import *
gSystem.Load('libRooFit')
import ROOT as root
root.gROOT.SetBatch(True)

import singleHelpers
from helpers import *
from xsec import *

#root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(1) #For MINUIT

canvas = root.TCanvas()

massVeryLowRange = [80,95]
massLowRange = [110,120]
massHighRange = [130,160]

maxMass = massHighRange[1]
minMass = massVeryLowRange[0]
dimuonMass = root.RooRealVar("dimuonMass","dimuonMass",minMass,maxMass)
dimuonMass.setRange("z",88,94)
dimuonMass.setRange("verylow",massVeryLowRange[0],massVeryLowRange[1])
dimuonMass.setRange("low",massLowRange[0],massLowRange[1])
dimuonMass.setRange("high",massHighRange[0],massHighRange[1])
dimuonMass.setRange("signal",massLowRange[1],massHighRange[0])
dimuonMass.setRange("lowsighigh",massLowRange[0],massHighRange[1])

channelName = "silly"
voitWidth = root.RooRealVar(channelName+"_voitWidth","voitWidth",2.4952)
voitWidth.setConstant(True)
voitmZ = root.RooRealVar(channelName+"_voitmZ","voitmZ",91,85,95)
voitSig = root.RooRealVar(channelName+"_voitSig","voitSig",1.5,0.0,30.0)
voitMmumu = root.RooVoigtian(channelName+"bak_voitMmumu","voitMmumu",dimuonMass,voitmZ,voitWidth,voitSig)
expParam = root.RooRealVar(channelName+"_expParam","expParam",0.01,-1,1)
expMm2Mmumu = root.RooGenericPdf(channelName+"bak_expMm2","exp(@0*@1*@1)*pow(@0,-2)",root.RooArgList(dimuonMass,expParam))
mixParam = root.RooRealVar(channelName+"_mixParam","mixParam",0.5,0,1)
pdfMmumu = root.RooAddPdf("bak","bak",root.RooArgList(voitMmumu,expMm2Mmumu),root.RooArgList(mixParam))

#####################################################################

filename = getDataStage2Directory()+"/DYJetsToLL_8TeV.root"
#filename = getDataStage2Directory()+"/ttbar_8TeV.root"
print filename
f = root.TFile(filename)
samples = ['ttbar_8TeV','DYJetsToLL_8TeV']
samples += ['ggHmumu125_8TeV','vbfHmumu125_8TeV']
scaleFactors = {}
#scaleFactors['ggHmumu125_8TeV'] = 100.
#scaleFactors['ttbar_8TeV'] = 10.
mDiMu = None
files = []
hstack = root.THStack("hstack","")
for sample in samples:
  f = root.TFile(getDataStage2Directory()+"/"+sample+".root")
  #mDiMuTmp = f.Get("Jets01FailPtG10BE/mDiMu")
  #mDiMuTmp = f.Get("Jets01PassPtG10BB/mDiMu")
  mDiMuTmp = f.Get("Jet2CutsGFPass/mDiMu")
  #mDiMuTmp = f.Get("Jet2CutsVBFPass/mDiMu")
  mDiMuTmp.Rebin(2)
  mDiMuTmp.Scale(xsec[sample]*1000.*lumiDict['8TeV']/nEventsMap[sample])
  mDiMuTmp.SetFillColor(colors[re.sub(r"_.*","",sample)])
  if scaleFactors.has_key(sample):
    mDiMuTmp.Scale(scaleFactors[sample])
  hstack.Add(mDiMuTmp)
  if mDiMu == None:
    mDiMu = mDiMuTmp.Clone("mDiMuHist")
    print "creating hist from ",sample
  else:
    mDiMu.Add(mDiMuTmp)
    print "adding hist from ",sample
  files.append(f)

#####################################################################

dimuonMassRooDataHist = root.RooDataHist("template","template",root.RooArgList(dimuonMass),mDiMu)

voitMmumu.fitTo(dimuonMassRooDataHist,root.RooFit.Range("z"),root.RooFit.SumW2Error(False),PRINTLEVEL)
voitmZ.setConstant(True)
voitSig.setConstant(True)

fr = pdfMmumu.fitTo(dimuonMassRooDataHist,root.RooFit.Range("lowsighigh"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
chi2 = pdfMmumu.createChi2(dimuonMassRooDataHist,root.RooFit.Range("lowsighigh"))

plotMmumu = dimuonMass.frame(root.RooFit.Range("lowsighigh"))

#############################################################################

dimuonMass.Print()
plotMmumu.Print()
dimuonMassRooDataHist.plotOn(plotMmumu)
errVisArg = root.RooFit.VisualizeError(fr,1,True)
errFillArg = root.RooFit.FillStyle(3001)
pdfMmumu.plotOn(plotMmumu,root.RooFit.Range("lowsighigh"),errVisArg,errFillArg,root.RooFit.FillColor(root.kCyan))
dimuonMassRooDataHist.plotOn(plotMmumu)
pdfMmumu.plotOn(plotMmumu,root.RooFit.Range("lowsighigh"))
pdfMmumu.plotOn(plotMmumu,root.RooFit.LineStyle(2))
pdfMmumu.plotOn(plotMmumu,root.RooFit.Components(channelName+"bak_expMm2"),root.RooFit.Range("lowsighigh"),errVisArg,errFillArg,root.RooFit.FillColor(root.kGreen-7))
pdfMmumu.plotOn(plotMmumu,root.RooFit.Components(channelName+"bak_expMm2"),root.RooFit.LineStyle(2),root.RooFit.LineColor(root.kGreen+1),root.RooFit.Range("lowsighigh"))
pdfMmumu.plotOn(plotMmumu,root.RooFit.Components(channelName+"bak_voitMmumu"),root.RooFit.Range("lowsighigh"),errVisArg,errFillArg,root.RooFit.FillColor(root.kRed-7))
pdfMmumu.plotOn(plotMmumu,root.RooFit.Components(channelName+"bak_voitMmumu"),root.RooFit.LineStyle(2),root.RooFit.LineColor(root.kRed+1),root.RooFit.Range("lowsighigh"))

#axisHist = root.TH2F("axisHist","",1,massLowRange[0],massHighRange[1],1,0,1200)
axisHist = root.TH2F("axisHist","",1,massLowRange[0],massHighRange[1],1,0,75)
axisHist.Draw()
hstack.Draw('histsame')
plotMmumu.Draw('same')
canvas.SaveAs("TestExpTimesVoigtPMm2.png")

###########################################################################

print("chi2/ndf: {}".format(plotMmumu.chiSquare(pdfMmumu.getParameters(dimuonMassRooDataHist).getSize())))
for i in [voitmZ,voitSig,expParam,mixParam]:
  print("{0}: {1:.3g} +/- {2:.3%}".format(i.GetName(),i.getVal(),abs(i.getError()/i.getVal())))
for i in [voitmZ,voitSig,expParam,mixParam]:
  i.Print()
#yay = raw_input("Press Enter to continue...")

