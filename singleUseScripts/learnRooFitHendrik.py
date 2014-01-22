#! /usr/bin/env python

from ROOT import gSystem

import sys
import os
import re
import math
from ROOT import *
gSystem.Load('libRooFit')
import ROOT as root

import singleHelpers
from helpers import *

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
mMuMu = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)
mMuMu.setRange("z",88,94)
mMuMu.setRange("verylow",massVeryLowRange[0],massVeryLowRange[1])
mMuMu.setRange("low",massLowRange[0],massLowRange[1])
mMuMu.setRange("high",massHighRange[0],massHighRange[1])
mMuMu.setRange("signal",massLowRange[1],massHighRange[0])
mMuMu.setRange("lowsighigh",massLowRange[0],massHighRange[1])

bwmZ = root.RooRealVar("bwmZ","bwmZ",85,95)
bwSig = root.RooRealVar("bwSig","bwSig",0.0,30.0)
expLambda = root.RooRealVar("expLambda","expLambda",-1e-03,1e-06,-1e-06)
mixParam = root.RooRealVar("mixParam","mixParam",0.5,0,1)

phoMmumu = root.RooGenericPdf("phoMmumu","pow(mMuMu,-2)",root.RooArgList(mMuMu))
expMmumu = root.RooExponential("expMmumu","expMmumu",mMuMu,expLambda)
bwMmumu = root.RooBreitWigner("bwMmumu","bwMmumu",mMuMu,bwmZ,bwSig)


sumMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",root.RooArgList(bwMmumu,phoMmumu),root.RooArgList(mixParam))
#pdfMmumu = root.RooProdPdf("pdfMmumu","pdfMmumu",root.RooArgList(sumMmumu,expMmumu))
#pdfMmumu = phoMmumu
#pdfMmumu = bwMmumu
#pdfMmumu = sumMmumu
pdfMmumu = root.RooProdPdf("pdfMmumu","pdfMmumu",root.RooArgList(bwMmumu,expMmumu))

#####################################################################

filename = getDataStage2Directory()+"/DYJetsToLL_8TeV.root"
print filename
f = root.TFile(filename)
#f = root.TFile("input/ttbar.root")

#mDiMu = f.Get("mDiMu")
#mDiMu = f.Get("IncPresel/mDiMu")
mDiMu = f.Get("Jets01PassPtG10BB/mDiMu")
mDiMu.Rebin(2)

#####################################################################

mMuMuRooDataHist = root.RooDataHist("template","template",root.RooArgList(mMuMu),mDiMu)

bwMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("z"),root.RooFit.SumW2Error(False),PRINTLEVEL)
bwmZ.setConstant(True)
bwSig.setConstant(True)

#plotMmumuZ = mMuMu.frame(root.RooFit.Range("z"))
#mMuMuRooDataHist.plotOn(plotMmumuZ)
#bwMmumu.plotOn(plotMmumuZ,root.RooFit.Range("z"))
##axisHistZ = TH2F("axisHistZ","",1,88,94,1,0,2000)
##axisHistZ.Draw()
##plotMmumuZ.Draw("same")
#plotMmumuZ.Draw()
#yay = raw_input("Press Enter to continue...")
#canvas.Clear()

pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("lowsighigh"),root.RooFit.SumW2Error(False),PRINTLEVEL)
chi2 = pdfMmumu.createChi2(mMuMuRooDataHist)

plotMmumu = mMuMu.frame(root.RooFit.Range("lowsighigh"))

drawRange = [80,150]
#mMuMu.setRange(*drawRange)
mMuMu.Print()
plotMmumu.Print()
#mMuMuRooDataHist.plotOn(plotMmumu,root.RooFit.Range(*drawRange))
mMuMuRooDataHist.plotOn(plotMmumu)
pdfMmumu.plotOn(plotMmumu,root.RooFit.Range("lowsighigh"))
#pdfMmumu.plotOn(plotMmumu,root.RooFit.LineStyle(2))
pdfMmumu.plotOn(plotMmumu,root.RooFit.Components("phoMmumu"),root.RooFit.LineStyle(2),root.RooFit.LineColor(root.kGreen+1),root.RooFit.Range("lowsighigh"))
pdfMmumu.plotOn(plotMmumu,root.RooFit.Components("bwMmumu"),root.RooFit.LineStyle(2),root.RooFit.LineColor(root.kRed+1),root.RooFit.Range("lowsighigh"))

plotMmumu.Draw()

print("chi2: {}".format(chi2.getVal()))
for i in [bwmZ,bwSig,expLambda,mixParam]:
  print("{0}: {1:.3g} +/- {2:.3%}".format(i.GetName(),i.getVal(),abs(i.getError()/i.getVal())))
yay = raw_input("Press Enter to continue...")

