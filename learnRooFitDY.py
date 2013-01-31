#! /usr/bin/env python

from ROOT import gSystem

import sys
import os
import re
import math
from ROOT import *
gSystem.Load('libRooFit')
import ROOT as root

from helpers import *

#####################################################################

f = root.TFile("input/freezeSample/DYJetsToLL_8TeV.root")
hist = f.Get("IncPresel/mDiMu")
#hist.Rebin(2)

channelName = ""

SIGNALFIT = [110.,140.]
controlRegionLow=[110,120]
controlRegionHigh=[130,160]

#####################################################################

#root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(1) #For MINUIT

canvas = root.TCanvas()

mMuMu = root.RooRealVar("mMuMu","mMuMu",controlRegionLow[0],controlRegionHigh[1])
mMuMu.setRange("low",controlRegionLow[0],controlRegionLow[1])
mMuMu.setRange("high",controlRegionHigh[0],controlRegionHigh[1])
mMuMu.setRange("signal",controlRegionLow[1],controlRegionHigh[0])
mMuMu.setRange("signalfit",SIGNALFIT[0],SIGNALFIT[1])

voitWidth = root.RooRealVar(channelName+"_voitWidth","voitWidth",2.4952)
voitmZ = root.RooRealVar(channelName+"_voitmZ","voitmZ",85,95)
voitSig = root.RooRealVar(channelName+"_voitSig","voitSig",0.0,30.0)
voitMmumu = root.RooVoigtian("bak_voitMmumu","voitMmumu",mMuMu,voitmZ,voitWidth,voitSig)

expParam = root.RooRealVar(channelName+"_expParam","expParam",-1,0)
expMmumu = root.RooExponential("bak_expMmumu","expMmumu",mMuMu,expParam)

mixParam = root.RooRealVar(channelName+"_mixParam","mixParam",0,1)

pdfMmumu = root.RooAddPdf("bak","bak",root.RooArgList(voitMmumu,expMmumu),root.RooArgList(mixParam))

# Just For Z-Peak Part

mMuMuZ = root.RooRealVar("mMuMu","mMuMu",88.0,94.0)
voitMmumuZ = root.RooVoigtian("bak_voitMmumuZ","voitMmumuZ",mMuMuZ,voitmZ,voitWidth,voitSig)
mMuMuZRooDataHist = root.RooDataHist("bak_TemplateZ","bak_TemplateZ",root.RooArgList(mMuMuZ),hist)

voitMmumuZ.fitTo(mMuMuZRooDataHist,root.RooFit.SumW2Error(False),PRINTLEVEL)
voitmZ.setConstant(True)
voitSig.setConstant(True)

# Back to everywhere else

mMuMuRooDataHist = root.RooDataHist("bak_Template","bak_Template",root.RooArgList(mMuMu),hist)

expMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("high"),root.RooFit.SumW2Error(False),PRINTLEVEL)
expParam.setConstant(True)

fr = pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
fr.SetName("bak"+"_fitResult")

#####################################################################


plotMmumu = mMuMu.frame(controlRegionLow[0],controlRegionHigh[1])

mMuMuRooDataHist.plotOn(plotMmumu)

pdfMmumu.plotOn(plotMmumu)
pdfMmumu.plotOn(plotMmumu,root.RooFit.Range("signal"),root.RooFit.LineStyle(2))

plotMmumu.Draw()

chi2ondf =  plotMmumu.chiSquare()
tlatex = root.TLatex()
tlatex.SetNDC()
tlatex.SetTextFont(root.gStyle.GetLabelFont())
tlatex.SetTextSize(0.05)
tlatex.SetTextAlign(22)
tlatex.DrawLatex(0.75,0.85,"#chi^{{2}}/NDF = {0:.2f}".format(chi2ondf))
canvas.SaveAs("testRes.png")

for i in [voitWidth,voitmZ,voitSig,expParam,mixParam]:
  print("{0}: {1:.3g} +/- {2:.3%}".format(i.GetName(),i.getVal(),abs(i.getError()/i.getVal())))
yay = raw_input("Press Enter to continue...")

