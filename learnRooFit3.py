#!/usr/bin/env python

from ROOT import gSystem

import sys
import os
import re
import math
from ROOT import *
gSystem.Load('libRooFit')
import ROOT as root

canvas = root.TCanvas()

"""
mMuMu = root.RooRealVar("mMuMu","mMuMu",110,160.0)
mMuMu.setRange("low",110,120)
#mMuMu.setRange("high",130,140)
mMuMu.setRange("high",130,160)
#mMuMu.setRange("low",110,124)
#mMuMu.setRange("high",126,140)
"""
mMuMu = root.RooRealVar("mMuMu","mMuMu",80,180.0)
mMuMu.setRange("z",88,94)
mMuMu.setRange("low",110,120)
mMuMu.setRange("high",130,180)

voitWidth = root.RooRealVar("voitWidth","voitWidth",2.4952)
voitmZ = root.RooRealVar("voitmZ","voitmZ",85,95)
voitSig = root.RooRealVar("voitSig","voitSig",0.0,30.0)
voitMmumu = root.RooVoigtian("voitMmumu","voitMmumu",mMuMu,voitmZ,voitWidth,voitSig)

expParam = root.RooRealVar("expParam","expParam",-1,0)
expMmumu = root.RooExponential("expMmumu","expMmumu",mMuMu,expParam)

expParam2 = root.RooRealVar("expParam2","expParam2",-1,0)
expMmumu2 = root.RooExponential("expMmumu2","expMmumu2",mMuMu,expParam2)

expParam3 = root.RooRealVar("expParam3","expParam3",-1,0)
expMmumu3 = root.RooExponential("expMmumu3","expMmumu3",mMuMu,expParam3)

mean2 = root.RooRealVar("mean2","mean2",0)
sig2 = root.RooRealVar("sig2","sig2",0.0,100.0)
gausMmumu = root.RooGaussian("gausMmumu","gausMmumu",mMuMu,mean2,sig2)

rooDoubleVoig = RooFFTConvPdf("rooDoubleVoig","rooDoubleVoig",mMuMu,voitMmumu,gausMmumu)

mixParam = root.RooRealVar("mixParam","mixParam",0,1)
mixParam2 = root.RooRealVar("mixParam2","mixParam2",0,1)
#pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",polyMmumu,expMmumu,mixParam)
#pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",expMmumu,expMmumu2,mixParam)
#pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",polyMmumu,polyMmumu2,mixParam)
#pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",expMmumu2,expMmumu3,mixParam)
#pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",root.RooArgList(voitMmumu,expMmumu2,expMmumu3),root.RooArgList(mixParam,mixParam2))
pdfMmumu = voitMmumu
pdfMmumu = rooDoubleVoig

f = root.TFile("input/DYJetsToLL.root")
#f = root.TFile("input/ttbar.root")

mDiMu = f.Get("mDiMu")
mDiMu.Rebin(2)

templateMmumu = root.RooDataHist("template","template",root.RooArgList(mMuMu),mDiMu)

sig2.setVal(0.1)
sig2.setConstant(True)
pdfMmumu.fitTo(templateMmumu,root.RooFit.Range("z"))
voitmZ.setConstant(True)
voitSig.setConstant(True)
sig2.setConstant(False)

#fitResult = pdfMmumu.fitTo(templateMmumu,root.RooFit.Range("low,high"),root.RooFit.Save(True))
fitResult = pdfMmumu.fitTo(templateMmumu,root.RooFit.Range("low,high"),root.RooFit.Save(True))
#pdfMmumu.fitTo(templateMmumu,root.RooFit.Range(130,160))
chi2 = pdfMmumu.createChi2(templateMmumu)

plotMmumu = mMuMu.frame()

templateMmumu.plotOn(plotMmumu)
pdfMmumu.plotOn(plotMmumu)
pdfMmumu.plotOn(plotMmumu,root.RooFit.Range(80,200),root.RooFit.LineStyle(2))
"""
expMmumu2.plotOn(plotMmumu,root.RooFit.LineColor(root.kGreen-1))
expMmumu3.plotOn(plotMmumu,root.RooFit.LineColor(root.kGreen-1))
expMmumu3.plotOn(plotMmumu,root.RooFit.LineColor(root.kGreen-1),root.RooFit.Range(80,200),root.RooFit.LineStyle(2))
expMmumu2.plotOn(plotMmumu,root.RooFit.LineColor(root.kGreen-1),root.RooFit.Range(80,200),root.RooFit.LineStyle(2))
expMmumu.plotOn(plotMmumu,root.RooFit.LineColor(root.kRed-1),root.RooFit.Range(80,200),root.RooFit.LineStyle(2))
"""
plotMmumu.Draw()

print fitResult
fitResult.Print()
print("chi2/ndf: {}".format(chi2.getVal()/(mDiMu.GetNbinsX()-1)))
yay = raw_input("Press Enter to continue...")
canvas.SaveAs("learnMmumu.png")

