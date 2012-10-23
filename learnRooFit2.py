#! /usr/bin/env python

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
mMuMu = root.RooRealVar("mMuMu","mMuMu",100,180.0)
mMuMu.setRange("low",115,120)
mMuMu.setRange("high",130,135)
mMuMu.setRange("signal",120,130)

expParam = root.RooRealVar("expParam","expParam",-1,0)
expMmumu = root.RooExponential("expMmumu","expMmumu",mMuMu,expParam)

expParam2 = root.RooRealVar("expParam2","expParam2",-1,0)
expMmumu2 = root.RooExponential("expMmumu2","expMmumu2",mMuMu,expParam2)

expParam3 = root.RooRealVar("expParam3","expParam3",-1,0)
expMmumu3 = root.RooExponential("expMmumu3","expMmumu3",mMuMu,expParam3)


a1 = root.RooRealVar("a1","a1",-5,5)
a2 = root.RooRealVar("a2","a2",-5,5)
a3 = root.RooRealVar("a3","a3",-5,5)
a4 = root.RooRealVar("a4","a4",-5,5)
a5 = root.RooRealVar("a5","a5",-5,5)
a6 = root.RooRealVar("a6","a6",-5,5)
a7 = root.RooRealVar("a7","a7",-5,5)
shift = root.RooRealVar("shift","shift",-110)
shiftedX = root.RooFormulaVar("shiftedX","shiftedX","mMuMu+shift",root.RooArgList(mMuMu,shift))

#polyMmumu = root.RooChebychev("polyMmumu","polyMmumu",mMuMu,root.RooArgList(a1,a2,a3))
polyMmumu = root.RooPolynomial("polyMmumu","polyMmumu",shiftedX,root.RooArgList(a1,a2))

b1 = root.RooRealVar("b1","b1",-50,50)
b2 = root.RooRealVar("b2","b2",-50,50)
b3 = root.RooRealVar("b3","b3",-50,50)
b4 = root.RooRealVar("b4","b4",-5,5)
b5 = root.RooRealVar("b5","b5",-5,5)
b6 = root.RooRealVar("b6","b6",-5,5)
b7 = root.RooRealVar("b7","b7",-5,5)
polyMmumu2 = root.RooChebychev("polyMmumu2","polyMmumu2",mMuMu,root.RooArgList(b1,b2,b3))
#polyMmumu2 = root.RooPolynomial("polyMmumu2","polyMmumu2",mMuMu,root.RooArgList(b1,b2,b3,b4,b5,b6,b7))

mixParam = root.RooRealVar("mixParam","mixParam",0,1)
mixParam2 = root.RooRealVar("mixParam2","mixParam2",0,1)
#pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",polyMmumu,expMmumu,mixParam)
#pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",expMmumu,expMmumu2,mixParam)
#pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",polyMmumu,polyMmumu2,mixParam)
#pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",expMmumu2,expMmumu3,mixParam)
#pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",root.RooArgList(expMmumu,expMmumu2,expMmumu3),root.RooArgList(mixParam,mixParam2))
pdfMmumu = polyMmumu

f = root.TFile("input/DYJetsToLL.root")
#f = root.TFile("input/ttbar.root")

mDiMu = f.Get("mDiMu")
#mDiMu.Rebin(2)

templateMmumu = root.RooDataHist("template","template",root.RooArgList(mMuMu),mDiMu)

fitResult = pdfMmumu.fitTo(templateMmumu,root.RooFit.Range("low,signal,high"),root.RooFit.Save(True))
chi2 = pdfMmumu.createChi2(templateMmumu)

plotMmumu = mMuMu.frame()

templateMmumu.plotOn(plotMmumu)
pdfMmumu.plotOn(plotMmumu)
a1val = a1.getVal()
a2val = a2.getVal()
a1Err = a1.getError()
a2Err = a2.getError()
a1.setVal(a1val+a1Err)
pdfMmumu.plotOn(plotMmumu,root.RooFit.LineColor(root.kGreen-1),root.RooFit.Range("low,signal,high"),root.RooFit.LineStyle(2))
a1.setVal(a1val-a1Err)
pdfMmumu.plotOn(plotMmumu,root.RooFit.LineColor(root.kGreen-1),root.RooFit.Range("low,signal,high"),root.RooFit.LineStyle(2))
plotMmumu.Draw()
a1.setVal(a1val)
a2.setVal(a2val+a2Err)
pdfMmumu.plotOn(plotMmumu,root.RooFit.LineColor(root.kRed-1),root.RooFit.Range("low,signal,high"),root.RooFit.LineStyle(2))
a2.setVal(a2val-a2Err)
pdfMmumu.plotOn(plotMmumu,root.RooFit.LineColor(root.kRed-1),root.RooFit.Range("low,signal,high"),root.RooFit.LineStyle(2))
plotMmumu.Draw()

print fitResult
fitResult.Print()
print("chi2/ndf: {}".format(chi2.getVal()/(mDiMu.GetNbinsX()-1)))
yay = raw_input("Press Enter to continue...")
canvas.SaveAs("learnMmumu.png")

