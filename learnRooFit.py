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

mMuMu = root.RooRealVar("mMuMu","mMuMu",110.0,150.0)
mMuMu.setRange("low",110,120)
mMuMu.setRange("high",130,145)
#mva = root.RooRealVar("mva","mva",-1,1)
mva = root.RooRealVar("mva","mva",-0.8,-0.1)

a0 = root.RooRealVar("a0","a0",-10.0,10.0)
a1 = root.RooRealVar("a1","a1",-10.0,10.0)
a2 = root.RooRealVar("a2","a3",-10.0,10.0)
a3 = root.RooRealVar("a3","a3",-10.0,10.0)
a4 = root.RooRealVar("a4","a4",-10.0,10.0)

polyArgsMmumu = root.RooArgList(a0,a1,a2,a3)#,a4)
polyMmumu = root.RooPolynomial("polyFunc","polyFunc",mMuMu,polyArgsMmumu)

f = root.TFile("input/DYJetsToLL.root")

mDiMu = f.Get("mDiMu")

templateMmumu = root.RooDataHist("template","template",root.RooArgList(mMuMu),mDiMu)

polyMmumu.fitTo(templateMmumu,root.RooFit.Range("low,high"))

###################

BDTHistMuonOnlyVMass = f.Get("BDTHistMuonOnlyVMass")
tmpHist = BDTHistMuonOnlyVMass.Clone("tmpHist")
tmpAxis = tmpHist.GetXaxis()
BDTHistLowControl = tmpHist.ProjectionY("_LowControl",tmpAxis.FindBin(110.0),tmpAxis.FindBin(115.0))
BDTHistHighControl = tmpHist.ProjectionY("_HighControl",tmpAxis.FindBin(125.0),tmpAxis.FindBin(150.0))

print("N Events in Low Control: {0:.1f}, in High Control: {1:.1f}".format(BDTHistLowControl.Integral(),BDTHistHighControl.Integral()))

BDTHist = BDTHistLowControl.Clone()
BDTHist.Add(BDTHistHighControl)

BDTHist.Rebin(BDTHist.GetNbinsX()/100)

templateMva = root.RooDataHist("templateMva","templateMva",root.RooArgList(mva),BDTHist)

pdfMva = root.RooHistPdf("pdfMva","pdfMva",root.RooArgSet(mva),templateMva)

pdfMva.fitTo(templateMva)

###################

plotMmumu = mMuMu.frame()

templateMmumu.plotOn(plotMmumu)
polyMmumu.plotOn(plotMmumu)
plotMmumu.Draw()
canvas.SaveAs("learnMmumu.png")

plotMva = mva.frame()

templateMva.plotOn(plotMva)
pdfMva.plotOn(plotMva)
plotMva.Draw()
canvas.SaveAs("learnMva.png")

###################

canvas.Clear()
canvas.Divide(2,2)

multipdf = root.RooProdPdf("multipdf","multipdf",RooArgList(polyMmumu,pdfMva))

mDiMuBinning = root.RooFit.Binning(50,110,150)
yvar = root.RooFit.YVar(mva,root.RooFit.Binning(50,-1,1))
print mDiMu
print mDiMuBinning
print yvar
multipdfHist = multipdf.createHistogram("hh_model_1",mMuMu,mDiMuBinning,yvar)

canvas.cd(3)

multipdfHist.GetXaxis().SetRangeUser(110,150)
multipdfHist.Draw("surf")
BDTHistMuonOnlyVMass.SetLineColor(kRed)
BDTHistMuonOnlyVMass.Rebin2D(BDTHistMuonOnlyVMass.GetNbinsX()/400,BDTHistMuonOnlyVMass.GetNbinsY()/50)
BDTHistMuonOnlyVMass.GetXaxis().SetRangeUser(110,150)
BDTHistMuonOnlyVMass.Draw("surf same")

#canvas.SaveAs("learn2D.png")

BDTHistMuonOnlyVMass.Draw("surf")
multipdfHist.GetXaxis().SetRangeUser(110,150)
multipdfHist.Draw("surf same")

canvas.cd(1)
BDTHistMuonOnlyVMass.Draw("surf")
multipdfHist.Draw("surf same")

canvas.cd(2)
BDTHistMuonOnlyVMass.Draw("colz")

canvas.cd(4)
multipdfHist.Draw("colz")

canvas.SaveAs("learn2D.png")

