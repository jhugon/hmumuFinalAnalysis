#! /usr/bin/env python
import sys
import os
import re
import math
from ROOT import *
import ROOT as root

root.gSystem.Load("libRooFit")

canvas = root.TCanvas()

mMuMu = root.RooRealVar("mMuMu","mMuMu",110.0,150.0)
#mva = root.RooRealVar("mva","mva",-1,1)
mva = root.RooRealVar("mva","mva",-0.8,-0.1)

a0 = root.RooRealVar("a0","a0",-10.0,10.0)
a1 = root.RooRealVar("a1","a1",-10.0,10.0)
a2 = root.RooRealVar("a2","a2",-10.0,10.0)
a3 = root.RooRealVar("a3","a3",-10.0,10.0)

polyArgsMmumu = root.RooArgList(a0,a1,a2)
polyMmumu = root.RooPolynomial("polyFunc","polyFunc",mMuMu,polyArgsMmumu)

f = root.TFile("input/DYJetsToLL.root")

mDiMu = f.Get("mDiMu")

templateMmumu = root.RooDataHist("template","template",root.RooArgList(mMuMu),mDiMu)

polyMmumu.fitTo(templateMmumu)

###################

BDTHistMuonOnly = f.Get("BDTHistMuonOnly")
BDTHistMuonOnly.Rebin(10)

templateMva = root.RooDataHist("templateMva","templateMva",root.RooArgList(mva),BDTHistMuonOnly)

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

BDTHistMuonOnlyVMass = f.Get("BDTHistMuonOnlyVMass")

multipdf = root.RooProdPdf("multipdf","multipdf",RooArgList(polyMmumu,pdfMva))
